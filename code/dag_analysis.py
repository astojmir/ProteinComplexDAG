#
# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
#
# Code author:  Aleksandar Stojmirovic
#

import os
import sys
import numpy as np
from collections import defaultdict
from collections import namedtuple
from collections import OrderedDict
from collections import Counter
from operator import itemgetter
from itertools import chain
from itertools import imap
import cPickle as pickle
import base64
from qmbpmn.common.utils.filesys import makedirs2

from utils import load_results
from utils import participation_ratio

from dictgraph import DictGraph

from dag_node import DAGNode
from dag_component import DAGComponent

import infflow
import enrichment


_rngf = lambda a, b: (lambda o: a<= o.pratio < b)


class DAGAnalysis(object):

    ssum_opts = [#('hgem', '-shgem -w1e-07'),
                 ('avgw', None),
                 ('wsum', None),
                 ]
    component_prefix = 'K'
    subcomponent_sep = '-'
    cut_prefix = 'T'

    _stats_funcs = [('TOTAL', lambda o: True),
                    ('MERGED', lambda o: o.name[0] == 'M'),
                    ('FIXED', lambda o: o.name[0] == 'F'),
                    ('pr(2.5)', _rngf(0.0, 2.5)),
                    ('pr(4.0)', _rngf(2.5, 4.0)),
                    ('pr(6.0)', _rngf(4.0, 6.0)),
                    ('pr(9.0)', _rngf(6.0, 9.0)),
                    ('pr(15.0)', _rngf(9.0, 15.0)),
                    ('pr(30.0)', _rngf(15.0, 30.0)),
                    ('LARGE', _rngf(30.0, 10000.0))]

    _nattr_funcs = \
    [(lambda obj: 'K' if obj.name[0] == obj.component_prefix \
      else obj.node_type + 'c', '{}', 'NodeType'),
     (lambda obj: len(obj.nodes) if hasattr(obj, 'nodes') else 1,
      '{}', 'NumNodes'),
     (lambda obj: obj.name if obj.name[0] == obj.component_prefix  else '',
      '{}', 'NodeLabel'),
     ]


    def __init__(self, dag_filename):

        self.dag_filename = os.path.abspath(dag_filename)
        self.nodes = None
        self.components = None

        self.enrichment_dir = None
        self.termdb_file = None
        self.saddlesum_cmd_path = None

        self.component_counter = None
        self.scale = None
        self.p2c_counts = None
        self._nr2complex = None

        data = load_results(self.dag_filename)
        self.__dict__.update(data)
        self._set_enrichment_params()

    def __getitem__(self, name):

        if name in self.nodes:
            return self.nodes[name]
        if name in self.components:
            return self.components[name]
        raise KeyError(name)

    def __contains__(self, name):
        return name in self.nodes or name in self.components

    @property
    def nr2complex(self):
        if not hasattr(self, '_nr2complex') or self._nr2complex is None:
            self._nr2complex = dict((nr_id, complex_id)
                                    for complex_id in self.complex2nr
                                    for nr_id in self.complex2nr[complex_id])
        return self._nr2complex

    @property
    def avg_ppi_path_length(self):

        G = infflow.get_pp_graph(self.pp_edges, self.min_pp_count, self.beta,
                                 self.adb)
        SPL = infflow.get_pp_solver(G)
        return infflow.get_avg_pp_path_length(SPL)[0]

    def get_original_components(self, min_size=2):
        return [cmpn for cmpn in self.components.itervalues()
                if len(cmpn.nodes) >= min_size]

    def get_active_componets(self, min_size=2):
        return self.get_original_components(min_size)


    # ************************************************************************
    # Construction
    # ************************************************************************

    def initial_processing(self, pratio_cutoff=2.5, pubmed_ids=None,
                           component_prefix='K'):

        self._filtered_small_complexes = list()
        self._filtered_unlabeled_complexes = list()

        self.component_prefix = component_prefix
        self._filter_small_complexes(pratio_cutoff)
        if pubmed_ids is not None:
            self._filter_unlabeled_complexes(pubmed_ids)
        self._apply_filters()
        self._init_nodes()
        self._init_components()

    def _filter_small_complexes(self, pratio_cutoff=2.5):

        self.pratio_cutoff = pratio_cutoff
        excluded_complexes = \
          [self.complexes[i] for i, W in enumerate(self.itm_weights) \
           if participation_ratio(W[0]) < pratio_cutoff]
        self._filtered_small_complexes = excluded_complexes
        self.excluded_complexes = excluded_complexes

    def _filter_unlabeled_complexes(self, pubmed_ids):
        """
        Remove all nodes that have no baits and the only pubmed id is in
        pubmed_ids.
        """

        unlabeled = list()
        selected_pubmeds = set(pubmed_ids)
        for complex_id in self.complexes:
            evd = self.adb[complex_id]
            num_baits = sum(1 for p, v in evd.proteins.iteritems()
                            if v[0] == evd.bait_code)
            remove = all([num_baits == 0,
                         len(evd.pubmed_ids) == 1,
                         evd.pubmed_ids[0] in selected_pubmeds])
            if remove:
                unlabeled.append(complex_id)

        self._filtered_unlabeled_complexes = unlabeled
        self.excluded_complexes += unlabeled

    def _apply_filters(self, verbose=True):

        G = self.graph
        self.filtered_graph = G.filtered_graph(self.excluded_complexes)
        if verbose:
            original = G.count_edges()
            retained = self.filtered_graph.count_edges()
            print "Reducing DAG"
            print "  * Removed for small pratio: %d" % \
                  len(self._filtered_small_complexes)
            print "  * Removed from selected pubmeds: %d" % \
                  len(self._filtered_unlabeled_complexes)
            print "  * Original edges: %d" % original
            print "  * Retained after filtering: %d" % retained

    def _init_nodes(self):

        if self.beta is not None:
            self.scale = self.min_pc_count / (1 - self.beta)
        else:
            self.scale = 1.0
        self.p2c_counts = defaultdict(int)
        for complex_id in self.complexes:
            for p in self.adb[complex_id].proteins:
                self.p2c_counts[p] += 1
        node_types = self.filtered_graph.classify_nodes(self.complexes,
                                                        self.excluded_complexes)
        self.nodes = dict()
        for complex_id in self.complexes:
            node = DAGNode(self, complex_id, node_types[complex_id],
                           enrich_config=self)
            self.nodes[complex_id] = node

    def _init_components(self):

        self.components = OrderedDict()
        self.component_counter = dict()
        component_list = self.filtered_graph.connected_components()
        component_list.reverse()
        for i, cmpn in enumerate(component_list):
            component_id = str.format('{}{:0=3d}', self.component_prefix, i)
            nodes = [self.nodes[complex_id] for complex_id in cmpn]
            self.components[component_id] = \
              DAGComponent(self, component_id, nodes, enrich_config=self)
            self.component_counter[component_id] = 0

    def _set_enrichment_params(self):
        res_dir, f = os.path.split(self.dag_filename)
        self.enrichment_dir = os.path.join(res_dir, f.split('.')[0] + '.enrich')

    # ************************************************************************
    # Load / Save
    # ************************************************************************

    @staticmethod
    def load(input_file):
        with open(input_file, 'rb') as fp:
            obj = pickle.load(fp)
        obj.dag_filename = os.path.abspath(input_file)
        obj._set_enrichment_params()
        return obj

    def save(self, output_file=None):

        if output_file is None:
            output_file = self.dag_filename
        else:
            self.dag_filename = output_file
        with open(output_file, 'wb') as fp:
            pickle.dump(self, fp, 2)

    # ************************************************************************
    # Enrichment - run all
    # ************************************************************************

    def run_enrichment(self, termdb_file, cmd_path=None, verbose=True):

        self.termdb_file = os.path.abspath(termdb_file)
        self.saddlesum_cmd_path = cmd_path
        objs = self.nodes.values() + self.components.values()
        for k, obj in enumerate(objs):
            obj.run_saddlesum()
            if verbose and (k+1) % 500 == 0:
                msg_fmt = 'Completed {:d} / {:d} saddlesum runs.'
                print msg_fmt.format(k+1, len(objs))

    # ************************************************************************
    # Non transitive links
    # ************************************************************************

    @property
    def non_transitive_relative_distances(self):
        return [item for cmpn in self.get_original_components() \
                     for item in cmpn.non_transitive_relative_distances]

    def get_essential_non_transitive_links(self):

        R = self.filtered_graph.transitive_reduction_graph()[0]
        ntrd0 = self.non_transitive_relative_distances
        ntrd1 = []
        for x, y, z, dxy, dyz, dxz in ntrd0:
            if y in R[x] and z in R[y]:
                ntrd1.append((x, y, z, dxy, dyz, dxz))
        ntrd1.sort(key=itemgetter(5))
        return ntrd1

    # ************************************************************************
    # Reports
    # ************************************************************************

    def report_components(self, components=None, fp=sys.stdout):

        if components is None:
            components = self.components.values()
        for cmpn in components:
            fp.write(cmpn.report_line(self.pp_edges))
            fp.write('\n')

    def report_nodes(self, nodes=None, fp=sys.stdout):

        if nodes is None:
            nodes = sorted(self.nodes.itervalues(),
                           key=lambda node: int(node.node_id[1:]))
        for node in nodes:
            fp.write(node.report_summary_line(self.beta, self.enrichment_dir,
                     [tag for tag, _ in self.ssum_opts]))
            fp.write('\n')

    def report_enrichment(self, obj_id, suffix, fp=sys.stdout):

        if obj_id in self.nodes:
            obj = self.nodes[obj_id]
        else:
            obj = self.components[obj_id]
        obj.report_enrichment(suffix, fp)

    def report_details(self, obj_id, scale=None, fp=sys.stdout):

        if scale is None:
            scale = self.scale
        if obj_id in self.nodes:
            obj = self.nodes[obj_id]
        else:
            obj = self.components[obj_id]
        obj.report_details(fp, scale=scale)

    def _collect_stats(self, node_types=None):

        if node_types is None:
            nodes = self.nodes.values()
        else:
            ntypes = set(node_types)
            nodes = [node for node in self.nodes.itervalues() \
                     if node.node_type in ntypes]
        row = [sum(imap(f, nodes)) for  _, f in self._stats_funcs]
        return row

    @staticmethod
    def _print_stats_rows(stats, fp):
        for head, row in stats:
            fp.write('{:<30.30}'.format(head))
            for x in row:
                fp.write(' {!s:>9.9}'.format(x))
            fp.write('\n')

    def report_stats(self, fp=sys.stdout):

        nr_sizes = map(len, self.complex2nr.itervalues())
        avg_pratio = sum(node.pratio for node in self.nodes.itervalues()) \
                     / len(self.nodes)

        header = [h for  h, _ in self._stats_funcs]
        stats = \
        [('Total raw evidences', [len(self.adb.raw_evidences)]),
         ('Filtered evidences', [len(self.adb.removed)]),
         ('Nonredundant complexes', [sum(nr_sizes)]),
         ('Total nodes', self._collect_stats()),
         ('Removed (small pratio)', self._collect_stats(['x'])),
         ('Singletons', self._collect_stats(['S'])),
         ('Connected nodes', self._collect_stats(['m', 'i', 'M'])),
         ('Minima', self._collect_stats(['m'])),
         ('Maxima', self._collect_stats(['M'])),
         ('Inner nodes', self._collect_stats(['i'])),
         ('Edges', [self.filtered_graph.count_edges()]),
         ('Clustering edges', [sum(n*(n-1)/2 for n in nr_sizes)]),
         ('Non-transitive edges',
          [len(self.non_transitive_relative_distances)]),
         ('Average pratio', ['{:9.2f}'.format(avg_pratio)]),
         ('Components', [len(self.get_original_components())]),
         ]

        self._print_stats_rows([['', header]], fp)
        fp.write('{0:{2}^{1}}\n'.format('', len(header)*10+30, '-'))
        self._print_stats_rows(stats, fp)

    # ************************************************************************
    # Component relationship graphs
    # ************************************************************************

    def write_component_graph_data(self, component_id, file_prefix,
                                   max_depth=100, min_component_size=20):

        edges_filename = file_prefix + '.edges.tab'
        nattr_filename = file_prefix + '.nattr.tab'

        S = self.components[component_id].subcomponent_graph(max_depth,
                                                             True,
                                                             min_component_size)
        with open(edges_filename, 'w') as fp:
            S.write_edges(fp)
        with open(nattr_filename, 'w') as fp:
            S.write_node_attrs(self, self._nattr_funcs, fp)

    def write_component_graph_data2(self, component_id, file_prefix,
                                    max_depth=100, min_component_size=20):

        edges_filename = file_prefix + '.edges.tab'
        nattr_filename = file_prefix + '.nattr.tab'

        S = self.components[component_id].subcomponent_graph(max_depth,
                                                             True,
                                                             min_component_size)
        U = S.U
        A = DictGraph()
        A.default_factory = lambda : defaultdict(int)
        # S should be bipartite graph component - complex - component
        for src in U.keys():
            if src[0] == self.component_prefix:
                for cut_node in U[src]:
                    for tgt in U[cut_node].keys():
                        if tgt[0] == self.component_prefix and src < tgt:
                            A[src][tgt] += 1

        with open(edges_filename, 'w') as fp:
            A.write_edges(fp)
        with open(nattr_filename, 'w') as fp:
            A.write_node_attrs(self, self._nattr_funcs, fp)

    def write_protein_graph_data(self, prot, file_prefix,
                                 show_filtered_nodes=False):

        edges_filename = file_prefix + '.edges.tab'
        nattr_filename = file_prefix + '.nattr.tab'
        nattr_funcs = DAGComponent.get_nattr_funcs([])

        shown_nodes = [node.node_id for node in self.nodes.itervalues()
                       if prot in node.evd.proteins]
        if show_filtered_nodes:
            G = self.graph.subgraph(shown_nodes)
        else:
            G = self.filtered_graph.subgraph(shown_nodes)

        R = G.transitive_reduction_graph()[0]
        with open(edges_filename, 'w') as fp:
            R.write_edges(fp)
        with open(nattr_filename, 'w') as fp:
            R.write_node_attrs(self, nattr_funcs, fp)


