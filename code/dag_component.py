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

import sys
import os
from collections import defaultdict
from collections import Counter
from operator import itemgetter
from enrichment import EnrichmentMixin
from dictgraph import DictGraph
from cmembers import Association


class DAGComponent(EnrichmentMixin):
    # Collection of nodes

    @staticmethod
    def get_nattr_funcs(*args, **kwargs):

        def _node_type_func(obj):
            return obj.node_type

        return  [(_node_type_func, '{}', 'NodeType'),
                 (lambda obj: obj.pratio, '{:.2f}', 'Pratio'),
                 ]

    def __init__(self, dag, component_id, nodes, **kwargs):

        EnrichmentMixin.__init__(self, component_id, **kwargs)

        self.dag = dag
        self.component_id = component_id
        self.nodes = nodes
        self._graph = None

        self.proteins = self._init_proteins()
        self.complexes = sorted(node.node_id for node in self.nodes)

    def __repr__(self):
        return str.format("DAGComponent({})", self.component_id)

    def __len__(self):
        return len(self.nodes)

    @property
    def title(self):
        return 'COMPONENT {}'.format(self.component_id)

    def _init_proteins(self):

        proteins = defaultdict(lambda : [Association.unknown_code, 0.0])
        wsum = 0.0
        for node in self.nodes:
            evdprots = node.evd.proteins
            weight = len(node.evd.evidences)  # This may change later
            for p in evdprots:
                proteins[p][0] = Association.max_code(proteins[p][0],
                                                      evdprots[p][0])
                proteins[p][1] += weight * evdprots[p][1]
            wsum += weight
        for p in proteins:
            proteins[p][1] /= wsum
        proteins.default_factory = None
        return proteins

    def find_node(self, node_id):
        found = [node for node in self.nodes if node.node_id == node_id]
        return found[0] if found else None


    @property
    def graph(self):
        if self._graph is None:
            self._graph = self.dag.filtered_graph.subgraph(self.complexes)
        return self._graph

    @property
    def non_transitive_relative_distances(self):

        G = self.graph
        E, num_added, added = G.transitive_closure_graph()
        added_dists = []
        for x, y, z in sorted(added):
            xnode = self.find_node(x)
            ynode = self.find_node(y)
            znode = self.find_node(z)
            if y in G[x]:
                dxy = G[x][y]
            else:
                dxy = xnode.rdist(ynode)
            if z in G[y]:
                dyz = G[y][z]
            else:
                dyz = ynode.rdist(znode)
            dxz = xnode.rdist(znode)
            added_dists.append((x, y, z, dxy, dyz, dxz))
        added_dists.sort(key=itemgetter(-1))
        return added_dists

    @property
    def raw_evidences(self):
        evidences = [evd2 for node in self.nodes \
                          for evd2 in node.evd.evidences]
        return evidences

    @property
    def pubmed_ids(self):
        pmids = set(pmid for node in self.nodes \
                         for pmid in node.evd.pubmed_ids)
        return sorted(pmids)

    @property
    def bait_evidences(self):
        bait_evidences = [evd2 for evd2 in self.raw_evidences \
                          if evd2.has_baits]
        return bait_evidences

    @property
    def distinct_baits(self):
        baits = sorted(p for p, v in self.proteins.iteritems() \
                       if v[0] == Association.bait_code)
        return baits

    @property
    def node_types(self):
        node_types = [node.node_type for node in self.nodes]
        return node_types

    def get_inner_pplinks(self):
        matches = set(pair for node in self.nodes \
                           for pair in node.get_inner_pplinks())
        return matches

    def _get_weights_wsum(self):
        # Produces maximal cover of the entire component
        weights = defaultdict(float)
        for node in self.nodes:
            for p, x in zip(node.itm_prots, node.itm_vals):
                if x > weights[p]:
                    weights[p] = x
        return weights.items()

    def _get_weights_avgw(self):
        # Uses 'plain' weights
        weights = [(p, v[1]) for p, v in self.proteins.iteritems()]
        return weights

    @property
    def base_name(self):
        return self.component_id.split(self.dag.subcomponent_sep)[0]

    def subcomponent_id(self, k):
        fmt = '{}{}{:0=2d}'
        return fmt.format(self.base_name, self.dag.subcomponent_sep, k)


    def report_data(self):

        matches = self.get_inner_pplinks()
        pulldown_pubmeds = set(self.pubmed_ids)
        pplinks_pubmeds = set(pmid for pair in matches \
                              for pmid in self.dag.pp_edges[pair])
        all_pubmeds = pulldown_pubmeds | pplinks_pubmeds
        node_type_counter = Counter(self.node_types)

        all_nodes = ', '.join(self.complexes)
        all_baits = ', '.join(self.distinct_baits)

        data = \
          [(self.component_id, '{}'),
           (len(self.nodes), '{:6d}'),
           (node_type_counter['m'], '{:6d}'),
           (node_type_counter['i'], '{:6d}'),
           (node_type_counter['M'], '{:6d}'),
           (len(self.non_transitive_relative_distances), '{:6d}'),
           (len(self.proteins), '{:6d}'),
           (len(pulldown_pubmeds), '{:6d}'),
           (len(self.raw_evidences), '{:6d}'),
           (len(self.distinct_baits), '{:6d}'),
           (len(matches), '{:6d}'),
           (len(pplinks_pubmeds), '{:6d}'),
           (len(all_pubmeds), '{:6d}'),
           (', '.join(self.complexes), '{}'),
           (', '.join(self.proteins), '{}'),
           (', '.join(self.distinct_baits), '{}'),
           ]
        return data

    def report_details(self, fp=sys.stdout, show_unused=True, scale=None):

        if show_unused:
            unused_proteins = set(p for node in self.nodes
                                  for p in node.evd.proteins
                                  if p not in node.itm_prots)
        else:
            unused_proteins = set()
        if scale is None:
            scale = self.dag.scale

        itm_weights = defaultdict(float)
        for node in self.nodes:
            for p, x in zip(node.itm_prots, node.itm_vals):
                if x > itm_weights[p]:
                    itm_weights[p] = x

        data = list()
        for p, x in itm_weights.iteritems():
            data.append((p, scale*x) + tuple(self.proteins[p]))
        data.sort(key=itemgetter(3,0), reverse=True)

        diam = self.graph.undirected_diameter()
        member_nodes = ', '.join(node.node_id for node in self.nodes)

        pmids = ', '.join(self.pubmed_ids)
        evidences = ', '.join(evd.complex_id for evd in self.raw_evidences)

        fp.write("**** {}\n".format(self.title))
        fp.write("undirected_diameter\t{:d}\n".format(int(diam)))
        fp.write("num_member_nodes\t{:d}\n".format(len(self.nodes)))
        fp.write("member_nodes\t{}\n".format(member_nodes))
        fp.write("num_proteins\t{:d}\n".format(len(self.proteins)))
        fp.write("#\n")
        fp.write("num_raw_evidences\t{:d}\n".format(len(self.raw_evidences)))
        fp.write("raw_evidences\t{}\n".format(evidences))
        fp.write("num_pubmed_ids\t{:d}\n".format(len(self.pubmed_ids)))
        fp.write("pubmed_ids\t{}\n".format(pmids))
        fp.write("%\n")

        for p, x, b, y in data:
            fp.write('{}\t{:.2f}\t{}\t{:.2f}\n'.format(p, y, b, x))
        for p in sorted(unused_proteins):
            b = self.evd.proteins[p][0]
            fp.write('{}\tNone\t{}\n'.format(p, b))
        fp.write('----\n')


    def write_node_data(self, output_dir, prefix=''):

        edges_file = os.path.join(output_dir,
                                  prefix + self.component_id + '.edges.tab')
        nattr_file = os.path.join(output_dir,
                                  prefix + self.component_id + '.nattr.tab')
        G = DictGraph(self.graph)
        R = G.transitive_reduction_graph()[0]
        with open(edges_file, 'w') as fp:
            R.write_edges(fp)
        with open(nattr_file, 'w') as fp:
            R.write_node_attrs(self.dag, self.get_nattr_funcs(), fp)
