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
import numpy as np
from operator import itemgetter
from utils import participation_ratio
from enrichment import EnrichmentMixin
from dictgraph import DictGraph


class DAGNode(EnrichmentMixin):

    def __init__(self, dag, complex_id, node_type, **kwargs):

        EnrichmentMixin.__init__(self, complex_id, **kwargs)

        self.dag = dag
        self.node_id = complex_id
        self.node_type = node_type

        k = self.dag.complexes.index(complex_id)

        self.itm_vals = self.dag.itm_weights[k][0]
        self.itm_prots = [self.dag.proteins[j]
                          for j in self.dag.itm_weights[k][1]]
        self.evd = self.dag.adb[complex_id]
        self.pratio = participation_ratio(self.itm_vals)
        self.plain_weight = self.evd.total_weight()

    def __repr__(self):
        return str.format("DAGNode({})", self.node_id)

    def __len__(self):
        return 1

    @property
    def complexes(self):
        return [self.node_id]

    @property
    def title(self):
        return 'NODE {}'.format(self.node_id)

    def rdist(self, other):

        xlst = zip(self.itm_prots, self.itm_vals)
        zmap = dict(zip(other.itm_prots, other.itm_vals))
        sxx = self.itm_vals.sum()
        sxz = sum(min(x, zmap[p]) for p, x in xlst if p in zmap)
        dxz = 1.0 - sxz / sxx
        return dxz

    def get_inner_pplinks(self):

        prots = self.evd.proteins.keys()
        pairs = set(frozenset((p1, p2)) for i, p1 in enumerate(prots) \
                                        for p2 in prots[i+1:])
        matches = set(pair for pair in pairs if pair in self.dag.pp_edges)
        return matches

    def get_components(self):
        """ Return all components that contain this node """

        found = [cmpn.name for cmpn in self.dag.components.itervalues()
                 if cmpn.find_node(self.node_id) is not None]
        return found

    def _extended_protein_name(self, prot):
        return '{}-{}'.format(self.name, prot)

    def save_pp_graph(self, file_prefix):

        edges_file = file_prefix + '.edges.tab'
        nodes_file = file_prefix + '.nodes.tab'
        scale = self.dag.scale
        pp_links = self.get_inner_pplinks()
        PP = DictGraph()
        for p1 in self.evd.proteins:
            q1 = self._extended_protein_name(p1)
            PP[q1] = dict()
        for p1, p2 in (tuple(edge) for edge in pp_links):
            q1 = self._extended_protein_name(p1)
            q2 = self._extended_protein_name(p2)
            PP[q1][q2] = 1.0

        with open(edges_file, 'w') as fp:
            PP.write_edges(fp, edge_type='pp')

        wghts = dict(zip(self.itm_prots, self.itm_vals))
        with open(nodes_file, 'w') as fp:
            fp.write('ID\tcanonicalName\tWeight\tWidth\n')
            for p1 in self.itm_prots:
                q1 = self._extended_protein_name(p1)
                line = '{}\t{}\t{:.6f}\t{:.6f}\n'.format(q1, p1, scale*wghts[p1],
                                                         np.sqrt(scale*wghts[p1]))
                fp.write(line)

    def _get_weights_wsum(self):

        weights = zip(self.itm_prots, self.itm_vals)
        return weights

    def _get_weights_avgw(self):
        # Uses 'plain' weights
        weights = [(p, v[1]) for p, v in self.evd.proteins.iteritems()]
        return weights

    def get_summary_data(self, beta):

        i = self.itm_vals.argmax()

        itm_mass = self.itm_vals.sum()
        max_weight_prot = self.itm_prots[i]
        max_weight_erole = self.evd.proteins[max_weight_prot][0]
        max_relative_weight = self.itm_vals[i] / itm_mass
        max_weight_p2c = self.dag.p2c_counts[max_weight_prot]

        num_pubmeds = len(set(self.evd.pubmed_ids))
        num_baits = sum(1 for p, v in self.evd.proteins.iteritems() \
                        if v[0] == self.evd.bait_code)

        beta0_vals = [v[1] / self.dag.p2c_counts[p] \
                       for p, v in self.evd.proteins.iteritems() \
                       if v[1] > 0]
        beta0_weight = sum(beta0_vals)
        pratio0 = participation_ratio(np.array(beta0_vals))
        b = 1.0 - beta
        itm_mass_growth = int(100 * (itm_mass / b / beta0_weight - 1.0))
        itm_ppi_infl = int(100 * (1.0 - (b * beta0_weight / itm_mass)))
        diff_pratio = self.pratio - pratio0
        rel_diff_pratio = abs(diff_pratio) / self.pratio

        prot_ = str.format('{}({})', max_weight_prot, max_weight_erole)
        data = ((self.node_id, '{}'),
                (num_pubmeds, '{:6d}'),
                (num_baits, '{:6d}'),
                (self.plain_weight, '{:6.2f}'),
                # (itm_mass, '{:8.4f}'),
                (prot_, '{:14.14s}'),
                (max_relative_weight, '{:6.2f}'),
                (max_weight_p2c, '{:5d}'),
                (self.pratio, '{:6.2f}'),
                (self.node_type, '{}'),
                (itm_ppi_infl, '{:6d}'),
                # (itm_mass_growth, '{:6d}'),
                # (diff_pratio, '{:6.2f}'),
                # (rel_diff_pratio, '{:6.2f}'),
                )
        return data

    def report_details(self, fp=sys.stdout, show_unused=True, scale=None):

        if show_unused:
            unused_proteins = set(self.evd.proteins) - set(self.itm_prots)
        else:
            unused_proteins = set()
        if scale is None:
            scale = self.dag.scale

        data = list()
        for x, p in zip(self.itm_vals, self.itm_prots):
            data.append((p, scale*x) + tuple(self.evd.proteins[p]))
        data.sort(key=itemgetter(1,0), reverse=True)

        pmids = ', '.join(self.evd.pubmed_ids)
        evidences = ', '.join(evd.complex_id for evd in self.evd.evidences)

        fp.write("**** {}\n".format(self.title))
        fp.write("node_type\t{}\n".format(self.node_type))
        fp.write("num_proteins\t{:d}\n".format(len(self.evd.proteins)))
        fp.write("evidence_weight\t{:.2f}\n".format(self.plain_weight))
        fp.write("itm_weight\t{:.2f}\n".format(scale * self.itm_vals.sum()))
        fp.write("itm_participation_ratio\t{:.2f}\n".format(self.pratio))
        fp.write("#\n")
        fp.write("raw_evidences\t{}\n".format(evidences))
        fp.write("pubmed_ids\t{}\n".format(pmids))
        fp.write("%\n")

        for p, x, b, y in data:
            fp.write('{}\t{:.2f}\t{}\t{:.2f}\n'.format(p, y, b, x))
        for p in sorted(unused_proteins):
            b = self.evd.proteins[p][0]
            fp.write('{}\tNone\t{}\n'.format(p, b))

        fp.write('----\n')
