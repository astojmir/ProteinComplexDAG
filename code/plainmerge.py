#! /usr/bin/env python
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
import getopt
import numpy as np
from qmbpmn.common.utils.filesys import makedirs2
from cmembers import AssociationDatabase
import itermerge


_HELP = \
"""
SYNOPSIS:

    %s [OPTIONS] label cmembers_file alpha
    %s -h|--help

OPTIONS:

    -h, --help                  Print this message and exit

""" % (__file__, __file__)

_EMSG = "Insufficient arguments.\n\n" + _HELP



class PlainAssocCluster(itermerge.AssocCluster):

    def compute_weights(self):

        self.x = self.v
        self.x_val = self.x[self.jix]
        self.s = self.x_val.sum()


class PlainAssocDataset(itermerge.AssocDataset):

    def __init__(self, cmembers_file, alpha):

        # Most attributes are unused
        self.num_merges = 0
        self.merge_order = []

        self.pplinks_file = None
        self.score_cutoff = None
        self.cmembers_file = cmembers_file

        self.alpha = alpha
        self.beta = None
        self.rtol = None
        self.weight_range = itermerge.ITM_WEIGHTS_MEMBERS

        self.min_pc_count = None
        self.min_pp_count = None
        self.shuffle = False

        self.clusters = None
        self.G = None
        self.proteins = None
        self.SPL = None

        self.p2c_counts = None

        self.adb = AssociationDatabase.load(cmembers_file)
        self.pp_edges = None
        self._initialize_clusters()

    def _get_solver(self):
        pass

    def _initialize_clusters(self):

        # Add proteins that are only connected to complexes
        all_adb_proteins = set()
        for cmplx in self.adb.associations():
            all_adb_proteins.update(cmplx.members())
        self.proteins = list(all_adb_proteins)

        protein2ix = dict((p,i) for i,p in enumerate(self.proteins))
        asc = self.adb.associations()
        m = len(asc)

        print "Initializing %d clusters (%d proteins)" % (m, len(self.proteins))

        self.clusters = []
        for j, evd in enumerate(asc):
            v = np.zeros(len(self.proteins), dtype='d')
            for p in evd.proteins:
                v[protein2ix[p]] = evd.proteins[p][1]
            self.clusters.append(PlainAssocCluster([evd.complex_id], v, j))

        all_clusters = list(enumerate(self.clusters))

        print "Computing initial weights"
        for i, clstr in all_clusters:
            clstr.compute_weights()

    def _merge(self, minval, min_i, min_j):

        assert min_j < min_i

        clstr1 = self.clusters[min_j]
        clstr2 = self.clusters[min_i]
        c1 = sum(len(self.adb[complex_id])
                 for complex_id in clstr1.complex_ids)
        c2 = sum(len(self.adb[complex_id])
                 for complex_id in clstr2.complex_ids)
        v = (c1*clstr1.v + c2*clstr2.v) / (c1 + c2)
        complex_ids = clstr1.complex_ids + clstr2.complex_ids
        affected_proteins = (clstr1.v > 0.0) & (clstr2.v > 0.0)

        merged_clstr = PlainAssocCluster(complex_ids, v, min_j)

        self.clusters[min_j] = merged_clstr
        del self.clusters[min_i]

        all_clusters = list(enumerate(self.clusters))
        merged_clstr.compute_weights()
        merged_clstr.set_flags(all_clusters)

        affected_clusters = [(min_j, merged_clstr)]
        for i, clstr in all_clusters:
            clstr.set_flags([(min_j, merged_clstr)])
            clstr.resize(min_i)

        self.num_merges += 1
        self.merge_order.append((minval, clstr1.complex_ids,
                                 clstr2.complex_ids))
        print "(%d) MERGE %6d/%6d %6d %.4f %s %s" % (self.num_merges,
                                                int(affected_proteins.sum()),
                                                int((v>0).sum()),
                                                len(affected_clusters),
                                                minval,
                                                str(clstr1.complex_ids),
                                                str(clstr2.complex_ids))
        return affected_clusters


if __name__ == '__main__':

    opts, args = getopt.getopt(sys.argv[1:], 'h', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.stdout.write(_HELP)
            sys.exit()

    if len(args) < 3:
        sys.stderr.write(_EMSG)
        sys.exit()

    data_file = args[0]
    cmembers_file = args[1]
    alpha = float(args[2])


    assert 0.0 <= alpha < 1.0

    ### ***** This could come handy later
    # initial_merge_order = None
    # if initial_merge_filename is not None:
    #     res = load_results(initial_merge_filename)
    #     initial_merge_order = res['merge_order']

    print "PLAIN run, alpha=%.2f" % alpha
    print "Merging similar complexes."
    A = PlainAssocDataset(cmembers_file, alpha)
    A.iterative_merge()
    A.save_results(data_file)
