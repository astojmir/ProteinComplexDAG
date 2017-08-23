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
import random
import cPickle as pickle
import numpy as np
from collections import defaultdict
from cmembers import AssociationDatabase
from dictgraph import DictGraph
import infflow


ITM_WEIGHTS_FULL = 0
ITM_WEIGHTS_MEMBERS = 1


class AssocCluster(object):

    rtol = 0.01
    weight_range = ITM_WEIGHTS_MEMBERS

    def __init__(self, complex_ids, v, self_index):

        self.complex_ids = complex_ids
        self.v = v     # Clustering (plain, average) weights

        self.x = None  # ITM weights
        self.s = None  # Self-similarity
        self.jix = None # Significant indices

        self.m = self_index  # Size of similarity array
        # Set if significant indices are totally disjoint
        self.flags = np.ones(self.m, dtype='b')
        # Similarities to other clusters
        self.sim = np.zeros(self.m, dtype='d')
        # Symmetric distances
        self.rho = np.zeros(self.m, dtype='d')


        self.member_mask = self.v > 0
        self.outside_mask = self.v == 0

        if self.weight_range == ITM_WEIGHTS_MEMBERS:
            self.jix = np.arange(len(self.v), dtype='int')[self.member_mask]

    def compute_itm_weights(self, SPL, pc):

        # This could be changed if we use v to weight pc links as well as cp
        # ones
        p = pc.copy()
        p[self.v == 0] = 0.0

        f = SPL.solve(p, True)
        h = SPL.solve(self.v, False)
        self.x = h * f

        if self.weight_range == ITM_WEIGHTS_FULL:
            self.s = self.x.sum()
            N = self.rtol * self.s

            # Find most significant indices
            x_ix0 = np.argsort(self.x)
            x0 = self.x[x_ix0]
            k = ((x0.cumsum() - N) >= 0.0).argmax()
            self.jix = x_ix0[k:]
            self.x_val = self.x[self.jix]

        else:
            self.x_val = self.x[self.member_mask]
            self.x[self.outside_mask] = 0.0
            self.s = self.x_val.sum()

    def set_flags(self, others):

        relevant_others = [(i, clstr) for i, clstr in others if i < self.m]
        for i, clstr in relevant_others:
            if not (set(self.jix) & set(clstr.jix)):
                self.flags[i] = False
                self.sim[i] = 0.0
                self.rho[i] = 1.0

    def update(self, others):

        relevant_others = [(i, clstr) for i, clstr in others if i < self.m]
        x_val = self.x_val
        y_val = np.empty(len(self.jix), dtype='d')
        for i, clstr in relevant_others:
            if self.flags[i]:
                y_val[:] = clstr.x[self.jix]
                np.minimum(x_val, y_val, y_val)
                sxy = y_val.sum()
                self.sim[i] = sxy
                self.rho[i] = 1.0 - (sxy / max(self.s, clstr.s))

    def get_minimum_rho(self, self_ix):

        min_j = np.argmin(self.rho[:self.m])
        minval = self.rho[min_j]
        return minval, min_j

    def resize(self, i):

        if i < self.m:
            self.m -= 1
            # This creates a copy
            self.sim = np.delete(self.sim, i)
            self.flags = np.delete(self.flags, i)
            self.rho = np.delete(self.rho, i)
            # Alternative
            # self.sim[i:self.m] = self.sim[i+1:self.m+1]

    def is_affected(self, protein_ixs):

        if self.jix is None:
            return True
        return np.any(protein_ixs[self.jix])


class AssocDataset(object):

    saved_attrs = ['pplinks_file',
                   'score_cutoff',
                   'cmembers_file',
                   'alpha',
                   'beta',
                   'rtol',
                   'weight_range',
                   'proteins',
                   'merge_order',
                   'min_pc_count',
                   'min_pp_count',
                   ]
    extra_attrs = ['adb',
                   'pp_edges',
                   ]

    def __init__(self, pplinks_file, cmembers_file, alpha, beta, rtol,
                 score_cutoff=0.21, weight_range=ITM_WEIGHTS_MEMBERS,
                 initial_merge_order=None, min_pc_count=1,
                 min_pp_count=1, shuffle=False):

        self.num_merges = 0
        self.merge_order = []

        self.pplinks_file = pplinks_file
        self.score_cutoff = score_cutoff
        self.cmembers_file = cmembers_file

        self.alpha = alpha
        self.beta = beta
        self.rtol = AssocCluster.rtol = rtol
        self.weight_range = AssocCluster.weight_range = weight_range

        self.min_pc_count = min_pc_count
        self.min_pp_count = min_pp_count
        self.shuffle = shuffle

        self.clusters = None
        self.G = None
        self.proteins = None
        self.SPL = None

        self.p2c_counts = None

        self.adb = AssociationDatabase.load(cmembers_file)
        self.pp_edges = infflow.load_pp_edges(self.pplinks_file,
                                              self.score_cutoff)
        self._set_minimum_counts()
        if self.beta >= 1.0:
            self._find_beta()
        self._get_pp_solver()
        self._initialize_clusters(initial_merge_order)

    def __getstate__(self):

        odict = self.__dict__.copy()
        odict['SPL'] = None
        return odict

    def __setstate__(self, odict):

        self.__dict__.update(odict)
        self._get_solver()

    @property
    def all_clusters(self):

        return list(enumerate(self.clusters))

    @property
    def pc(self):

        # pc = (1.0 - self.beta) / self.p2c_counts
        pc = np.zeros(len(self.p2c_counts), dtype='d')
        pc[:] = 1.0 - self.beta
        _lrg_counts = self.p2c_counts > self.min_pc_count
        _sml_counts = self.p2c_counts <= self.min_pc_count
        pc[_lrg_counts] /= self.p2c_counts[_lrg_counts]
        pc[_sml_counts] /= self.min_pc_count
        return pc

    def _find_beta(self):

        print 'Solving for beta - path length = %.2f' % self.beta
        self._avg_ppi_path_length = self.beta
        self.beta = infflow.find_beta(self.pp_edges, self.min_pp_count,
                                      self.beta, self.adb, verbose=True)

    def _set_minimum_counts(self):

        if self.min_pc_count == 0:
            self.min_pc_count = float(self.adb.get_weighted_median_pc_degree())
            print "Set min_pc_count to %.1f" % self.min_pc_count
        if self.min_pp_count == 0:
            self.min_pp_count = float(infflow.get_median_pp_degree(self.pp_edges))
            print "Set min_pp_count to %.1f" % self.min_pp_count

    def _get_pp_solver(self):

        print "Computing graph"
        self.G = infflow.get_pp_graph(self.pp_edges, self.min_pp_count,
                                      self.beta, self.adb)
        self.proteins = map(str, self.G.nodes)
        self.SPL = infflow.get_pp_solver(self.G)

    def _initialize_clusters(self, initial_merge_order=None):

        if self.shuffle:
            # Random shuffling of proteins - testing only
            print "Random shuffling of proteins."
            random.shuffle(self.proteins)

        protein2ix = dict((p,i) for i,p in enumerate(self.proteins))
        self.p2c_counts = np.zeros(len(self.proteins), dtype='int')
        asc = self.adb.associations()
        m = len(asc)

        fmt = "Initializing {:d} clusters ({:d} proteins)"
        print fmt.format(m, len(self.proteins))

        self.clusters = []
        for j, evd in enumerate(asc):
            v = np.zeros(len(self.proteins), dtype='d')
            for p in evd.proteins:
                if evd.proteins[p][1] > 0.0:
                    v[protein2ix[p]] = evd.proteins[p][1]
                self.p2c_counts[protein2ix[p]] += 1
            self.clusters.append(AssocCluster([evd.complex_id], v, j))

        if initial_merge_order is not None:
            self.replay_merge(initial_merge_order)

        print "Computing initial ITM weights"
        for clstr in self.clusters:
            clstr.compute_itm_weights(self.SPL, self.pc)

    def _compute_initial_distances(self):

        print "Setting flags"
        for i, clstr in self.all_clusters:
            clstr.set_flags(self.all_clusters)
            if (i+1) % 1000 == 0:
                print (i+1)
        print (i+1)
        print "Computing initial distances"
        for i, clstr in self.all_clusters:
            clstr.update(self.all_clusters)
            if (i+1) % 1000 == 0:
                print (i+1)
        print (i+1)

    def iterative_merge(self):

        affected_clusters = None
        self._compute_initial_distances()
        minval, min_i, min_j = self._find_smallest_rho(self.all_clusters)

        while minval <= self.alpha:
            affected_clusters = self._merge(minval, min_i, min_j)
            self._update_clusters(affected_clusters)
            minval, min_i, min_j = self._find_smallest_rho(affected_clusters)

    def _find_index(self, complex_ids):

        i = None
        for j, clstr in enumerate(self.clusters):
            if complex_ids == clstr.complex_ids:
                i = j
                break
        assert i is not None
        return i

    def replay_merge(self, initial_merge_order):

        for minval, complex_ids1, complex_ids2 in initial_merge_order:
            if minval > self.alpha:
                break
            min_j = self._find_index(complex_ids1)
            min_i = self._find_index(complex_ids2)
            self._merge(minval, min_i, min_j, update=False)

    def _merge(self, minval, min_i, min_j, update=True):

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
        self.p2c_counts[affected_proteins] -= 1

        merged_clstr = AssocCluster(complex_ids, v, min_j)

        self.clusters[min_j] = merged_clstr
        del self.clusters[min_i]

        if update:
            merged_clstr.compute_itm_weights(self.SPL, self.pc)

            # Note that distances, similarities and flags are only set to
            # indices less than clstr.m (where clstr is a cluster). Since new
            # cluster is placed at min_j position (and min_j < min_i),
            # computing flags here will not affect clusters where m > min_j,
            # which may be resized due to deletion of min_i.
            merged_clstr.set_flags(self.all_clusters)

        affected_clusters = []
        for i, clstr in self.all_clusters:
            if update:
                clstr.set_flags([(min_j, merged_clstr)])
            clstr.resize(min_i)
            if update and clstr.is_affected(affected_proteins):
                affected_clusters.append((i, clstr))
                clstr.compute_itm_weights(self.SPL, self.pc)

        self.num_merges += 1
        self.merge_order.append((minval, clstr1.complex_ids,
                                 clstr2.complex_ids))

        fmt = "({:d}) MERGE {:d}/{:d} {:6d} {:.4f} {} {}"
        print fmt.format(self.num_merges, int(affected_proteins.sum()),
                         int((v>0).sum()), len(affected_clusters), minval,
                         str(clstr1.complex_ids), str(clstr2.complex_ids))

        return affected_clusters

    def _update_clusters(self, affected_clusters):

        affected_indices = set(i for i, _ in affected_clusters)
        for i, clstr in self.all_clusters[1:]:
            if i in affected_indices:
                clstr.update(self.all_clusters)
            else:
                clstr.update(affected_clusters)

    def _find_smallest_rho(self, affected_clusters):

        minval = 1.0
        min_i = 0
        min_j = 0
        for i, clstr in self.all_clusters[1:]:
            minval2, min_j2 = clstr.get_minimum_rho(i)
            if minval2 < minval:
                minval = minval2
                min_i = i
                min_j = min_j2
        return minval, min_i, min_j

    def save_results(self, output_file):

        m = len(self.clusters)
        itm_weights = []
        complexes = []
        complex2nr = {}
        for k, clstr in enumerate(self.clusters):
            itm_weights.append((clstr.x_val, clstr.jix))
            complex_id = self.adb.merge_clusters(clstr.complex_ids)
            complexes.append(complex_id)
            complex2nr[complex_id] = clstr.complex_ids

        G = DictGraph()
        for i in xrange(1, m):
            cmplx1 = complexes[i]
            clstr1 = self.clusters[i]
            for j in xrange(i):
                if clstr1.flags[j]:
                    cmplx2 = complexes[j]
                    clstr2 = self.clusters[j]
                    dist1 = 1.0 - (clstr1.sim[j] / clstr1.s)
                    dist2 = 1.0 - (clstr1.sim[j] / clstr2.s)
                    assert not (dist1 <= self.alpha and dist2 <= self.alpha)
                    if dist1 <= self.alpha:
                        G[cmplx1][cmplx2] = dist1
                    elif dist2 <= self.alpha:
                        G[cmplx2][cmplx1] = dist2

        output = {'complexes': complexes,
                  'df': 1.0,
                  'graph': G,
                  'itm_weights': itm_weights,
                  'complex2nr': complex2nr,
                  }
        output.update((attr, getattr(self, attr)) for attr in self.saved_attrs)
        output.update((attr, getattr(self, attr)) for attr in self.extra_attrs)
        with open(output_file, 'wb') as fp:
            pickle.dump(output, fp, 2)

    def save_initial_weights(self, output_file):

        X = np.empty((len(self.clusters), len(self.proteins)), dtype='d')
        complexes = []
        itm_weights = []
        for k, clstr in enumerate(self.clusters):
            itm_weights.append((clstr.x_val, clstr.jix))
            complex_id = self.adb.merge_clusters(clstr.complex_ids)
            complexes.append(complex_id)

        output = {'complexes': complexes,
                  'itm_weights': itm_weights,
                  }
        output.update((attr, getattr(self, attr)) for attr in self.saved_attrs)
        with open(output_file, 'wb') as fp:
            pickle.dump(output, fp, 2)
