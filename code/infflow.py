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

import numpy as np
from collections import defaultdict
from pplinks import PPILinkRow
from qmbpmn.common.graph.digraph import DirectedGraph
from qmbpmn.ITMProbe.core.laplacian import BasicLaplacian
from qmbpmn.common.utils.newton import rootfind_newton


def load_pp_edges(pplinks_file, score_cutoff):

    pp_edges = dict()
    with open(pplinks_file, 'r') as fp:
        for line in fp:
            row = PPILinkRow.fromline(line)
            if row.edgetype == 'D':
                continue
            p1, p2 = row.proteins[:2]
            if row.score >= score_cutoff and p1 != p2:
                pp_edges[frozenset(row.proteins)] = row.pmids
    return pp_edges


def get_pp_graph(pp_edges, min_pp_count, beta, adb=None):

    G = DirectedGraph()
    # Create pp graph
    p_edges = defaultdict(list)
    for edge in pp_edges:
        p1, p2 = tuple(edge)
        p_edges[p1].append(p2)
        p_edges[p2].append(p1)

    for p1 in p_edges:
        pp_deg = max(len(p_edges[p1]), min_pp_count)
        weight = beta / pp_deg
        for p2 in p_edges[p1]:
            G.insert_edge(p1, p2, weight)

    # Add proteins that are only connected to complexes
    if adb is not None:
        all_adb_proteins = set()
        for cmplx in adb.associations():
            all_adb_proteins.update(cmplx.members())
        for p in all_adb_proteins:
            G.insert_node(p)

    return G


def get_pp_solver(G):

    W = G.weighted_adjacency_matrix()
    W.row_weights[:] = 1.0
    df_mask = W.get_df_mask(1.0, None, 1.0, None)
    SPL = BasicLaplacian(W, df_mask)
    return SPL


def get_avg_pp_path_length(SPL):

    p = np.ones(SPL.L.shape[0], dtype='d')
    Gp = SPL.solve(p, True)
    GGp = SPL.solve(Gp, True)
    return Gp.mean(), (GGp - Gp).mean()


def find_beta(pp_edges, min_pp_count, target_avg_path_length, adb=None,
              verbose=False):

    assert target_avg_path_length >= 1.0

    def _root_func(x0):
        G = get_pp_graph(pp_edges, min_pp_count, x0, adb)
        SPL = get_pp_solver(G)
        t, beta_difft = get_avg_pp_path_length(SPL)
        fval = t - target_avg_path_length
        fpval = beta_difft / x0
        if verbose:
            print ', '.join(map(lambda y: '%.3f' % y, [x0, fval, fpval]))
        return fval, fpval, SPL

    a = 0.0     # lower bound for beta
    b = 0.9999  # upper bound for beta
    x0 = 0.8    # starting guess
    x, _, _, _ =  rootfind_newton(_root_func, x0, a, b)
    return x


def get_median_pp_degree(pp_edges):

    p2p_counts = defaultdict(int)
    for edge in pp_edges:
        p1, p2 = tuple(edge)
        p2p_counts[p1] += 1
        p2p_counts[p2] += 1
    v = sorted(p2p_counts.itervalues())
    return v[len(v) // 2]
