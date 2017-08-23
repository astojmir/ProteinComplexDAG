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

import cPickle as pickle
from operator import itemgetter


def load_results(results_file):

    with open(results_file, 'rb') as fp:
        obj_data = pickle.load(fp)
    return obj_data


def participation_ratio(x):

    s = x.sum()
    if s == 0.0:
        return 0.0
    return s ** 2 / (x**2).sum()


def show_edge_details(dag, edge):

    a, b = [dag.nr2complex[p] for p in edge]
    dag.report_details(a)
    dag.report_details(b)
    print '********* DISTANCES ***********'
    print '{:.4f}, {:.4f}'.format(dag[a].rdist(dag[b]),
                                  dag[b].rdist(dag[a]))


def venn3(A, B, C, verbose=False):

    res = [('A & B & C', sorted(A & B & C)),
           ('A - B - C', sorted(A - B - C)),
           ('A & B - C', sorted(A & B - C)),
           ('A & C - B', sorted(A & C - B)),
           ('B - A - C', sorted(B - A - C)),
           ('B & C - A', sorted(B & C - A)),
           ('C - A - B', sorted(C - A - B)),
           ]
    total = 1.0 * sum(len(x[1]) for x in res)
    if verbose:
        for lbl, X in res:
            print '|{}| = {:d} ({:.4f})'.format(lbl, len(X), len(X) / total)
    return res


def characterize_non_transitive_links(dag, max_tol=0.2):

    ntrd_all = dag.get_essential_non_transitive_links()
    ntrd_large = [item for item in ntrd_all if item[5] > max_tol]
    small_nodes = sorted(set(item[0] for item in ntrd_large))

    smdata = []
    for node_id in small_nodes:
        node = dag[node_id]
        smdata.append((node_id, len(node.evd.proteins), node.pratio))
    smdata.sort(key=itemgetter(1), reverse=True)
    return smdata, ntrd_all
