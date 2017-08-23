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

import os
import sys
import getopt
from utils import load_results
import itermerge


_HELP = \
"""
SYNOPSIS:

    %s [OPTIONS] label cmembers_file pplinks_file alpha beta [score_cutoff]
    %s -h|--help

OPTIONS:

    -h, --help                  Print this message and exit
    -p <min_pp_count>
    -c <min_pc_count>
    -m <filename>               Read initial merge order from filename

""" % (__file__, __file__)

_EMSG = "Insufficient arguments.\n\n" + _HELP


if __name__ == '__main__':

    initial_merge_filename = None
    initial_merge_order = None
    min_pp_count = 1
    min_pc_count = 1
    opts, args = getopt.getopt(sys.argv[1:], 'hm:p:c:', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.stdout.write(_HELP)
            sys.exit()
        elif o in ("-m"):
            initial_merge_filename = a
        elif o in ("-p", ):
            min_pp_count = float(a)
        elif o in ("-c", ):
            min_pc_count = float(a)

    weight_range = itermerge.ITM_WEIGHTS_MEMBERS
    score_cutoff = 0.21
    rtol = 0.01

    if len(args) < 5:
        sys.stderr.write(_EMSG)
        sys.exit()

    data_file = args[0]
    cmembers_file = args[1]
    pplinks_file = args[2]
    alpha = float(args[3])
    beta = float(args[4])

    if len(args) >= 6:
        score_cutoff = float(args[5])

    assert 0.0 <= alpha < 1.0
    assert 0.0 <= beta

    if initial_merge_filename is not None:
        print "Warning - all parameters reused except alpha."
        print "Warning - beta parameter only used for label."
        res = load_results(initial_merge_filename)
        initial_merge_order = res['merge_order']
        cmembers_file = res['cmembers_file']
        pplinks_file = res['pplinks_file']
        beta = res['beta']
        score_cutoff = res['score_cutoff']
        min_pp_count = res['min_pp_count']
        min_pc_count = res['min_pc_count']
        rtol = res['rtol']

    print "ITM (simplified) run, alpha=%.2f, beta=%.2f" % (alpha, beta)
    print "Merging similar complexes."
    A = itermerge.AssocDataset(pplinks_file, cmembers_file, alpha, beta,
                               rtol, score_cutoff, weight_range,
                               initial_merge_order, min_pc_count,
                               min_pp_count, False)
    A.iterative_merge()
    A.save_results(data_file)
