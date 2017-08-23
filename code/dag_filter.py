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
import argparse
from dag_analysis import DAGAnalysis


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process constructed DAGs.')

    parser.add_argument('results_file',
                        help='pickled results file')

    parser.add_argument('output_file',
                        help='output file')

    parser.add_argument('-p',
                        default=2.5,
                        type=float,
                        metavar='pratio_cutoff',
                        dest='pratio_cutoff',
                        help='exclude nodes with lower participation ratio'
                        ' (default: 2.5)')

    parser.add_argument('-u',
                        default='',
                        metavar='pubmed_list',
                        dest='unlabeled_pubmeds',
                        help='exclude nodes without baits arising solely from'
                             ' given publications')

    parser.add_argument('-k',
                        default='K',
                        dest='component_prefix',
                        help='set component prefix (single char, default: K)')


    args = parser.parse_args()

    dag = DAGAnalysis(args.results_file)

    pubmed_ids = None
    if args.unlabeled_pubmeds:
        pubmed_ids = map(str.strip, args.unlabeled_pubmeds.split(','))

    dag.initial_processing(args.pratio_cutoff,
                           pubmed_ids=pubmed_ids,
                           component_prefix=args.component_prefix,
                           )
    dag.save(args.output_file)
