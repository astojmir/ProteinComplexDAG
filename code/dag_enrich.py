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
import getopt
from dag_analysis import DAGAnalysis


_HELP = \
"""
Run enrichment analysis for DAG.

SYNOPSIS:

    %s [OPTIONS] dag_file termdb_file
    %s -h|--help

OPTIONS:

    -h, --help                  Print this message and exit
    -x <path-to-saddlesum>
    -v                          Verbose

""" % (__file__, __file__)

_EMSG = "Insufficient arguments.\n\n" + _HELP


if __name__ == '__main__':

    cmd_path = 'saddlesum'
    verbose = False
    opts, args = getopt.getopt(sys.argv[1:], 'x:hv', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.stdout.write(_HELP)
            sys.exit()
        elif o in ("-x", ):
            cmd_path = a
        elif o in ("-v", ):
            verbose = True

    if len(args) < 2:
        sys.stderr.write(_EMSG)
        sys.exit()

    dag_file = args[0]
    termdb_file = args[1]
    dag = DAGAnalysis.load(dag_file)
    dag.run_enrichment(termdb_file, cmd_path, verbose)
