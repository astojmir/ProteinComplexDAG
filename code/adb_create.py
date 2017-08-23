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
from cmembers import initial_import


_HELP = \
"""
Extract complexes and their members from ppiTrim output files.

SYNOPSIS:

    %s [OPTIONS] ppiTrim_file output_file
    %s -h|--help

ARGUMENTS:

    input_file                  File in ppiTrim MITAB 2.6 format.
    output_file                 Name of the output file

OPTIONS:

    -h, --help                  Print this message and exit
    -p <pmids_file>             File with Pubmed IDs to filter
    -n <complexes_per_protein>  Minimum complexes per protein. Non-bait
                                proteins present in fewer complexes are
                                removed.

OUTPUT:

    Dataset of complexes.
""" % (__file__, __file__)

_EMSG = "Insufficient arguments.\n\n" + _HELP


if __name__ == '__main__':

    pmids_file = None
    min_complexes_per_protein = None

    opts, args = getopt.getopt(sys.argv[1:], 'hn:p:', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.stdout.write(_HELP)
            sys.exit()
        elif o in ("-n", ):
            min_complexes_per_protein = int(a)
        elif o in ("-p", ):
            pmids_file = a

    if len(args) < 2:
        sys.stderr.write(_EMSG)
        sys.exit()

    ppiTrim_file = args[0]
    output_file = args[1]

    initial_import(ppiTrim_file, output_file, pmids_file,
                   min_complexes_per_protein)
