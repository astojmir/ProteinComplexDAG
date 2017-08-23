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
from collections import defaultdict
from operator import itemgetter
from cmembers import AssociationDatabase
from pplinks import PPILinkRow

_HELP = \
"""
Write a report based on association database.

SYNOPSIS:

    %s [OPTIONS] report_name input_file output_file
    %s -h|--help

ARGUMENTS:

    report_name                 Name of report function
    input_file                  Association database file
    output_file                 Output file

OPTIONS:

    -h, --help                  Print this message and exit

OUTPUT:

   Specified report
""" % (__file__, __file__)

_EMSG = "Insufficient arguments.\n\n" + _HELP


def tab(cmembers_file, output_file):

    adb = AssociationDatabase.load(cmembers_file)
    with open(output_file, 'w') as fp:
        adb.write_tab(fp)

def tab2(cmembers_file, output_file):

    adb = AssociationDatabase.load(cmembers_file)
    with open(output_file, 'w') as fp:
        adb.write_tab2(fp)

def details(cmembers_file, output_file):

    adb = AssociationDatabase.load(cmembers_file)
    with open(output_file, 'w') as fp:
        adb.write_details(fp)

def pubmed(cmembers_file, stats_file):

    adb = AssociationDatabase.load(cmembers_file)
    with open(stats_file, 'w') as fp:
        adb.write_pubmed_stats(fp)

def complex(cmembers_file, output_file, complex_id):

    adb = AssociationDatabase.load(cmembers_file)
    if output_file == '-':
        fp = sys.stdout
    else:
        fp = open(output_file, 'w')

    cls = adb[complex_id]
    cls.write_line2(fp)

def p2c(cmembers_file, output_file):

    p2c_counts = defaultdict(int)
    adb = AssociationDatabase.load(cmembers_file)
    for evd in adb.associations():
        for p in evd.proteins:
            if evd.proteins[p][1] > 0.0:
                p2c_counts[p] += 1
    items = sorted(p2c_counts.iteritems(), key=itemgetter(1,0))
    with open(output_file, 'w') as fp:
        for i, row in enumerate(items):
            fp.write('{}\t{}\t{}\n'.format(i, *row))

def p2p(pplinks_file, output_file, score_cutoff=0.21):

    score_cutoff = float(score_cutoff)
    p2p_counts = defaultdict(int)
    with open(pplinks_file, 'r') as fp:
        for line in fp:
            row = PPILinkRow.fromline(line)
            if row.edgetype == 'D':
                continue
            p1, p2 = row.proteins[:2]
            if row.score >= score_cutoff and p1 != p2:
                p2p_counts[p1] += 1
                p2p_counts[p2] += 1

    items = sorted(p2p_counts.iteritems(), key=itemgetter(1,0))
    with open(output_file, 'w') as fp:
        for i, row in enumerate(items):
            fp.write('{}\t{}\t{}\n'.format(i, *row))




if __name__ == '__main__':

    opts, args = getopt.getopt(sys.argv[1:], 'h', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.stdout.write(_HELP)
            sys.exit()

    if len(args) < 3:
        sys.stderr.write(_EMSG)
        sys.exit()

    func_name = args[0]
    input_file = args[1]
    output_file = args[2]

    f = globals()[func_name]
    f(input_file, output_file, *args[3:])
