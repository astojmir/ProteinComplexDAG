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
from ppiTrim import obo
from ppiTrim.parser import parse_mitab_file
from ppiTrim.parser import full_mitab_iterator
from collections import defaultdict
from collections import namedtuple
from operator import attrgetter

_HELP = \
"""
Extract direct interactions from ppiTrim output files.

SYNOPSIS:

    %s [OPTIONS] ppiTrim_file obo_file output_file
    %s -h|--help

ARGUMENTS:

    input_file                  File in ppiTrim MITAB 2.6 format.

OPTIONS:

    -h, --help                  Print this message and exit

OUTPUT:

    List of direct interactions printed to stdout.
""" % (__file__, __file__)

_EMSG = "Insufficient arguments.\n\n" + _HELP


# TO DO:
# - option for excluded pmids

PPIRecordData = namedtuple('PPIRecordData', ['score', 'ppiTrim_id', 'code',
                                             'conflicts', 'interaction_type',
                                             'pmid'])

class PPILinkRow(namedtuple('PPILinkRow', ['edgetype', 'itype', 'proteins',
                                           'score', 'codes', 'sources',
                                           'pmids'])):
    __slots__ = ()

    @classmethod
    def fromline(cls, line):

        fields = line.strip().split('\t')
        edgetype = fields[0]
        itype = fields[1]
        if itype == '-':
            itype = None
        proteins = tuple(fields[2].split('|'))
        if len(proteins) == 1:
            proteins = (proteins[0], proteins[0])
        score = float(fields[3])
        codes = fields[4].split('|')
        sources = fields[5].split('|')
        if len(fields) > 6:
            pmids = fields[6].split('|')
        else:
            pmids = []

        data = (edgetype, itype, proteins, score, codes, sources, pmids)
        return cls._make(data)

    def toline(self):

        row = [self.edgetype,
               self.itype if self.itype is not None else '-',
               '|'.join(self.proteins),
               '%.4f' % self.score,
               '|'.join(self.codes),
               '|'.join(self.sources),
               '|'.join(self.pmids),
               ]
        line = '%s\n' % '\t'.join(row)
        return line


class DirectInteractomeConstructor(object):

    def __init__(self, ppiTrim_file, ontology_file, hpr_cutoff=50,
                 score_exponent=1.6):

        self.ppiTrim_file = ppiTrim_file
        self.ontology_file = ontology_file
        self.hpr_cutoff = hpr_cutoff
        self.score_exponent = score_exponent
        self.hpr_map = self._calculate_hpr()
        self.obo_fp = None
        self.ontology = None
        self.enz_reaction_term = None
        self.direct_interaction_term = None
        self.term_cmp_cache = {}

    def __enter__(self):

        self.obo_fp = open(self.ontology_file, 'r')
        self.ontology = obo.OBOntology(self.obo_fp)
        self.enz_reaction_term = self.ontology.get_term('MI:0414')
        self.direct_interaction_term = self.ontology.get_term('MI:0407')
        return self

    def __exit__(self, type, value, traceback):

        self.obo_fp.close()

    def group_undirected(self):

        undirected_groups = defaultdict(list)

        with open(self.ppiTrim_file, 'rb') as input_fp:
            scanner = parse_mitab_file(input_fp, full_mitab_iterator)
            for intr, _ in scanner:
                network_edgetype = intr.edgetype
                conflicts = self._get_conflicts(intr)

                if network_edgetype == 'X':
                    proteins, score, code = self._score_undirected(intr)
                    if code is None:
                        continue
                    data = PPIRecordData(score, str(intr.checksum), code,
                                         conflicts,
                                         str(intr.interaction_type.term_id),
                                         self._get_pmid(intr))
                    undirected_groups[proteins].append(data)
        return undirected_groups

    def print_undirected(self, undirected_groups, fp=sys.stdout):

        undirected_keys = sorted((sorted(key), key) for key in undirected_groups)
        for proteins, key in undirected_keys:
            if len(proteins) == 1:
                proteins.append(proteins[0])
            src_list = undirected_groups[key]
            score, sources, codes, pmids = self._compute_total_score(src_list)
            if score is not None:
                row = PPILinkRow('X', None, proteins, score, codes, sources, pmids)
                fp.write(row.toline())

    def _calculate_hpr(self):
        """For each pubmed ID, count the number of interactions it reports """

        hpr = defaultdict(int)
        with open(self.ppiTrim_file, 'rb') as input_fp:
            scanner = parse_mitab_file(input_fp, full_mitab_iterator)

            # We assume a single pubmed ID for each interaction
            # 0 does not count
            for intr, _ in scanner:
                pmid = self._get_pmid(intr)
                if pmid != 0:
                    hpr[pmid] += 1
        return hpr

    @staticmethod
    def _get_maxsources(intr):
        maxsources = int([item.acc for item in intr.confidence.ids
                                       if item.db=='maxsources'][0])
        return maxsources

    @staticmethod
    def _get_pmid(intr):
        pmid = 0
        if len(intr.publications.ids):
            pmid = int(intr.publications.ids[0].acc)
        return pmid

    @staticmethod
    def _get_conflicts(intr):
        any_conflicts = [item.acc for item in intr.confidence.ids
                                  if item.db=='conflicts']
        if any_conflicts:
            conflicts = any_conflicts[0].split(',')
        else:
            conflicts = None
        return conflicts

    def _score_undirected(self, intr):
        """Score an undirected interaction"""

        # Here we only consider direct interactions, i.e. NOT pulldowns etc.

        pmid = self._get_pmid(intr)
        maxsources = self._get_maxsources(intr)

        # Score throughput of publication
        shpr = 1.0 if self.hpr_map[pmid] <= self.hpr_cutoff else 0.5
        chpr = 'L' if self.hpr_map[pmid] <= self.hpr_cutoff else 'H'

        # Score number of evidences of this type (usually just 1 but could be more)
        ssrc = maxsources
        csrc = str(maxsources)

        # Score quality - only directed interaction (or its subterms) counts
        term = self.ontology.get_term(intr.interaction_type.term_id)
        if term.compare_to(self.direct_interaction_term, self.term_cmp_cache) >= 0:
            sqlt = 1.0
            cqlt = 'D'
            score = ssrc * shpr * sqlt
            code = chpr + cqlt + csrc
        else:
            score = 0.0
            code = None

        return frozenset(p.alias.ids[0].acc for p in intr.interactors), score, code

    def _compute_total_score(self, src_list):
        """Compute a total score for an interaction"""
        # Deals with conflicts
        idmap = dict()
        direct_interaction = False

        evidences = sorted(src_list, reverse=True, key=attrgetter('score'))
        sources = [evd.ppiTrim_id for evd in evidences]
        codes = [evd.code for evd in evidences]
        pmids = [str(evd.pmid) for evd in evidences]

        for evd in evidences:
            idmap[evd.ppiTrim_id] = evd.score
            if evd.code[1] == 'D':
                direct_interaction = True

        src_wght = 0.0
        for evd in evidences:
            if evd.ppiTrim_id not in idmap:
                continue
            final_evd_score = evd.score

            if evd.conflicts is not None:
                conflict_score = max(idmap.get(k, -1.0) for k in evd.conflicts)
                final_evd_score = max(final_evd_score, conflict_score)
                for k in evd.conflicts:
                    full_conflict_id = 'ppiTrim:%s' % k
                    if full_conflict_id in idmap:
                        del idmap[full_conflict_id]
            src_wght += final_evd_score

        total_score = 1.0 - self.score_exponent ** -(src_wght)
        if not direct_interaction:
            total_score = None

        return total_score, sources, codes, pmids


def main(ppiTrim_file, obo_file, output_file):

    with open(output_file, 'wb') as fp:
        with DirectInteractomeConstructor(ppiTrim_file, obo_file) as cons:
            undirected_groups = cons.group_undirected()
            cons.print_undirected(undirected_groups, fp)


if __name__ == '__main__':

    opts, args = getopt.getopt(sys.argv[1:], 'h', ['help'])
    for o, a in opts:
        if o in ("-h", "--help"):
            sys.stdout.write(_HELP)
            sys.exit()

    if len(args) < 3:
        sys.stderr.write(_EMSG)
        sys.exit()

    ppiTrim_file = args[0]
    obo_file = args[1]
    output_file = args[2]

    main(ppiTrim_file, obo_file, output_file)

