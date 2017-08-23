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
import os
from collections import defaultdict
from collections import namedtuple
from qmbpmn.web.SaddleSum.results import get_saddlesum_results
from qmbpmn.web.SaddleSum.results import parse_saddlesum_tab_output
from qmbpmn.web.SaddleSum.results import _tab_sections
from qmbpmn.common.utils.filesys import makedirs2
from dictgraph import DictGraph


INCLUDED_SECTIONS = ('biological_process',
                     'cellular_component',
                     'molecular_function')

SaddleSumRow = namedtuple('SaddleSumRow', ['term_id', 'name', 'num_hits',
                                           'score', 'Evalue'])


def process_saddlesum_output(output_file, short_file):
    """
    Filters SaddleSum results (TAB format -F 'tab') to obtain a
    digest that contains only the most significant terms. It does so by taking
    each chain in the ontology DAG, from its minimum to its maximum and picking
    within it a single term with maximal significance. It then reports the set
    of all such maximal significance terms.
    """

    with open(output_file, 'rb') as fp:
        output = fp.read()

    section_data = parse_saddlesum_tab_output(output)
    namespaces, relationships, node_props = section_data[3:6]

    significant_terms = {}
    for ns, data in namespaces:
        for term_id, name, num_hits, score, Evalue, term_url in data:
            significant_terms[term_id] = float(Evalue)

    G = defaultdict(dict)
    target_nodes = set()
    for o1, _, o2 in relationships:
        G[o1][o2] = True
        target_nodes.add(o2)
    minima = set(G.keys()) - target_nodes

    assert all(term_id in significant_terms for term_id in minima)

    retained_terms = set()
    def _traverse(u, best_node):
        bval = significant_terms[best_node]
        uval = significant_terms.get(u, 1000)
        if uval < bval:
            best_node = u
        if u not in G:
            retained_terms.add(best_node)
            return
        for v in G[u]:
            _traverse(v, best_node)

    for u in minima:
        _traverse(u, u)

    with open(short_file, 'w') as fp:
        for ns, data in namespaces:
            fp.write('#\n# {}\n#\n'.format(ns))
            for _row in data:
                row = SaddleSumRow._make(_row[:5])
                if row.term_id in retained_terms:
                    fp.write('\t'.join(row))
                    fp.write('\n')


class EnrichmentMixin(object):

    full_fmt = '{}.{}.{}'
    short_fmt = '{}.{}.txt'
    output_format = 'tab'

    def __init__(self, name, enrich_config, **kwargs):

        super(EnrichmentMixin, self ).__init__()
        self.name = name
        self.config = enrich_config

    def get_raw_weights(self, suffix):

        weight_func = getattr(self, '_get_weights_' + suffix)
        weights = weight_func()
        used_prots = set(p for p, _ in weights)
        all_prots = self.config.proteins
        lines = list()
        for p, x in weights:
            lines.append(str.format('{}\t{:.5g}\n', p, x))
        for p in all_prots:
            if p not in used_prots:
                lines.append(str.format('{}\t{:.5g}\n', p, 0.0))
        return ''.join(lines)

    def _parse_saddlesum_results(self, included_sections=INCLUDED_SECTIONS,
                                 abbrv=True):

        output_dir = self.config.enrichment_dir
        raw_output_dir = os.path.join(output_dir, "tab")
        section2ix = dict((key, i) for i, key in enumerate(included_sections))
        results = dict()

        for suffix, _ in self.config.ssum_opts:

            if abbrv:
                input_name = self.short_fmt.format(self.name, suffix)
                input_file = os.path.join(output_dir, input_name)
            else:
                input_name = self.full_fmt.format(self.name, suffix,
                                                  self.output_format)
                input_file = os.path.join(raw_output_dir, input_name)

            sections = [None] * len(section2ix)
            with open(input_file, 'r') as fp:
                saddlesum_output = fp.read()
                not_found = False
            for title, data in _tab_sections(saddlesum_output):
                if title in section2ix:
                    lines = [SaddleSumRow._make(row[:5]) for row in data]
                    sections[section2ix[title]] = lines
            results[suffix] = sections
        return results

    def run_saddlesum(self):

        output_dir = self.config.enrichment_dir
        termdb_file = self.config.termdb_file
        cmd_path = self.config.saddlesum_cmd_path

        makedirs2(output_dir)
        raw_output_dir = os.path.join(output_dir, "tab")
        makedirs2(raw_output_dir)

        for suffix, ssum_opts in self.config.ssum_opts:

            full_name = self.full_fmt.format(self.name, suffix, self.output_format)
            output_file = os.path.join(raw_output_dir, full_name)
            short_name = self.short_fmt.format(self.name, suffix)
            short_file = os.path.join(output_dir, short_name)

            if not os.path.exists(short_file):
                if not os.path.exists(output_file):
                    opts = [] if ssum_opts is None else ssum_opts.split()
                    raw_weights = self.get_raw_weights(suffix)
                    output, _ = get_saddlesum_results(opts, raw_weights,
                                                      termdb_file,
                                                      self.output_format,
                                                      cmd_path)
                    with open(output_file, 'w') as fp:
                        fp.write(output)
                process_saddlesum_output(output_file, short_file)

    def get_enrichment_results(self, abbrv=True):

        self.run_saddlesum()
        ssum_results = self._parse_saddlesum_results(abbrv=abbrv)
        return ssum_results

    def get_enrichment_data(self):

        ssum_results = self.get_enrichment_results()
        suffixes = [suffix for suffix, _ in self.config.ssum_opts]

        data = []
        term_sets = []
        for suffix in suffixes:
            terms = set()
            for sect in ssum_results[suffix]:
                if sect is None:
                    sect_ = '-'
                else:
                    descs = ['{0.term_id}({0.name})'.format(row) \
                             for row in sect]
                    sect_ = '|'.join(descs)
                    terms |= set(row.term_id for row in sect)
                data.append((sect_, '{}'))
            data.append(('*****', '{}'))
            term_sets.append(terms)

        if len(term_sets) == 2: # two enrichment approaches to compare
            diff01 = len(term_sets[0] - term_sets[1])
            joint = len(term_sets[0] & term_sets[1])

            if len(term_sets[0]) == 0 and len(term_sets[1]) == 0:
                ecode = 'missing'
            elif term_sets[0] == term_sets[1]:
                ecode = 'coincide'
            elif len(term_sets[0]) > 0:
                if diff01 == 0:
                    ecode = 'full'
                elif joint > 0:
                    ecode = 'partial'
                else:
                    ecode = 'disjoint'
            else:
                ecode = 'none'
            data.append((ecode, '{}'))

        return tuple(data)

    def report_enrichment(self, suffix, fp=sys.stdout, scale=None):

        if scale is None:
            scale = self.dag.scale

        sections = self.get_enrichment_results()[suffix]
        sheadings = [ns.replace('_', ' ').title() for ns in INCLUDED_SECTIONS]
        header = ['Term ID', 'Name', 'Associations', 'Score', 'E-value']
        rowfmt = '{:<12s} {:<40s} {:>12s} {:>8s} {:>8s}\n'

        fp.write("**** {}\n".format(self.title))
        for stitle, section in zip(sheadings, sections):
            if section is not None:
                fp.write(str.format('^^^ {} ^^^\n', stitle))
                for row in section:
                    if suffix == 'wsum':
                        score = float(row.score)
                        _score = '{:.4f}'.format(scale * score)
                        row = row[:3] + (_score, ) + row[4:]
                    fp.write(rowfmt.format(*row))
                fp.write('#\n')
        fp.write('----\n')

    def find_links(self, other):

        G = self.dag.filtered_graph
        H = DictGraph()
        other_complexes = set(other.complexes)

        for cmplx1 in self.complexes:
            if cmplx1 in G:
                nbhd = [(cmplx2, rdist) for cmplx2, rdist in G[cmplx1].items()
                        if cmplx2 in other_complexes]
                if nbhd:
                    H[cmplx1] = dict(nbhd)
        return H
