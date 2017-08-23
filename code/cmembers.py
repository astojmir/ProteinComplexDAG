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
import sys
import copy
import numpy as np
from operator import attrgetter
from collections import defaultdict
from collections import Counter
from itertools import chain
from ppiTrim.parser import parse_mitab_file
from ppiTrim.parser import full_mitab_iterator
from ppiTrim.parser import partial_mitab_iterator


def _get_eqv_group(intr):
    grp = [intr.checksum.acc]
    any_conflicts = [item.acc for item in intr.confidence.ids
                              if item.db=='conflicts']
    if any_conflicts:
        grp += any_conflicts[0].split(',')
    return grp


def _offsets_by_bait_evidence_consumer(scanner):
    # We are grouping all evidences 'in conflict' together.
    # However, we must be careful because 'in conflict' (which really means 'we
    # believe it is eqivalent to' is not actually an equivalence relation (it
    # is not transitive).


    evdgrp_ixs = {} # maps ppiTrim_id -> index in evdgrp_list
    evdgrp_list = [] # contains lists of offsets

    for intr, lines in scanner:
        if not intr.is_complex() or len(intr.interactors) < 4:
            continue
        line_offsets = lines[0]

        eqvgrp = _get_eqv_group(intr)

        eqixs = set(evdgrp_ixs[key] for key in eqvgrp if key in evdgrp_ixs)
        assert len(eqixs) < 2
        if len(eqixs) == 0:
            i = len(evdgrp_list)
            evdgrp_list.append([])
            for key in eqvgrp:
                evdgrp_ixs[key] = i
        else:
            i = eqixs.pop()
            if len(eqixs) > 0:
                for key in eqvgrp:
                    j = evdgrp_ixs[key]
                    if j != i:
                        evdgrp_list[i] += evdgrp_list[j]
                        evdgrp_list[j] = None
                        evdgrp_ixs[key] = i
        evdgrp_list[i].extend(line_offsets)

    return [offsets for offsets in evdgrp_list if offsets is not None]


class Association(object):

    bait_code = 'B'
    prey_code = 'P'
    unknown_code = 'U'

    @classmethod
    def max_code(cls, code1, code2):

        if code1 == cls.bait_code or code2 == cls.bait_code:
            return cls.bait_code
        if code1 == cls.prey_code or code2 == cls.prey_code:
            return cls.prey_code
        return cls.unknown_code

    def write_line(self, fp):

        proteins = self.proteins
        has_baits = any(v[0] == self.bait_code for v in proteins.itervalues())

        fields = (self.complex_id,
                  '|'.join(self.pubmed_ids),
                  '|'.join(self.ppiTrim_ids),
                  '|'.join(self.detection_method_term_ids),
                  '|'.join(self.edgetype_codes),
                  str(self.maxsources),
                 )
        fp.write('\t'.join(fields))

        prot_data = [(pp, vv[0], vv[1]) for (pp, vv) in proteins.iteritems()]
        prot_data.sort(key=lambda item: (item[1], -item[2], item[0]))

        for p, code, w in prot_data:
            fp.write('\t{}|{}|{:.4f}'.format(p, code, w))
        fp.write('\n')

    def write_line2(self, fp):

        proteins = self.proteins
        has_baits = any(v[0] == self.bait_code for v in proteins.itervalues())

        fields = (self.complex_id,
                  str(len(set(self.pubmed_ids))),
                  '|'.join(sorted(set(self.pubmed_ids))),
                  str(len(self)),
                  '|'.join(evd.complex_id for evd in self.evidences),
                  str(self.total_weight()),
                  str(self.total_weight() / len(proteins)),
                 )
        fp.write('\t'.join(fields))

        prot_data = [(pp, vv[0], vv[1]) for (pp, vv) in proteins.iteritems()]
        prot_data.sort(key=lambda item: (item[1], -item[2], item[0]))

        for p, code, w in prot_data:
            fp.write('\t{}|{}|{:.4f}'.format(p, code, w))
        fp.write('\n')

    def members(self):
        return set(self.proteins.iterkeys())

    def member_weights(self):
        return set((p, v[1]) for p, v in self.proteins.iteritems())

    def total_weight(self):
        return sum(v[1] for v in self.proteins.itervalues())

    def report_details(self, fp=sys.stdout):

        ppiTrim_ids = ', '.join(self.ppiTrim_ids)
        pmids = ', '.join(self.pubmed_ids)
        dmethods = ', '.join(self.detection_method_term_ids)
        edgetypes = ', '.join(self.edgetype_codes)
        prot_data = [(p, v[0], v[1]) for (p, v) in self.proteins.iteritems()]
        prot_data.sort(key=lambda item: (item[1], -item[2], item[0]))

        fp.write("**** COMPLEX {}\n".format(self.complex_id))
        fp.write("num_proteins\t{:d}\n".format(len(self.proteins)))
        fp.write("evidence_weight\t{:.2f}\n".format(self.total_weight()))
        fp.write("#\n")
        fp.write("ppiTrim_ids\t{}\n".format(ppiTrim_ids))
        fp.write("pubmed_ids\t{}\n".format(pmids))
        fp.write("detection_methods\t{}\n".format(dmethods))
        fp.write("edgetype_codes\t{}\n".format(edgetypes))
        fp.write("%\n")

        for p, b, y in prot_data:
            fp.write('{}\t{:.2f}\t{}\n'.format(p, y, b))
        fp.write('----\n')


class AssociationEvidence(Association):

    def __init__(self):

        self._id = None
        self.complex_id = None
        self.pubmed_ids = None
        self.ppiTrim_ids = []
        self.detection_method_term_ids = []
        self.edgetype_codes = []
        self.maxsources = 0
        self.weight = 1.0
        self.proteins = defaultdict(lambda : [self.unknown_code, 0.0])
        self.has_baits = False

    def __len__(self):

        return 1

    def get_erole_code(self, p):

        if p.experimental_role[0].term_id == 'MI:0496':
            return self.bait_code
        if p.experimental_role[0].term_id == 'MI:0498':
            return self.prey_code
        return self.unknown_code

    def insert_complex(self, intr, bait_indices=[], weight=1.0):

        n = len(self.ppiTrim_ids)
        if self.pubmed_ids is None:
            self.pubmed_ids = [intr.publications.ids[0].acc]
        else:
            assert self.pubmed_ids[0] == intr.publications.ids[0].acc
        self.ppiTrim_ids.append(str(intr.checksum))
        self.detection_method_term_ids.append(intr.detection_method.term_id)
        self.edgetype_codes.append(str(intr.edgetype))
        _maxsources = [int(cnfid.acc) for cnfid in intr.confidence.ids
                       if cnfid.db == 'maxsources']
        if not _maxsources:
            _maxsources = 1
        else:
            _maxsources = _maxsources[0]
        self.maxsources += _maxsources

        for i, p in enumerate(intr.interactors):
            pname = p.alias.ids[0].acc
            w = (n * self.proteins[pname][1] + weight) / (n + 1)
            self.proteins[pname][1] = w
            if i in bait_indices:
                self.proteins[pname][0] = self.bait_code
                self.has_baits = True
            else:
                self.proteins[pname][0] = self.get_erole_code(p)

    @property
    def evidences(self):
        return [ self ]


class AssociationCluster(Association):


    def __init__(self, evidences):

        self._id = None
        self.complex_id = None
        self.evidences = evidences

    def __len__(self):

        return len(self.evidences)

    def add(self, evd):

        self.evidences.append(evd)

    @property
    def pubmed_ids(self):
        return [evd.pubmed_ids[0] for evd in self.evidences]

    @property
    def ppiTrim_ids(self):
        return [pid for evd in self.evidences for pid in evd.ppiTrim_ids]

    @property
    def detection_method_term_ids(self):
        return [tid for evd in self.evidences \
                    for tid in evd.detection_method_term_ids]

    @property
    def edgetype_codes(self):
        return [ecd for evd in self.evidences for ecd in evd.edgetype_codes]

    @property
    def maxsources(self):
        return sum(evd.maxsources for evd in self.evidences)

    @property
    def proteins(self):

        proteins = defaultdict(lambda : [self.unknown_code, 0.0])
        wsum = 0.0
        for evd in self.evidences:
            for pname in evd.proteins:
                proteins[pname][0] = self.max_code(proteins[pname][0],
                                                   evd.proteins[pname][0])
                proteins[pname][1] += evd.weight * evd.proteins[pname][1]
            wsum += evd.weight
        for pname in proteins:
            proteins[pname][1] /= wsum
        return proteins


class AssociationDatabase(object):

    def __init__(self, id_prefix='C', id_counter=0, merge_prefix='M',
                 fix_prefix='F'):

        self.id_prefix = id_prefix
        self.id_counter = id_counter
        self.merge_prefix = merge_prefix
        self.fix_prefix = fix_prefix

        self.raw_evidences = []
        self.clusters = {}
        self.singletons = {}
        self.removed = {}

    def __len__(self):

        return len(self.singletons) + len(self.clusters)

    def __getitem__(self, complex_id):

        if complex_id in self.singletons:
            return self.singletons[complex_id]
        if complex_id in self.clusters:
            return self.clusters[complex_id]
        return self.removed[complex_id][0]

    def remove_complex(self, complex_id, code):
        if complex_id in self.singletons:
            dct = self.singletons
        else:
            dct = self.clusters
        self.removed[complex_id] = (dct.pop(complex_id), code)

    def valid_complexes_iterator(self):

        for cls in self.clusters.itervalues():
            yield cls
        for evd in self.singletons.itervalues():
            yield evd

    @staticmethod
    def load(input_file):
        with open(input_file, 'rb') as fp:
            obj = pickle.load(fp)
        return obj

    def save(self, output_file):

        with open(output_file, 'wb') as fp:
            pickle.dump(self, fp, 2)

    def set_complex_id(self, evd):

        evd.complex_id = '{0.id_prefix}{0.id_counter:0=5d}'.format(self)
        evd._id = self.id_counter
        self.id_counter += 1

    def set_merged_id(self, cls):

        assert len(cls.evidences) > 1
        evd = min(cls.evidences, key=attrgetter('_id'))
        cls.complex_id = '{0.merge_prefix}{1._id:0=5d}'.format(self, evd)
        cls._id = evd._id

    def set_fixed_id(self, evd, old_complex_id):

        assert old_complex_id[0] == self.id_prefix
        evd.complex_id = '{0.fix_prefix}{1._id:0=5d}'.format(self, evd)

    def get_id_pattern(self):
        prefixes = [self.id_prefix,
                    self.merge_prefix,
                    self.fix__prefix]
        return '[{}][0-9]{{5}}'.format(''.join(prefixes))

    @staticmethod
    def _get_bait_indices(intr):

        baits = [i for i, p in enumerate(intr.interactors) \
                   if p.experimental_role[0].term_id == 'MI:0496']
        return baits

    def insert_evidences(self, ppiTrim_file):

        with open(ppiTrim_file, 'rt') as input_fp:
            scanner = parse_mitab_file(input_fp, full_mitab_iterator)
            evdgrp_list = _offsets_by_bait_evidence_consumer(scanner)

            for offsets in evdgrp_list:
                scanner = parse_mitab_file(input_fp, partial_mitab_iterator,
                                           (offsets,))
                evd = AssociationEvidence()
                j = 0
                for intr, _ in scanner:
                    evd.insert_complex(intr)
                    j += 1
                self.set_complex_id(evd)
                self.raw_evidences.append(evd)
                self.singletons[evd.complex_id] = evd

        # Must remove lambda as default_factory to be able to pickle
        for evd in self.raw_evidences:
            evd.proteins.default_factory = None

    def get_pubmed_stats(self):

        counts = defaultdict(int)
        nonlabeled = defaultdict(int)
        nl_proteins = defaultdict(int)
        baits = defaultdict(int)
        preys = defaultdict(int)

        for evd in self.raw_evidences:
            pmid = evd.pubmed_ids[0]
            counts[pmid] += 1
            num_baits = sum(1 for p in evd.proteins \
                            if evd.proteins[p][0] == evd.bait_code)
            baits[pmid] += num_baits
            num_preys = sum(1 for p in evd.proteins \
                            if evd.proteins[p][0] == evd.prey_code)
            preys[pmid] += num_preys
            if not num_baits and not num_preys:
                nonlabeled[pmid] += 1


        pmids = sorted(counts.keys(), key=lambda pmid: counts[pmid])

        stats = []
        for pmid in pmids:
            stats.append((pmid,
                          counts[pmid],
                          nonlabeled[pmid],
                          baits[pmid],
                          preys[pmid]))
        return stats

    def write_pubmed_stats(self, fp):

        lines = self.get_pubmed_stats()
        for line in lines:
            fp.write('\t'.join(map(str, line)))
            fp.write('\n')

    def associations(self):

        asc = sorted(chain(self.clusters.itervalues(),
                           self.singletons.itervalues()),
                     key=attrgetter('_id'))
        return asc

    def write_tab(self, fp):

        for evd in self.associations():
            evd.write_line(fp)

    def write_tab2(self, fp):

        for evd in self.associations():
            evd.write_line2(fp)

    def write_details(self, fp):

        for evd in self.raw_evidences:
            evd.report_details(fp)

    def merge_clusters(self, complex_ids):

        evidences = []
        if len(complex_ids) == 1:
            cls = self[complex_ids[0]]
        else:
            for clid in complex_ids:
                if clid in self.singletons:
                    evidences.append(self.singletons[clid])
                    del self.singletons[clid]
                elif clid in self.clusters:
                    evidences.extend(self.clusters[clid].evidences)
                    del self.clusters[clid]
                else:
                    msg = "Cluster ID {} does not exist.".format(clid)
                    raise RuntimeError(msg)
            cls = AssociationCluster(evidences)
            self.set_merged_id(cls)
            self.clusters[cls.complex_id] = cls
        return cls.complex_id

    def merge_all_identical_entries(self):

        assert len(self.clusters) == 0
        unvisited = set(self.singletons.iterkeys())
        while len(unvisited):
            clid1 = unvisited.pop()
            M1 = self.singletons[clid1].member_weights()
            complex_ids = [clid1]
            for clid2 in unvisited:
                M2 = self.singletons[clid2].member_weights()
                if M1 == M2:
                    complex_ids.append(clid2)
            if len(complex_ids) > 1:
                self.merge_clusters(complex_ids)
                for clid in complex_ids[1:]:
                    unvisited.remove(clid)

    def fix_rare_proteins(self, min_complexes_per_protein=4):

        assert len(self.clusters) == 0
        p2c_counts = defaultdict(int)
        for evd in self.singletons.itervalues():
            for prot in evd.proteins:
                p2c_counts[prot] += 1

        removed_proteins = set(prot for prot in p2c_counts \
                               if p2c_counts[prot] < min_complexes_per_protein)
        affected_complexes = 0
        fixed_instances = 0
        fixed_complexes = 0
        removed_baits = 0
        removed_small = 0
        for complex_id, evd in self.singletons.items():
            fixed_prots = set()
            rare_bait = False
            for prot in evd.proteins:
                erole_code = evd.proteins[prot][0]
                if prot in removed_proteins:
                    fixed_prots.add(prot)
                    if erole_code == evd.bait_code:
                        rare_bait = True

            if rare_bait:
                affected_complexes += 1
                self.remove_complex(complex_id, 'B')
                removed_baits += 1
                fixed_instances += len(fixed_prots)

            elif fixed_prots:
                affected_complexes += 1
                # Must make proper deep copy of proteins
                new_proteins = {}
                for prot in evd.proteins:
                    erole_code, weight = evd.proteins[prot]
                    if prot in fixed_prots:
                        weight = 0.0
                    new_proteins[prot] = [erole_code, weight]
                total_weight = sum(v[1] for v in new_proteins.itervalues())

                if total_weight < 4:
                    self.remove_complex(complex_id, 'S')
                    removed_small += 1
                else:
                    # Shallow copy of evd - shared everything but complex_id and
                    # proteins
                    new_evd = copy.copy(evd)
                    new_evd.proteins = new_proteins
                    self.set_fixed_id(new_evd, complex_id)
                    self.singletons[new_evd.complex_id] = new_evd
                    self.remove_complex(complex_id, 'F')
                    fixed_complexes += 1
                fixed_instances += len(fixed_prots)



        msg = ("Found {} proteins with less than {} instances.\n"
               "  Affected complexes: {}\n"
               "  Complexes removed due to rare bait: {}\n"
               "  Complexes removed as being too small after fixing: {}\n"
               "  Fixed complexes: {}\n"
               "  Total removed instances: {}"
               )

        print msg.format(len(removed_proteins), min_complexes_per_protein,
                         affected_complexes, removed_baits, removed_small,
                         fixed_complexes, fixed_instances)

    def filter_pubmed_ids(self, pmids):

        assert len(self.clusters) == 0
        rejected = 0
        for pmid in pmids:
            for complex_id, evd in self.singletons.items():
                if pmid in evd.pubmed_ids:
                    self.removed[complex_id] = self.singletons.pop(complex_id)
                    rejected += 1
        print "Pubmed ID filter: removed {} complexes.".format(rejected)

    def filter_inconsistent_clusters(self, avg_weight_cutoff=0.0):

        rejected = 0
        for complex_id, cls in self.clusters.items():
            w = sum(v[1] for v in cls.proteins.itervalues())
            w /= len(cls.proteins)
            if w < avg_weight_cutoff:
                print (complex_id,
                       [evd.complex_id for evd in cls.evidences],
                       w)
                self.removed[complex_id] = self.clusters.pop(complex_id)
                rejected += 1
        print "Consistency filter: removed {} complexes.".format(rejected)

    def get_avg_complex_size(self):

        w = sum(len(cls.proteins) for cls in self.clusters.itervalues())
        w += sum(len(evd.proteins) for evd in self.singletons.itervalues())
        n = len(self.clusters) + len(self.singletons)
        return 1.0 * w / n

    def get_p2c_cum_mass(self):

        p2c_counts = defaultdict(int)
        for evd in self.associations():
            for p in evd.proteins:
                if evd.proteins[p][1] > 0.0:
                    p2c_counts[p] += 1
        count_freqs = Counter(p2c_counts.itervalues())
        counts = np.array(sorted(count_freqs.iterkeys()), dtype='d')
        mass = np.array([k*count_freqs[k] for k in counts], dtype='d')
        cum_mass = mass.cumsum()
        rel_cum_mass = cum_mass / cum_mass[-1]
        return counts, cum_mass, rel_cum_mass

    def get_weighted_median_pc_degree(self):

        x, y, z = self.get_p2c_cum_mass()
        cutoff = x[z<=0.5][-1]
        return cutoff

    def plain_weights(self):

        # Enumerate all proteins
        complexes = []
        used_prots = set()

        associations = self.associations()
        for cls in associations:
            complexes.append(cls.complex_id)
            for p in cls.proteins:
                used_prots.add(p)

        proteins = sorted(used_prots)
        protein2ix = dict((p, i) for i, p in enumerate(proteins))
        X = np.zeros((len(complexes), len(proteins)), dtype='d')

        for i, cls in enumerate(associations):
            for p in cls.proteins:
                j = protein2ix[p]
                X[i, j] = cls.proteins[p][1]

        return proteins, complexes, X


def initial_import(ppiTrim_file, output_file, pmids_file=None,
                   min_complexes_per_protein=None):

    if pmids_file is None:
        pmids = []
    else:
        with open(pmids_file, 'rt') as fp:
            pmids = [line.strip() for line in fp]

    adb = AssociationDatabase()
    adb.insert_evidences(ppiTrim_file)
    adb.filter_pubmed_ids(pmids)
    if min_complexes_per_protein is not None:
        adb.fix_rare_proteins(min_complexes_per_protein)
    adb.merge_all_identical_entries()
    adb.save(output_file)
