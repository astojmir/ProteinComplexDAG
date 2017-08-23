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
from collections import defaultdict
from operator import itemgetter
from unweighted import single_source_shortest_path_length

class DictGraph(defaultdict):

    def __init__(self, initdict=None, writeable=True):

        defaultdict.__init__(self, dict)
        if initdict is not None:
            self._subgraph_copy(self, initdict)

        self._writeable = writeable
        if not self._writeable:
            self.default_factory = None

        # Cached items
        self._T = None
        self._U = None

    def __reduce__(self):

        state0 = super(DictGraph, self).__reduce__()
        state1 = (state0[0],     # Callable to create initial object
                  (None, ),      # Arguments to callable
                  self.__dict__, # Object state
                  state0[3],     # List iterator (empty)
                  state0[4],     # Dict iterator (set by parent)
                  )
        return state1

    @staticmethod
    def _get_all_nodes(initdict):

        all_nodes = set(initdict.iterkeys())
        for node1 in initdict:
            all_nodes |= set(initdict[node1].iterkeys())
        return all_nodes

    @property
    def nodes(self):

        return self._get_all_nodes(self)

    def count_edges(self):

        num_edges = sum(len(nbhd) for nbhd in self.itervalues())
        return num_edges

    def update(self, other):

        for node in other:
            self[node].update(other[node])

    @staticmethod
    def _subgraph_copy(dctgraph, other, nodes=None):

        H, G = dctgraph, other
        if nodes is None:
            _nodes = set(DictGraph._get_all_nodes(G))
        else:
            _nodes = set(nodes)

        for node1 in _nodes:
            nbhd = dict((node2, val) for node2, val in G[node1].iteritems() \
                        if node2 in _nodes)
            H[node1] = nbhd
        return H

    def subgraph(self, nodes):

        H = self.__class__()
        return self._subgraph_copy(H, self, nodes)

    def filtered_graph(self, excluded_nodes):

        included_nodes = self.nodes - set(excluded_nodes)
        return self.subgraph(included_nodes)

    def undirected_graph(self):

        U = self.__class__()
        for node1 in self:
            for node2, _ in self[node1].iteritems():
                U[node1][node2] = 1.0
                U[node2][node1] = 1.0
        return U

    @property
    def U(self):

        if self._U is None:
            self._U = self.undirected_graph()
        return self._U

    def transposed_graph(self):

        T = self.__class__()
        for node1 in self:
            for node2, val in self[node1].iteritems():
                T[node2][node1] = val
        return T

    @property
    def T(self):

        if self._T is None:
            self._T = self.transposed_graph()
        return self._T

    def out_degrees(self):

        _out_degrees = sorted(( (node, len(self[node])) for node in self ),
                              key=itemgetter(1))
        return _out_degrees

    def in_degrees(self):

        return self.T.out_degrees()

    def classify_nodes(self, all_complexes, excluded_complexes=None):

        srcs = set()
        dsts = set()
        excl = set() if excluded_complexes is None else set(excluded_complexes)

        for node1 in self:
            if len(self[node1]):
                srcs.add(node1)
                for node2 in self[node1].keys():
                    dsts.add(node2)

        minima = srcs - dsts - excl
        maxima = dsts - srcs - excl
        inner = (srcs & dsts) - excl
        singletons = set(all_complexes) - (srcs | dsts) - excl

        nt = {}
        nt.update((cmplx, 'x') for cmplx in excl)
        nt.update((cmplx, 'm') for cmplx in minima)
        nt.update((cmplx, 'M') for cmplx in maxima)
        nt.update((cmplx, 'i') for cmplx in inner)
        nt.update((cmplx, 'S') for cmplx in singletons)
        return nt

    def transitive_reduction_graph(self):
        """ Transitive reduction """
        # Assumes self is a DAG

        R = self.__class__(self)
        nodes = R.keys()
        deleted = set()
        for i, x in enumerate(nodes):
            x_nbhd = R[x].keys()
            for y in x_nbhd:
                if y in R:
                    for z in x_nbhd:
                        if z in R[y]:
                            deleted.add((x,z))

        for x, z in deleted:
            del R[x][z]
        return R, len(deleted)

    def transitive_closure_graph(self, new_link_val=1.0):
        """ Transitive closure """
        # Assumes self is a DAG

        E = self.__class__(self)
        num_added = 0
        nodes = E.keys()
        added = set()
        for y in nodes:
            for x in nodes:
                if y in E[x]:
                    for z in E[y]:
                        if z not in E[x]:
                            E[x][z] = new_link_val
                            added.add((x, y, z))
                            num_added += 1
        return E, num_added, added

    def is_transitively_closed(self):

        _, num_added, _ = self.transitive_closure_graph()
        return num_added == 0

    def undirected_diameter(self, cutoff=None):

        diam = 0
        U = self.U
        for src in U:
            paths = single_source_shortest_path_length(U, src, cutoff)
            diam = max(paths.itervalues())
            if cutoff is not None and diam >= cutoff:
                break
        return diam

    def connected_components(self):

        U = self.U
        nodes = set(U.keys())
        components = list()

        while nodes:
            w = nodes.pop()
            _visited = set([w])
            _unvisited = set(U[w].iterkeys())

            while _unvisited:
                u = _unvisited.pop()
                _visited.add(u)
                nodes.remove(u)
                if u in U:
                    _unvisited.update(v for v in U[u].keys() \
                                      if v not in _visited)
            components.append(sorted(_visited))
        components.sort(key=lambda _lst: (len(_lst), ''.join(_lst)))
        return components

    def get_component_data(self):

        components = self.connected_components()
        component_map = dict((node, cmpn) \
                             for cmpn in components \
                             for node in cmpn)
        closed = [self.subgraph(C).is_transitively_closed() \
                  for C in components]
        non_closed = [C for C, _cld in zip(components, closed) if not _cld]
        return components, component_map, closed, non_closed

    def forward_closure(self, sources):

        _visited = set([src for src in sources])
        _linked_sources = [src for src in sources if src in self]
        _unvisited = set(node for src in _linked_sources \
                              for node in self[src].iterkeys())

        while _unvisited:
            u = _unvisited.pop()
            _visited.add(u)
            if u in self:
                _unvisited.update(v for v in self[u].keys() if v not in _visited)
        return _visited

    def backward_closure(self, destinations):

        return self.T.forward_closure(destinations)

    def bf_components(self):
        """ Backward-forward closures of maxima """

        nt = classify_nodes(self, [])
        maxima = [node for node, code in nt.iteritems() if code == 'M']

        bf_closures = []
        for dst in maxima:
            bclosure = self.backward_closure([dst])
            bf_closures.append(self.forward_closure(bclosure))

        components = list()
        for i, clr1 in enumerate(bf_closures):
            match = False
            for clr2 in bf_closures[i+1:]:
                if clr1 == clr2:
                    match = True
                    break
            if not match:
                components.append(list(clr1))
        components.sort(key=lambda _lst: (len(_lst), ''.join(_lst)))
        return components

    def write_edges(self, fp=sys.stdout, edge_type='partOf', node_suffix=''):

        header = ['Node1', 'EdgeType', 'Node2', 'RelDist']
        for a in self:
            if len(self[a]):
                for b in self[a]:
                    rdist = self[a][b]
                    aa = a + node_suffix
                    bb = b + node_suffix
                    fields = [aa, edge_type, bb, '{:.4g}'.format(rdist)]
                    fp.write('\t'.join(fields))
                    fp.write('\n')
            else:
                fp.write(a + node_suffix)
                fp.write('\n')

    def write_node_attrs(self, container, attr_funcs, fp=sys.stdout):

        header = ['ID'] + [title  for _, _, title in attr_funcs]
        fp.write('\t'.join(header))
        fp.write('\n')
        for node_id in sorted(self.nodes):
            if node_id in container:
                obj = container[node_id]
                fp.write(node_id)
                fp.write('\t')
                for attr_func, fmt, _ in attr_funcs:
                    attr = attr_func(obj)
                    if attr is not None:
                        fp.write(fmt.format(attr))
                    fp.write('\t')
                fp.write('\n')

    def expand_nodes(self, complex2nr):

        H = self.__class__()
        for cmplx1 in complex2nr:
            mapped_ids = sorted(complex2nr[cmplx1])
            for i, nrc11 in enumerate(mapped_ids):
                for nrc12 in mapped_ids[i+1:]:
                    H[nrc11][nrc12] = 0.0
                    H[nrc12][nrc11] = 0.0

            if not len(self[cmplx1]):
                continue
            for cmplx2 in self[cmplx1].keys():
                for nrc1 in complex2nr[cmplx1]:
                    for nrc2 in complex2nr[cmplx2]:
                        H[nrc1][nrc2] = self[cmplx1][cmplx2]
        return H

    def get_edge_set(self):

        es = set((a, b) for a in self for b in self[a])
        return es

    def compare_edges(self, other, selected_nodes):

        num_ematches = 0
        emiss12 = []
        emiss21 = []

        for cmplx in selected_nodes:
            nb1 = self.get(cmplx, dict())
            nb2 = other.get(cmplx, dict())

            es1 = set((cmplx, cmplx1) for cmplx1 in nb1.keys() if nb1[cmplx1])
            es2 = set((cmplx, cmplx2) for cmplx2 in nb2.keys() if nb2[cmplx2])

            num_ematches += len(es1 & es2)

            emiss12 += list(es1 - es2)
            emiss21 += list(es2 - es1)

        emiss12.sort()
        emiss21.sort()

        res = (num_ematches,
               emiss12,
               emiss21,
               )
        return res
