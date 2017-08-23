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
from dag_component import DAGComponent

_HELP = \
"""
Write a report based on association database.

SYNOPSIS:

    %s [OPTIONS] report_name input_file output_file
    %s -h|--help

ARGUMENTS:

    report_name                 Name of report function
    input_file                  DAG file
    output_file                 Output file

OPTIONS:

    -h, --help                  Print this message and exit

OUTPUT:

   Specified report
""" % (__file__, __file__)


_EMSG = "Insufficient arguments.\n\n" + _HELP


def node_details(dag, fp):

    for node_id in sorted(dag.nodes, key=lambda k: k[1:]):
        node = dag.nodes[node_id]
        node.report_details(fp)


def original_component_details(dag, fp):

    components = dag.get_original_components()
    components.sort(key=lambda cmpn: len(cmpn.nodes), reverse=True)
    for cmpn in components:
        cmpn.report_details(fp)


def node_enrichment(dag, fp, suffix):

    for node_id in sorted(dag.nodes, key=lambda k: k[1:]):
        node = dag.nodes[node_id]
        node.report_enrichment(suffix, fp)


def original_component_enrichment(dag, fp, suffix):

    components = dag.get_original_components()
    components.sort(key=lambda cmpn: len(cmpn.nodes), reverse=True)
    for cmpn in components:
        cmpn.report_enrichment(suffix, fp)


def stats(dag, fp):

    dag.report_stats(fp)


def summary(dag, fp, name):
    """ Summary as a row in a LaTeX table."""

    _nodes = dag.nodes.itervalues
    nodes_tot = len(dag.nodes)
    nodes_x = sum(1 for node in _nodes() if node.node_type == 'x')
    nodes_S = sum(1 for node in _nodes() if node.node_type == 'S')
    nodes_M = sum(1 for node in _nodes() if node.node_type == 'M')
    nodes_i = sum(1 for node in _nodes() if node.node_type == 'i')
    nodes_m = sum(1 for node in _nodes() if node.node_type == 'm')
    nodes_con = nodes_M + nodes_i + nodes_m

    R = dag.filtered_graph.transitive_reduction_graph()[0]
    edges_tot = dag.filtered_graph.count_edges()
    edges_reduced = R.count_edges()
    edges_nt_total = len(dag.non_transitive_relative_distances)
    edges_nt_essential = len(dag.get_essential_non_transitive_links())

    row = [nodes_tot, nodes_x, nodes_S, nodes_M, nodes_i, nodes_m, nodes_con,
           edges_tot, edges_reduced, edges_nt_total]

    fp.write(name.replace('_', '-'))
    fp.write(' & ')
    fp.write(' & '.join(map(str, row)))
    fp.write(' \\tabularnewline\n')


def _get_cmpn_func(components):

    def f(node):
        found = [cmpn.name for cmpn in components
                 if cmpn.find_node(node.node_id) is not None]
        if found:
            return found[0]
        return "ISOLATED"
    return f


def final_dag_data(dag, fp, file_prefix, use_original_components=True):

    edges_filename = file_prefix + '.edges.tab'
    nattr_filename = file_prefix + '.nattr.tab'
    shown_nodes = [node.node_id for node in dag.nodes.itervalues()]
    nattr_funcs = DAGComponent.get_nattr_funcs([])

    if use_original_components:
        components = dag.get_original_components()
    else:
        components = dag.get_active_components()
    nattr_funcs.append((_get_cmpn_func(components), '{}', 'Component'))

    S = dag.filtered_graph.subgraph(shown_nodes)
    with open(edges_filename, 'w') as fp2:
        S.write_edges(fp2)
    with open(nattr_filename, 'w') as fp2:
        S.write_node_attrs(dag, nattr_funcs, fp2)


def expanded_dag_data(dag, fp):

    S = dag.graph.expand_nodes(dag.complex2nr)
    S.write_edges(fp)


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

    dag = DAGAnalysis.load(input_file)
    f = globals()[func_name]

    if output_file != '-':
        fp = open(output_file, 'w')
    else:
        fp = sys.stdout

    f(dag, fp, *args[3:])
