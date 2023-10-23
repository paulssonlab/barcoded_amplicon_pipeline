import click
import networkx as nx
from cytoolz import compose
from gfapy import Gfa

from paulssonlab.util.cli import split_delimited_list

# SEE: https://stackoverflow.com/a/76464205
filter_gfa_options = compose(
    click.option(
        "-i",
        "--include",
        default=[],
        multiple=True,
        show_default=True,
        callback=split_delimited_list,
        help="Include only these segments.",
    ),
    click.option(
        "-I",
        "--include-prefix",
        default=[],
        multiple=True,
        show_default=True,
        callback=split_delimited_list,
        help="Include only the segments starting with this prefix.",
    ),
    click.option(
        "-x",
        "--exclude",
        default=[],
        multiple=True,
        show_default=True,
        callback=split_delimited_list,
        help="Exclude these segments.",
    ),
    click.option(
        "-X",
        "--exclude-prefix",
        default=[],
        multiple=True,
        show_default=True,
        callback=split_delimited_list,
        help="Exclude segments starting with this prefix.",
    ),
)


def _sign_to_angle_bracket(s, reverse=False):
    if reverse:
        if s == "-":
            return ">"
        elif s == "+":
            return "<"
        else:
            raise NotImplementedError
    else:
        if s == "-":
            return "<"
        elif s == "+":
            return ">"
        else:
            raise NotImplementedError


def gfa_to_dag(gfa):
    graph = nx.DiGraph()
    for edge in gfa.edges:
        graph.add_edge(
            f"{_sign_to_angle_bracket(edge.from_orient)}{edge.from_name}",
            f"{_sign_to_angle_bracket(edge.to_orient)}{edge.to_name}",
        )
        graph.add_edge(
            f"{_sign_to_angle_bracket(edge.from_orient, reverse=True)}{edge.from_name}",
            f"{_sign_to_angle_bracket(edge.to_orient, reverse=True)}{edge.to_name}",
        )
    return graph


def dag_endpoints(graph, wccs=None):
    sources = []
    sinks = []
    if wccs is None:
        wccs = nx.weakly_connected_components(graph)
    for wcc in wccs:
        subgraph = nx.subgraph(graph, wcc)
        sources.extend(
            node
            for node in subgraph.nodes
            if subgraph.in_degree[node] != 0 and subgraph.out_degree[node] == 0
        )
        sinks.extend(
            node
            for node in subgraph.nodes
            if subgraph.in_degree[node] == 0 and subgraph.out_degree[node] != 0
        )
    return sources, sinks


def gfa_endpoints(gfa):
    graph = gfa_to_dag(gfa)
    return dag_endpoints(graph)


def dag_forward_segments(graph, wccs=None):
    if wccs is None:
        wccs = nx.weakly_connected_components(graph)
    return list(max(wccs, key=lambda wcc: sum(s[0] == ">" for s in wcc)))


def gfa_forward_segments(gfa):
    graph = gfa_to_dag(gfa)
    return dag_forward_segments(graph)


def _exclude_gfa_line(line, segments):
    if line.RECORD_TYPE == "S":
        return line.name in segments
    else:
        fields = line.REFERENCE_FIELDS
        if line.NAME_FIELD:
            fields = [line.NAME_FIELD, *fields]
        return any(getattr(line, field) in segments for field in fields)


def filter_gfa(gfa, include=[], include_prefix=[], exclude=[], exclude_prefix=[]):
    if isinstance(include, str):
        include = [include]
    if isinstance(include_prefix, str):
        include = [include_prefix]
    if isinstance(exclude, str):
        include = [exclude]
    if isinstance(exclude_prefix, str):
        include = [exclude_prefix]
    include = set(include)
    exclude = set(exclude)
    segments_to_delete = []
    for segment in gfa.segments:
        delete = False
        if include or include_prefix:
            delete = True
            if segment.name in include:
                delete = False
            if any(segment.name.startswith(prefix) for prefix in include_prefix):
                delete = False
        if segment.name in exclude:
            delete = True
        if any(segment.name.startswith(prefix) for prefix in exclude_prefix):
            delete = True
        if delete:
            segments_to_delete.append(segment)
    segments_to_delete = set(segments_to_delete)
    lines = [
        str(line)
        for line in gfa.lines
        if not _exclude_gfa_line(line, segments_to_delete)
    ]
    return Gfa(lines)
