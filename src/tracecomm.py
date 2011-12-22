###########################################################################
# Copyright (C) 2011 Minh Van Nguyen <mvngu.name@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# http://www.gnu.org/licenses/
###########################################################################

# NOTE: This script MUST be run from Sage.
#
# Construct a digraph from data on community matches.  Use that digraph to
# determine traces of all dynamic communities.  Each trace is a path from a
# community birth to death, or up to the latest snapshot year.

from sage.all_cmdline import *
# Use optparse because Sage 4.7.1 has Python 2.6.  From Python 2.7 onwards,
# we could use argparse.
import optparse
import os
import sys

###################################
# helper functions
###################################

def add_missing_matches(E, M, V, R):
    """
    Update the edge set with fronts missing observations.  As we go along,
    we also update the vertex set.

    INPUT:

    - E -- edge set given as a dictionary of dictionaries.
    - M -- fronts with missing observations, given as a dictionary of lists.
    - V -- a vertex set to update.
    - R -- reverse mapping from vertex ID to year/community index.
     """
    for u in M:
        # add one front with missing observation
        i = 1
        v = u + "," + str(i)
        pchange = M[u][0][u]
        update_edge_set(E, V, R, (u, v, pchange), M)
        # add the remaining fronts
        for _ in M[u][1:]:
            v = u + "," + str(i)
            w = u + "," + str(i + 1)
            update_edge_set(E, V, R, (v, w, pchange), M)
            i += 1
        # has indirect continuation after missing observations?
        # take care of resumption of observations
        for k in E[u].keys():
            # community at y has observation at y + i
            if not k.startswith(u):
                v = u + "," + str(i)
                w = k
                pchange = E[u][k]
                del E[u][k]
                update_edge_set(E, V, R, (v, w, pchange), M)

def birth(year, dirname):
    """
    Get all the birth of dynamic communities in the given year.  The birth
    or beginning of a dynamic community is a community in year y without
    a matching to some community in previous years.

    INPUT:

    - year -- a year for which data are available.  If there are
      community births in the given year, we return a list of all such
      communities.  Otherwise, we return the empty list.
    - dirname -- path to directory with match data.
    """
    f = open(os.path.join(dirname, "birth", "birth-%d.dat" % year), "r")
    n = int(f.readline().strip())  # first line is number of births
    # zero births
    if n < 1:
        f.close()
        return []
    # at least one birth
    B = []
    for line in f:
        yi = line.strip()  # year/community index
        B.append(yi)
    f.close()
    return B

def extract_vertices_pchange(match):
    """
    Extract a pair of nodes from the given match data and the corresponding
    percentage change in community sizes.  The match data follows the format:

        year1 community_index1 year2 community_index2 pchange

    The year year1 is assumed to be the same as or earlier than year2.  The
    match data format above means that the community with index
    community_index1 at snapshot year year1 matches the community with index
    community_index2 at snapshot year year2.  The percentage change in
    community size from year1 to year2 is given by pchange.  The first node is
    taken to be the string "year1 community_index1" and the second node is the
    string "year2 community_index2".

    INPUT:

    - match -- a line of match data.
    """
    D = match.split()
    u = " ".join([D[0], D[1]])
    v = " ".join([D[2], D[3]])
    p = D[4]
    return (u, v, p)

def first_year(Year, dirname):
    """
    The first year with matching data.  We only return the index of that
    year from within the given list of years.

    INPUT:

    - Year -- list of years in nondecreasing order.
    - dirname -- path to directory with match data.
    """
    i = 0
    while not os.path.isfile(os.path.join(dirname, "match-%d.dat" % Year[i])):
        i += 1
        if i >= len(Year):
            raise ValueError("Expected matching data, but received none")
    return i

def get_years(fname):
    """
    Get all the snapshot years from the given file.  The resulting list of
    years will be sorted in nondecreasing order.

    INPUT:

    - fname -- file with all the snapshot years, one per line.

    OUTPUT:

    List of all snapshot years, sorted in nondecreasing order.
    """
    Year = []
    f = open(fname, "r")
    for line in f:
        Year.append(int(line.strip()))
    f.close()
    return sorted(Year)

def trace_community(v, G):
    """
    Trace all dynamic communities from the given start node.

    INPUT:

    - v -- the node from which to start tracing dynamic communities.  This
      must be a node in the given digraph.
    - G -- a directed graph representation of community matches.
    """
    # Get all shortest paths starting from v.  This will also give us
    # partial paths, i.e. paths that don't end with leaves.  The result is
    # a dictionary where keys are end points of paths and corresponding
    # values are lists of nodes in the paths from v to the keys.
    D = G.shortest_paths(v)
    # Prune partial paths, leaving only paths from v to leaves.  Such full
    # paths are traces of dynamic communities.
    for u in D.keys():
        if G.out_degree(u) > 0:
            del D[u]
    return D.values()

def trace_dynamic_communities(Year, G, dirname):
    """
    Trace all dynamic communities using the given digraph.  Each community
    is traced from its birth to its death, or possibly up to the latest
    year for which we have data.

    INPUT:

    - Year -- list of years for which data are available.  The list is
      assumed to be in nondecreasing order.
    - G -- a digraph representation of community matches.  We use this
      digraph representation to trace the evolution of each dynamic
      community.
    - dirname -- path to directory with match data.
    """
    assert(isinstance(G, DiGraph))
    # Get all community births.  Each of these is the start of some
    # dynamic community.
    i = first_year(Year, dirname)
    B = []
    for y in Year[i - 1:]:
        map(B.append, birth(y, dirname))
    # get full traces of all dynamic communities
    D = []
    for v in B:
        map(D.append, trace_community(v, G))
    return D

def update_edge_set(E, V, R, match, M):
    """
    Update the given edge set with the provided match data.

    WARNING:  This function will modify its arguments.

    INPUT:

    - E -- an edge set to update.  The edge set is assumed to be implemented
      as a dictionary of dictionaries.
    - V -- a vertex set to update.
    - R -- reverse mapping from vertex ID to year/community index.
    - match -- a line of match data read from a match file; or a tuple of
      match data.  In the case of tuple, we assume that all elements are
      distinct.
    - M -- dictionary of matches denoting fronts with missing observations.
    """
    u = None
    v = None
    p = None
    if isinstance(match, str):
        u, v, p = extract_vertices_pchange(match)
        # match denoting front with missing observation
        if u == v:
            if u in M:
                M[u].append({v:p})
            else:
                M.setdefault(u, [{v:p}])
            return
    elif isinstance(match, tuple):
        u, v, p = match
    else:
        raise ValueError("Invalid match data")
    # front with observation
    if u in E:
        E[u].setdefault(v, p)
    else:
        E.setdefault(u, {v:p})
        update_vertex_set(V, R, u)
    update_vertex_set(V, R, v)

def update_vertex_set(V, R, v):
    """
    Update the given vertex set with the provided vertex.

    WARNING:  This function will modify its arguments.

    INPUT:

    - V -- a vertex set to update.
    - R -- reverse mapping from vertex ID to year/community index.
    - v -- a node.
    """
    if v not in V:
        n = len(V)
        V.setdefault(v, n)
        R.setdefault(n, v)

def write_dynamic_communities(D, dirname, n):
    """
    Output to a file all traces of dynamic communities.  We are only
    interested in the full trace of each dynamic community, i.e. the trace
    from birth to death, or up to the latest snapshot year.  The full trace
    of each dynamic community is written on one line.

    INPUT:

    - D -- list of full traces of dynamic communities.  Each full trace is
      given as a list.
    - dirname -- path to directory with match data.
    - n -- the maximum length of any trace.  The length of a trace is the
      number of communities in it.  The maximum length of any trace is the
      number of years for which data are available.
    """
    if not os.path.isdir(os.path.join(dirname, "trace")):
        os.makedirs(os.path.join(dirname, "trace"))
    for i in range(1, n + 1):
        Di = filter(lambda x: len(x) == i, D)  # all traces of a given length
        f = open(os.path.join(dirname, "trace", "trace-%d.dat" % i), "w")
        f.write("%d\n" % len(Di))  # number of traces with given length
        for trace in Di:
            s = "->".join(trace)
            f.write(s + "\n")
        f.close()

def write_edges(E, V, dirname):
    """
    Write to a file the given edge list.

    INPUT:

    - E -- edge set.
    - V -- vertex set; this is the mapping from year/community index to
      vertex ID.
    - dirname -- path to directory with match data.
    """
    f = open(os.path.join(dirname, "edgelist.dat"), "w")
    for u in sorted(E.keys()):
        for v in sorted(E[u].keys()):
            f.write("%d %d %s\n" % (V[u], V[v], E[u][v]))
    f.close()

def write_vertices(V, dirname):
    """
    Write to a file the given mapping from year/community index to vertex ID.

    INPUT:

    - V -- vertex set.
    - dirname -- path to directory with match data.
    """
    f = open(os.path.join(dirname, "vertex-id.dat"), "w")
    for key in sorted(V.keys()):
        f.write("%s,%d\n" % (key, V[key]))
    f.close()

###################################
# script starts here
###################################

# setup parser for command line options
s = "Track events in life-cycle of dynamic communities.\n"
s += "Usage: %prog arg1 arg2 ..."
parser = optparse.OptionParser(usage=s)
parser.add_option("--size", metavar="integer", required=True, type=int,
                  help="minimum community size")
parser.add_option("--theta", metavar="float", required=True, type=float,
                  help="matching threshold")
parser.add_option("--delta", metavar="integer", required=True, type=int,
                  help="maximum allowable consecutive missing observations")
parser.add_option("--year", metavar="file",
                  help="file containing all snapshot years, one per line")
parser.add_option("--eventdir", metavar="path",
                  help="directory to read/write event data")
options, _ = parser.parse_args()

# get command line arguments & sanity checks
if ((options.size is None) or (options.theta is None)
    or (options.delta is None) or (options.year is None)
    or (options.eventdir is None)):
    raise optparse.OptionValueError(
        "All options must be used. Use -h for help.")
size = options.size
theta = options.theta
delta = options.delta
yearfname = options.year
eventdir = options.eventdir
if size < 1:
    raise ValueError("Invalid minimum size")
if theta < 0.0:
    raise ValueError("Invalid theta")
if delta < 0:
    raise ValueError("Invalid delta")

dirname = os.path.join(eventdir, "size-%d" % size,
                       "theta-%s_delta-%d" % (str(theta), delta))
if not os.path.isdir(dirname):
    raise ValueError("No such directory: %s" % dirname)
Year = get_years(yearfname)

# determine the first snapshot year having matching data
i = first_year(Year, dirname)
# mapping of year/community index to vertex ID
V = {}
# reverse mapping from vertex ID to year/community index
R = {}
# Mapping of matching data to edge.  Each edge follows the format (u,v;p),
# where each of u and v is a string of year/community index and p is the
# label for the edge and denotes the percentage change in community size.
E = {}
# matches denoting fronts with mising observations
M = {}
# all community births without matches
B = []
# build vertex and edge sets
for y in range(i, len(Year)):
    # process all match data
    f = open(os.path.join(dirname, "match-%d.dat" % Year[y]), "r")
    for line in f:
        match = line.strip()
        update_edge_set(E, V, R, match, M)
    f.close()
    # Process all community births.  In some cases, it is possible that a
    # new dynamic community is born at snapshot year y and subsequently has
    # no observations at all.  It is hence considered dead from y + 1
    # onwards.  A case is when delta = 0.
    for yi in birth(Year[y - 1], dirname):
        if yi not in V:
            B.append(yi)
            update_vertex_set(V, R, yi)
# Process all community births in the latest snapshot year.  In some cases,
# it is possible that a new dynamic community is born at snapshot year y and
# subsequently has no observations at all.  It is hence considered dead from
# y + 1 onwards.  A case is when delta = 0.
for yi in birth(Year[-1], dirname):
    if yi not in V:
        B.append(yi)
        update_vertex_set(V, R, yi)
add_missing_matches(E, M, V, R)
write_vertices(V, dirname)
write_edges(E, V, dirname)
# construct digraph
G = DiGraph(E, multiedges=False)
for yi in B:
    if yi not in G:
        G.add_vertex(yi)
# trace all dynamic communities
D = trace_dynamic_communities(Year, G, dirname)
write_dynamic_communities(D, dirname, len(Year))
