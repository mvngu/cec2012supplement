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

# This is a C extension for Python.  You should not use this directly.  Write
# a Python script to import specific functions from this module.

import csv
import os
import sys

###########################################################################
# Sanitize names in bibliographic datasets.
###########################################################################

cdef clean_hyphen(line):
    """
    Remove spaces before or after hyphens.

    WARNING:

    You should first run the function insert_period() before
    invoking this function.

    INPUT:

    - line -- a line of text from a dataset.

    OUTPUT:

    The same line but with spaces immediately before or after each hyphen
    removed.
    """
    L = line.strip()
    A = [u"A", u"B", u"C", u"D", u"E", u"F", u"G", u"H", u"I", u"J",
         u"K", u"L", u"M", u"N", u"O", u"P", u"Q", u"R", u"S", u"T",
         u"U", u"V", u"W", u"X", u"Y", u"Z"]
    for a in A:
        # space comes before hyphen, e.g. " -A"
        old_pattern = u" -" + a
        new_pattern = u"-" + a
        if old_pattern in L:
            L = L.replace(old_pattern, new_pattern)
        # space comes after hyphen & following initial, e.g. "A.- "
        old_pattern = a + u".- "
        new_pattern = a + u".-"
        if old_pattern in L:
            L = L.replace(old_pattern, new_pattern)
        # space comes after hyphen & following lower case letter, e.g. "a- "
        old_pattern = a.lower() + u"- "
        new_pattern = a.lower() + u"-"
        if old_pattern in L:
            L = L.replace(old_pattern, new_pattern)
    return L

def clean_names(infile):
    """
    Name cleansing for the given file.

    WARNING: This function modifies the file you feed it.

    INPUT:

    - infile -- the file on which to perform name cleansing.
    """
    outfile = infile
    tmpfile = "%s.tmp" % outfile
    fin = open(infile, "r")
    fout = open(tmpfile, "w")
    for line in fin:
        L = insert_period(line)
        L = clean_hyphen(L)
        fout.write("%s\n" % L)
    fin.close()
    fout.close()
    os.system("mv %s %s" % (tmpfile, outfile))

cdef insert_period(line):
    """
    Insert periods after initials.

    INPUT:

    - line -- a line of text from a dataset.

    OUTPUT:

    The same line but where each initial is immediately followed by a period.
    """
    L = line.strip()
    A = [u"A", u"B", u"C", u"D", u"E", u"F", u"G", u"H", u"I", u"J",
         u"K", u"L", u"M", u"N", u"O", u"P", u"Q", u"R", u"S", u"T",
         u"U", u"V", u"W", u"X", u"Y", u"Z"]
    # period after single initial
    for a in A:
        # initial of first name
        old_pattern = u"," + a + u" "
        new_pattern = u"," + a + u". "
        if old_pattern in L:
            L = L.replace(old_pattern, new_pattern)
        # middle initial
        old_pattern = u" " + a + u" "
        new_pattern = u" " + a + u". "
        if old_pattern in L:
            L = L.replace(old_pattern, new_pattern)
    # period after double initial, e.g. A-B
    for a in A:
        for b in A:
            # periods for both & following comma, e.g. ",A-B"
            old_pattern = u"," + a + u"-" + b + u" "
            new_pattern = u"," + a + u".-" + b + u". "
            if old_pattern in L:
                L = L.replace(old_pattern, new_pattern)
            # periods for both, e.g. A-B
            old_pattern = u" " + a + u"-" + b + u" "
            new_pattern = u" " + a + u".-" + b + u". "
            if old_pattern in L:
                L = L.replace(old_pattern, new_pattern)
            # period for only first initial, e.g. A-B.
            old_pattern = u"," + a + u"-" + b + u". "
            new_pattern = u"," + a + u".-" + b + u". "
            if old_pattern in L:
                L = L.replace(old_pattern, new_pattern)
            # period for only second initial, e.g. A.-B
            old_pattern = u"," + a + u".-" + b + u" "
            new_pattern = u"," + a + u".-" + b + u". "
            if old_pattern in L:
                L = L.replace(old_pattern, new_pattern)
    return L

###########################################################################
# This is used to determine snapshots of an evolving network.  We consider
# yearly snapshots.  Work out all the authors in a yearly dataset.  We use
# names as provided and consider each unique name as a unique author. Then we
# sort the list of unique author names in alphabetical order.  This is
# followed by a mapping of each author name to a unique nonnegative integer.
###########################################################################

def build_snapshots(yearfname, datadir, edgelistdir):
    """
    Construct all the snapshots of an evolving network.  Each network snapshot
    is an edgelist, where nodes are author IDs.  This function also writes to
    file a mapping of author names to nonnegative integers.

    INPUT:

    - yearfname -- file containing all snapshot years, one per line.
    - datadir -- directory with data files, one for each year.  The data file
      for a year will be used to construct the edge list for that year.
    - edgelistdir -- directory to output edge list for each snapshot.  The
      mapping of author names to integer IDs will also be written to a file
      under this directory.
    """
    # mapping of each author to a unique nonnegative integer
    ID = {}
    # All unique edges.  This mean we don't have multiple edges, but we might
    # have self-loops.  When you process this edge list, ignore all self-loops.
    E = {}
    # get snapshot graph of network for each time step
    for year in get_years(yearfname):
        fname = os.path.join(datadir, "%d.dat" % year)
        A = get_names(fname)
        handle_names(A)
        to_integer(A, ID)
        edge_list(fname,
                  os.path.join(edgelistdir, "edgelist-%d.dat" % year),
                  ID, E)
    # write out mapping of author names to nonnegative integers
    A = sorted(ID.keys())
    f = open(os.path.join(edgelistdir, "author-id.dat"), "w")
    for a in A:
        s = a + "," + str(ID[a]) + "\n"
        f.write(s)
    f.close()

def edge_list(fdataset, fedgelist, D, E):
    """
    Create edges from authors of papers.  If a paper has n authors, then each
    of those n authors is connected to the other n - 1 authors via an edge.
    If a paper has only one author, then we create a self-loop.  However, note
    that when reading the resulting edge list, we should ignore all self-loops.

    WARNING:

    We will modify the argument E.

    INPUT:

    - fdataset -- the file containing the dataset to read.

    - fedgelist -- the file to which we write the edge list.

    - D -- dictionary of IDs.  Each ID is a mapping of an author to a
      nonnegative integer.  We use each such integer as the node name.

    - E -- edge list; use a dictionary to maintain unique edges.  We will
      extend this edge list with *new* edges only.

    OUTPUT:

    An edge list representing a coauthorship graph.  This edge list is written
    to the file fedgelist.
    """
    cdef i
    cdef j
    fin = open(fdataset, "r")
    reader = csv.reader(fin)
    for line in reader:
        # get the authors for this paper and their corresponding IDs
        A = map(lambda x: x.strip(), line[1:])
        handle_names(A)
        ID = map(lambda x: D[x], A)
        ID = sorted(ID)
        # lone author; self-loop
        if len(ID) == 1:
            a = str(ID[0]) + " " + str(ID[0])
            E.setdefault(a, True)
        # multiple authors; create a complete graph on those number of authors
        # for i <- 0,...,|A|-2
        #     for j <- i+1,...,|A|-1
        #         add edge A[i], A[j]
        elif len(ID) > 1:
            for i in range(len(ID) - 1):
                for j in range(i + 1, len(ID)):
                    a = str(ID[i]) + " " + str(ID[j])
                    E.setdefault(a, True)
        else:
            raise ValueError("Should have at least one author, but got none.")
    fin.close()
    fout = open(fedgelist, "w")
    for e in E:
        fout.write(e + "\n")
    fout.close()

# NOTE: This function is stolen from Sage 4.7.1, straight from the file
# sage-4.7.1/local/lib/python2.6/site-packages/sage/misc/flatten.py
# We just Cythonize the Python version.
cdef flatten(in_list, ltypes=(list, tuple), long max_level=sys.maxint):
    """
    Flattens a nested list.

    INPUT:

    - in_list -- a list or tuple
    - ltypes -- optional list of particular types to flatten
    - max_level -- the maximum level to flatten

    OUTPUT:

    a flat list of the entries of in_list

    EXAMPLES:

        sage: flatten([[1,1],[1],2])
        [1, 1, 1, 2]
        sage: flatten([[1,2,3], (4,5), [[[1],[2]]]])
        [1, 2, 3, 4, 5, 1, 2]
        sage: flatten([[1,2,3], (4,5), [[[1],[2]]]],max_level=1)
        [1, 2, 3, 4, 5, [[1], [2]]]
        sage: flatten([[[3],[]]],max_level=0)
        [[[3], []]]
        sage: flatten([[[3],[]]],max_level=1)
        [[3], []]
        sage: flatten([[[3],[]]],max_level=2)
        [3]

    In the following example, the vector isn't flattened because
    it is not given in the ltypes input.

        sage: flatten((['Hi',2,vector(QQ,[1,2,3])],(4,5,6)))
        ['Hi', 2, (1, 2, 3), 4, 5, 6]

    We give the vector type and then even the vector gets flattened:

        sage: flatten((['Hi',2,vector(QQ,[1,2,3])], (4,5,6)), ltypes=(list, tuple,sage.modules.vector_rational_dense.Vector_rational_dense))
        ['Hi', 2, 1, 2, 3, 4, 5, 6]

    We flatten a finite field.

        sage: flatten(GF(5))
        [0, 1, 2, 3, 4]
        sage: flatten([GF(5)])
        [Finite Field of size 5]
        sage: flatten([GF(5)], ltypes = (list, tuple, sage.rings.finite_rings.finite_field_prime_modn.FiniteField_prime_modn))
        [0, 1, 2, 3, 4]

    Degenerate cases:

        sage: flatten([[],[]])
        []
        sage: flatten([[[]]])
        []
    """
    cdef long current_level
    cdef long index
    cdef long len_v
    cdef long old_level
    index = 0
    current_level = 0
    len_v = -1
    new_list = [x for x in in_list]
    level_list = [0] * len(in_list)
    while index < len(new_list):
        while (isinstance(new_list[index], ltypes)
               and current_level < max_level):
            v = list(new_list[index])
            len_v = len(v)
            new_list[index : index + 1] = v
            old_level = level_list[index]
            level_list[index : index + 1] = [0] * len_v
            if len_v != 0:
                current_level += 1
                level_list[index + len_v - 1] = old_level + 1
            else:
                current_level -= old_level
                index -= 1
                break
        # If len_v==0, then index points to a previous element, so we
        # don't need to do anything.
        if len_v != 0:
            current_level -= level_list[index]
        index += 1
    return new_list

def get_names(fname):
    """
    Get all author names.  This would likely result in duplicates.  The given
    file is assumed to be in CSV format, where each line follows the
    specification:

        yyyy,author1,author2,...

    INPUT:

    - fname -- the file name containing the dataset to read.

    OUTPUT:

    List of all author names from dataset in fname.
    """
    A = []
    with open(fname, "rb") as f:
        reader = csv.reader(f)
        for row in reader:
            a = map(lambda x: x.strip(), row[1:])
            A.append(a)
    return flatten(A)

cdef handle_names(names):
    """
    Some basic name handling for author names on a preprint.  Add more fine
    grained rules here to weed out possibly duplicate names.  For example, if
    an author name is given as any of the following forms

        Full Name, jr
        Full Name, Jr
        Full Name Jr

    or any other variations, we convert the name to follow the form

        Full Name jr

    WARNING:

    This function will modify its argument.

    INPUT:

    - names -- list of author names on a preprint.

    OUTPUT:

    Ensure each name conforms to some rules for the purposes of our study.
    """
    cdef int i
    # the string "Jr", "jr", "jr.", or "Jr." appears as separate list element
    while (("jr" in names)
           or ("Jr" in names)
           or ("jr." in names)
           or ("Jr." in names)):
        i = -1
        try:
            i = names.index("jr")
        except:
            pass
        try:
            i = names.index("Jr")
        except:
            pass
        try:
            i = names.index("jr.")
        except:
            pass
        try:
            i = names.index("Jr.")
        except:
            pass
        if i == -1:
            raise ValueError("Problem with name handling.")
        names[i - 1] += " jr"
        del names[i]
    # the string "Jr.", "jr.", or "Jr" appears at end of a name
    for i in range(len(names)):
        if ("Jr." in names[i]) or ("jr." in names[i]) or ("Jr" in names[i]):
            a = names[i].split()
            a[-1] = "jr"
            names[i] = "".join(a)

cdef to_integer(A, D):
    """
    Map unique author names to unique nonnegative integers.

    INPUT:

    - A -- list of author names.

    - D -- a mapping to extend.  This is essentially a dictionary.  New author
      names are added to this mapping.

    OUTPUT:

    A mapping of unique author names to nonnegative integers.  The input list
    A is first removed of all duplicates.  The resulting list is then sorted
    in alphabetical order.  Based on this sorted list, we map author names to
    integers and extend the dictionary D with new unique authors.
    """
    L = sorted(list(set(A)))
    for author in L:
        D.setdefault(author, len(D))

###########################################################################
# Communities with a specified minimum size.
###########################################################################

def comm_minsize(yearfname, int size, comdir, dirname):
    """
    Let y be the first year in which we have communities with the given minimum
    size minsize.  For some t > y, it is possible that the year t does not have
    any communities with minimum size minsize.  In that case, we create the
    file comm-t.dat and have zero as the very first line of that file.  The
    file would contain nothing else.  The variable minyear stores the earliest
    year in which we have communities with size at least minsize.

    INPUT:

    - yearfname -- file containing all snapshot years, one per line.
    - size -- minimum community size.
    - comdir -- directory with communities per snapshot year.
    - dirname -- directory to write all communities with the given minimum
      size.
    """
    cdef int minyear
    cdef int year
    minyear = -1
    for year in get_years(yearfname):
        C = read_communities(year, comdir, size)
        if (len(C) > 0) and (minyear < 0):
            minyear = year
        write_communities(C, year, minyear, dirname)

cdef write_communities(comm, int year, int minyear, dirname):
    """
    Write the given list of communities to a file.  The file name follows the
    format "comm-yyyy.dat", where "yyyy" is a four-digit year.  It is
    written under the directory named by dirname.  The very first line of
    the output file lists the number of communities.  Each subsequent line
    lists all the node IDs belonging a community.

    INPUT:

    - comm -- list of communities.  Each element is itself a list of
      community IDs.  All IDs belonging to a community will be written on
      one line.  If there are no communities in the given snapshot year, we
      do nothing.
    - year -- the snapshot year.
    - minyear -- the earliest year in which we have communities with a given
      minimum size.
    - dirname -- name of the directory under which to output the given
      communities.
    """
    if len(comm) < 1:
        # Haven't set the minimum year in which we have communities with
        # given minimum size.
        if minyear < 0:
            return
        # we have set it
        assert(minyear < year)
    fname = os.path.join(dirname, "comm-%d.dat" % year)
    f = open(fname, "w")
    f.write(str(len(comm)) + "\n")
    for C in comm:
        s = " ".join(map(str, C))
        f.write(s + "\n")
    f.close()

###########################################################################
# Community matches with direct/indirect continuation.  Missing observations
# are allowed, depending on the command line argument.
###########################################################################

cdef get_years(fname):
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

cdef double jaccard(front, step):
    """
    The Jaccard coefficient between the given two communities.  The Jaccard
    coefficient measures the degree to which two communities are similar.

    INPUT:

    - front -- a front community.
    - step -- a step community.
    """
    cdef double m
    cdef double n
    A = set(front)
    B = set(step)
    n = float(len(A.intersection(B)))
    m = float(len(A.union(B)))
    return n / m

def match(yearfname, dirname, comdir, double theta, int delta):
    """
    Match step communities in consecutive time steps.

    INPUT:

    - yearfname -- file containing all snapshot years, one per line.
    - dirname -- directory under which to write match data.
    - comdir -- directory under which community structures are found.
    - theta -- the matching threshold.
    - delta -- maximum allowable number of missing observations.
    """
    cdef int a
    cdef int cindex
    cdef int fi
    cdef int i
    cdef int mi
    cdef int y
    cdef int yindex
    cdef double p
    cdef double sim
    # Let the time step be t.  The communities at time t are referred to as the
    # step communities.  Those at time t - 1 are the front communities.  Any
    # communities at time < t - 1 without observations are called missing
    # communities.  Missing communities are lumped together with the fronts at
    # time t.  We allow communities to have missing observations for a number
    # of consecutive time steps.  If the number of missing observations exceed
    # a fixed missing threshold, then we consider the corresponding missing
    # community as dead.
    # key := "year community_index"
    # value := no. missing observations
    Missing = {}
    Year = get_years(yearfname)
    i = first_year(Year, comdir)
    for y in range(i + 1, len(Year)):
        Front = read_communities(Year[y - 1], comdir)
        # no front communities to work with; move on to next snapshot year
        if len(Front) == 0:
            continue
        Step = read_communities(Year[y], comdir)
        mfile = open(os.path.join(dirname, "match-%d.dat" % Year[y]), "w")
        # try to match each front community with each step community
        for i, Fi in enumerate(Front):
            fi = 0  # assume front F_i has zero matches so far
            for a, Ca in enumerate(Step):
                sim = jaccard(Fi, Ca)
                if sim > theta:
                    fi += 1
                    p = pchange(Fi, Ca)
                    mfile.write("%d %d %d %d %.10f\n" %
                                (Year[y - 1], i, Year[y], a, p))
            # keep track of fronts with missing observations
            if fi == 0:
                # Start with zero missing count because we will re-match this
                # missing community with each of the step communities.
                key = "%d %d" % (Year[y - 1], i)
                Missing.setdefault(key, 0)
        # try to match each missing community with each step community
        for M in Missing.keys():
            yindex, cindex = tuple(map(int, M.split()))
            Mi = read_one_community(yindex, cindex, comdir)
            mi = 0  # assume missing community M_i has zero matches so far
            for a, Ca in enumerate(Step):
                sim = jaccard(Mi, Ca)
                if sim > theta:
                    mi += 1
                    p = pchange(Mi, Ca)
                    mfile.write("%d %d %d %d %.10f\n" %
                                (yindex, cindex, Year[y], a, p))
            # keep track of missing communities not yet considered dead
            if mi == 0:
                Missing[M] += 1
                if Missing[M] <= delta:
                    p = 0.0
                    mfile.write("%d %d %d %d %.10f\n" %
                                (yindex, cindex, yindex, cindex, p))
            # missing community M_i matches some step community C_a
            # so M_i no longer has missing observations
            else:
                del Missing[M]
        mfile.close()
        # remove any dead communities
        for M in Missing.keys():
            if Missing[M] > delta:
                del Missing[M]

cdef double pchange(front, step):
    """
    The percentage change in community sizes for the given two communities.

    INPUT:

    - front -- a front community.
    - step -- a step community.
    """
    cdef double m
    cdef double n
    m = float(len(front))
    n = float(len(step))
    return (n - m) / m

cdef read_communities(int year, dirname, int size=1):
    """
    Read all communities from a snapshot with the given year.  The communities
    of a snapshot are stored in a file.  The very first line of the file lists
    the number of communities.  Each subsequent line lists the node IDs making
    up a community.

    INPUT:

    - year -- the snapshot year.
    - dirname -- the directory under which communities are stored.
    - size -- (default: 1) the minimum community size.  The default value of
      1 means that we accept communities of sizes >= 1.  If we want all
      communities of a given minimum size, then set 'size' to some other
      positive integer.
    """
    cdef int ncomm
    fname = os.path.join(dirname, "comm-%d.dat" % year)
    f = open(fname, "r")
    ncomm = int(f.readline())
    Comm = []
    for line in f:
        C = line.strip().split()
        if len(C) >= size:
            Comm.append(C)
    f.close()
    if size == 1:
        assert(ncomm == len(Comm))
    return Comm

cdef read_one_community(int year, int i, dirname):
    """
    Get one community from a snapshot with the given year.

    INPUT:

    - year -- the snapshot year.
    - i -- the community index.  In general, community index starts from zero.
    - dirname -- the directory under which communities are stored.
    """
    cdef int ncomm
    assert(i >= 0)
    fname = os.path.join(dirname, "comm-%d.dat" % year)
    f = open(fname, "r")
    ncomm = int(f.readline().strip())  # first line is number of communities
    assert(i < ncomm)
    for _ in range(i):                 # ignore the first i - 1 communities
        f.readline()
    line = f.readline()
    f.close()
    return line.strip().split()

###########################################################################
# Track key events in the life-cycle of dynamic communities.  Here, we
# assume that we already have data on community matches.
###########################################################################

cdef int count_communities(int year, dirname):
    """
    The number of communities in the given snapshot year.

    INPUT:

    - year -- the snapshot year.
    - dirname -- name of the directory from which to read file storing
      communities.
    """
    cdef int ncomm
    f = None
    ncomm = 0
    fname = os.path.join(dirname, "comm-%d.dat" % year)
    try:
        f = open(fname, "r")
        ncomm = int(f.readline().strip())
        f.close()
    except IOError:
        # file doesn't exist
        ncomm = 0
    return ncomm

cdef dead_fronts(fname, D, int delta):
    """
    Get all the dead fronts from trace data in the given file.  A dead front
    should follow the following format:

        yyyy ID,n

    Here, yyyy is the year, ID is the community ID, and n = delta is the
    number of missing observations.  Together, "yyyy ID" means a community
    identified as ID in year yyyy.  The sum yyyy + n means that community
    "yyyy ID" died during the year yyyy + n.  It should be noted that we
    must have n = delta.

    WARNING: This function modifies its arguments.

    INPUT:

    - fname -- path to a file with trace data.  A trace is merely a record of
      community matches that together is the timeline of some dynamic
      community.
    - D -- a dictionary storing dead fronts corresponding to each year.  Thus
      D[y] is a list of all dynamic communities that died during year y.
    - delta -- the maximum number of missing observations allowed.
    """
    cdef int n
    cdef int year
    f = open(fname, "r")
    f.readline()  # ignore first line, which is #dynamic communities
    for line in f:
        front = line.strip().split("->")[-1]
        if "," not in front:
            continue
        commID, nstr = front.split(",")
        n = int(nstr)
        if n < delta:
            continue
        ystr, _ = commID.split(" ")
        year = int(ystr)
        D[year + n].append(front)
    f.close()

cdef empty_list_of_lists(int n):
    """
    A list of n empty lists.  That is, we want a list of n elements, where
    each element is itself an empty list.

    INPUT:

    - n -- how many empty lists we want.
    """
    return [[] for _ in range(n)]

cdef etype_comm(fname):
    """
    Get all the communities from the given file that satisfy a required
    life-cycle event.  The event might be contraction, expansion, merge,
    resume, split, stable, or stagnant.

    INPUT:

    - fname -- path to a file with data for a required event.
    """
    cdef int IDA
    cdef int IDB
    cdef int nevent
    cdef int yearA
    cdef int yearB
    f = None
    try:
        f = open(fname, "r")
    except IOError:
        return {}
    # number of occurrences of the required event
    nevent = int(f.readline().strip())
    if nevent < 1:
        return {}
    Event = {}
    # mental model is
    # Event = {commB0: [commA0, commA1, ...], commB1: [...], ...}
    # If the event is contraction, then community (yearA,commA) contracts
    # into community (yearB,commB).
    # If the event is expansion, then community (yearA,commA) expands into
    # community (yearB,commB).
    # If the event is merge, then community (yearA,commA) merge into
    # community (yearB,commB).
    # If the event is resume, then community (yearA,commA) has a resume
    # observation in community (yearB,commB).
    # If the event is split, then community (yearA,commA) split into
    # community (yearB,commB).
    # If the event is stable, then community (yearA,commA) has stable change
    # in community (yearB,commB).
    # If the event is stagnant, then community (yearA,commA) has stagnant
    # change in community (yearB,commB).
    yearA = 0  # index where first year is located
    IDA = 1    # index where ID of first community is located
    yearB = 2  # index where second year is located
    IDB = 3    # index where ID of second community is located
    for line in f:
        dat = line.strip().split(" ")[:-1]  # get everything except last elem
        commA = "%s %s" % (dat[yearA], dat[IDA])
        commB = "%s %s" % (dat[yearB], dat[IDB])
        if commB in Event:
            Event[commB].append(commA)
        else:
            Event.setdefault(commB, [commA])
    f.close()
    return Event

def events_except_death(yearfname, dirname, comdir, double gamma,
                        double kappa):
    """
    Obtain data for all events, except death.  We don't consider death here.

    INPUT:

    - yearfname -- file containing all snapshot years, one per line.
    - dirname -- directory under which to write event data.
    - comdir -- directory under which can be found communities with a given
      minimum size.
    - gamma -- the growth threshold.
    - kappa -- the contraction threshold.
    """
    cdef int Cs
    cdef int Ct
    cdef int i
    cdef int k
    cdef int ncomm
    cdef int nFront
    cdef int nStep
    cdef int s
    cdef int t
    cdef int y
    Year = get_years(yearfname)
    # determine the first snapshot year having communities
    i = first_year(Year, comdir)
    # key events for the very first snapshot with a community
    ncomm = count_communities(Year[i], comdir)
    Birth = [str(Year[i]) + " " + str(k) for k in range(ncomm)]
    write_event(Year[i], dirname, "birth", Birth)
    write_event(Year[i], dirname, "expand", [])
    write_event(Year[i], dirname, "contract", [])
    write_event(Year[i], dirname, "split", [])
    write_event(Year[i], dirname, "merge", [])
    write_event(Year[i], dirname, "stable", [])
    write_event(Year[i], dirname, "stagnant", [])
    write_event(Year[i], dirname, "direct", [])
    write_event(Year[i], dirname, "missing", [])
    write_event(Year[i], dirname, "resume", [])

    # key events for subsequent snapshots
    for y in range(i + 1, len(Year)):
        # prepare for action
        nFront = count_communities(Year[y - 1], comdir)
        nStep = count_communities(Year[y], comdir)
        C = empty_list_of_lists(nStep)   # [[], [], ..., []]
        F = empty_list_of_lists(nFront)  # [[], [], ..., []]
        FR = {}  # all fronts that are communities with resuming observations
        Birth = []
        Expand = []
        Contract = []
        Split = []
        Merge = []
        Stable = []
        Stagnant = []
        Direct = [] # direct continuation: F_i matches exactly C_a & vice versa
        Missing = []  # missing observations
        Resume = []   # communities w/o observations having observations again
        Match = {}    # mapping of matches to corresponding percentage changes
        # get all matches and update counts of matches between fronts and steps
        f = open(os.path.join(dirname, "match-%d.dat" % Year[y]), "r")
        for line in f:
            match = line.strip()
            if is_missing_observation(match):
                update_missing(match, Missing)
                continue
            update_stagnant(match, Stagnant)
            update_stable(gamma, kappa, match, Stable)
            update_expand(gamma, match, Expand)
            update_contract(kappa, match, Contract)
            s, Cs, t, Ct, p = year_index_pchange(match)
            assert(Year[y] == t)
            if is_resuming_observation(match):
                update_resuming(match, Resume)
                key = (s, Cs)
                if key in FR:
                    FR[key].append(Ct)
                else:
                    FR.setdefault(key, [Ct])
                C[Ct].append((s, Cs))
            else:
                assert(Year[y - 1] == s)
                # Use tuples because a front or a resuming community can match
                # a step community.  A 2-tuple allows us to keep track of the
                # snapshot year and community index for each match.
                C[Ct].append((s, Cs))
                F[Cs].append(Ct)
            key = "%d %d %d %d" % (s, Cs, t, Ct)
            Match.setdefault(key, p)
        f.close()
        # determine birth, split, merge; don't consider death here
        update_events(Year[y], F, FR, C, Match, Birth, Split, Merge, Direct)
        write_event(Year[y], dirname, "birth", Birth)
        write_event(Year[y], dirname, "expand", Expand)
        write_event(Year[y], dirname, "contract", Contract)
        write_event(Year[y], dirname, "split", Split)
        write_event(Year[y], dirname, "merge", Merge)
        write_event(Year[y], dirname, "stable", Stable)
        write_event(Year[y], dirname, "stagnant", Stagnant)
        write_event(Year[y], dirname, "direct", Direct)
        write_event(Year[y], dirname, "missing", Missing)
        write_event(Year[y], dirname, "resume", Resume)

    # post-processing and clean-ups of event data
    Event = ["contract", "expand", "merge", "resume", "split", "stable",
             "stagnant"]
    for etype in Event:
        for y in Year:
            fname = os.path.join(dirname, etype, "%s-%d.dat" % (etype, y))
            write_etype_comm(fname, etype_comm(fname))

def events_only_death(yearfname, dirname, int delta):
    """
    Obtain data for death only.

    INPUT:

    - yearfname -- file containing all snapshot years, one per line.
    - dirname -- directory under which to write event data.
    - delta -- maximum allowable consecutive missing observations.
    """
    cdef int i
    cdef int maxlen
    Year = get_years(yearfname)
    Death = {}
    for y in Year:
        Death.setdefault(y, [])
    maxlen = len(Year)
    # get all dead fronts during each year and write to file
    for i in range(1, maxlen + 1):
        fname = os.path.join(dirname, "trace", "trace-%d.dat" % i)
        dead_fronts(fname, Death, delta)
    write_death(Death, dirname)

cdef int first_year(Year, dirname):
    """
    The first snapshot year having communities.  Depending on the minimum
    required size of each community, some snapshot years might be neglected
    because each of the communities in that year has size less than the
    required minimum.  We only return the index of that year from within the
    given list of years.

    INPUT:

    - Year -- list of all the snapshot years, sorted in increasing order.
    - dirname -- name of the directory from which to read file storing
      communities.
    """
    cdef int i
    cdef int ncomm
    i = 0
    ncomm = count_communities(Year[i], dirname)
    while ncomm < 1:
        i += 1
        if i >= len(Year):
            raise ValueError("Expected community data, but received none")
        ncomm = count_communities(Year[i], dirname)
    return i

cdef bint is_missing_observation(match):
    """
    Whether the given match is a case of some dynamic community missing an
    observation in some snapshot.

    INPUT:

    - match -- a community match between a front and a step.  The match
      follows this format:

          year1 comm_index1 year2 comm_index2 pchange

      The key "year1" is the snapshot year for a front community and
      "comm_index1" is a community index.  Similarly, "year2" is the snapshot
      year for a step community and "comm_index2" is a community index.  The
      key "pchange" is the percentage change in community size.  For a missing
      observation, the match follows the format:

          year1 comm_index1 year1 comm_index1 0.0
    """
    cdef int Cs
    cdef int Ct
    cdef int s
    cdef int t
    s, Cs, t, Ct, _ = year_index_pchange(match)
    if s == t:
        if Cs == Ct:
            return True
        else:
            raise ValueError("Expected a missing observation: %s" % match)
    return False

cdef bint is_resuming_observation(match):
    """
    Whether the given match is a case of some dynamic community that has
    been missing observations in the previous snapshot, but now has an
    observation in the current snapshot.

    INPUT:

    - match -- a community match between a front and a step.  The match
      follows this format:

          year1 comm_index1 year2 comm_index2 pchange

      The key "year1" is the snapshot year for a front community and
      "comm_index1" is a community index.  Similarly, "year2" is the snapshot
      year for a step community and "comm_index2" is a community index.  The
      key "pchange" is the percentage change in community size.  For a
      resuming observation, the absolute value of the difference between year1
      and year2 must be > 1.
    """
    cdef int s
    cdef int t
    s, _, t, _, _ = year_index_pchange(match)
    if abs(s - t) > 1:
        return True
    return False

cdef update_contract(double kappa, match, S):
    """
    Update the list that maintains all communities with contraction in size.

    INPUT:

    - kappa -- contraction threshold.
    - match -- a community match between a front and a step.  The match
      follows this format:

          year1 comm_index1 year2 comm_index2 pchange

      The key "year1" is the snapshot year for a front community and
      "comm_index1" is a community index.  Similarly, "year2" is the snapshot
      year for a step community and "comm_index2" is a community index.  The
      key "pchange" is the percentage change in community size.
    - S -- list of contracting communities.  This list will be updated if the
      given match results in a community with decreased size.
    """
    cdef double p
    p = float(match.split()[-1])
    if (p < 0.0) and (abs(p) > kappa):
        S.append(match)

cdef update_events(int y, F, FR, C, Match, Birth, Split, Merge, Direct):
    """
    Determine the community birth, split, and merge for the given
    snapshot year.  Don't consider death here.

    INPUT:

    - y -- current snapshot year.
    - F -- list of fronts with corresponding matching step communities.
    - FR -- dictionary of communities with resuming observations.
    - C -- list of steps with corresponding matching front communities.
    - Match -- dictionary mapping matches to corresponding percentage changes
      in community sizes.
    - Birth -- this will store all the new communities in this snapshot year.
    - Split -- this will store all community splits in this snapshot year.
    - Merge -- this will store all community merges in this snapshot year.
    - Direct -- this will store all communities with direct continuation.
      Given a front F_i and a step C_a, if F_i matches only C_a and vice
      versa, then we say that there is a direct continuation from F_i to C_a.
    """
    cdef int a
    cdef int cindex
    cdef int i
    cdef int index
    cdef int year
    # iterate over each front community
    for i, f in enumerate(F):
        # NOTE: Don't consider community death here.  Do it after we have
        # constructed a graph of dynamic communities.
        #
        # F_i matches exactly C_a, which in turn matches exactly F_i
        # direct continuation
        if len(f) == 1:
            a = f[0]
            if len(C[a]) == 1:
                match = "%d %d %d %d" % (y - 1, i, y, a)
                pchange = Match[match]
                Direct.append("%s %s" % (match, pchange))
        # F_i matches multiple step communities: community split
        elif len(f) > 1:
            for a in f:
                match = "%d %d %d %d" % (y - 1, i, y, a)
                pchange = Match[match]
                Split.append("%s %s" % (match, pchange))
    # iterate over each community with resuming observations
    for year, index in FR:
        # NOTE: Don't consider community death here.
        #
        # M_i matches exactly C_a, which in turn matches exactly M_i
        # direct continuation
        if len(FR[(year, index)]) == 1:
            a = FR[(year, index)][0]
            if len(C[a]) == 1:
                match = "%d %d %d %d" % (year, index, y, a)
                pchange = Match[match]
                Direct.append("%s %s" % (match, pchange))
        # M_i matches multiple step communities: community split
        elif len(FR[(year, index)]) > 1:
            for a in FR[(year, index)]:
                match = "%d %d %d %d" % (year, index, y, a)
                pchange = Match[match]
                Split.append("%s %s" % (match, pchange))
    # iterate over each step community
    for a, c in enumerate(C):
        # step C_a does not match any front communities: community birth
        if len(c) == 0:
            Birth.append("%d %d" % (y, a))
        # NOTE: We don't consider the case where C_a matches exactly F_i and
        # vice versa.  This has already been considered in the loop for front
        # communities.
        #
        # C_a matches multiple front/missing communities: community merge
        elif len(c) > 1:
            for year, cindex in c:
                match = "%d %d %d %d" % (year, cindex, y, a)
                pchange = Match[match]
                Merge.append("%s %s" % (match, pchange))

cdef update_expand(double gamma, match, S):
    """
    Update the list that maintains all communities with expansion in size.

    INPUT:

    - gamma -- growth threshold.
    - match -- a community match between a front and a step.  The match
      follows this format:

          year1 comm_index1 year2 comm_index2 pchange

      The key "year1" is the snapshot year for a front community and
      "comm_index1" is a community index.  Similarly, "year2" is the snapshot
      year for a step community and "comm_index2" is a community index.  The
      key "pchange" is the percentage change in community size.
    - S -- list of expanding communities.  This list will be updated if the
      given match results in a community with increased size.
    """
    cdef double p
    p = float(match.split()[-1])
    if p > gamma:
        S.append(match)

cdef update_missing(match, S):
    """
    Update the list that maintains all communities with missing observations.

    INPUT:

    - match -- a community match between a front and a step.  The match
      follows this format:

          year1 comm_index1 year2 comm_index2 pchange

      The key "year1" is the snapshot year for a front community and
      "comm_index1" is a community index.  Similarly, "year2" is the snapshot
      year for a step community and "comm_index2" is a community index.  The
      key "pchange" is the percentage change in community size.  For a missing
      observation, the match follows the format:

          year1 comm_index1 year1 comm_index1 0.0
    - S -- list of communities with missing observations.
    """
    cdef int Cs
    cdef int Ct
    cdef int s
    cdef int t
    s, Cs, t, Ct, _ = year_index_pchange(match)
    assert(s == t)
    assert(Cs == Ct)
    S.append(match)

cdef update_resuming(match, S):
    """
    Update the list that maintains all communities with resuming observations.

    INPUT:

    - match -- a community match between a front and a step.  The match
      follows this format:

          year1 comm_index1 year2 comm_index2 pchange

      The key "year1" is the snapshot year for a front community and
      "comm_index1" is a community index.  Similarly, "year2" is the snapshot
      year for a step community and "comm_index2" is a community index.  The
      key "pchange" is the percentage change in community size.  For a resuming
      observation, the absolute value of the difference between year1 and
      year2 must be > 1.
    - S -- list of communities with resuming observations.
    """
    cdef int s
    cdef int t
    s, _, t, _, _ = year_index_pchange(match)
    assert(abs(s - t) > 1)
    S.append(match)

cdef update_stable(double gamma, double kappa, match, S):
    """
    Update the list that maintains all communities with stable community sizes.
    A stable community cannot be a stagnant community.  A community is stable
    when it has changed in size.  If growth has occurred, the percentage
    size increase does not exceed the growth threshold gamma.  Where
    contraction has occurred, the percentage size decrease does not exceed
    the contraction threshold kappa.

    INPUT:

    - gamma -- growth threshold.
    - kappa -- contraction threshold.
    - match -- a community match between a front and a step.  The match
      follows this format:

          year1 comm_index1 year2 comm_index2 pchange

      The key "year1" is the snapshot year for a front community and
      "comm_index1" is a community index.  Similarly, "year2" is the snapshot
      year for a step community and "comm_index2" is a community index.  The
      key "pchange" is the percentage change in community size.
    - S -- list of stable communities.  This list will be updated if the
      given match results in a community with a relatively stable size.
    """
    cdef double p
    p = float(match.split()[-1])
    if p == 0.0:
        return
    if p > gamma:
        return
    if (p < 0.0) and (abs(p) > kappa):
        return
    assert((0.0 < p <= gamma) or (0.0 < abs(p) <= kappa))
    S.append(match)

cdef update_stagnant(match, S):
    """
    Update the list that maintains all communities with no increase or
    decrease in size.

    INPUT:

    - match -- a community match between a front and a step.  The match
      follows this format:

          year1 comm_index1 year2 comm_index2 pchange

      The key "year1" is the snapshot year for a front community and
      "comm_index1" is a community index.  Similarly, "year2" is the snapshot
      year for a step community and "comm_index2" is a community index.  The
      key "pchange" is the percentage change in community size.  If pchange is
      zero, then we say that the corresponding match results in a stagnant
      community.
    - S -- list of stable communities.  This list will be updated if the given
      match results in stable community size.
    """
    cdef double p
    p = float(match.split()[-1])
    if p == 0.0:
        S.append(match)

cdef write_death(D, dirname):
    """
    Output to a file data on community deaths.

    INPUT:

    - D -- dictionary of dynamic communities that are dead.  The key is a
      year.  The value is a list of all fronts of dynamic communities that die
      during that year.
    - dirname -- path to directory with match data.
    """
    cdef int year
    if not os.path.isdir(os.path.join(dirname, "death")):
        os.makedirs(os.path.join(dirname, "death"))
    # write all dead communities for each year to a file
    # very first line is the number of dead communities in that year
    for year in D:
        s = "%d\n" % len(D[year])
        for front in sorted(D[year]):
            s += "%s\n" % front
        fname = os.path.join(dirname, "death", "death-%d.dat" % year)
        f = open(fname, "w")
        f.write(s)
        f.close()

cdef write_etype_comm(fname, D):
    """
    Output to the given file data on a particular life-cycle event.

    INPUT:

    - fname -- file to which we write.
    - D -- dictionary of a particular life-cycle event.  The key is a
      dynamic community.  The value is a list of all dynamic communities that
      relates to the given dynamic community through the event. Here are some
      examples.

      * Contraction: The key is a contracted dynamic community.  The value is
        a list of all dynamic communities that contract into that given
        dynamic community.
      * Expansion: The key is an expanded dynamic community.  The value is a
        list of all dynamic communities that expand into that given dynamic
        community.
      * Merge: The key is a merged dynamic community.  The value is a list of
        all dynamic communities that merged in that given dynamic community.
      * Resume: The key is a resume dynamic community.  The value is a list of
        all dynamic communities that have resume observation.
      * Split: The key is a split dynamic community.  The value is a list of
        all dynamic communities that split in that given dynamic community.
      * Stable: The key is a dynamic community with stable size.  The value is
        a list of all dynamic communities whose sizes change into the dynamic
        community referenced by the key.
      * Stagnant: The key is a dynamic community with stagnant size.  The value
        is a list of all dynamic communities whose sizes change into the
        dynamic community referenced by the key.
    """
    # very first line is the count
    f = open(fname, "w")
    f.write("%d\n" % len(D))
    for commB in sorted(D):
        contract = ",".join(sorted(D[commB]))
        f.write("%s<-%s\n" % (commB, contract))
    f.close()

cdef write_event(int year, dirname, etype, E):
    """
    Write to a file all the community indices satisfying a particular
    community key event.  The very first line of the output file is an
    integer counting the number of communities satisfying the required event.

    INPUT:

    - year -- the snapshot year.
    - dirname -- name of directory under which to write file.
    - etype -- the type of community event.  Supported key events are:
      birth, death, expand, contract, split, merge, stable, stagnant, direct,
      missing, resume.
    - E -- list of all community indices whose corresponding communities
      satisfy the key event under consideration.
    """
    events = ("birth", "death", "expand", "contract", "split", "merge",
              "stable", "stagnant", "direct", "missing", "resume")
    if etype not in events:
        raise ValueError("Unsupported key community event: %s" % etype)
    try:
        os.makedirs(os.path.join(dirname, etype))
    except OSError:
        # directory already exists
        pass
    fname = os.path.join(dirname, etype, "%s-%s.dat" % (etype, str(year)))
    f = open(fname, "w")
    f.write("%s\n" % str(len(E)))
    for e in E:
        f.write("%s\n" % str(e))
    f.close()

cdef year_index_pchange(match):
    """
    Extract the years, community indices, and percentage change from the
    given match.  The match follows this format:

        year1 comm_index1 year2 comm_index2 pchange

    The key "year1" is the snapshot year for a front community and
    "comm_index1" is a community index.  Similarly, "year2" is the snapshot
    year for a step community and "comm_index2" is a community index.  The
    key "pchange" is the percentage change in community size.

    INPUT:

    - match -- a community match between a front and a step.
    """
    cdef int Cs
    cdef int Ct
    cdef int s
    cdef int t
    M = match.split()
    s = int(M[0])   # year of front community
    Cs = int(M[1])  # index of front community
    t = int(M[2])   # year of step community
    Ct = int(M[3])  # index of step community
    p = M[4]
    return (s, Cs, t, Ct, p)
