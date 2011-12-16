###########################################################################
# Copyright (C) 2011 Minh Van Nguyen <nguyenminh2@gmail.com>
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

# NOTE: This script could be run from Python or Sage.

import os

# This is used to extract publication metadata from the genetic programming
# dataset.  The metadata we are interested in for each publication include
# information about the year of publication, the list of authors, and where
# relevant the list of editors.  The publication metadata is represented in
# the following CSV format
#
# yyyy,Author1,Author2,...[,editor,Editor1,Editor2,...]
#
# where the items in square brackets are optional.  For an edited volume, it
# is possible for the list of authors to be blank and we only have a list of
# editors.  In that case, we use ",," to denote empty authors list.  For
# example:
#
# yyyy,,editor,Editor1,Editor2,...

####################
# helper functions
####################

def comma_separated(arr):
    """
    Concatenate each element of the given array into one string.  The
    resulting string has comma separated values.  Essentially, a comma
    should precede each value.

    INPUT:

    - arr -- an array of values to concatenate into a CSV string.
    """
    A = map(lambda x: "," + x, arr)
    return "".join(A)

def extract_names(s, name):
    """
    Extract a list of author/editor names from the given string.  The
    resulting names are appended to the given array.

    INPUT:

    - s -- a string containing author/editor names.  This string should be in
      the refer format for bibliography.  Thus it should start with "%A " to
      indicate an authors list, or start with "%E " to indicate an editors
      list.
    - name -- an array of author/editor names.  The results will be appended
      to this array.
    """
    line = ""
    if s.startswith("%A "):
        line = s.split("%A ")[-1].strip()
    elif s.startswith("%E "):
        line = s.split("%E ")[-1].strip()
    else:
        raise ValueError
    a = []
    # multiple authors/editors
    if " and " in line:
        a = line.split(" and ")
    # single author/editor
    else:
        a = [line]
    a = map(lambda x: x.strip(), a)
    for e in a:
        name.append(e)

def extract_year(s, year):
    """
    Extract the publication year from the given string.

    INPUT:

    - s -- a string containing the publication year of some ...er...
      publication.  This string should be in the refer format for
      bibliography.  Thus it should start with "%D " to indicate the
      year of publication.
    - year -- an array of publication years.  Each publication should ideally
      have exactly one publication year.  The extracted year will be stored
      in this array.
    """
    if not s.startswith("%D "):
        raise ValueError
    yr = s.split("%D ")[-1].strip()
    year.append(yr)

def write_metadata(author, editor, year, fname):
    """
    Write to a file the metadata for a publication.  The metadata we are
    interested in for each publication include information about the year of
    publication, the list of authors, and where relevant the list of editors.
    The publication metadata is represented in the following CSV format

        yyyy,Author1,Author2,...[,editor,Editor1,Editor2,...]

    where the items in square brackets are optional.  For an edited volume, it
    is possible for the list of authors to be blank and we only have a list of
    editors.  Each publication must have an associated, and only one,
    publication year; an error is raised otherwise.  Furthermore, a
    publication must have at least one author.  If the publication has no
    authors, then it must have at least one editor.  An error is raised if the
    publication has no authors nor editors.

    INPUT:

    - author -- array of author names; can be an empty array.
    - editor -- array of editor names; can be empty.
    - year -- an array of publication years; cannot be empty.
    - fname -- the name of the file to which we write the metadata.
    """
    # no publication year or multiple publication years
    if len(year) != 1:
        raise ValueError("Multiple publication years.")
    # neither authors nor editors provided
    if (len(author) < 1) and (len(editor) < 1):
        raise ValueError("No author nor editor names provided.")
    # now write stuff to file
    # first, include both authors and editors
    s = year[0]
    if len(author) > 0:
        s += comma_separated(author)
    if len(editor) > 0:
        s += comma_separated(editor)
    f = open(os.path.join("..", "editors", fname), "a")
    f.write(s + "\n")
    f.close()
    # next, consider the author names but exclude the editors
    if len(author) > 0:
        s = year[0]
        s += comma_separated(author)
        f = open(os.path.join("..", "noeditors", fname), "a")
        f.write(s + "\n")
        f.close()

####################
# script starts here
####################

# Extract metadata from bibliography maintained in the refer format.  Use
# the bibliography that was last updated on 2011/08/24.
PREFIX = ".."
infname = "gp-bib.dat"
outfname = "gp-bib.csv"
author = []  # the authors of a publication
editor = []  # the editors of a publication
year = []    # the publication year
# the case of including both authors and editors
f = open(os.path.join(PREFIX, "editors", outfname), "w")
f.close()
# now only author names; exclude the editors
f = open(os.path.join("noeditors", outfname), "w")
f.close()

f = open(os.path.join(PREFIX, infname), "r")
for line in f:
    ln = line.strip()
    # Nonempty line.  This indicates that we are processing the first line
    # of a publication.  In general, it a nonempty string here means that
    # we are reading a data line corresponding to some publication.
    if ln:
        if ln.startswith("%A "):
            extract_names(ln, author)
        elif ln.startswith("%E "):
            extract_names(ln, editor)
        elif ln.startswith("%D "):
            extract_year(ln, year)
    # Empty line.  We've finished processing data for a publication.
    else:
        write_metadata(author, editor, year, outfname)
        author = []
        editor = []
        year = []
f.close()
