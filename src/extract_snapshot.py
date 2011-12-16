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

# NOTE: This script should be run from Python.
#
# This is used to determine snapshots of an evolving network.  We consider
# yearly snapshots.  Work out all the authors in a yearly dataset.  We use
# names as provided and consider each unique name as a unique author. Then we
# sort the list of unique author names in alphabetical order.  This is
# followed by a mapping of each author name to a unique nonnegative integer.

# third-party library
from dycom import build_snapshots
# standard library
import argparse

####################
# script starts here
####################

# setup parser for command line options
s = "Extract yearly snapshots from bibliographic database."
parser = argparse.ArgumentParser(description=s)
parser.add_argument("--year", metavar="file", required=True,
                    help="file containing all snapshot years, one per line")
parser.add_argument("--datadir", metavar="dir", required=True,
                    help="directory with data files, one for each year")
parser.add_argument("--edgelistdir", metavar="dir", required=True,
                    help="directory to output edge list for each snapshot")
args = parser.parse_args()

# get the command line arguments
yearfname = args.year
datadir = args.datadir
edgelistdir = args.edgelistdir

build_snapshots(yearfname, datadir, edgelistdir)
