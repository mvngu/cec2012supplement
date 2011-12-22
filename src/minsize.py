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

# NOTE: This script should be run from Python.
#
# Communities with a specified minimum size.

# third-party library
from dycom import comm_minsize
# standard library
import argparse
import os

###################################
# script starts here
###################################

# setup parser for command line options
s = "Extract communities with a specified minimum size.\n"
parser = argparse.ArgumentParser(description=s)
parser.add_argument("--size", metavar="integer", required=True, type=int,
                    help="minimum community size")
parser.add_argument("--year", metavar="file", required=True,
                    help="file containing all snapshot years, one per line")
parser.add_argument("--comdir", metavar="path", required=True,
                    help="directory with communities per snapshot year")
parser.add_argument("--outdir", metavar="path", required=True,
                    help="directory to write all communities with minsize")
args = parser.parse_args()

# get the command line arguments
size = args.size
yearfname = args.year
comdir = args.comdir
outdir = args.outdir
# sanity check
if size < 1:
    raise ValueError("Invalid minimum size")

dirname = os.path.join(outdir, "size-%d" % size)
if not os.path.isdir(dirname):
    os.makedirs(dirname)
comm_minsize(yearfname, size, comdir, dirname)
