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
# Community matches with direct/indirect continuation.  Missing observations
# are allowed, depending on the command line argument.

# third party library
from dycom import match
# built-in library
import argparse
import os

###################################
# script starts here
###################################

# setup parser for command line options
s = "Community matches with direct/indirect continuation.\n"
parser = argparse.ArgumentParser(description=s)
parser.add_argument("--size", metavar="integer", required=True, type=int,
                    help="minimum community size")
parser.add_argument("--theta", metavar="float", required=True, type=float,
                    help="matching threshold")
parser.add_argument("--delta", metavar="integer", required=True, type=int,
                    help="maximum allowable consecutive missing observations")
parser.add_argument("--year", metavar="file", required=True,
                    help="file containing all snapshot years, one per line")
parser.add_argument("--eventdir", metavar="path", required=True,
                    help="directory to read/write event data")
args = parser.parse_args()

# get command line arguments & sanity checks
size = args.size
theta = args.theta
delta = args.delta
yearfname = args.year
eventdir = args.eventdir
if size < 1:
    raise ValueError("Invalid minimum size")
if theta < 0.0:
    raise ValueError("Invalid theta")
if delta < 0:
    raise ValueError("Invalid delta")

dirname = os.path.join(eventdir, "size-%d" % size,
                       "theta-%s_delta-%d" % (str(theta), delta))
if not os.path.isdir(dirname):
    os.makedirs(dirname)
comdir = os.path.join(eventdir, "size-%d" % size)
match(yearfname, dirname, comdir, theta, delta)
