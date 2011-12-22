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
# Track key events in the life-cycle of dynamic communities.  This script
# assumes that we already have data on community matches.

# third-party library
from dycom import events_except_death
from dycom import events_only_death
# standard library
import argparse
import os

###################################
# script starts here
###################################

# setup parser for command line options
s = "Track events in life-cycle of dynamic communities.\n"
parser = argparse.ArgumentParser(description=s)
parser.add_argument("--deathonly", metavar="boolean", required=True,
                    help="[True|False] process only data on death")
parser.add_argument("--size", metavar="integer", required=True, type=int,
                    help="minimum community size")
parser.add_argument("--theta", metavar="float", required=True, type=float,
                    help="matching threshold")
parser.add_argument("--delta", metavar="integer", required=True, type=int,
                    help="maximum allowable consecutive missing observations")
parser.add_argument("--gamma", metavar="float", required=True, type=float,
                    help="growth threshold")
parser.add_argument("--kappa", metavar="float", required=True, type=float,
                    help="contraction threshold")
parser.add_argument("--year", metavar="file", required=True,
                    help="file containing all snapshot years, one per line")
parser.add_argument("--eventdir", metavar="path", required=True,
                    help="directory to read/write event data")
args = parser.parse_args()

# get command line arguments & sanity checks
deathonly = args.deathonly
size = args.size
theta = args.theta
delta = args.delta
gamma = args.gamma
kappa = args.kappa
yearfname = args.year
eventdir = args.eventdir
if size < 1:
    raise ValueError("Invalid minimum size")
if theta < 0.0:
    raise ValueError("Invalid theta")
if delta < 0:
    raise ValueError("Invalid delta")
if gamma < 0.0:
    raise ValueError("Invalid gamma")
if kappa < 0.0:
    raise ValueError("Invalid kappa")

dirname = os.path.join(eventdir, "size-%d" % size,
                       "theta-%s_delta-%d" % (str(theta), delta))
if not os.path.isdir(dirname):
    raise ValueError("No such directory: %s" % dirname)
comdir = os.path.join(eventdir, "size-%d" % size)

if deathonly in ("True", "true"):
    events_only_death(yearfname, dirname, delta)
elif deathonly in ("False", "false"):
    events_except_death(yearfname, dirname, comdir, gamma, kappa)
else:
    raise ValueError("Invalid value for boolean deathonly")
