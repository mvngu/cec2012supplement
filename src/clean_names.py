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
# WARNING: This script modifies the file you feed it.
#
# Sanitize names in bibliographic datasets.

# third-party library
from dycom import clean_names
# standard library
import argparse

####################
# script starts here
####################

# setup parser for command line options
s = "Sanitize names.\n"
s += "WARNING: This script modifies the file you feed it."
parser = argparse.ArgumentParser(description=s)
parser.add_argument("file", help="the file you want to clean up names")
args = parser.parse_args()

infile = args.file
clean_names(infile)
