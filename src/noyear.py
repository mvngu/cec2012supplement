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

# NOTE: This script should be run from Sage, not from Python.

# Find any publications without possibly corresponding years.

fname = "gp-bib.dat"
f = open(fname, "r")
hastitle = False
hasyear = False
nline = 0
for line in f:
    ln = line.strip()
    nline += 1
    # nonempty line
    if ln:
        if ln.startswith("%T "):
            hastitle = True
            continue
        if ln.startswith("%D "):
            hasyear = True
            continue
    # empty line
    else:
        if hastitle and hasyear:
            hastitle = False
            hasyear = False
            continue
        else:
            hastitle = False
            hasyear = False
            print(nline)
            continue
f.close()
