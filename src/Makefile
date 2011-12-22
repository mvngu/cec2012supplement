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

FILE =
MELB = /home/mvngu/usr/bin
SAGE = /scratch/mvngu/usr/bin

all:
	gcc $(FILE).c -Wall -W -I/home/mvngu/usr/include/igraph -L/home/mvngu/usr/lib -ligraph -I/home/mvngu/usr/include -largtable2 -o $(FILE)

cython-melb:
	$(MELB)/python setup.py build_ext --inplace

cython-sage:
	$(SAGE)/python setup.py build_ext --inplace

geom:
	gcc $(FILE).c -Wall -W -I/scratch/mvngu/usr/include/igraph -L/scratch/mvngu/usr/lib -ligraph -o $(FILE)

mungerabah:
	gcc $(FILE).c -Wall -W -I/home/pgrad/minguyen/usr/include/igraph -L/home/pgrad/minguyen/usr/lib -ligraph -o $(FILE)

sage:
	gcc $(FILE).c -Wall -W -I/scratch/mvngu/usr/include/igraph -L/scratch/mvngu/usr/lib -ligraph -I/scratch/mvngu/usr/include -largtable2 -o $(FILE)

clean:
	rm -rf *~
	rm -rf *.so
	rm -rf build/
	rm -rf dycom.c
	rm -rf $(FILE)
	rm -rf valgrind.log

valgrind:
	valgrind --verbose --leak-check=full --show-reachable=yes --leak-resolution=high --log-file=valgrind.log --num-callers=25 --track-origins=yes ./$(FILE)
