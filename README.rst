Genetic Programming Bibliography
================================

Downloaded on 25th August 2011 from

http://www.cs.bham.ac.uk/~wbl/biblio/

by Minh Van Nguyen <mvngu.name@gmail.com>.  The bibliography is maintained
and managed by William Langdon, Steven Gustafson, and John Koza.

We consider two versions of this dataset:

* First version: consider both authors and editors.  See the directory
  ``editors/`` for this version.

* Second version: only consider authors and exclude editors.  See the
  directory ``noeditors/`` for this other version of the GP dataset.


Datasets
--------

Here's an explanation of the dataset files we used:

* ``gp-bib.txt`` -- The original dataset downloaded from the above URL.

* ``gp-bib.ref`` -- The same as gp-bib.txt, but the bibliography is
  structured according to ref format.

* ``gp-bib.dat`` -- The same as ``gp-bib.ref``, but containing only
  fields and values of interest for our investigation.  We parse the
  file ``gp-bib.dat`` using ``src/noyear.py`` and
  ``src/extract_dataset.py``.  The script ``src/noyear.py`` looks for
  any publications without a corresponding publication year.  The
  script ``src/extract_dataset.py`` extracts the relevant publication
  metadata and writes the result in a CSV file.  Two versions are
  written to files.  The first version is ``editors/gp-bib.csv``,
  which contains both editors and authors.  The other version is
  ``noeditors/gp-bib.csv``, which excludes editors and retains
  authors.


Source code
-----------

The directory ``src/`` contains source code of programs to preprocess
the GP dataset.  It also has programs to track dynamic communities
based on a community life-cycle.


License
-------

All source code are licensed under the terms of the GNU General Public
License, either version 2 or (at your option) any later version.  For
more information, see http://www.gnu.org/licenses/
