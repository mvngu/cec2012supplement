/**************************************************************************
 * Copyright (C) 2011 Minh Van Nguyen <mvngu.name@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * http://www.gnu.org/licenses/
 *************************************************************************/

/* Basic network statistics.  These are the summary statistics for a
 * network.
 */

#define _GNU_SOURCE  /* asprintf */

#include <argtable2.h>
#include <assert.h>
#include <ctype.h>
#include <igraph.h>
#include <stdio.h>

void avg_component_size(const igraph_vector_t *csize,
                        igraph_real_t *res);
void get_years(const char *infile,
               igraph_vector_t *year);
void read_communities(FILE *f,
                      igraph_vector_ptr_t *C);
void size_largest_comp(const igraph_vector_t *csize,
                       const igraph_integer_t nvert,
                       igraph_real_t *res);
void vector_of_vectors_init(igraph_vector_ptr_t *V,
                            const igraph_integer_t size);
void write_header(const char *fname);
void write_results(const char *fname,
                   const igraph_integer_t year,
                   const igraph_real_t assortativity,
                   const igraph_real_t clustcoeff,
                   const igraph_real_t charpathlen,
                   const igraph_real_t avg_compsize,
                   const igraph_integer_t ncomp,
                   const igraph_real_t bigcomp_size);

/* The average component size.  The size of a connected component is the
 * number of vertices in it.  Let V be a vector of component sizes, where V[i]
 * is the size of component i.  Let G be the graph under consideration, having
 * n vertices.  Define S = sum_{i=0}^k V[i] to be the sum of all component
 * sizes.  Then the average component size is A = S / |V|.
 *
 * - csize -- vector of component sizes.
 * - res -- the result will be stored here.
 */
void avg_component_size(const igraph_vector_t *csize,
                        igraph_real_t *res) {
  *res = (igraph_real_t)igraph_vector_sum(csize);
  *res = *res / (igraph_real_t)igraph_vector_size(csize);
}

/* Get all the snapshot years from the given file.  Sort the list of snapshot
 * years in nondecreasing order.
 *
 * - infile -- file containing all the snapshot years.
 * - year -- vector of snapshot years.  The result will be stored here.  This
 *   should be initialized beforehand.
 */
void get_years(const char *infile,
               igraph_vector_t *year) {
  FILE *f;
  long int y;  /* a snapshot year */
  int c;       /* a character */

  f = fopen(infile, "r");

  /* skip all white spaces */
  do {
    c = getc(f);
  } while (isspace(c));
  ungetc(c, f);

  /* read all snapshot years, one per line */
  while (!feof(f)) {
    /* get all the nodes belonging to a community */
    fscanf(f, "%li", &y);
    igraph_vector_push_back(year, (igraph_integer_t)y);
    /* skip all white spaces */
    do {
      c = getc(f);
    } while (isspace(c));
    ungetc(c, f);
  }

  fclose(f);
  igraph_vector_sort(year);
}

/* Read in the collection of communities for a network snapshot.
 *
 * - f -- a stream from which to read the communities.
 * - C -- an initialized vector of vectors.  The communities for a snapshot
 *   will be stored here.  Each community will be stored as a vector of
 *   integers, where each integer denotes a vertex ID.
 */
void read_communities(FILE *f,
                      igraph_vector_ptr_t *C) {
  long int i;      /* generic index */
  long int ncomm;  /* # of communities */
  long int node;   /* a node belonging to a community */
  int c;
  igraph_vector_t *p;

  /* skip all white spaces */
  do {
    c = getc(f);
  } while (isspace(c));
  ungetc(c, f);
  /* The very first line of the given file contains only the number of */
  /* communities.  Read in this number and initialize the community vector */
  /* C to have that many vectors. */
  fscanf(f, "%li", &ncomm);
  vector_of_vectors_init(C, (igraph_integer_t)ncomm);

  /* skip all white spaces */
  do {
    c = getc(f);
  } while (isspace(c));
  ungetc(c, f);

  /* Index i is now the community index.  We start from community with */
  /* index 0.  Each newline character read indicates that we are to */
  /* increment the community index by one. */
  i = 0;
  while (!feof(f)) {
    /* get all the nodes belonging to a community */
    fscanf(f, "%li", &node);
    p = (igraph_vector_t *)VECTOR(*C)[i];
    igraph_vector_push_back(p, (igraph_integer_t)node);
    /* skip all white spaces */
    do {
      c = getc(f);
      /* All the nodes belonging to a community are assumed to be on one */
      /* line.  If we encounter a newline character, then we know that */
      /* we have read all the nodes of a community. */
      if (c == '\n') {
        i++;
      }
    } while (isspace(c));
    ungetc(c, f);
  }
}

/* Relative size of the largest component.
 *
 * - csize -- vector of component sizes.
 * - nvert -- the number of vertices.  The size of the largest component will
 *   be normalized using this integer.
 * - res -- the result will be stored here.
 */
void size_largest_comp(const igraph_vector_t *csize,
                       const igraph_integer_t nvert,
                       igraph_real_t *res) {
  igraph_integer_t i;  /* generic index */
  igraph_integer_t m;  /* index of component with largest size */

  m = 0;
  for (i = 0; i < igraph_vector_size(csize); i++) {
    if (VECTOR(*csize)[m] < VECTOR(*csize)[i]) {
      m = i;
    }
  }
  *res = (igraph_real_t)VECTOR(*csize)[m] / (igraph_real_t)nvert;
}

/* Initialize a vector of vectors.  This is just a vector, each of whose
 * elements is a pointer to a vector.
 *
 * - V -- pointer to an initialized vector of vectors.  The result will be
 *   stored here.
 * - size -- the number of elements in V.
 */
void vector_of_vectors_init(igraph_vector_ptr_t *V,
                            const igraph_integer_t size) {
  igraph_integer_t i;
  igraph_vector_t *p;

  for (i = 0; i < size; i++) {
    p = igraph_Calloc(1, igraph_vector_t);
    igraph_vector_init(p, 0);
    igraph_vector_ptr_push_back(V, p);
  }
}

/* Write header for the results file.
 *
 * - fname -- the name of the file to which we write the header.
 */
void write_header(const char *fname) {
  FILE *file;

  file = fopen(fname, "w");
  fprintf(file, "year,deg assort,avg clustering coeff,avg path length,avg comp size,#component,rel size big comp\n");
  fclose(file);
}

/* Write results to a file.  Results are written according to the following
 * CSV file format:
 *
 * time,assortativity,clustering coefficient, characteristic path length
 *
 * - fname -- the name of the file to which we write results.
 * - year -- the current time step, usually the year.
 * - assortativity -- the degree assortativity.
 * - clustcoeff -- the average clustering coefficient.
 * - charpathlen -- the average path length, otherwise known as the
 *   characteristic path length.
 * - avg_compsize -- the average component size.
 * - ncomp -- the number of components.
 * - bigcomp_size -- the relative size of the largest component.
 */
void write_results(const char *fname,
                   const igraph_integer_t year,
                   const igraph_real_t assortativity,
                   const igraph_real_t clustcoeff,
                   const igraph_real_t charpathlen,
                   const igraph_real_t avg_compsize,
                   const igraph_integer_t ncomp,
                   const igraph_real_t bigcomp_size) {
  FILE *file;

  file = fopen(fname, "a");
  fprintf(file, "%li,%lf,%lf,%lf,%lf,%li,%lf\n",
          (long int)year, (double)assortativity, (double)clustcoeff,
          (double)charpathlen, (double)avg_compsize, (long int)ncomp,
          (double)bigcomp_size);
  fclose(file);
}

/* Evolution of various statistics of a network.  Currently we consider the
 * following statistics:
 *
 * - degree assortativity (Newman)
 * - average clustering coefficient (Watts-Strogatz)
 * - average path length
 * - average component size
 * - number of components
 * - relative size of the largest component
 */
int main(int argc,
         char **argv) {
  FILE *infile;                 /* file to read from */
  char *infname;
  igraph_vector_t year;         /* all the snapshot years */
  igraph_t G;                   /* a graph */
  igraph_integer_t i;           /* generic index */
  igraph_integer_t ncomp;       /* number of components */
  igraph_real_t a;              /* degree assortativity */
  igraph_real_t c;              /* average clustering coefficient */
  igraph_real_t p;              /* average path length */
  igraph_real_t avg_comp_size;  /* average component size */
  igraph_real_t bigcomp_size;   /* size of largest component */
  igraph_vector_t csize;        /* vector of component sizes */

  /* setup the table of command line options */
  struct arg_lit *help = arg_lit0(NULL, "help", "print this help and exit");
  struct arg_file *fyear = arg_file0(NULL, "year", NULL,
                           "file containing all snapshot years, one per line");
  struct arg_str *edgelistdir = arg_strn(NULL, "edgelistdir", "path", 0,
                                    argc + 2,
                                    "directory from which to read edge lists");
  struct arg_file *fstats = arg_file0(NULL, "outfile", NULL,
                                      "file to which we write statistics");
  struct arg_end *end = arg_end(20);
  void *argtable[] = {help, fyear, edgelistdir, fstats, end};

  /* parse the command line as defined by argtable[] */
  arg_parse(argc, argv, argtable);

  /* print usage information when --help is passed in */
  if (help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Summary statistics for network snapshots.\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }
  /* special case: all command line arguments, except --help, must be used */
  if ((fyear->count < 1) || (edgelistdir->count < 1) || (fstats->count < 1)) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Summary statistics for network snapshots.\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    /* deallocate each non-null entry in argtable[] */
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }

  write_header(fstats->filename[0]);
  igraph_vector_init(&year, 0);
  get_years(fyear->filename[0], &year);

  /* iterate over each snapshot and get various statistics on it */
  for (i = 0; i < igraph_vector_size(&year); i++) {
    printf("%li\n", (long int)VECTOR(year)[i]);
    fflush(stdout);
    /* Construct the snapshot graph from its edge list.  The graph is */
    /* a simple graph, i.e. no self-loops nor multiple edges. */
    asprintf(&infname, "%s/edgelist-%li.dat", edgelistdir->sval[0],
             (long int)VECTOR(year)[i]);
    infile = fopen(infname, "r");
    igraph_read_graph_edgelist(&G, infile, 0, IGRAPH_UNDIRECTED);
    fclose(infile);
    free(infname);
    igraph_simplify(&G, /*no multiple edges*/ 1, /*no self-loops*/ 1,
                    /*edge_comb*/ 0);

    /* average degree assortativity */
    igraph_assortativity_degree(&G, &a, 1);
    /* Average clustering coefficient following Watts and Strogatz.  For */
    /* a vertex with degree < 2, we set its average clustering coefficient */
    /* to zero. */
    igraph_transitivity_avglocal_undirected(&G, &c, IGRAPH_TRANSITIVITY_ZERO);
    /* Average path length, otherwise known as the characteristic path */
    /* length.  In case the graph is unconnected, we use the average of the */
    /* geodesics of all connected components. */
    igraph_average_path_length(&G, &p, IGRAPH_UNDIRECTED, /*unconn*/ 1);

    /* average component size; number of components; size of largest comp */
    igraph_vector_init(&csize, 0);
    igraph_clusters(&G, /*membership*/ NULL, &csize, &ncomp, IGRAPH_WEAK);
    assert(ncomp == igraph_vector_size(&csize));
    avg_component_size(&csize, &avg_comp_size);
    size_largest_comp(&csize, igraph_vcount(&G), &bigcomp_size);

    /* output results to a file */
    write_results(fstats->filename[0], VECTOR(year)[i], a, c, p,
                  avg_comp_size, ncomp, bigcomp_size);

    /* clean up */
    igraph_destroy(&G);
    igraph_vector_destroy(&csize);
  }

  igraph_vector_destroy(&year);
  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

  return EXIT_SUCCESS;
}
