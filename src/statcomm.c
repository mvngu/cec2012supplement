/**************************************************************************
 * Copyright (C) 2011 Minh Van Nguyen <nguyenminh2@gmail.com>
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

/* Statistics on communities. */

#define _GNU_SOURCE  /* asprintf */

#include <argtable2.h>
#include <ctype.h>
#include <igraph.h>
#include <stdio.h>

void avg_community_statistics(const igraph_vector_ptr_t *C,
                              const igraph_t *G,
                              const igraph_integer_t min_commsize,
                              igraph_integer_t *ncomm,
                              igraph_real_t *avg_comm_plength,
                              igraph_real_t *avg_comm_coeff,
                              igraph_real_t *avg_comm_assort,
                              igraph_real_t *avg_comm_size);
void get_years(const char *infile,
               igraph_vector_t *year);
void read_communities(FILE *f,
                      igraph_vector_ptr_t *C);
void vector_of_vectors_init(igraph_vector_ptr_t *V,
                            const igraph_integer_t size);
void write_header(const char *fname);
void write_results(const char *fname,
                   const igraph_integer_t time,
                   const igraph_integer_t ncomm,
                   const igraph_real_t avg_comm_plength,
                   const igraph_real_t avg_comm_coeff,
                   const igraph_real_t avg_comm_assort,
                   const igraph_real_t avg_comm_size);

/* Various average community statistics.  For each community, we find various
 * statistics on it.  Then we average all such statistics.  The following
 * statistics are considered:
 *
 * - average community characteristic path length
 * - average community clustering coefficient
 * - average community degree assortativity
 * - average community size
 *
 * - C -- a vector of vectors; each element of the vector is itself a
 *   vector.  Each such vector stores the node IDs belonging to one community.
 * - G -- the snapshot graph.
 * - min_commsize -- the minimum size of each community.  Must be at least 1.
 *   We will only consider all communities with the given minimum number of
 *   nodes.
 * - ncomm -- the total number of communities with the given minimum size.
 *   The result will be stored here.
 * - avg_comm_plength -- the average community characteristic path length
 *   will be stored here.
 * - avg_comm_coeff -- the average community clustering coefficient will be
 *   stored here.
 * - avg_comm_assort -- the average community degree assortativity will be
 *   stored here.
 * - avg_commsize -- the average community size will be stored here.
 */
void avg_community_statistics(const igraph_vector_ptr_t *C,
                              const igraph_t *G,
                              const igraph_integer_t min_commsize,
                              igraph_integer_t *ncomm,
                              igraph_real_t *avg_comm_plength,
                              igraph_real_t *avg_comm_coeff,
                              igraph_real_t *avg_comm_assort,
                              igraph_real_t *avg_comm_size) {
  igraph_t g;  /* subgraph of G; graph representing community C[i] */
  igraph_integer_t i;
  igraph_integer_t commsize;  /* the size of a community */
  igraph_vs_t vs;
  igraph_real_t plength;      /* characteristic path length of a community */
  igraph_real_t coeff;        /* clustering coefficient of a community */
  igraph_real_t assort;       /* community degree assortativity */
  igraph_real_t sum_plength;  /* cumulative sum of char. path length */
  igraph_real_t sum_coeff;    /* cumulative sum of clustering coefficient */
  igraph_real_t sum_assort;   /* cumulative sum of degree assortativity */
  igraph_integer_t sum_size;  /* cumulative sum of community size */
  igraph_integer_t nplength;  /* # char. path length values we have summed */
  igraph_integer_t ncoeff;    /* # clustering coeff. values we have summed */
  igraph_integer_t nassort;   /* # assortativity values we have summed */
  igraph_integer_t nsize;     /* # community size values we have summed */

  sum_plength = 0.0;
  sum_coeff = 0.0;
  sum_assort = 0.0;
  sum_size = 0;
  nplength = 0;
  ncoeff = 0;
  nassort = 0;
  nsize = 0;
  for (i = 0; i < igraph_vector_ptr_size(C); i++) {
    /* neglect community with size < min_commsize */
    commsize = igraph_vector_size((igraph_vector_t *)VECTOR(*C)[i]);
    if (commsize < min_commsize) {
      continue;
    }

    igraph_vs_vector(&vs, (igraph_vector_t *)VECTOR(*C)[i]);
    igraph_induced_subgraph(G, &g, vs, IGRAPH_SUBGRAPH_AUTO);

    /* Average path length, otherwise known as the characteristic path */
    /* length.  In case the subgraph is unconnected, we use the average of */
    /* the geodesics of all connected components. */
    igraph_average_path_length(&g, &plength, IGRAPH_UNDIRECTED, /*unconn*/1);
    sum_plength += plength;
    nplength++;

    /* Average clustering coefficient following Watts and Strogatz.  For */
    /* a vertex with degree < 2, we set its average clustering coefficient */
    /* to zero. */
    igraph_transitivity_avglocal_undirected(&g, &coeff,
                                            IGRAPH_TRANSITIVITY_ZERO);
    sum_coeff += coeff;
    ncoeff++;

    /* average community degree assortativity */
    igraph_assortativity_degree(&g, &assort, 1);
    /* only consider assortativity value whenever this is defined */
    if (!isnan(assort)) {
      sum_assort += assort;
      nassort++;
    }

    /* average community size */
    sum_size += commsize;
    nsize++;

    /* clean up */
    igraph_vs_destroy(&vs);
    igraph_destroy(&g);
  }

  *avg_comm_plength = sum_plength / (igraph_real_t)nplength;
  *avg_comm_coeff = sum_coeff / (igraph_real_t)ncoeff;
  *avg_comm_assort = sum_assort / (igraph_real_t)nassort;
  *avg_comm_size = sum_size / (igraph_real_t)nsize;
  *ncomm = nsize;
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
  fprintf(file, "year,#comm,deg assort,avg clustering coeff,avg path length,avg comm size\n");
  fclose(file);
}

/* Write results to a file.  Results are written according to the following
 * CSV file format:
 *
 * year,assortativity,clustering coefficient, characteristic path length,community size
 *
 * - fname -- the name of the file to which we write results.
 * - time -- the current time step.
 * - assortativity -- average community degree assortativity.
 * - clustcoeff -- average community clustering coefficient.
 * - charpathlen -- average community characteristic path length.
 * - avg_commsize -- average community size.
 */
void write_results(const char *fname,
                   const igraph_integer_t time,
                   const igraph_integer_t ncomm,
                   const igraph_real_t avg_comm_plength,
                   const igraph_real_t avg_comm_coeff,
                   const igraph_real_t avg_comm_assort,
                   const igraph_real_t avg_comm_size) {
  FILE *file;

  file = fopen(fname, "a");
  fprintf(file, "%li,%li,%lf,%lf,%lf,%lf\n",
          (long int)time, (long int)ncomm, (double)avg_comm_assort,
          (double)avg_comm_coeff, (double)avg_comm_plength,
          (double)avg_comm_size);
  fclose(file);
}

/* Evolution of various statistics of a network.  Here we consider statistics
 * on communities, not the whole network.
 */
int main(int argc,
         char **argv) {
  FILE *infile;                    /* file to read from */
  char *infname;
  char *outfname;                  /* name of file to write to */
  igraph_vector_t year;            /* all the snapshot years */
  igraph_t G;                      /* a graph */
  igraph_integer_t i;              /* generic index */
  igraph_integer_t min_commsize;   /* minimum #nodes for each community */
  igraph_integer_t ncomm;          /* # communities of given minimum size */
  igraph_real_t avg_comm_assort;   /* avg community degree assortativity */
  igraph_real_t avg_comm_coeff;    /* avg community clustering coeff */
  igraph_real_t avg_comm_plength;  /* avg community path length */
  igraph_real_t avg_comm_size;     /* avg community size */
  igraph_vector_ptr_t C;           /* vector of communities */

  /* setup the table of command line options */
  struct arg_lit *help = arg_lit0(NULL, "help", "print this help and exit");
  struct arg_str *minsize = arg_strn(NULL, "minsize", "integer", 0,
                                     argc + 2, "minimum community size");
  struct arg_file *fyear = arg_file0(NULL, "year", NULL,
                           "file containing all snapshot years, one per line");
  struct arg_str *edgelistdir = arg_strn(NULL, "edgelistdir", "path", 0,
                                    argc + 2,
                                    "directory from which to read edge lists");
  struct arg_str *comdir = arg_strn(NULL, "comdir", "path", 0, argc + 2,
                               "directory with files of community structures");
  struct arg_str *statdir = arg_strn(NULL, "statdir", "path", 0, argc + 2,
                                     "directory to which we write statistics");
  struct arg_end *end = arg_end(20);
  void *argtable[] = {help, minsize, fyear, edgelistdir, comdir, statdir, end};

  /* parse the command line as defined by argtable[] */
  arg_parse(argc, argv, argtable);

  /* print usage information when --help is passed in */
  if (help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Statistics on communities for a network snapshot.\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }
  /* special case: all command line arguments, except --help, must be used */
  if ((minsize->count < 1) || (fyear->count < 1) || (edgelistdir->count < 1)
      || (comdir->count < 1) || (statdir->count < 1)) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Statistics on communities for a network snapshot.\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }

  min_commsize = (igraph_integer_t)atoi(minsize->sval[0]);
  asprintf(&outfname, "%s/statistics-%ld.csv", statdir->sval[0],
           (long int)min_commsize);
  write_header(outfname);
  igraph_vector_init(&year, 0);
  get_years(fyear->filename[0], &year);

  /* iterate over each snapshot and get various statistics on it */
  for (i = 0; i < igraph_vector_size(&year); i++) {
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

    /* get statistics on this snapshot graph */
    igraph_vector_ptr_init(&C, 0);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&C, igraph_vector_destroy);
    asprintf(&infname, "%s/comm-%li.dat", comdir->sval[0],
             (long int)VECTOR(year)[i]);
    infile = fopen(infname, "r");
    read_communities(infile, &C);
    fclose(infile);
    free(infname);
    avg_community_statistics(&C, &G, min_commsize, &ncomm, &avg_comm_plength,
                             &avg_comm_coeff, &avg_comm_assort,
                             &avg_comm_size);

    /* output results to a file */
    write_results(outfname, (igraph_integer_t)VECTOR(year)[i], ncomm,
                  avg_comm_plength, avg_comm_coeff, avg_comm_assort,
                  avg_comm_size);

    /* clean up */
    igraph_destroy(&G);
    igraph_vector_ptr_destroy_all(&C);
  }

  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
  free(outfname);
  igraph_vector_destroy(&year);

  return 0;
}
