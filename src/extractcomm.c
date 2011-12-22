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

/* Extract community structures from network snapshots.  This program
 * supports some algorithms for community detection.
 */

#define _GNU_SOURCE  /* asprintf */

#include <argtable2.h>
#include <ctype.h>
#include <igraph.h>
#include <stdio.h>
#include <string.h>

void extract_communities_blondel(igraph_vector_ptr_t *C,
                                 const igraph_t *G);
void extract_communities_clauset(igraph_vector_ptr_t *C,
                                 const igraph_t *G);
void extract_communities_newman(igraph_vector_ptr_t *C,
                                const igraph_t *G);
void extract_communities_raghavan(igraph_vector_ptr_t *C,
                                  const igraph_t *G);
void get_years(const char *infile,
               igraph_vector_t *year);
void vector_of_vectors_init(igraph_vector_ptr_t *V,
                            igraph_integer_t size);

/* Get all step communities from a given snapshot graph.  Get community
 * structure via multilevel optimization of modularity as described in
 * (Blondel et al. 2008).
 *
 * - C -- pointer to an initialized object.  The step communities will be
 *   extracted from the graph G and stored here.
 * - G -- a snapshot graph from which to extract communities.
 *
 * References
 *
 * (Blondel et al. 2008) V. D. Blondel, J.-L. Guillaume, R. Lambiotte, and
 *     E. Lefebvre. Fast unfolding of communities in large networks. Journal
 *     of Statistical Mechanics, 2008(October):P10008, 2008.
 */
void extract_communities_blondel(igraph_vector_ptr_t *C,
                                 const igraph_t *G) {
  igraph_vector_t membership;
  igraph_integer_t a;      /* community index */
  igraph_integer_t i;      /* generic index */
  igraph_integer_t ncomm;  /* number of communities */
  igraph_vector_t *p;

  igraph_vector_init(&membership, 0);
  igraph_community_multilevel(G, /*weights*/ NULL, &membership,
                              /*memberships*/ NULL, /*modularity*/ NULL);

  /* Find out how many communities we are dealing with.  Allocate enough */
  /* space to track that many communities.  The index of C is therefore */
  /* the index of some community.  Note that indices for communities */
  /* start from zero.  If there is only one community, then the index */
  /* of that community is 0.  Hence we add one to the largest community */
  /* index. */
  ncomm = (igraph_integer_t)igraph_vector_max(&membership) + 1;
  vector_of_vectors_init(C, ncomm);
  /* Place all vertices belonging to one community into a vector. */
  for (i = 0; i < igraph_vector_size(&membership); i++) {
    a = VECTOR(membership)[i];
    p = (igraph_vector_t *)VECTOR(*C)[a];
    igraph_vector_push_back(p, i);
  }

  igraph_vector_destroy(&membership);
}

/* Get all step communities from a given snapshot graph.  Get community
 * structure via greedy optimization of modularity.  We use the method of
 * (Clauset et al. 2004) with some optimization by (Wakita & Tsurumi 2007).
 * Then map the community structure to the membership vector.  Use the
 * membership vector to determine collection of nodes belonging to each
 * community.  Such collections of nodes are treated as step communities for
 * time t.
 *
 * - C -- pointer to an initialized object.  The step communities will be
 *   extracted from the graph G and stored here.
 * - G -- a snapshot graph from which to extract communities.
 *
 * References
 *
 * (Clauset et al. 2004) Aaron Clauset, M. E. J. Newman, and Cristopher Moore.
 *     Finding community structure in very large networks. Physical Review E,
 *     70(6):066111, 2004.
 * (Wakita & Tsurumi 2007) Ken Wakita and Toshiyuki Tsurumi. Finding community
 *     structure in mega-scale social networks (extended abstract). In
 *     Carey L. Williamson, Mary Ellen Zurko, Peter F. Patel-Schneider, and
 *     Prashant J. Shenoy, eds, WWW, pp.1275--1276. ACM, 2007.
 */
void extract_communities_clauset(igraph_vector_ptr_t *C,
                                 const igraph_t *G) {
  igraph_matrix_t merges;
  igraph_vector_t membership;
  igraph_vector_t modularity;
  igraph_integer_t a;      /* community index */
  igraph_integer_t i;      /* generic index */
  igraph_integer_t ncomm;  /* number of communities */
  igraph_vector_t *p;

  igraph_vector_init(&membership, 0);
  igraph_vector_init(&modularity, 0);
  igraph_matrix_init(&merges, 0, 0);
  igraph_community_fastgreedy(G, 0, &merges, &modularity, &membership);

  /* Find out how many communities we are dealing with.  Allocate enough */
  /* space to track that many communities.  The index of C is therefore */
  /* the index of some community.  Note that indices for communities */
  /* start from zero.  If there is only one community, then the index */
  /* of that community is 0.  Hence we add one to the largest community */
  /* index. */
  ncomm = (igraph_integer_t)igraph_vector_max(&membership) + 1;
  vector_of_vectors_init(C, ncomm);
  /* place all vertices belonging to one community into a vector */
  for (i = 0; i < igraph_vector_size(&membership); i++) {
    a = VECTOR(membership)[i];
    p = (igraph_vector_t *)VECTOR(*C)[a];
    igraph_vector_push_back(p, i);
  }

  igraph_matrix_destroy(&merges);
  igraph_vector_destroy(&membership);
  igraph_vector_destroy(&modularity);
}

/* Get all step communities from a given snapshot graph.  Get community
 * structure via Newman's method of leading eigenvector.  Then map the
 * community structure to the membership vector.  Use the membership vector to
 * determine collection of nodes belonging to each community.   Such
 * collections of nodes are treated as step communities for time t.
 *
 * - C -- pointer to an initialized object.  The step communities will be
 *   extracted from the graph G and stored here.
 * - G -- a snapshot graph from which to extract communities.  The graph is
 *   assumed to be simple and undirected.
 *
 * Reference
 *
 * (Newman 2006) M. E. J. Newman. Finding community structure in networks
 *     using the eigenvectors of matrices. Physical Review E, 74(3):036104,
 *     2006.
 */
void extract_communities_newman(igraph_vector_ptr_t *C,
                                const igraph_t *G) {
  igraph_matrix_t merges;
  igraph_vector_t membership;
  igraph_integer_t a;      /* community index */
  igraph_integer_t i;      /* generic index */
  igraph_integer_t nstep;
  igraph_integer_t ncomm;  /* number of communities */
  igraph_vector_t *p;
  igraph_arpack_options_t options;
  int ret;

  igraph_set_error_handler(igraph_error_handler_ignore);

  /* extract community structures */
  igraph_vector_init(&membership, 0);
  igraph_matrix_init(&merges, 0, 0);
  /* While the function returns an error, invoke it with a smaller number */
  /* of steps.  First, try with the maximum number of steps. */
  nstep = igraph_vcount(G) + 1;
  ret = 1;
  while (ret) {
    nstep--;
    igraph_vector_destroy(&membership);
    igraph_matrix_destroy(&merges);
    igraph_vector_init(&membership, 0);
    igraph_matrix_init(&merges, 0, 0);
    igraph_arpack_options_init(&options);
    ret = igraph_community_leading_eigenvector(G, &merges, &membership,
              /*n steps*/ nstep, &options, /*modularity*/ NULL, /*start*/ 0,
              /*eigenvalues*/ NULL, /*eigenvectors*/ NULL, /*history*/ NULL,
              /*callback*/ NULL, /*callback extrac*/ NULL);
  }

  /* Find out how many communities we are dealing with.  Allocate enough */
  /* space to track that many communities.  The index of C is therefore */
  /* the index of some community.  Note that indices for communities */
  /* start from zero.  If there is only one community, then the index */
  /* of that community is 0.  Hence we add one to the largest community */
  /* index. */
  ncomm = (igraph_integer_t)igraph_vector_max(&membership) + 1;
  vector_of_vectors_init(C, ncomm);
  /* Place all vertices belonging to one community into a vector. */
  for (i = 0; i < igraph_vector_size(&membership); i++) {
    a = VECTOR(membership)[i];
    p = (igraph_vector_t *)VECTOR(*C)[a];
    igraph_vector_push_back(p, i);
  }

  /* clean up */
  igraph_matrix_destroy(&merges);
  igraph_vector_destroy(&membership);
}

/* Get all step communities from a given snapshot graph.  Get community
 * structure based on label propagation as described in (Raghavan et al. 2007).
 *
 * - C -- pointer to an initialized object.  The step communities will be
 *   extracted from the graph G and stored here.
 * - G -- a snapshot graph from which to extract communities.
 *
 * References
 *
 * (Raghavan et al. 2007) U. N. Raghavan, R. Albert, and S. Kumara. Near
 *     linear time algorithm to detect community structures in large-scale
 *     networks. Physical Review E, 76(3):036106, 2007.
 */
void extract_communities_raghavan(igraph_vector_ptr_t *C,
                                  const igraph_t *G) {
  igraph_vector_t membership;
  igraph_integer_t a;      /* community index */
  igraph_integer_t i;      /* generic index */
  igraph_integer_t ncomm;  /* number of communities */
  igraph_vector_t *p;

  /* singleton graph: the vertex is its own community */
  if (igraph_vcount(G) == 1) {
    ncomm = 1;
    vector_of_vectors_init(C, ncomm);
    p = (igraph_vector_t *)VECTOR(*C)[0];
    igraph_vector_push_back(p, 0);
    return;
  }

  igraph_vector_init(&membership, 0);
  igraph_community_label_propagation(G, &membership, /*weights*/ NULL,
                                     /*initial states*/ NULL, /*fixed*/ NULL,
                                    /*modularity*/ NULL);

  /* Find out how many communities we are dealing with.  Allocate enough */
  /* space to track that many communities.  The index of C is therefore */
  /* the index of some community.  Note that indices for communities */
  /* start from zero.  If there is only one community, then the index */
  /* of that community is 0.  Hence we add one to the largest community */
  /* index. */
  ncomm = (igraph_integer_t)igraph_vector_max(&membership) + 1;
  vector_of_vectors_init(C, ncomm);
  /* Place all vertices belonging to one community into a vector. */
  for (i = 0; i < igraph_vector_size(&membership); i++) {
    a = VECTOR(membership)[i];
    /* community index == -1, so we place vertex i in its own community */
    if (a < 0) {
      p = igraph_Calloc(1, igraph_vector_t);
      igraph_vector_init(p, 0);
      igraph_vector_push_back(p, i);
      igraph_vector_ptr_push_back(C, p);
    } else {
      p = (igraph_vector_t *)VECTOR(*C)[a];
      igraph_vector_push_back(p, i);
    }
  }

  igraph_vector_destroy(&membership);
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

/* Read in a graph and output all the communities of that graph to a file.
 * The format of communities is as follows:
 *
 * v_0 v_1 ... v_k
 *
 * where each v_j refers to a node in a community.  All v_j on a line refer
 * to all nodes in a community.  The very first line of the file is the number
 * of communities.
 */
int main(int argc,
         char **argv) {
  FILE *file;
  char *fname;
  igraph_vector_t year;   /* all the snapshot years */
  igraph_integer_t i;     /* generic index */
  igraph_integer_t j;     /* generic index */
  igraph_integer_t k;     /* generic index */
  igraph_t G;             /* network snapshot graph */
  igraph_vector_ptr_t C;  /* vector of communities */
  igraph_vector_t *v;

  /* setup the table of command line options */
  struct arg_lit *help = arg_lit0(NULL, "help", "print this help and exit");
  struct arg_file *fyear = arg_file0(NULL, "year", NULL,
                           "file containing all snapshot years, one per line");
  struct arg_str *edgelistdir = arg_strn(NULL, "edgelistdir", "path", 0,
                                    argc + 2,
                                    "directory from which to read edge lists");
  struct arg_str *algo = arg_strn(NULL, "algo", "string", 0, argc + 2,
          "community detection algorithm: blondel, clauset, newman, raghavan");
  struct arg_str *comdir = arg_strn(NULL, "comdir", "path", 0, argc + 2,
                                 "directory under which to write communities");
  struct arg_end *end = arg_end(20);
  void *argtable[] = {help, fyear, edgelistdir, algo, comdir, end};

  /* parse the command line as defined by argtable[] */
  arg_parse(argc, argv, argtable);

  /* print usage information when --help is passed in */
  if (help->count > 0) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Extract communities from network snapshots.\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }
  /* special case: all command line arguments, except --help, must be used */
  if ((fyear->count < 1) || (edgelistdir->count < 1) || (algo->count < 1)
      || (comdir->count < 1)) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Extract communities from network snapshots.\n");
    arg_print_glossary(stdout, argtable, "  %-20s %s\n");
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }
  /* special case: must use a supported algorithm for community detection */
  if ((strcmp(algo->sval[0], "blondel") != 0)
      && (strcmp(algo->sval[0], "clauset") != 0)
      && (strcmp(algo->sval[0], "newman") != 0)
      && (strcmp(algo->sval[0], "raghavan") != 0)) {
    printf("Usage: %s", argv[0]);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Extract communities from network snapshots.\n");
    arg_print_glossary(stdout, argtable, "  %-20s %s\n");
    arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
    exit(EXIT_FAILURE);
  }

  igraph_vector_init(&year, 0);
  get_years(fyear->filename[0], &year);

  /* Construct a graph G from its edge list.  Extract all communities from */
  /* G and output those communities to a file. */
  for (i = 0; i < igraph_vector_size(&year); i++) {
    printf("%li\n", (long int)VECTOR(year)[i]);
    fflush(stdout);
    /* reconstruct graph from its edge list */
    asprintf(&fname, "%s/edgelist-%li.dat", edgelistdir->sval[0],
             (long int)VECTOR(year)[i]);
    file = fopen(fname, "r");
    igraph_read_graph_edgelist(&G, file, /*nvert*/ 0, IGRAPH_UNDIRECTED);
    fclose(file);
    free(fname);
    igraph_simplify(&G, /*no multiple edges*/ 1, /*no self-loops*/ 1,
                    /*edge_comb*/ 0);

    /* extract communities and output those communities to a file */
    igraph_vector_ptr_init(&C, 0);
    IGRAPH_VECTOR_PTR_SET_ITEM_DESTRUCTOR(&C, igraph_vector_destroy);
    if (strcmp(algo->sval[0], "blondel") == 0) {
      extract_communities_blondel(&C, &G);
    } else if (strcmp(algo->sval[0], "clauset") == 0) {
      extract_communities_clauset(&C, &G);
    } else if (strcmp(algo->sval[0], "newman") == 0) {
      extract_communities_newman(&C, &G);
    } else if (strcmp(algo->sval[0], "raghavan") == 0) {
      extract_communities_raghavan(&C, &G);
    }
    asprintf(&fname, "%s/comm-%li.dat", comdir->sval[0],
             (long int)VECTOR(year)[i]);
    file = fopen(fname, "w");
    /* write the number of communities as the very first line */
    fprintf(file, "%li\n", (long int)igraph_vector_ptr_size(&C));
    /* Write the nodes of each community on one line.  Don't include */
    /* trailing white spaces. */
    for (j = 0; j < igraph_vector_ptr_size(&C); j++) {
      v = (igraph_vector_t *)VECTOR(C)[j];
      for (k = 0; k < (igraph_vector_size(v) - 1); k++) {
        fprintf(file, "%li ", (long int)VECTOR(*v)[k]);
      }
      k = igraph_vector_size(v) - 1;
      fprintf(file, "%li\n", (long int)VECTOR(*v)[k]);
    }
    fclose(file);
    free(fname);

    /* clean ups */
    igraph_destroy(&G);
    igraph_vector_ptr_destroy_all(&C);
  }

  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));
  igraph_vector_destroy(&year);

  return 0;
}
