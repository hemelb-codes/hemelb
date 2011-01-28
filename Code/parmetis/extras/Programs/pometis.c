/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * pometis.c
 *
 * This is the entry point of the ILUT
 *
 * Started 10/19/95
 * George
 *
 * $Id: parmetis.c,v 1.5 2003/07/30 21:18:54 karypis Exp $
 *
 */

#include <parmetisbin.h>

/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  int i, j, npes, mype, optype, nparts, adptf, options[10];
  idxtype *order, *sizes;
  GraphType graph;
  MPI_Comm comm;
  int numflag=0, wgtflag=0, ndims=3, edgecut;
  int mtype, rtype, nseps, seed, dbglvl;
  float ubfrac;

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (argc != 5) {
    if (mype == 0) {
      printf("Usage: %s <graph-file> <op-type> <seed> <dbglvl>\n", argv[0]);
      printf("  op-type: 1=ParNodeND_V3, 2=SerNodeND\n");
    }

    MPI_Finalize();
    exit(0);
  }

  optype = atoi(argv[2]);

  if (mype == 0) 
    printf("reading file: %s\n", argv[1]);
  ParallelReadGraph(&graph, argv[1], comm);

  order = idxsmalloc(graph.nvtxs, mype%nparts, "main: order");
  sizes = idxmalloc(2*npes, "main: sizes");

  switch (optype) {
    case 1: 
      options[0] = 1;
      options[PMV3_OPTION_SEED]   = atoi(argv[3]);
      options[PMV3_OPTION_DBGLVL] = atoi(argv[4]);
      ParMETIS_V3_NodeND(graph.vtxdist, graph.xadj, graph.adjncy, &numflag, 
          options, order, sizes, &comm);
      break;
    case 2: 
      ParMETIS_SerialNodeND(graph.vtxdist, graph.xadj, graph.adjncy, &numflag, 
          options, order, sizes, &comm);
      break;
    default:
      if (mype == 0) 
        printf("Uknown optype of %d\n", optype);
      MPI_Finalize();
      exit(0);
  }
  MALLOC_CHECK(NULL);

  WriteOVector(argv[1], graph.vtxdist, order, comm);

  /* print the partition sizes and the separators */
  if (mype == 0) {
    for (i=0; i<2*npes-1; i++)
      printf("%6d ", sizes[i]);
    printf("\n");
  }


  GKfree(&order, &sizes, &graph.vtxdist, &graph.xadj, &graph.adjncy, &graph.vwgt, LTERM);

  MPI_Comm_free(&comm);

  MPI_Finalize();

  return 0;
}

