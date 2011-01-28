/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * main.c
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
  idxtype *part, *sizes;
  GraphType graph;
  float ipc2redist, *xyz, *tpwgts, ubvec[MAXNCON];
  MPI_Comm comm;
  int numflag=0, wgtflag=0, ndims=3, edgecut;

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (argc != 8) {
    if (mype == 0)
      printf("Usage: %s <graph-file> <op-type> <nparts> <adapth-factor> <ipc2redist> <dbglvl> <seed>\n", argv[0]);

    MPI_Finalize();
    exit(0);
  }

  optype = atoi(argv[2]);
  nparts = atoi(argv[3]);
  adptf = atoi(argv[4]);
  ipc2redist = atof(argv[5]);

  options[0] = 1;
  options[PMV3_OPTION_DBGLVL] = atoi(argv[6]);
  options[PMV3_OPTION_SEED] = atoi(argv[7]);

  if (mype == 0) printf("reading file: %s\n", argv[1]);
  ParallelReadGraph(&graph, argv[1], comm);

  /* Remove the edges for testing */
  /*idxset(graph.vtxdist[mype+1]-graph.vtxdist[mype]+1, 0, graph.xadj); */

  sset(graph.ncon, 1.05, ubvec);
  tpwgts = fmalloc(nparts, "tpwgts");
  sset(nparts*graph.ncon, 1.0/(float)nparts, tpwgts);

  /*
  ChangeToFortranNumbering(graph.vtxdist, graph.xadj, graph.adjncy, mype, npes); 
  numflag = 1;

  nvtxs = graph.vtxdist[mype+1]-graph.vtxdist[mype];
  nedges = graph.xadj[nvtxs];
  printf("%d %d\n", idxsum(nvtxs, graph.xadj), idxsum(nedges, graph.adjncy));
  */


  xyz = NULL;
  if (optype >= 20) 
    xyz = ReadTestCoordinates(&graph, argv[1], ndims, comm);
  if (mype == 0) printf("finished reading file: %s\n", argv[1]);
  
  part = idxsmalloc(graph.nvtxs, mype%nparts, "main: part");
  sizes = idxmalloc(2*npes, "main: sizes");

  switch (optype) {
    case 1: 
      ParMETIS_V3_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag, 
               &numflag, &graph.ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
      WritePVector(argv[1], graph.vtxdist, part, MPI_COMM_WORLD); 
      break;
    case 2:
      options[PMV3_OPTION_PSR] = COUPLED;
      ParMETIS_V3_RefineKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, 
               &wgtflag, &numflag, &graph.ncon, &nparts, tpwgts, ubvec, options, 
	       &edgecut, part, &comm);
      break;
    case 3:
      options[PMV3_OPTION_PSR] = COUPLED;
      graph.vwgt = idxsmalloc(graph.nvtxs, 1, "main: vwgt");
      if (npes > 1) {
        AdaptGraph(&graph, adptf, comm);
      }
      else {
        wgtflag = 0;
        ParMETIS_V3_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag, 
                 &numflag, &graph.ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);

        printf("Initial partitioning with edgecut of %d\n", edgecut);
        for (i=0; i<graph.ncon; i++) {
          for (j=0; j<graph.nvtxs; j++) {
            if (part[j] == i)
              graph.vwgt[j*graph.ncon+i] = adptf; 
            else
              graph.vwgt[j*graph.ncon+i] = 1; 
          }
        }
      }
/*
wgtflag = 0;
ParMETIS_V3_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag, 
&numflag, &graph.ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
if (mype == 0) printf("Initial partitioning with edgecut of %d\n", edgecut);
*/
      wgtflag = 3;
      ParMETIS_V3_AdaptiveRepart(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, NULL, 
               graph.adjwgt, &wgtflag, &numflag, &graph.ncon, &nparts, tpwgts, ubvec, 
	       &ipc2redist, options, &edgecut, part, &comm);
      break;
    case 4: 
      ParMETIS_V3_NodeND(graph.vtxdist, graph.xadj, graph.adjncy, &numflag, options, part, sizes, &comm);
      /* WriteOVector(argv[1], graph.vtxdist, part, comm);   */
      MALLOC_CHECK(NULL);
      break;
    case 5: 
      ParMETIS_SerialNodeND(graph.vtxdist, graph.xadj, graph.adjncy, &numflag, options, part, sizes, &comm);
      /* WriteOVector(argv[1], graph.vtxdist, part, comm);  */ 
      printf("%d %d %d %d %d %d %d\n", sizes[0], sizes[1], sizes[2], sizes[3], sizes[4], sizes[5], sizes[6]);
      break;
    case 11: 
      /* TestAdaptiveMETIS(graph.vtxdist, graph.xadj, graph.adjncy, part, options, adptf, comm); */
      break;
    case 20: 
      ParMETIS_V3_PartGeomKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag, 
        &numflag, &ndims, xyz, &graph.ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
      break;
    case 21: 
      ParMETIS_V3_PartGeom(graph.vtxdist, &ndims, xyz, part, &comm);
      break;
    case 22: 
      /* ParMETIS_PartGeomRefine(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, &wgtflag, &numflag, &ndims, xyz, options, &edgecut, part, &comm); */
      break;
  }

  /* printf("%d %d\n", idxsum(nvtxs, graph.xadj), idxsum(nedges, graph.adjncy)); */

  GKfree(&part, &sizes, &graph.vtxdist, &graph.xadj, &graph.adjncy, &graph.vwgt, &xyz, LTERM);

  MPI_Comm_free(&comm);

  MPI_Finalize();

  return 0;
}


/*************************************************************************
* This function changes the numbering to be from 1 instead of 0
**************************************************************************/
void ChangeToFortranNumbering(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, int mype, int npes)
{
  int i, nvtxs, nedges;

  nvtxs = vtxdist[mype+1]-vtxdist[mype];
  nedges = xadj[nvtxs];

  for (i=0; i<npes+1; i++)
    vtxdist[i]++;
  for (i=0; i<nvtxs+1; i++)
    xadj[i]++;
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  return;
}
