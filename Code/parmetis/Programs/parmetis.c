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
  idxtype *part=NULL, *sizes=NULL;
  GraphType graph;
  float ipc2redist, *xyz=NULL, *tpwgts=NULL, ubvec[MAXNCON];
  MPI_Comm comm;
  int numflag=0, wgtflag=0, ndims, edgecut;

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

  optype     = atoi(argv[2]);
  nparts     = atoi(argv[3]);
  adptf      = atoi(argv[4]);
  ipc2redist = atof(argv[5]);

  options[0] = 1;
  options[PMV3_OPTION_DBGLVL] = atoi(argv[6]);
  options[PMV3_OPTION_SEED]   = atoi(argv[7]);

  if (mype == 0) 
    printf("reading file: %s\n", argv[1]);
  ParallelReadGraph(&graph, argv[1], comm);

  /* Remove the edges for testing */
  /*idxset(graph.vtxdist[mype+1]-graph.vtxdist[mype]+1, 0, graph.xadj); */

  sset(graph.ncon, 1.05, ubvec);
  tpwgts = fmalloc(nparts*graph.ncon, "tpwgts");
  sset(nparts*graph.ncon, 1.0/(float)nparts, tpwgts);

  /*
  ChangeToFortranNumbering(graph.vtxdist, graph.xadj, graph.adjncy, mype, npes); 
  numflag = 1;

  nvtxs = graph.vtxdist[mype+1]-graph.vtxdist[mype];
  nedges = graph.xadj[nvtxs];
  printf("%d %d\n", idxsum(nvtxs, graph.xadj), idxsum(nedges, graph.adjncy));
  */


  if (optype >= 20)  
    xyz = ReadTestCoordinates(&graph, argv[1], &ndims, comm);

  if (mype == 0) 
    printf("finished reading file: %s\n", argv[1]);
  
  part  = idxsmalloc(graph.nvtxs, mype%nparts, "main: part");
  sizes = idxmalloc(2*npes, "main: sizes");

  switch (optype) {
    case 1: 
      wgtflag = 3;
      ParMETIS_V3_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
          graph.adjwgt, &wgtflag, &numflag, &graph.ncon, &nparts, tpwgts, ubvec, 
          options, &edgecut, part, &comm);
      WritePVector(argv[1], graph.vtxdist, part, MPI_COMM_WORLD); 
      break;
    case 2:
      wgtflag = 3;
      options[PMV3_OPTION_PSR] = PARMETIS_PSR_COUPLED;
      ParMETIS_V3_RefineKway(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
          graph.adjwgt, &wgtflag, &numflag, &graph.ncon, &nparts, tpwgts, ubvec, 
          options, &edgecut, part, &comm);
      break;
    case 3:
      options[PMV3_OPTION_PSR] = PARMETIS_PSR_COUPLED;
      graph.vwgt = idxsmalloc(graph.nvtxs, 1, "main: vwgt");
      if (npes > 1) {
        AdaptGraph(&graph, adptf, comm);
      }
      else {
        wgtflag = 3;
        ParMETIS_V3_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
            graph.adjwgt, &wgtflag, &numflag, &graph.ncon, &nparts, tpwgts, 
            ubvec, options, &edgecut, part, &comm);

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

      wgtflag = 3;
      ParMETIS_V3_AdaptiveRepart(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
          NULL, graph.adjwgt, &wgtflag, &numflag, &graph.ncon, &nparts, tpwgts, ubvec, 
	  &ipc2redist, options, &edgecut, part, &comm);
      break;
    case 4: 
      ParMETIS_V3_NodeND(graph.vtxdist, graph.xadj, graph.adjncy, &numflag, options, 
          part, sizes, &comm);
      /* WriteOVector(argv[1], graph.vtxdist, part, comm);   */
      break;

    case 5: 
      ParMETIS_SerialNodeND(graph.vtxdist, graph.xadj, graph.adjncy, &numflag, options, 
          part, sizes, &comm);
      /* WriteOVector(argv[1], graph.vtxdist, part, comm);  */ 
      printf("%" IDXTYPE_FORMAT " %" IDXTYPE_FORMAT " %" IDXTYPE_FORMAT " %" IDXTYPE_FORMAT " %" IDXTYPE_FORMAT " %" IDXTYPE_FORMAT " %" IDXTYPE_FORMAT "\n", sizes[0], sizes[1], sizes[2], sizes[3], sizes[4], sizes[5], sizes[6]);
      break;
    case 11: 
      /* TestAdaptiveMETIS(graph.vtxdist, graph.xadj, graph.adjncy, part, options, adptf, comm); */
      break;
    case 20: 
      wgtflag = 3;
      ParMETIS_V3_PartGeomKway(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
          graph.adjwgt, &wgtflag, &numflag, &ndims, xyz, &graph.ncon, &nparts, 
          tpwgts, ubvec, options, &edgecut, part, &comm);
      break;
    case 21: 
      ParMETIS_V3_PartGeom(graph.vtxdist, &ndims, xyz, part, &comm);
      break;
  }

  /* printf("%d %d\n", idxsum(nvtxs, graph.xadj), idxsum(nedges, graph.adjncy)); */

  GKfree((void **)&part, &sizes, &tpwgts, &graph.vtxdist, &graph.xadj, &graph.adjncy, 
         &graph.vwgt, &graph.adjwgt, &xyz, LTERM);

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
