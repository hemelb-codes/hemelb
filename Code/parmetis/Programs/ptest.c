/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * main.c
 * 
 * This file contains code for testing teh adaptive partitioning routines
 *
 * Started 5/19/97
 * George
 *
 * $Id: ptest.c,v 1.3 2003/07/22 21:47:20 karypis Exp $
 *
 */

#include <parmetisbin.h>

#define NCON    5

/*************************************************************************/
/*! Entry point of the testing routine */
/*************************************************************************/
int main(int argc, char *argv[])
{
  int mype, npes;
  MPI_Comm comm;

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (argc != 2 && argc != 3) {
    if (mype == 0)
      printf("Usage: %s <graph-file> [coord-file]\n", argv[0]);

    MPI_Finalize();
    exit(0);
  }

  TestParMetis_GPart(argv[1], (argc == 3 ? argv[2] : NULL), comm); 

  MPI_Comm_free(&comm);

  MPI_Finalize();

  return 0;
}



/***********************************************************************************/
/*! This function tests the various graph partitioning and ordering routines */
/***********************************************************************************/
void TestParMetis_GPart(char *filename, char *xyzfile, MPI_Comm comm)
{
  int ncon, nparts, npes, mype, opt2, realcut;
  GraphType graph, mgraph;
  idxtype *part, *mpart, *savepart, *order, *sizes;
  int numflag=0, wgtflag=0, options[10], edgecut, ndims;
  float ipc2redist, *xyz=NULL, *tpwgts = NULL, ubvec[MAXNCON];

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  ParallelReadGraph(&graph, filename, comm);
  if (xyzfile)
    xyz = ReadTestCoordinates(&graph, xyzfile, &ndims, comm);
  MPI_Barrier(comm);

  part   = idxmalloc(graph.nvtxs, "TestParMetis_V3: part");
  tpwgts = fmalloc(MAXNCON*npes*2, "TestParMetis_V3: tpwgts");
  sset(MAXNCON, 1.05, ubvec);

  graph.vwgt = idxsmalloc(graph.nvtxs*5, 1, "TestParMetis_GPart: vwgt");


  /*======================================================================
  / ParMETIS_V3_PartKway
  /=======================================================================*/
  options[0] = 1;
  options[1] = 3;
  options[2] = 1;
  wgtflag = 2;
  numflag = 0;
  edgecut = 0;

  for (nparts=2*npes; nparts>=npes/2 && nparts > 0; nparts = nparts/2) {
    for (ncon=1; ncon<=NCON; ncon++) {
      if (ncon > 1 && nparts > 1)
        Mc_AdaptGraph(&graph, part, ncon, nparts, comm);
      else
        idxset(graph.nvtxs, 1, graph.vwgt);

      if (mype == 0)
        printf("\nTesting ParMETIS_V3_PartKway with ncon: %d, nparts: %d\n", ncon, nparts);

      sset(nparts*ncon, 1.0/(float)nparts, tpwgts);
      ParMETIS_V3_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
          NULL, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, 
          &edgecut, part, &comm);

      realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
      if (mype == 0) {
        printf("ParMETIS_V3_PartKway reported a cut of %d [%s:%d]\n", edgecut,
            (edgecut == realcut ? "OK" : "ERROR"), realcut);
      }

      if (mype == 0)
        printf("\nTesting ParMETIS_V3_RefineKway with ncon: %d, nparts: %d\n", ncon, nparts);

      options[3] = PARMETIS_PSR_UNCOUPLED;
      ParMETIS_V3_RefineKway(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
          NULL, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, 
          &edgecut, part, &comm);

      realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
      if (mype == 0) {
        printf("ParMETIS_V3_RefineKway reported a cut of %d [%s:%d]\n", edgecut,
            (edgecut == realcut ? "OK" : "ERROR"), realcut);
      }
    }
  }


  /*======================================================================
  / ParMETIS_V3_PartGeomKway 
  /=======================================================================*/
  if (xyzfile != NULL) {
    options[0] = 1;
    options[1] = 3;
    options[2] = 1;
    wgtflag = 2;
    numflag = 0;

    for (nparts=2*npes; nparts>=npes/2 && nparts > 0; nparts = nparts/2) {
      for (ncon=1; ncon<=NCON; ncon++) {
        if (ncon > 1)
          Mc_AdaptGraph(&graph, part, ncon, nparts, comm);
        else
          idxset(graph.nvtxs, 1, graph.vwgt);
  
        if (mype == 0)
          printf("\nTesting ParMETIS_V3_PartGeomKway with ncon: %d, nparts: %d\n", ncon, nparts);
  
        sset(nparts*ncon, 1.0/(float)nparts, tpwgts);
        ParMETIS_V3_PartGeomKway(graph.vtxdist, graph.xadj, graph.adjncy, graph.vwgt, 
            NULL, &wgtflag, &numflag, &ndims, xyz, &ncon, &nparts, tpwgts, ubvec, 
            options, &edgecut, part, &comm);
  
        realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
        if (mype == 0) 
          printf("ParMETIS_V3_PartGeomKway reported a cut of %d [%s:%d]\n", edgecut,
              (edgecut == realcut ? "OK" : "ERROR"), realcut);
      }
    }
  }



  /*======================================================================
  / ParMETIS_V3_PartGeom 
  /=======================================================================*/
  if (xyz != NULL) {
    wgtflag = 0;
    numflag = 0;
    if (mype == 0)
      printf("\nTesting ParMETIS_V3_PartGeom\n");

      ParMETIS_V3_PartGeom(graph.vtxdist, &ndims, xyz, part, &comm); 

    realcut = ComputeRealCut(graph.vtxdist, part, filename, comm);
    if (mype == 0) 
      printf("ParMETIS_V3_PartGeom reported a cut of %d\n", realcut);
  }


  /*======================================================================
  / Coupled ParMETIS_V3_RefineKway 
  /=======================================================================*/
  options[0] = 1;
  options[1] = 3;
  options[2] = 1;
  options[3] = PARMETIS_PSR_COUPLED;
  nparts = npes;
  wgtflag = 0;
  numflag = 0;
  ncon = 1;
  sset(nparts*ncon, 1.0/(float)nparts, tpwgts);

  if (mype == 0)
    printf("\nTesting coupled ParMETIS_V3_RefineKway with default options (before move)\n");

  ParMETIS_V3_RefineKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, 
      &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, 
      part, &comm);





  /* Compute a good partition and move the graph. Do so quietly! */
  options[0] = 0;
  nparts = npes;
  wgtflag = 0;
  numflag = 0;
  ncon = 1;
  sset(nparts*ncon, 1.0/(float)nparts, tpwgts);
  ParMETIS_V3_PartKway(graph.vtxdist, graph.xadj, graph.adjncy, NULL, NULL, 
      &wgtflag, &numflag, &ncon, &npes, tpwgts, ubvec, options, &edgecut, 
      part, &comm);
  TestMoveGraph(&graph, &mgraph, part, comm);
  GKfree((void **)&(graph.vwgt), LTERM);
  mpart    = idxsmalloc(mgraph.nvtxs, mype, "TestParMetis_V3: mpart");
  savepart = idxmalloc(mgraph.nvtxs, "TestParMetis_V3: savepart");



  /*======================================================================
  / Coupled ParMETIS_V3_RefineKway 
  /=======================================================================*/
  options[0] = 1;
  options[1] = 3;
  options[2] = 1;
  options[3] = PARMETIS_PSR_COUPLED;
  nparts  = npes;
  wgtflag = 0;
  numflag = 0;

  for (ncon=1; ncon<=NCON; ncon++) {
    if (mype == 0)
      printf("\nTesting coupled ParMETIS_V3_RefineKway with ncon: %d, nparts: %d\n", ncon, nparts);

    sset(nparts*ncon, 1.0/(float)nparts, tpwgts);
      ParMETIS_V3_RefineKway(mgraph.vtxdist, mgraph.xadj, mgraph.adjncy, NULL, NULL, 
          &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, 
          mpart, &comm);

    realcut = ComputeRealCutFromMoved(graph.vtxdist, mgraph.vtxdist, part, mpart, 
                  filename, comm);
    if (mype == 0) 
      printf("ParMETIS_V3_RefineKway reported a cut of %d [%s:%d]\n", edgecut,
          (edgecut == realcut ? "OK" : "ERROR"), realcut);
  }


  /*======================================================================
  / ParMETIS_V3_AdaptiveRepart
  /=======================================================================*/
  mgraph.vwgt  = idxsmalloc(mgraph.nvtxs*NCON, 1, "TestParMetis_V3: mgraph.vwgt");
  mgraph.vsize = idxsmalloc(mgraph.nvtxs, 1, "TestParMetis_V3: mgraph.vsize");
  AdaptGraph(&mgraph, 4, comm); 
  options[0] = 1;
  options[1] = 7;
  options[2] = 1;
  options[3] = PARMETIS_PSR_COUPLED;
  wgtflag = 2;
  numflag = 0;

  for (nparts=2*npes; nparts>=npes/2; nparts = nparts/2) {
    options[0] = 0;
    ncon    = 1;
    wgtflag = 0;
    sset(nparts*ncon, 1.0/(float)nparts, tpwgts);
    ParMETIS_V3_PartKway(mgraph.vtxdist, mgraph.xadj, mgraph.adjncy, NULL, NULL, 
        &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, 
        savepart, &comm);

    options[0] = 1;
    wgtflag    = 2;

    for (ncon=1; ncon<=NCON; ncon++) {
      sset(nparts*ncon, 1.0/(float)nparts, tpwgts);

      if (ncon > 1)
        Mc_AdaptGraph(&mgraph, savepart, ncon, nparts, comm);
      else
        AdaptGraph(&mgraph, 4, comm); 

      for (ipc2redist=1000.0; ipc2redist>=0.001; ipc2redist/=1000.0) {
        idxcopy(mgraph.nvtxs, savepart, mpart);

        if (mype == 0)
          printf("\nTesting ParMETIS_V3_AdaptiveRepart with ipc2redist: %.3f, ncon: %d, nparts: %d\n", 
              ipc2redist, ncon, nparts);

        ParMETIS_V3_AdaptiveRepart(mgraph.vtxdist, mgraph.xadj, mgraph.adjncy, 
            mgraph.vwgt, mgraph.vsize, NULL, &wgtflag, &numflag, &ncon, &nparts, 
            tpwgts, ubvec, &ipc2redist, options, &edgecut, mpart, &comm);

        realcut = ComputeRealCutFromMoved(graph.vtxdist, mgraph.vtxdist, part, mpart, 
                      filename, comm);
        if (mype == 0) 
          printf("ParMETIS_V3_AdaptiveRepart reported a cut of %d [%s:%d]\n", edgecut,
              (edgecut == realcut ? "OK" : "ERROR"), realcut);
      }
    }
  }

  GKfree((void **)&tpwgts, &part, &mpart, &savepart, &xyz, &mgraph.vwgt, &mgraph.vsize, LTERM);
}



/******************************************************************************/
/*! This function takes a partition vector that is distributed and reads in
    the original graph and computes the edgecut */
/******************************************************************************/
int ComputeRealCut(idxtype *vtxdist, idxtype *part, char *filename, MPI_Comm comm)
{
  int i, j, nvtxs, mype, npes, cut;
  idxtype *xadj, *adjncy, *gpart;
  MPI_Status status;

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (mype != 0) {
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_DATATYPE, 0, 1, comm);
  }
  else {  /* Processor 0 does all the rest */
    gpart = idxmalloc(vtxdist[npes], "ComputeRealCut: gpart");
    idxcopy(vtxdist[1], part, gpart);

    for (i=1; i<npes; i++) 
      MPI_Recv((void *)(gpart+vtxdist[i]), vtxdist[i+1]-vtxdist[i], IDX_DATATYPE, i, 1, comm, &status);

    ReadMetisGraph(filename, &nvtxs, &xadj, &adjncy);

    /* OK, now compute the cut */
    for (cut=0, i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (gpart[i] != gpart[adjncy[j]])
          cut++;
      }
    }
    cut = cut/2;

    GKfree((void **)&gpart, &xadj, &adjncy, LTERM);

    return cut;
  }
  return 0;
}


/******************************************************************************/
/*! This function takes a partition vector of the original graph and that of the
    moved graph and computes the cut of the original graph based on the moved
    graph */
/*******************************************************************************/
int ComputeRealCutFromMoved(idxtype *vtxdist, idxtype *mvtxdist, idxtype *part, 
    idxtype *mpart, char *filename, MPI_Comm comm)
{
  int i, j, nvtxs, mype, npes, cut;
  idxtype *xadj, *adjncy, *gpart, *gmpart, *perm, *sizes;
  MPI_Status status;


  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  if (mype != 0) {
    MPI_Send((void *)part, vtxdist[mype+1]-vtxdist[mype], IDX_DATATYPE, 0, 1, comm);
    MPI_Send((void *)mpart, mvtxdist[mype+1]-mvtxdist[mype], IDX_DATATYPE, 0, 1, comm);
  }
  else {  /* Processor 0 does all the rest */
    gpart = idxmalloc(vtxdist[npes], "ComputeRealCut: gpart");
    idxcopy(vtxdist[1], part, gpart);
    gmpart = idxmalloc(mvtxdist[npes], "ComputeRealCut: gmpart");
    idxcopy(mvtxdist[1], mpart, gmpart);

    for (i=1; i<npes; i++) {
      MPI_Recv((void *)(gpart+vtxdist[i]), vtxdist[i+1]-vtxdist[i], IDX_DATATYPE, i, 1, comm, &status);
      MPI_Recv((void *)(gmpart+mvtxdist[i]), mvtxdist[i+1]-mvtxdist[i], IDX_DATATYPE, i, 1, comm, &status);
    }

    /* OK, now go and reconstruct the permutation to go from the graph to mgraph */
    perm  = idxmalloc(vtxdist[npes], "ComputeRealCut: perm");
    sizes = idxsmalloc(npes+1, 0, "ComputeRealCut: sizes");

    for (i=0; i<vtxdist[npes]; i++)
      sizes[gpart[i]]++;
    MAKECSR(i, npes, sizes);
    for (i=0; i<vtxdist[npes]; i++)
      perm[i] = sizes[gpart[i]]++;

    /* Ok, now read the graph from the file */
    ReadMetisGraph(filename, &nvtxs, &xadj, &adjncy);

    /* OK, now compute the cut */
    for (cut=0, i=0; i<nvtxs; i++) {
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (gmpart[perm[i]] != gmpart[perm[adjncy[j]]])
          cut++;
      }
    }
    cut = cut/2;

    GKfree((void **)&gpart, &gmpart, &perm, &sizes, &xadj, &adjncy, LTERM);

    return cut;
  }

  return 0;
}



/******************************************************************************
* This function takes a graph and its partition vector and creates a new
* graph corresponding to the one after the movement
*******************************************************************************/
void TestMoveGraph(GraphType *ograph, GraphType *omgraph, idxtype *part, MPI_Comm comm)
{
  int npes, mype;
  CtrlType ctrl;
  WorkSpaceType wspace;
  GraphType *graph, *mgraph;
  int options[5] = {0, 0, 1, 0, 0};

  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);

  SetUpCtrl(&ctrl, npes, 0, comm); 
  ctrl.CoarsenTo = 1;  /* Needed by SetUpGraph, otherwise we can FP errors */
  graph = SetUpGraph(&ctrl, ograph->vtxdist, ograph->xadj, NULL, ograph->adjncy, NULL, 0);
  AllocateWSpace(&ctrl, graph, &wspace);

  SetUp(&ctrl, graph, &wspace);
  graph->where = part;
  graph->ncon = 1;
  mgraph = Mc_MoveGraph(&ctrl, graph, &wspace);

  omgraph->gnvtxs = mgraph->gnvtxs;
  omgraph->nvtxs = mgraph->nvtxs;
  omgraph->nedges = mgraph->nedges;
  omgraph->vtxdist = mgraph->vtxdist;
  omgraph->xadj = mgraph->xadj;
  omgraph->adjncy = mgraph->adjncy;
  mgraph->vtxdist = NULL;
  mgraph->xadj = NULL;
  mgraph->adjncy = NULL;
  FreeGraph(mgraph);

  graph->where = NULL;
  FreeInitialGraphAndRemap(graph, 0, 1);
  FreeWSpace(&wspace);
}  


/*****************************************************************************
*  This function sets up a graph data structure for partitioning
*****************************************************************************/
GraphType *SetUpGraph(CtrlType *ctrl, idxtype *vtxdist, idxtype *xadj,
   idxtype *vwgt, idxtype *adjncy, idxtype *adjwgt, int wgtflag)
{
  int mywgtflag;

  mywgtflag = wgtflag;
  return Mc_SetUpGraph(ctrl, 1, vtxdist, xadj, vwgt, adjncy, adjwgt, &mywgtflag);
}


