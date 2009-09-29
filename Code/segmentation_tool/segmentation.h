#ifndef SEGMENTATION
#define SEGMENTATION

#include "config.h"
#include "timing.h"
#include "editing.h"
#ifdef MESH
#include "rt.h"
#endif


#define A           0
#define B           1
#define C           2
#ifndef MESH
#define ITERS_MAX   10
#endif


#ifndef MESH
void segFromVoxelCoordsToSiteCoords (short int voxel[3], short int site[3], Vis *vis);
void segFromSiteCoordsToVoxelCoords (short int site[3], short int voxel[3], Vis *vis);
#else
void segFromMeshCoordsToSiteCoords (double pos[3], short int site[3], Vis *vis);
void segFromSiteCoordsToMeshCoords (short int site[3], double pos[3], Vis *vis);
#endif
void segFromSiteCoordsToGridCoords (short int site[3], int *b_id, int *s_id, Vis *vis);
Site *segSitePtr (short int site[3], Vis *vis);
void segSetSiteData (Site *site_p, unsigned int config, unsigned int label);
Site *segStackSitePtr (Vis *vis);
int segIsSuperficialSite (short int site[3], Vis *vis);
int segRayVsDisc (BoundaryTriangle *t_p, double x1[3], double x2[3], double t_max, double *t);
#ifndef MESH
double segInterpolatedGrey (short int site[3], Vis *vis);
int segEstimateNormal (short int site[3], double nor[3], Vis *vis);
void segEstimateDiameter (short int site[3], double nor[3], double *diameter, Vis *vis);
#else
int segEstimateExtrema (Hit *first_hit, Hit *second_hit, Vis *vis);
void segEstimateDiameter (double *diameter, Vis *vis);
#endif
#ifndef MESH
int segCreateOptimisedTriangle (unsigned int site_type, short int site[3], Vis *vis);
#else
int segCreateOptimisedTriangle (unsigned int site_type, Vis *vis);
void segCalculateBoundarySiteData (unsigned int site_cfg, short int site[3], double boundary_nor[3],
				   double *boundary_dist, double wall_nor[3], double *wall_dist,
				   double cut_dist[14], Vis *vis);
#endif
#ifndef MESH
int segIsSegmentIntercepted (short int site[3], int dir, Site *site_p, Vis *vis);
#else
int segIsSegmentIntercepted (short int site[3], int dir, Site *site_p, Hit *hit, Vis *vis);
#endif
int segSegmentation (Vis *vis);
#ifndef MESH
int segUpdateSegmentation (Vis *vis);
#endif
void segSetBoundaryConfigurations (Vis *vis);


#endif // SEGMENTATION
