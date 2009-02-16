#ifndef SEGMENTATION
#define SEGMENTATION

#include "config.h"
#include "math.h"
#include "timing.h"
#include "editing.h"


#define A           0
#define B           1
#define C           2
#define ITERS_MAX   10


void segFromVoxelCoordsToSiteCoords (short int voxel[3], short int site[3], Vis *vis);
void segFromSiteCoordsToVoxelCoords (short int site[3], short int voxel[3], Vis *vis);
void segFromSiteCoordsToGridCoords (short int site[3], int *b_id, int *s_id, Vis *vis);
Site *segSitePtr (short int site[3], Vis *vis);
void segSetSiteData (Site *site_p, unsigned int config, unsigned int label);
Site *segStackSitePtr (Vis *vis);
int segIsSuperficialSite (short int site[3], Vis *vis);
int segRayVsDisc (Triangle *t_p, float x1[3], float x2[3], float t_max, float *t);
float segInterpolatedGrey (short int site[3], Vis *vis);
int segEstimateNormal (short int site[3], float nor[], Vis *vis);
void segEstimateDiameter (short int site[3], float nor[3], float *diameter, Vis *vis);
int segCreateOptimisedTriangle (unsigned int site_type, short int site[3], Vis *vis);
int segIsSegmentIntercepted (short int site[3], int dir, Site *site_p, Vis *vis);
int segSegmentation (Vis *vis);
int segUpdateSegmentation (Vis *vis);
void segSetBoundaryConfigurations (Vis *vis);


#endif // SEGMENTATION
