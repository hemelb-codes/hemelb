#ifndef SEGMENTATION
#define SEGMENTATION

#include "config.h"
#include "timing.h"
#include "editing.h"
#include "rt.h"

#define A           0
#define B           1
#define C           2

void segFromMeshCoordsToSiteCoords(double pos[3], site_index site[3], Vis *vis);
void segFromSiteCoordsToMeshCoords(site_index site[3], double pos[3], Vis *vis);
void segFromSiteCoordsToGridCoords(site_index site[3],
                                   block_id *b_id,
                                   site_id *s_id,
                                   Vis *vis);
Site *segSitePtr(site_index site[3], Vis *vis);
void segSetSiteData(Site *site_p, config_data config, unsigned int label);
Site *segStackSitePtr(Vis *vis);
int segIsSuperficialSite(site_index site[3], Vis *vis);
int segRayVsDisc(BoundaryTriangle *t_p,
                 double x1[3],
                 double x2[3],
                 double t_max,
                 double *t);
int segEstimateExtrema(Hit *first_hit, Hit *second_hit, Vis *vis);
void segEstimateDiameter(double *diameter, Vis *vis);
triangle_id segCreateOptimisedTriangle(site_index, Vis *vis);
void segCalculateBoundarySiteData(config_data site_cfg,
                                  site_index site[3],
                                  double boundary_nor[3],
                                  double *boundary_dist,
                                  double wall_nor[3],
                                  double *wall_dist,
                                  double cut_dist[14],
                                  Vis *vis);
int segIsSegmentIntercepted(site_index site[3],
                            int dir,
                            Site *site_p,
                            Hit *hit,
                            Vis *vis);
int segSegmentation(Vis *vis);
void segSetBoundaryConfigurations(Vis *vis);

#endif // SEGMENTATION
