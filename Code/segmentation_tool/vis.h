#ifndef VIS
#define VIS

#define GL_GLEXT_PROTOTYPES


#include "usage.h"
#include "segmentation.h"
#include "io.h"
#include "menu.h"
#ifdef MESH
#include "rt.h"
#endif


void visColorPalette (int iters, int res_factor, double col[3]);
ScreenVoxel *visScreenVoxelPtr (short int x[2], Vis *vis);
void visOpenWindow (int pixels_x, int pixels_y);
void visProjection (Vis *vis);
void visProject (double px1[3], double px2[3]);
void visAntiProject (double px1[3], double px2[3]);
void visCalculateSceneCenter (Vis *vis);
void visCreateCubeDisplayList (Vis *vis);
void visDeleteCubeDisplayList (void);
void visVisualiseString (double r, double g, double b, int x, int y, char *string, void *font);
void visVisualiseTrianglePars (int b_id, int t_id, Vis *vis);
#ifndef MESH
void visVisualiseSiteData (short int site[3], Vis *vis);
#else
void visVisualiseSiteData (Vis *vis);
#endif
void visVisualiseVisData (Vis *vis);
void visVisualiseTriangles (Vis *vis);
void visVisualiseDiscs (Vis *vis);
void visVisualiseActiveBoundaryVoxel (Vis *vis);
#ifdef MESH
void visVisualiseHitData (Hit *first_hit, Hit *second_hit, Vis *vis);
void visVisualiseMesh (Vis *vis);
#else
void visVisualiseSelectedSlice (Vis *vis);
void visUpdateSegmentation (Vis *vis);
#endif
void visVisualiseFluidSitesWithPoints (Vis *vis);
void visVisualiseFluidSitesWithCubes (Vis *vis);
#ifdef MESH
void visVisualiseInterceptions (Vis *vis);
#endif
void visVisualiseSystem (Vis *vis);
void GLUTCALLBACK Visualise (void);
void visInitBoundaries (Vis *vis);
void visEndBoundaries (Vis *vis);
void visInit (int argc, char *argv[], Vis *vis);
void visEnd (Vis *vis);

#endif // VIS
