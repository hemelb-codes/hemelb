#ifndef VIS
#define VIS

#define GL_GLEXT_PROTOTYPES


#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include "usage.h"
#include "segmentation.h"
#include "io.h"
#include "menu.h"


void visColorPalette (int iters, int res_factor, float col[3]);
ScreenVoxel *visScreenVoxelPtr (short int x[2], Vis *vis);
void visOpenWindow (int pixels_x, int pixels_y);
void visProjection (Vis *vis);
void visProject (float px1[3], float px2[3]);
void visAntiProject (float px1[3], float px2[3]);
void visRescaleViewpoint (float scale, Vis *vis);
void visCalculateSceneCenter (Vis *vis);
void visCreateCubeDisplayList (void);
void visDeleteCubeDisplayList (void);
void visVisualiseString (float r, float g, float b, int x, int y, char *string, void *font);
void visVisualiseTrianglePars (int b_id, int t_id, Vis *vis);
void visVisualiseSiteData (short int site[3], Vis *vis);
void visVisualiseVisData (Vis *vis);
void visVisualiseTriangles (Vis *vis);
void visVisualiseDiscs (Vis *vis);
void visVisualiseActiveBoundaryVoxel (Vis *vis);
void visVisualiseSelectedSlice (Vis *vis);
void visVisualiseFluidSitesWithPoints (Vis *vis);
void visVisualiseFluidSitesWithCubes (Vis *vis);
void visVisualiseSystem (Vis *vis);
void GLUTCALLBACK Visualise (void);
void visInitBoundaries (Vis *vis);
void visEndBoundaries (Vis *vis);
void visInit (int argc, char *argv[], Vis *vis);
void visEnd (Vis *vis);

#endif // VIS
