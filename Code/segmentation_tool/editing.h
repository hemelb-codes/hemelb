#ifndef EDITING
#define EDITING

#include "config.h"
#include "math.h"
#include "vis.h"
#include "segmentation.h"

#define EPSILON   1.0e-30
#define VERTEX    0
#define SITE      1


void editTriangleCenter (Triangle *t_p, float c[3]);
void editTriangleNormal (Triangle *t_p, float nor[3]);
float editTriangleRadius (Triangle *t_p);
void editCalculateTriangleData (Triangle *t_p);
void editDeleteTriangle (int b_id, int t_id, Vis *vis);
void editInvertTriangleNormal (int b_id, int t_id, Vis *vis);
void editChangeTrianglePars (int b_id, int t_id, float dp_avg, float dp_amp, float dp_phs, Vis *vis);
void editRescaleSystem (Vis *vis);
void editRescaleTriangles (float scale, Vis *vis);
void editProjectBoundariesToScreenVoxels (Vis *vis);
void editMoveTriangleVertexWithMouse (int x0[3], Vis *vis);
void editRotateTriangleWithMouse (float dx0, float dy0, int b_id, int t_id, Vis *vis);
void editScaleTriangleWithMouse (float scaling_factor, int b_id, int t_id, Vis *vis);
void editRotateViewpointWithMouse (float dx0, float dy0, Vis *vis);
void editMouseFunction (int button, int state, int x, int y, Vis *vis);
void editMotionFunction (int x, int y, Vis *vis);
void editPassiveMotionFunction (int x, int y, Vis *vis);
void GLUTCALLBACK MouseFunction (int button, int state, int x, int y);
void GLUTCALLBACK MotionFunction (int x, int y);
void GLUTCALLBACK PassiveMotionFunction (int x, int y);
void GLUTCALLBACK Reshape (GLsizei w, GLsizei h);

#endif // EDITING
