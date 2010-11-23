#ifndef EDITING
#define EDITING

#include "config.h"
#include "vis.h"
#include "segmentation.h"

#define VERTEX    0
#define SITE      1

void editTriangleCenter(BoundaryTriangle *t_p, double c[3]);
void editTriangleNormal(BoundaryTriangle *t_p, double nor[3]);
double editTriangleRadius(BoundaryTriangle *t_p);
void editCalculateTriangleData(BoundaryTriangle *t_p);
void editDeleteTriangle(block_id b_id, triangle_id t_id, Vis *vis);
void editInvertTriangleNormal(block_id b_id, triangle_id t_id, Vis *vis);
void editChangeTrianglePars(block_id b_id,
                            triangle_id t_id,
                            double dp_avg,
                            double dp_amp,
                            double dp_phs,
                            Vis *vis);
void editRescaleGrid(Vis *vis);
void editProjectBoundariesToScreenVoxels(Vis *vis);
void editMoveTriangleVertexWithMouse(int x0[3], Vis *vis);
void editRotateTriangleWithMouse(double dx0,
                                 double dy0,
                                 block_id b_id,
                                 triangle_id t_id,
                                 Vis *vis);
void editScaleTriangleWithMouse(double scaling_factor,
                                block_id b_id,
                                triangle_id t_id,
                                Vis *vis);
void editRotateViewpointWithMouse(double dx0, double dy0, Vis *vis);
void editMouseFunction(int button, int state, int x, int y, Vis *vis);
void editMotionFunction(int x, int y, Vis *vis);
void editPassiveMotionFunction(int x, int y, Vis *vis);
void GLUTCALLBACK MouseFunction(int button, int state, int x, int y);
void GLUTCALLBACK MotionFunction(int x, int y);
void GLUTCALLBACK PassiveMotionFunction(int x, int y);
void GLUTCALLBACK Reshape(GLsizei w, GLsizei h);

#endif // EDITING
