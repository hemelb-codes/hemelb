#ifdef MESH
#ifndef RAYTRACING
#define RAYTRACING

#include "config.h"
#include "segmentation.h"


int rtTriangleVsRay (MeshTriangle *triangle, Ray *ray, Hit *hit);
void rtAABBvsRay (double aabb[], Ray *ray);
void rtInitRayTracing (Mesh *mesh);
int rtVoxelVsRay (Voxel *voxel, Ray *ray, Hit *hit, Mesh *mesh);
int rtTraceRay (Ray *ray, Hit *hit, Mesh *mesh);
int rtTracePrimaryRay (Hit *first_hit, Vis *vis);
int rtTraceSecondaryRay (Hit *first_hit, Hit *second_hit, Vis *vis);
void rtEndRayTracing (Mesh *mesh);


#endif // RAYTRACING
#endif // MESH
