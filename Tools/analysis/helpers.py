import math

def space(grid_coords,voxel_size,origin):
  return [grid_coord*voxel_size+origin for grid_coord in grid_coords] 

def aligned_cylindrical(grid_coords,voxel_size,origin):
  cartesian=space(grid_coords,voxel_size,origin)
  return [[pow(pow(x,2)+pow(y,2),0.5),math.atan2(y,x),z] for x,y,z in cartesian]

eval_helpers=globals()
