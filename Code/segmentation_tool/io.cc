#include <cstdio>
#include <rpc/rpc.h>

#include "io_dicom.h"
#include "io.h"

#ifdef USE_TIFFLIB
#include <tiffio.h>
#include <cstdlib>
#include <cstring>
#endif // USE_TIFFLIB


#ifdef MESH
unsigned int PRESSURE_EDGE_MASK = 1U << 31U;
#endif


#ifndef MESH
/*
int ioGetDir (string dir, vector<string> &files)
{
  DIR *dp;
  
  struct dirent *dirp;
  
  
  if ((dp = opendir(dir.c_str())) == NULL)
    {
      cout << "Error(" << errno << ") opening " << dir << endl;
      return errno;
    }
  
  while ((dirp = readdir(dp)) != NULL)
    {
      files.push_back(string(dirp->d_name));
    }
  closedir(dp);
  
  return 0;
}



void ioGetFileNames (Vis *vis)
{
  string dir = string(vis->input_path);
  
  ioGetDir (dir, vis->file_list);
}



void ioReadFirstSlice (Vis *vis)
{
  char ppm_type[16];
  
  string file_name;
  
  ifstream input_file;
  
  
  file_name = string(vis->input_path) + vis->file_list[2];
  
  input_file.open(file_name.c_str());
  
  input_file >> ppm_type;
  input_file >> vis->input_voxels[0];
  input_file >> vis->input_voxels[1];
  input_file >> vis->pixel_depth;
  
  input_file.close();
}


void ioReadSlice (int slice_id, Vis *vis)
{
  int voxel[3];
  
  char ppm_type[16];
  
  string file_name;
  
  ifstream input_file;
  
  voxel[2] = slice_id;
  
  file_name = string(vis->input_path) + vis->file_list[ slice_id+2 ];
  
  input_file.open(file_name.c_str());
  
  input_file >> ppm_type;
  input_file >> vis->input_voxels[0];
  input_file >> vis->input_voxels[1];
  input_file >> vis->pixel_depth;
  
  for (voxel[1] = 0; voxel[1] < vis->input_voxels[1]; voxel[1]++)
    {
      for (voxel[0] = 0; voxel[0] < vis->input_voxels[0]; voxel[0]++)
	{
	  int temp;
	  
	  input_file >> temp;
	  vis->voxel[ VoxelId(voxel,vis->input_voxels) ] = temp;
	  
	  vis->grey_min = fmin(vis->grey_min, temp);
	  vis->grey_max = fmax(vis->grey_max, temp);
	}
    }

  input_file.close();

}
*/


void ioReadConfig (Vis *vis)
{
  int l;
  
  // ioGetFileNames (vis);
  
  // vis->input_voxels[2] = vis->file_list.size() - 2;
  
  DICOM_read DICOM_input = DICOM_read(vis->input_path);
  
  DICOM_input.DICOM_volume_size(vis->input_voxels[0], vis->input_voxels[1], vis->input_voxels[2]);
  
  vis->output_voxels[0] = vis->input_voxels[0] * vis->res_factor;
  vis->output_voxels[1] = vis->input_voxels[1] * vis->res_factor;
  vis->output_voxels[2] = (int)(vis->input_voxels[2] * (vis->slice_size/vis->pixel_size) *
				(double)vis->output_voxels[0]/(double)vis->input_voxels[0]);
  
  for (l = 0; l < 3; l++)
    vis->blocks[l] = vis->output_voxels[l] >> SHIFT;
  
  for (l = 0; l < 3; l++)
    if ((vis->blocks[l] << SHIFT) < vis->output_voxels[l]) ++vis->blocks[l];
  
  for (l = 0; l < 3; l++)
    vis->sites[l] = vis->blocks[l] * BLOCK_SIZE;
  
  for (l = 0; l < 3; l++)
    vis->dim[l] = vis->output_voxels[l] / vis->res_factor;
  
  for (l = 0; l < 3; l++)
    vis->half_dim[l] = 0.5 * vis->dim[l];
  
  vis->system_size = fmax(vis->dim[0], fmax(vis->dim[1], vis->dim[2]));
  
  vis->tot_blocks = vis->blocks[0] * vis->blocks[1] * vis->blocks[2];
  
  vis->voxel = new short int[vis->input_voxels[0] * vis->input_voxels[1]
      * vis->input_voxels[2]];
  
  vis->block = new Block[vis->tot_blocks];
  
  vis->stack_sites_max = SITES_PER_BLOCK * 100000;
  
  vis->stack_site = new Site[vis->stack_sites_max];
  /*
  for (i = 0; i < vis->input_voxels[2]; i++)
    {
      ioReadSlice (i, vis);
    }
  */
  DICOM_input.DICOM_copy_data(vis->voxel, vis->input_voxels[0] * vis->input_voxels[1] * vis->input_voxels[2]);
  
  vis->grey_min = +1e9;
  vis->grey_max = -1e9;
  
  for (l = 0; l < vis->input_voxels[0] * vis->input_voxels[1] * vis->input_voxels[2]; l++)
    {
      if (vis->voxel[l] < -30000 || vis->voxel[l] > 30000)
	{
	  continue;
	}
      if (vis->voxel[l] < vis->grey_min) vis->grey_min = vis->voxel[l];
      if (vis->voxel[l] > vis->grey_max) vis->grey_max = vis->voxel[l];
    }
  vis->coord[A] = new Coord[COORD_BUFFER_SIZE_A];
  vis->coord[B] = new Coord[COORD_BUFFER_SIZE_B];
  vis->coord[C] = new Coord[COORD_BUFFER_SIZE_C];
  
  vis->stress_type = VON_MISES_STRESS;
}
#else // MESH


void ioReadConfig (Vis *vis)
{
  /* Read triangles from the STL file, into vis->mesh. The triangles
   * are scaled up by 1% (presumably to prevent leaks?) and translated
   * such that the centre of mass (assuming equal triangle mass) is at
   * the origin.
   *  
   * It then does some set up for the voxels that I don't yet
   * understand.
   *  
   * An STL file looks like:
   *  
   * solid ascii
   *   facet normal -0.518317 -0.799163 -0.304445
   *     outer loop
   *       vertex 36.9201 20.5062 0.419999
   *       vertex 36.8551 20.5519 0.410923
   *       vertex 36.8745 20.5816 0.299708
   *     endloop
   *   endfacet
   *   ... more facets ...
   * endsolid
   */
  double min_x[3], MaxXValue[3];	// Bounding box of mesh?
  double dx[3];			// Unknown
  
  int l, m, n;			// Indices
  
  char key_word[16];
  
  ifstream input_file;
  
  
  vis->mesh.triangles_max = 10000;
  vis->mesh.triangle = new MeshTriangle[vis->mesh.triangles_max];
  vis->mesh.triangles = 0;	// This counts the number of triangles
  
  for (l = 0; l < 3; l++)
    {
      min_x[l] =  1.0e+30;
      MaxXValue[l] = -1.0e+30;
      vis->mesh.centre[l] =  0.0;
    }
  input_file.open(vis->input_file);
  
  // Scan to the first facet
  while (strcmp(key_word, "facet"))
    input_file >> key_word;
  
  // Until the end of the STL
  while (strcmp(key_word, "endsolid"))
    {
      // If we've run out of triangle space, alloc some more
      if (vis->mesh.triangles == vis->mesh.triangles_max)
	{
	  vis->mesh.triangles_max *= 2;
	  vis->mesh.triangle = (MeshTriangle *)realloc(vis->mesh.triangle,
						       sizeof(MeshTriangle) * vis->mesh.triangles_max);
	}
      // Scan to the normal entry
      while (strcmp(key_word, "normal"))
	input_file >> key_word;
      
      input_file >> vis->mesh.triangle[ vis->mesh.triangles ].nor[0];
      input_file >> vis->mesh.triangle[ vis->mesh.triangles ].nor[1];
      input_file >> vis->mesh.triangle[ vis->mesh.triangles ].nor[2];
      
      // Normalise the normal
      for (l = 0; l < 3; l++)
	vis->mesh.triangle[ vis->mesh.triangles ].nor[l] /=
	  sqrt(ScalarProd (vis->mesh.triangle[ vis->mesh.triangles ].nor, vis->mesh.triangle[ vis->mesh.triangles ].nor));
      
      // Scan to loop
      while (strcmp(key_word, "loop"))
	input_file >> key_word;
      
      for (m = 0; m < 3; m++)	// Must be a triangle
	{
	  input_file >> key_word; // Ditch "vertex", the read the coords
	  input_file >> vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[0];
	  input_file >> vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[1];
	  input_file >> vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[2];
	}
      float triangle_ctr_x[3];	// This is the centroid, strictly
      
      for (l = 0; l < 3; l++)
	triangle_ctr_x[l] = 0.0;
      
      for (m = 0; m < 3; m++)	// vertices
	for (l = 0; l < 3; l++)	// dimensions
	  {
	    triangle_ctr_x[l] += vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l];
	  }
      for (l = 0; l < 3; l++)
	triangle_ctr_x[l] /= 3.0;
      // Got the centroid

      for (m = 0; m < 3; m++) // vertices
	{
	for (l = 0; l < 3; l++)	// dimensions
	  {
	    // I think this makes the triangle 1% larger, keeping the same centroid
	    vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l] = triangle_ctr_x[l] +
	      1.01 * (vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l] - triangle_ctr_x[l]);
	  }
	}
      
      for (m = 0; m < 3; m++)
	{
	  for (l = 0; l < 3; l++)
	    {
	      // Update the min & max
	      min_x[l] = fmin(min_x[l], vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l]);
	      MaxXValue[l] = fmax(MaxXValue[l], vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l]);
	      // Increment running total of positions.
	      vis->mesh.centre[l] += vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l];
	    }
	}
      ++vis->mesh.triangles;
      
      input_file >> key_word;	// endloop
      input_file >> key_word;	// endfacet
      input_file >> key_word;	// facet
    }
  input_file.close();
  
  // Compute average from running tot.
  for (l = 0; l < 3; l++)
    vis->mesh.centre[l] /= vis->mesh.triangles * 3;
  
  vis->mesh.voxel_size = 0.0;
  
  for (n = 0; n < vis->mesh.triangles; n++) // for each tri
    {
      // Translate such that vis->mesh.centre is the origin
      for (m = 0; m < 3; m++)	// vertex
	for (l = 0; l < 3; l++)	// dimen
	  vis->mesh.triangle[n].v[m].pos[l] -= vis->mesh.centre[l];
      
      // Add perimeter of tri to vis->mesh.voxel_size
      for (l = 0; l < 3; l++)
	dx[l] = vis->mesh.triangle[n].v[1].pos[l] - vis->mesh.triangle[n].v[0].pos[l];
      
      vis->mesh.voxel_size += sqrtf(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      
      for (l = 0; l < 3; l++)
	dx[l] = vis->mesh.triangle[n].v[2].pos[l] - vis->mesh.triangle[n].v[0].pos[l];
      
      vis->mesh.voxel_size += sqrtf(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      
      for (l = 0; l < 3; l++)
	dx[l] = vis->mesh.triangle[n].v[2].pos[l] - vis->mesh.triangle[n].v[1].pos[l];
      
      vis->mesh.voxel_size += sqrtf(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
    }
  // Now the average side length of a tri
  vis->mesh.voxel_size /= (vis->mesh.triangles * 3);
    
  // Set number of voxels in along each axis to 1.5 * extent of
  // surface, over the size of a voxel.
  for (l = 0; l < 3; l++)
    vis->mesh.voxels[l] = (int)ceil(1.5 * (MaxXValue[l] - min_x[l]) / vis->mesh.voxel_size);
  
  // This is the total number of voxels
  vis->mesh.voxels[3] = vis->mesh.voxels[0] * vis->mesh.voxels[1] * vis->mesh.voxels[2];
  
  vis->mesh.voxel = new Voxel [vis->mesh.voxels[3]];
  
  // Dimensions of voxels in STL units
  for (l = 0; l < 3; l++)
    vis->mesh.dim[l] = vis->mesh.voxel_size * vis->mesh.voxels[l];

  for (l = 0; l < 3; l++)
    vis->mesh.half_dim[l] = 0.5 * vis->mesh.dim[l];
       
  for (l = 0; l < 3; l++)
    vis->input_voxels[l] = vis->mesh.voxels[l];
  
  for (l = 0; l < 3; l++)
    vis->output_voxels[l] = vis->input_voxels[l] * vis->res_factor;
  
  // Number of superblocks
  for (l = 0; l < 3; l++)
    vis->blocks[l] = vis->output_voxels[l] >> SHIFT;
  
  // Deal with the fact that bit shifting right is implicitly integer
  // division by a power of 2 (i.e. truncates)
  for (l = 0; l < 3; l++) {
    if ((vis->blocks[l] << SHIFT) < vis->output_voxels[l])
      ++vis->blocks[l];
  }
  
  for (l = 0; l < 3; l++)
    vis->sites[l] = vis->blocks[l] * BLOCK_SIZE;
  
  for (l = 0; l < 3; l++)
    vis->dim[l] = vis->mesh.dim[l];
  
  for (l = 0; l < 3; l++)
    vis->half_dim[l] = 0.5 * vis->dim[l];
  
  vis->system_size = fmax(vis->dim[0], fmax(vis->dim[1], vis->dim[2]));
  
  vis->tot_blocks = vis->blocks[0] * vis->blocks[1] * vis->blocks[2];
  
  vis->block = new Block[vis->tot_blocks];
  
  for (n = 0; n < vis->tot_blocks; n++)
    vis->block[n].site = NULL;
  
  vis->stack_sites_max = SITES_PER_BLOCK * 100000;
  
  vis->stack_site = new Site[vis->stack_sites_max];
  
  vis->coord[A] = new Coord [COORD_BUFFER_SIZE_A];
  vis->coord[B] = new Coord [COORD_BUFFER_SIZE_B];
  vis->coord[C] = new Coord [COORD_BUFFER_SIZE_C];
  
  vis->stress_type = SHEAR_STRESS;
}
#endif // MESH


void ioReadCheckpoint (Vis *vis)
{
  FILE *system_config;
  XDR xdr_config;  
  
  
  system_config = fopen (vis->checkpoint, "r");
  xdrstdio_create (&xdr_config, system_config, XDR_DECODE);
#ifndef MESH
  xdr_short  (&xdr_config, &vis->selected_voxel[0]);
  xdr_short  (&xdr_config, &vis->selected_voxel[1]);
  xdr_short  (&xdr_config, &vis->selected_voxel[2]);
  xdr_double (&xdr_config, &vis->selected_grey);
  xdr_int    (&xdr_config, &vis->res_factor);
#else
  xdr_double (&xdr_config, &vis->seed_pos[0]);
  xdr_double (&xdr_config, &vis->seed_pos[1]);
  xdr_double (&xdr_config, &vis->seed_pos[2]);
  xdr_int    (&xdr_config, &vis->res_factor);
  
  segFromMeshCoordsToSiteCoords (vis->seed_pos, vis->seed_site, vis);
#endif
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      xdr_int (&xdr_config, &vis->boundary[n].triangles);
      
      for (int m = 0; m < vis->boundary[n].triangles; m++)
	{
	  for (int k = 0; k < 3; k++)
	    for (int l = 0; l < 3; l++)
	      xdr_double (&xdr_config, &vis->boundary[n].triangle[m].v[k].pos[l]);
	  
	  editCalculateTriangleData (&vis->boundary[n].triangle[m]);
	  
	  if (n == INLET_BOUNDARY || n == OUTLET_BOUNDARY)
	    {
	      xdr_double (&xdr_config, &vis->boundary[n].triangle[m].pressure_avg);
	      xdr_double (&xdr_config, &vis->boundary[n].triangle[m].pressure_amp);
	      xdr_double (&xdr_config, &vis->boundary[n].triangle[m].pressure_phs);
	    }
	  if (n == INLET_BOUNDARY)
	    {
	      xdr_int (&xdr_config, &vis->boundary[n].triangle[m].normal_sign);
	    }
	}
    }
  xdr_destroy (&xdr_config);
}


void ioWriteCheckpoint (Vis *vis)
{
  FILE *system_config;
  XDR xdr_config;
  
  
  system_config = fopen (vis->checkpoint, "w");
  xdrstdio_create (&xdr_config, system_config, XDR_ENCODE);
#ifndef MESH
  xdr_short  (&xdr_config, &vis->selected_voxel[0]);
  xdr_short  (&xdr_config, &vis->selected_voxel[1]);
  xdr_short  (&xdr_config, &vis->selected_voxel[2]);
  xdr_double (&xdr_config, &vis->selected_grey);
  xdr_int    (&xdr_config, &vis->res_factor);
#else
  xdr_double (&xdr_config, &vis->seed_pos[0]);
  xdr_double (&xdr_config, &vis->seed_pos[1]);
  xdr_double (&xdr_config, &vis->seed_pos[2]);
  xdr_int    (&xdr_config, &vis->res_factor);
#endif
  for (unsigned int n = 0; n < BOUNDARIES; n++)
    {
      xdr_int (&xdr_config, &vis->boundary[n].triangles);
      
      for (int m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  for (int k = 0; k < 3; k++)
	    for (int l = 0; l < 3; l++)
	      xdr_double (&xdr_config, &vis->boundary[n].triangle[m].v[k].pos[l]);
	  
	  if (n == INLET_BOUNDARY || n == OUTLET_BOUNDARY)
	    {
	      xdr_double (&xdr_config, &vis->boundary[n].triangle[m].pressure_avg);
	      xdr_double (&xdr_config, &vis->boundary[n].triangle[m].pressure_amp);
	      xdr_double (&xdr_config, &vis->boundary[n].triangle[m].pressure_phs);
	    }
	  if (n == INLET_BOUNDARY)
	    {
	      xdr_int (&xdr_config, &vis->boundary[n].triangle[m].normal_sign);
	    }
	}
    }
  xdr_destroy (&xdr_config);
}


void ioWriteConfig (Vis *vis) {
  /* Write the config data for HemeLB.  Format (from reverse
   * engineering HemeLB code) is as follows.
   *
   * All values encoded using XDR format. Uses int, double and u_int.
   * 
   * System parameters:
   *   double stress_type
   *   int blocks_x
   *   int blocks_y
   *   int blocks_z
   *   int block_size 
   *
   * For each block (all blocks_x * blocks_y * blocks_z of them): 
   *
   *   int flag (indicates presence of non-solid sites in the block)
   *   
   *   If flag == 0 go to next block
   *
   *   Otherwise for each site in the block (all block_size^3):
   *
   *     u_int site_data -- this is a bit field which indicates site
   *     type (OR with SITE_TYPE_MASK to get bits zero and one; 00 =
   *     solid, 01 = fluid, 10 = inlet, 11 = outlet) or edgeness (set
   *     bit with PRESSURE_EDGE_MASK)
   * 
   *     If solid or simple fluid, go to next site
   *     
   *     If inlet or outlet (irrespective of edge state) {
   *       double boundary_normal[3]
   *       double boundary_dist
   *     }
   *
   *     If edge bit set {
   *       double wall_normal[3]
   *       double wall_dist
   *     }
   *     
   *     double mDistanceToWall[14]
   */
  FILE *system_config;
  XDR xdr_config;
  
#ifdef MESH
  double boundary_nor[3], boundary_dist;
  double wall_nor[3], wall_dist, cut_dist[14];
#endif
  double stress_type;
  
  int i, j, k;
  int m, n;
  int are_all_solid_sites;
  int flag;

  unsigned int *site_type;
  short int site[3];
  
  Block *block_p;
  
  
  system_config = fopen (vis->output_config, "w");
  xdrstdio_create (&xdr_config, system_config, XDR_ENCODE);
  
#ifndef MESH
  stress_type = VON_MISES_STRESS;
#else
  stress_type = SHEAR_STRESS;
#endif
  xdr_double (&xdr_config, &stress_type);
  xdr_int    (&xdr_config, &vis->blocks[0]);
  xdr_int    (&xdr_config, &vis->blocks[1]);
  xdr_int    (&xdr_config, &vis->blocks[2]);
  
  m = BLOCK_SIZE;
  xdr_int    (&xdr_config, &m);
  
  n = -1;
  
  for (i = 0; i < vis->blocks[0] * BLOCK_SIZE; i+=BLOCK_SIZE) {
    for (j = 0; j < vis->blocks[1] * BLOCK_SIZE; j+=BLOCK_SIZE) {
      for (k = 0; k < vis->blocks[2] * BLOCK_SIZE; k+=BLOCK_SIZE) {
	n++;
	
	/* flag == 1 indicates that the current macro block contains
	 * non-solid sites (i.e. ones where we have to something).
	 */
	block_p = &vis->block[n];
	if (block_p->site == NULL) {
	  /* If block pointer is null, there are no fluid sites.
	   */
	  flag = 0;
	} else {
	  /* If not, we must check the whole block.*/
	  are_all_solid_sites = 1;
	  for (m = 0; m < SITES_PER_BLOCK; m++) {
	    if (block_p->site[m].cfg != SOLID_TYPE) {
	      are_all_solid_sites = 0;
	      break;
	    }
	  }	      
	  if (are_all_solid_sites) {
	    flag = 0;
	  } else {
	    flag = 1;
	  }
	}
	
	xdr_int (&xdr_config, &flag);
	
	// No fluid site, skip.
	if (flag == 0)
	  continue;
	
	m = -1;
	for (site[0] = i; site[0] < i + BLOCK_SIZE; site[0]++) {
	  for (site[1] = j; site[1] < j + BLOCK_SIZE; site[1]++) {
	    for (site[2] = k; site[2] < k + BLOCK_SIZE; site[2]++){ 
	      m++;
	      
#ifdef MESH
	      site_type = &(block_p->site[m].cfg);
	      
	      if (*site_type == SOLID_TYPE) {
		// Site is solid, no further data needed
		xdr_u_int (&xdr_config, site_type);
		
		
	      } else if (*site_type == FLUID_TYPE) {
		// Site is fluid, no further data needed
		xdr_u_int (&xdr_config, site_type);
		
		
	      } else {
		// Site is complicated!
		segCalculateBoundarySiteData (*site_type, site,
					      boundary_nor, &boundary_dist,
					      wall_nor, &wall_dist, cut_dist, vis);
		
		// Note this alters site[m].cfg, which is reason for
		// deferring its output. (Set PRESSURE_EDGE_MASK
		// bit(s) to on)
		if (wall_dist < sqrt(3.0)) {
		  *site_type |= PRESSURE_EDGE_MASK;
		}
		xdr_u_int (&xdr_config, site_type);
		
		if (((*site_type & SITE_TYPE_MASK) == INLET_TYPE) ||
		    ((*site_type & SITE_TYPE_MASK) == OUTLET_TYPE)) {
		  for (int l = 0; l < 3; l++)
		    xdr_double (&xdr_config, &boundary_nor[l]);
		  
		  xdr_double (&xdr_config, &boundary_dist);
		}
		
		if ((*site_type & SITE_TYPE_MASK) == FLUID_TYPE ||
		    (*site_type & PRESSURE_EDGE_MASK)) {
		  for (int l = 0; l < 3; l++)
		    xdr_double (&xdr_config, &wall_nor[l]);
		  
		  xdr_double (&xdr_config, &wall_dist);
		}
		for (int l = 0; l < 14; l++)
		  xdr_double (&xdr_config, &cut_dist[l]);
	      }
#else // no MESH
	      xdr_u_int (&xdr_config, &block_p->site[m].cfg);
#endif // MESH
	    } // site[2]
	  }   // site[1]
	}     // site[0]
      }	// k
    }   // j
  }     // i
  xdr_destroy (&xdr_config);
}

void ioWritePars (Vis *vis)
{
  FILE *pars = fopen (vis->output_pars, "w");
  
  double nor[3], pos[3];
  
  int l, n;
  
  BoundaryTriangle *t_p;
  
  
  fprintf (pars, "%i\n", vis->boundary[ INLET_BOUNDARY ].triangles);
  
  for (n = 0; n < vis->boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      t_p = &vis->boundary[ INLET_BOUNDARY ].triangle[n];
      
      editTriangleNormal (t_p, nor);
      
      if (t_p->normal_sign == 1)
	{
	  for (l = 0; l < 3; l++)
	    nor[l] = -nor[l];
	}
      fprintf (pars, "%le %le %le\n", t_p->pressure_avg, t_p->pressure_amp, t_p->pressure_phs);
    }
  
  fprintf (pars, "%i\n", vis->boundary[ OUTLET_BOUNDARY ].triangles);
  
  for (n = 0; n < vis->boundary[ OUTLET_BOUNDARY ].triangles; n++)
    {
      t_p = &vis->boundary[ OUTLET_BOUNDARY ].triangle[n];
      
      editTriangleNormal (t_p, nor);
      
      fprintf (pars, "%le %le %le\n", t_p->pressure_avg, t_p->pressure_amp, t_p->pressure_phs);
    }
  for (n = 0; n < vis->boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      t_p = &vis->boundary[ INLET_BOUNDARY ].triangle[n];
      
      editTriangleNormal (t_p, nor);
      
      if (t_p->normal_sign == 1)
	{
	  for (l = 0; l < 3; l++)
	    nor[l] = -nor[l];
	}
      fprintf (pars, "%e %e %e\n", nor[0], nor[1], nor[2]);
    }
  for (n = 0; n < vis->boundary[ OUTLET_BOUNDARY ].triangles; n++)
    {
      t_p = &vis->boundary[ OUTLET_BOUNDARY ].triangle[n];
      
      editTriangleNormal (t_p, nor);
      
      fprintf (pars, "%e %e %e\n", nor[0], nor[1], nor[2]);
    }
 for (n = 0; n < vis->boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      t_p = &vis->boundary[ INLET_BOUNDARY ].triangle[n];
      
      editTriangleCenter (t_p, pos);
      
      fprintf (pars, "%e %e %e\n", pos[0], pos[1], pos[2]);
    }
  for (n = 0; n < vis->boundary[ OUTLET_BOUNDARY ].triangles; n++)
    {
      t_p = &vis->boundary[ OUTLET_BOUNDARY ].triangle[n];
      
      editTriangleCenter (t_p, pos);
      
      fprintf (pars, "%e %e %e\n", pos[0], pos[1], pos[2]);
    }
  fclose (pars);
}


void ioWriteCoords(Vis *vis) {
  // Write the required information to translate between output
  // coordinate system and the input STL file.
  FILE *coords = fopen (vis->output_coords, "w");
  if (coords == NULL) {
    printf("Cannot open output file for coordinates: '%s'\n", vis->output_coords);
    return;
  }
  
  fprintf(coords, "voxel_size = %e\n",
	  vis->mesh.voxel_size);
  fprintf(coords, "res_factor = %d\n",
	  vis->res_factor);
  fprintf(coords, "voxels =  %d, %d, %d\n",
	  vis->mesh.voxels[0], vis->mesh.voxels[1], vis->mesh.voxels[2]);
  fprintf(coords, "dim =  %e, %e, %e\n", 
	  vis->mesh.dim[0], vis->mesh.dim[1], vis->mesh.dim[2]);
  fprintf(coords, "half_dim = %e, %e, %e\n", 
	  vis->mesh.half_dim[0], vis->mesh.half_dim[1], vis->mesh.half_dim[2]);
  fprintf(coords, "mesh_centre = %e, %e, %e\n",
	  vis->mesh.centre[0], vis->mesh.centre[1], vis->mesh.centre[2]);
  fprintf(coords, "seed_pos = %e, %e, %e\n",
	  vis->seed_pos[0], vis->seed_pos[1], vis->seed_pos[2]);
  fprintf(coords, "seed_site = %d, %d, %d\n",
	  vis->seed_site[0], vis->seed_site[1], vis->seed_site[2]);

  fclose(coords);
}

#ifdef USE_TIFFLIB

void ioSaveWindowImage (char *file_name) {
	
	TIFF *tif = TIFFOpen("snapshot.tif", "w");
	
	int pix_x, pix_y;
	int i, j;
	
	char *image_data = NULL;  
	char *image_data_p = NULL;
	char *row_data;
	
	glReadBuffer (GL_FRONT);
	
	pix_x = mScreen.pixels[0];
	pix_y = mScreen.pixels[1];
	
	image_data = new char [pix_x * pix_y * 3];
	
	row_data = new char[pix_x * 3];
	
	image_data_p = image_data;
	
	for (j = 0; j < pix_y; j++) {
		glReadPixels (0, j, pix_x, 1, GL_RGB, GL_UNSIGNED_BYTE, row_data);
		for (i = pix_x-1; i >=0; i--) {
			*image_data_p++ = row_data[ i*3+2 ];
			*image_data_p++ = row_data[ i*3+1 ];
			*image_data_p++ = row_data[ i*3 ];
		}
    }
	
	delete[] row_data;
	
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, pix_x);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, pix_y);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_BOTRIGHT);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);	
	
	TIFFWriteEncodedStrip(tif, 0, image_data, pix_x*pix_y*3);

	delete[] image_data;
	
	TIFFClose(tif);
}

#else // USE_TIFFLIB

void ioSaveWindowImage (char *file_name)
{
  FILE *ppm_image_file = fopen (file_name, "wb");
  
  int pix_x, pix_y;
  int i, j;
  
  unsigned char *image_data = NULL;  
  unsigned char *image_data_p = NULL;
  unsigned char *row_data;
  
  glReadBuffer (GL_FRONT);
  
  pix_x = screen.pixels[0];
  pix_y = screen.pixels[1];
  
  image_data = new unsigned char[pix_x * pix_y * 3];
  
  row_data = new unsigned char[pix_x * 3];
  
  image_data_p = image_data;

  for (j = 0; j < pix_y; j++)
    {
      glReadPixels (0, j, pix_x, 1, GL_RGB, GL_UNSIGNED_BYTE, row_data);
      
      for (i = 0; i < pix_x; i++)
	{
	  *image_data_p++ = row_data[ i*3   ];
	  *image_data_p++ = row_data[ i*3+1 ];
	  *image_data_p++ = row_data[ i*3+2 ];
	}
    }
  delete[] row_data;
  
  fprintf (ppm_image_file, "P6\n%i %i\n255\n", pix_x, pix_y);
  
  for (j = pix_y - 1; j >= 0; j--)
    {
      fwrite (image_data + j * pix_x * 3, 1, pix_x * 3, ppm_image_file);
    }
  
  delete[] image_data;
  
  fclose (ppm_image_file);
}

#endif // USE_TIFFLIB



