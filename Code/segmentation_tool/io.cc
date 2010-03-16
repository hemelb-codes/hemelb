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
  
  vis->voxel = (short int *)malloc(sizeof(short int) *
				   vis->input_voxels[0] *
				   vis->input_voxels[1] *
				   vis->input_voxels[2]);
  
  vis->block = (Block *)malloc(sizeof(Block) * vis->tot_blocks);
  
  vis->stack_sites_max = SITES_PER_BLOCK * 100000;
  
  vis->stack_site = (Site *)malloc(sizeof(Site) * vis->stack_sites_max);
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
  vis->coord[A] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_A);
  vis->coord[B] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_B);
  vis->coord[C] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_C);
  
  vis->stress_type = VON_MISES_STRESS;
}
#else // MESH


void ioReadConfig (Vis *vis)
{
  double min_x[3], max_x[3];
  double ctr_x[3];
  double dx[3];
  
  int l, m, n;
  
  char key_word[16];
  
  ifstream input_file;
  
  
  vis->mesh.triangles_max = 10000;
  vis->mesh.triangle = (MeshTriangle *)malloc(sizeof(MeshTriangle) * vis->mesh.triangles_max);
  vis->mesh.triangles = 0;
  
  for (l = 0; l < 3; l++)
    {
      min_x[l] =  1.0e+30;
      max_x[l] = -1.0e+30;
      ctr_x[l] =  0.0;
    }
  input_file.open(vis->input_file);
  
  while (strcmp(key_word, "facet"))
    input_file >> key_word;
  
  while (strcmp(key_word, "endsolid"))
    {
      if (vis->mesh.triangles == vis->mesh.triangles_max)
	{
	  vis->mesh.triangles_max *= 2;
	  vis->mesh.triangle = (MeshTriangle *)realloc(vis->mesh.triangle,
						       sizeof(MeshTriangle) * vis->mesh.triangles_max);
	}
      while (strcmp(key_word, "normal"))
	input_file >> key_word;
      
      input_file >> vis->mesh.triangle[ vis->mesh.triangles ].nor[0];
      input_file >> vis->mesh.triangle[ vis->mesh.triangles ].nor[1];
      input_file >> vis->mesh.triangle[ vis->mesh.triangles ].nor[2];
      
      for (l = 0; l < 3; l++)
	vis->mesh.triangle[ vis->mesh.triangles ].nor[l] /=
	  sqrt(ScalarProd (vis->mesh.triangle[ vis->mesh.triangles ].nor, vis->mesh.triangle[ vis->mesh.triangles ].nor));
      
      while (strcmp(key_word, "loop"))
	input_file >> key_word;
      
      for (m = 0; m < 3; m++)
	{
	  input_file >> key_word;
	  input_file >> vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[0];
	  input_file >> vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[1];
	  input_file >> vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[2];
	}
      float triangle_ctr_x[3];
      
      for (l = 0; l < 3; l++)
	triangle_ctr_x[l] = 0.0;
      
      for (m = 0; m < 3; m++)
	for (l = 0; l < 3; l++)
	  {
	    triangle_ctr_x[l] += vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l];
	  }
      for (l = 0; l < 3; l++)
	triangle_ctr_x[l] /= 3.0;
      
      for (m = 0; m < 3; m++)
	for (l = 0; l < 3; l++)
	  {
	    vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l] = triangle_ctr_x[l] +
	      1.01 * (vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l] - triangle_ctr_x[l]);
	  }
      for (m = 0; m < 3; m++)
	{
	  for (l = 0; l < 3; l++)
	    {
	      min_x[l] = fmin(min_x[l], vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l]);
	      max_x[l] = fmax(max_x[l], vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l]);
	      ctr_x[l] += vis->mesh.triangle[ vis->mesh.triangles ].v[m].pos[l];
	    }
	}
      ++vis->mesh.triangles;
      
      input_file >> key_word;
      input_file >> key_word;
      input_file >> key_word;
    }
  input_file.close();
  
  for (l = 0; l < 3; l++)
    ctr_x[l] /= vis->mesh.triangles * 3;
  
  vis->mesh.voxel_size = 0.0;
  
  for (n = 0; n < vis->mesh.triangles; n++)
    {
      for (m = 0; m < 3; m++)
	for (l = 0; l < 3; l++)
	  vis->mesh.triangle[n].v[m].pos[l] -= ctr_x[l];
      
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
  vis->mesh.voxel_size /= (vis->mesh.triangles * 3);
  
  for (l = 0; l < 3; l++)
    vis->mesh.voxels[l] = (int)ceil(1.5 * (max_x[l] - min_x[l]) / vis->mesh.voxel_size);
  
  vis->mesh.voxels[3] = vis->mesh.voxels[0] * vis->mesh.voxels[1] * vis->mesh.voxels[2];
  
  vis->mesh.voxel = (Voxel *)malloc(sizeof(Voxel) * vis->mesh.voxels[3]);
  
  for (l = 0; l < 3; l++)
    vis->mesh.dim[l] = vis->mesh.voxel_size * vis->mesh.voxels[l];
  
  for (l = 0; l < 3; l++)
    vis->mesh.half_dim[l] = 0.5 * vis->mesh.dim[l];
      
  for (l = 0; l < 3; l++)
    vis->input_voxels[l] = vis->mesh.voxels[l];
  
  for (l = 0; l < 3; l++)
    vis->output_voxels[l] = vis->input_voxels[l] * vis->res_factor;
  
  for (l = 0; l < 3; l++)
    vis->blocks[l] = vis->output_voxels[l] >> SHIFT;
  
  for (l = 0; l < 3; l++)
    if ((vis->blocks[l] << SHIFT) < vis->output_voxels[l]) ++vis->blocks[l];
  
  for (l = 0; l < 3; l++)
    vis->sites[l] = vis->blocks[l] * BLOCK_SIZE;
  
  for (l = 0; l < 3; l++)
    vis->dim[l] = vis->mesh.dim[l];
  
  for (l = 0; l < 3; l++)
    vis->half_dim[l] = 0.5 * vis->dim[l];
  
  vis->system_size = fmax(vis->dim[0], fmax(vis->dim[1], vis->dim[2]));
  
  vis->tot_blocks = vis->blocks[0] * vis->blocks[1] * vis->blocks[2];
  
  vis->block = (Block *)malloc(sizeof(Block) * vis->tot_blocks);
  
  for (n = 0; n < vis->tot_blocks; n++)
    vis->block[n].site = NULL;
  
  vis->stack_sites_max = SITES_PER_BLOCK * 100000;
  
  vis->stack_site = (Site *)malloc(sizeof(Site) * vis->stack_sites_max);
  
  vis->coord[A] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_A);
  vis->coord[B] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_B);
  vis->coord[C] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_C);
  
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


void ioWriteConfig (Vis *vis)
{
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
  
  for (i = 0; i < vis->blocks[0] * BLOCK_SIZE; i+=BLOCK_SIZE)
    for (j = 0; j < vis->blocks[1] * BLOCK_SIZE; j+=BLOCK_SIZE)
      for (k = 0; k < vis->blocks[2] * BLOCK_SIZE; k+=BLOCK_SIZE)
	{
	  block_p = &vis->block[++n];
	  
	  flag = 0;
	  
	  if (block_p->site == NULL)
	    {
	      xdr_int (&xdr_config, &flag);
	      continue;
	    }
	  are_all_solid_sites = 1;
	  
	  for (m = 0; m < SITES_PER_BLOCK; m++)
	    {
	      if (block_p->site[m].cfg != SOLID_TYPE)
		{
		  are_all_solid_sites = 0;
		  break;
		}
	    }	      
	  if (are_all_solid_sites)
	    {
	      xdr_int (&xdr_config, &flag);
	      continue;
	    }
	  flag = 1;
	  xdr_int (&xdr_config, &flag);
	  
	  m = -1;
	  
	  for (site[0] = i; site[0] < i + BLOCK_SIZE; site[0]++)
	    for (site[1] = j; site[1] < j + BLOCK_SIZE; site[1]++)
	      for (site[2] = k; site[2] < k + BLOCK_SIZE; site[2]++)
		{
#ifndef MESH
		  xdr_u_int (&xdr_config, &block_p->site[++m].cfg);
#else
		  if (block_p->site[++m].cfg == SOLID_TYPE)
		    {
                      continue;
                    }
                  if (block_p->site[m].cfg != FLUID_TYPE)
		    {
		      segCalculateBoundarySiteData (block_p->site[m].cfg, site,
						    boundary_nor, &boundary_dist,
						    wall_nor, &wall_dist, cut_dist, vis);
		      
		      if (wall_dist < sqrt(3.0))
			{
			  block_p->site[m].cfg |= PRESSURE_EDGE_MASK;
			}
		    }
		  xdr_u_int (&xdr_config, &block_p->site[m].cfg);
		  
		  if (block_p->site[m].cfg == FLUID_TYPE)
		    {
                      continue;
                    }
		  if (((block_p->site[m].cfg & SITE_TYPE_MASK) & INLET_TYPE) ||
                      ((block_p->site[m].cfg & SITE_TYPE_MASK) & OUTLET_TYPE))
	            {
	              for (int l = 0; l < 3; l++)
			xdr_double (&xdr_config, &boundary_nor[l]);
	              
	              xdr_double (&xdr_config, &boundary_dist);
	            }
		  if ((block_p->site[m].cfg & SITE_TYPE_MASK) == FLUID_TYPE ||
                      (block_p->site[m].cfg & PRESSURE_EDGE_MASK))
	            {
	              for (int l = 0; l < 3; l++)
			xdr_double (&xdr_config, &wall_nor[l]);
		      
		      xdr_double (&xdr_config, &wall_dist);
		    }
		  for (int l = 0; l < 14; l++)
	            xdr_double (&xdr_config, &cut_dist[l]);
#endif // MESH
		}
	}
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


#ifdef USE_TIFFLIB

void ioSaveWindowImage (char *file_name) {
	
	TIFF *tif = TIFFOpen("snapshot.tif", "w");
	
	int pix_x, pix_y;
	int i, j;
	
	char *image_data = NULL;  
	char *image_data_p = NULL;
	char *row_data;
	
	glReadBuffer (GL_FRONT);
	
	pix_x = screen.pixels[0];
	pix_y = screen.pixels[1];
	
	image_data = (char *)malloc(sizeof(char) * pix_x * pix_y * 3);
	
	row_data = (char *)malloc(sizeof(char) * pix_x * 3);
	
	image_data_p = image_data;
	
	for (j = 0; j < pix_y; j++) {
		glReadPixels (0, j, pix_x, 1, GL_RGB, GL_UNSIGNED_BYTE, row_data);
		for (i = pix_x-1; i >=0; i--) {
			*image_data_p++ = row_data[ i*3+2 ];
			*image_data_p++ = row_data[ i*3+1 ];
			*image_data_p++ = row_data[ i*3 ];
		}
    }
	
	free((char *)row_data);
	
	TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, pix_x);
	TIFFSetField(tif, TIFFTAG_IMAGELENGTH, pix_y);
	TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
	TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_BOTRIGHT);
	TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);	
	
	TIFFWriteEncodedStrip(tif, 0, image_data, pix_x*pix_y*3);

	free((char *)image_data);
	
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
  
  image_data = (unsigned char *)malloc(sizeof(unsigned char) * pix_x * pix_y * 3);
  
  row_data = (unsigned char *)malloc(sizeof(unsigned char) * pix_x * 3);
  
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
  free((unsigned char *)row_data);
  
  fprintf (ppm_image_file, "P6\n%i %i\n255\n", pix_x, pix_y);
  
  for (j = pix_y - 1; j >= 0; j--)
    {
      fwrite (image_data + j * pix_x * 3, 1, pix_x * 3, ppm_image_file);
    }
  
  free((unsigned char *)image_data);
  
  fclose (ppm_image_file);
}

#endif // USE_TIFFLIB



