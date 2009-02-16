#include "io_dicom.h"

#include "io.h"

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
*/

/*
void ioGetFileNames (Vis *vis)
{
  string dir = string(vis->input_path);
  
  ioGetDir (dir, vis->file_list);
}
*/


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
	  
	  vis->grey_min = fminf(vis->grey_min, temp);
	  vis->grey_max = fmaxf(vis->grey_max, temp);
	}
    }

  input_file.close();

}


void ioReadConfig (Vis *vis)
{

  int i, l;
  
  //ioGetFileNames (vis);
  
  //vis->input_voxels[2] = vis->file_list.size() - 2;

  DICOM_read DICOM_input = DICOM_read(vis->input_path);

  DICOM_input.DICOM_volume_size(vis->input_voxels[0], vis->input_voxels[1], vis->input_voxels[2]);

  vis->output_voxels[0] = vis->input_voxels[0] * vis->res_factor;
  vis->output_voxels[1] = vis->input_voxels[1] * vis->res_factor;
  vis->output_voxels[2] = (int)(vis->input_voxels[2] * (vis->slice_size/vis->pixel_size) *
				(float)vis->output_voxels[0]/(float)vis->input_voxels[0]);
  
  vis->scale[0] = vis->res_factor;
  vis->scale[1] = vis->res_factor;
  vis->scale[2] = vis->res_factor * vis->slice_size / vis->pixel_size;
  
  for (l = 0; l < 3; l++)
    vis->inv_scale[l] = 1.F / vis->scale[l];
  
  for (l = 0; l < 3; l++)
    vis->blocks[l] = vis->output_voxels[l] >> SHIFT;
  
  for (l = 0; l < 3; l++)
    if ((vis->blocks[l] << SHIFT) < vis->output_voxels[l]) ++vis->blocks[l];
  
  for (l = 0; l < 3; l++)
    vis->sites[l] = vis->blocks[l] * BLOCK_SIZE;
  
  for (l = 0; l < 3; l++)
    vis->dim[l] = vis->output_voxels[l];
  
  for (l = 0; l < 3; l++)
    vis->half_dim[l] = 0.5F * vis->dim[l];
  
  vis->system_size = fmaxf(vis->dim[0], fmaxf(vis->dim[1], vis->dim[2]));
  
  vis->tot_blocks = vis->blocks[0] * vis->blocks[1] * vis->blocks[2];
  
  //vis->voxel = (unsigned short int *)malloc(sizeof(unsigned short int) *

  vis->voxel = (signed short int *)malloc(sizeof(signed short int) *
					    vis->input_voxels[0] *
					    vis->input_voxels[1] *
					    vis->input_voxels[2]);
  
  vis->block = (Block *)malloc(sizeof(Block) * vis->tot_blocks);
  
  vis->stack_sites_max = SITES_PER_BLOCK * 100000;
  
  vis->stack_site = (Site *)malloc(sizeof(Site) * vis->stack_sites_max);
  
  vis->grey_max = -1.0e9;
  vis->grey_min = +1.0e9;
  
/*  for (i = 0; i < vis->input_voxels[2]; i++)
    {
      ioReadSlice (i, vis);
    } */

  DICOM_input.DICOM_copy_data(vis->voxel, vis->input_voxels[0] * vis->input_voxels[1] * vis->input_voxels[2]);

  for(l = 0; l<vis->input_voxels[0] * vis->input_voxels[1] * vis->input_voxels[2]; l++) {
    if( vis->voxel[l] > vis->grey_max ) vis->grey_max = vis->voxel[l];
    if( vis->voxel[l] < vis->grey_min ) vis->grey_min = vis->voxel[l];
  }

  printf("grey min %0.2f\n", vis->grey_min);
  printf("grey max %0.2f\n", vis->grey_max);

  vis->grey_min = -2000.0;
  vis->grey_max = 6000.0;


  vis->coord[A] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_A);
  vis->coord[B] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_B);
  vis->coord[C] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_C);

}


void ioWriteCheckpoint (Vis *vis)
{
  FILE *system_config;
  XDR xdr_config;
  
  int k, l, m, n;
  
  unsigned int data;
  
  
  system_config = fopen (vis->checkpoint, "w");
  xdrstdio_create (&xdr_config, system_config, XDR_ENCODE);
  
  xdr_float  (&xdr_config, &vis->slice_size);
  xdr_float  (&xdr_config, &vis->pixel_size);
  
  float res_factor = vis->res_factor;
  
  xdr_float  (&xdr_config, &res_factor);
  
  xdr_int    (&xdr_config, &vis->input_voxels[0]);
  xdr_int    (&xdr_config, &vis->input_voxels[1]);
  xdr_int    (&xdr_config, &vis->input_voxels[2]);
  
  xdr_short  (&xdr_config, &vis->selected_voxel[0]);
  xdr_short  (&xdr_config, &vis->selected_voxel[1]);
  xdr_short  (&xdr_config, &vis->selected_voxel[2]);
  xdr_float  (&xdr_config, &vis->selected_grey);
  
  xdr_int    (&xdr_config, &vis->blocks[0]);
  xdr_int    (&xdr_config, &vis->blocks[1]);
  xdr_int    (&xdr_config, &vis->blocks[2]);
  
  for (n = 0; n < vis->input_voxels[0]*vis->input_voxels[1]*vis->input_voxels[2]; n++)
    {
      data = vis->voxel[n];
      
      xdr_u_int (&xdr_config, &data);
    }
  
  for (n = 0; n < BOUNDARIES; n++)
    {
      xdr_int (&xdr_config, &vis->boundary[n].triangles);
      
      for (m = 0; m < vis->boundary[ n ].triangles; m++)
	{
	  for (k = 0; k < 3; k++)
	    for (l = 0; l < 3; l++)
	      xdr_float (&xdr_config, &vis->boundary[n].triangle[m].v[k].pos[l]);
	  
	  if (n == INLET_BOUNDARY || n == OUTLET_BOUNDARY)
	    {
	      xdr_float (&xdr_config, &vis->boundary[n].triangle[m].pressure_avg);
	      xdr_float (&xdr_config, &vis->boundary[n].triangle[m].pressure_amp);
	      xdr_float (&xdr_config, &vis->boundary[n].triangle[m].pressure_phs);
	    }
	}
    }
  xdr_destroy (&xdr_config);
}


void ioReadCheckpoint (Vis *vis)
{
  FILE *system_config;
  XDR xdr_config;
  
  float res_factor;
  
  int k, l, m, n;
  
  unsigned int data;
  
  
  system_config = fopen (vis->checkpoint, "r");
  xdrstdio_create (&xdr_config, system_config, XDR_DECODE);
  
  xdr_float  (&xdr_config, &vis->slice_size);
  xdr_float  (&xdr_config, &vis->pixel_size);
  xdr_float  (&xdr_config, &res_factor);
  
  vis->res_factor = res_factor;
  
  xdr_int    (&xdr_config, &vis->input_voxels[0]);
  xdr_int    (&xdr_config, &vis->input_voxels[1]);
  xdr_int    (&xdr_config, &vis->input_voxels[2]);
  
  xdr_short  (&xdr_config, &vis->selected_voxel[0]);
  xdr_short  (&xdr_config, &vis->selected_voxel[1]);
  xdr_short  (&xdr_config, &vis->selected_voxel[2]);
  xdr_float  (&xdr_config, &vis->selected_grey);
  
  xdr_int    (&xdr_config, &vis->blocks[0]);
  xdr_int    (&xdr_config, &vis->blocks[1]);
  xdr_int    (&xdr_config, &vis->blocks[2]);
  
  vis->output_voxels[0] = vis->input_voxels[0] * vis->res_factor;
  vis->output_voxels[1] = vis->input_voxels[1] * vis->res_factor;
  vis->output_voxels[2] = (int)(vis->input_voxels[2] * (vis->slice_size/vis->pixel_size) *
				(float)vis->output_voxels[0]/(float)vis->input_voxels[0]);
  
  for (l = 0; l < 3; l++)
    vis->scale[l] = (float)vis->output_voxels[l] / (float)vis->input_voxels[l];
  
  for (l = 0; l < 3; l++)
    vis->inv_scale[l] = 1.F / vis->scale[l];
  
  for (l = 0; l < 3; l++)
    vis->sites[l] = vis->blocks[l] * BLOCK_SIZE;
  
  for (l = 0; l < 3; l++)
    vis->dim[l] = vis->output_voxels[l];
  
  for (l = 0; l < 3; l++)
    vis->half_dim[l] = 0.5F * vis->dim[l];
  
  vis->system_size = fmaxf(vis->dim[0], fmaxf(vis->dim[1], vis->dim[2]));
  
  vis->tot_blocks = vis->blocks[0] * vis->blocks[1] * vis->blocks[2];
  
//  vis->voxel = (unsigned short int *)malloc(sizeof(unsigned short int) *
//					    vis->input_voxels[0] * vis->input_voxels[1] * vis->input_voxels[2]);

  vis->voxel = (signed short int *)malloc(sizeof(signed short int) *
					    vis->input_voxels[0] * vis->input_voxels[1] * vis->input_voxels[2]);
  
  vis->grey_max = -1.0e9;
  vis->grey_min = +1.0e9;
  
  for (n = 0; n < vis->input_voxels[0] * vis->input_voxels[1] * vis->input_voxels[2]; n++)
    {
      xdr_u_int (&xdr_config, &data);
      
      vis->voxel[n] = data;
      
      vis->grey_min = fminf(vis->grey_min, data);
      vis->grey_max = fmaxf(vis->grey_max, data);
    }
  
  vis->block = (Block *)malloc(sizeof(Block) * vis->tot_blocks);
  
  vis->stack_sites_max = SITES_PER_BLOCK * 100000;
  
  vis->stack_site = (Site *)malloc(sizeof(Site) * vis->stack_sites_max);
  
  vis->screen_voxel = (ScreenVoxel *)malloc(sizeof(ScreenVoxel) * vis->screen_voxels * vis->screen_voxels);
  
  vis->coord[A] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_A);
  vis->coord[B] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_B);
  vis->coord[C] = (Coord *)malloc(sizeof(Coord) * COORD_BUFFER_SIZE_C);
  
  vis->boundary[ INLET_BOUNDARY  ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U<<BOUNDARY_ID_BITS));
  vis->boundary[ OUTLET_BOUNDARY ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U<<BOUNDARY_ID_BITS));
  vis->boundary[ WALL_BOUNDARY   ].triangle = (Triangle *)malloc(sizeof(Triangle) * (1U<<BOUNDARY_ID_BITS));
  
  for (n = 0; n < BOUNDARIES; n++)
    {
      xdr_int (&xdr_config, &vis->boundary[n].triangles);
      
      for (m = 0; m < vis->boundary[n].triangles; m++)
	{
	  for (k = 0; k < 3; k++)
	    for (l = 0; l < 3; l++)
	      xdr_float (&xdr_config, &vis->boundary[n].triangle[m].v[k].pos[l]);
	  
	  editCalculateTriangleData (&vis->boundary[n].triangle[m]);
	  
	  if (n == INLET_BOUNDARY || n == OUTLET_BOUNDARY)
	    {
	      xdr_float (&xdr_config, &vis->boundary[n].triangle[m].pressure_avg);
	      xdr_float (&xdr_config, &vis->boundary[n].triangle[m].pressure_amp);
	      xdr_float (&xdr_config, &vis->boundary[n].triangle[m].pressure_phs);
	    }
	}
    }
  xdr_destroy (&xdr_config);
}


void ioWriteConfig (Vis *vis)
{
  FILE *system_config;
  XDR xdr_config;

  int i, j, k;
  int m, n;
  int are_all_solid_sites;
  int flag;
  
  Block *block_p;
  
  
  system_config = fopen (vis->output_config, "w");
  xdrstdio_create (&xdr_config, system_config, XDR_ENCODE);
  
  double dummy = 1.0;
  
  xdr_double (&xdr_config, &dummy);
  xdr_int    (&xdr_config, &vis->blocks[0]);
  xdr_int    (&xdr_config, &vis->blocks[1]);
  xdr_int    (&xdr_config, &vis->blocks[2]);
  
  m = BLOCK_SIZE;
  xdr_int    (&xdr_config, &m);
  
  n = -1;
  
  for (i = 0; i < vis->blocks[0]; i++)
    for (j = 0; j < vis->blocks[1]; j++)
      for (k = 0; k < vis->blocks[2]; k++)
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
	  
	  for (m = 0; m < SITES_PER_BLOCK; m++)
	    {
	      xdr_u_int (&xdr_config, &block_p->site[m].cfg);
	    }
	}
  xdr_destroy (&xdr_config);
}


void ioWritePars (Vis *vis)
{
  FILE *pars = fopen (vis->output_pars, "w");
  
  float nor[3];
  
  int l, n;
  
  Triangle *t_p;
  
  
  fprintf (pars, "%i\n", vis->boundary[ INLET_BOUNDARY ].triangles);
  
  for (n = 0; n < vis->boundary[ INLET_BOUNDARY ].triangles; n++)
    {
      fprintf (pars, "%f %f %f\n",
	       vis->boundary[ INLET_BOUNDARY ].triangle[n].pressure_avg,
	       vis->boundary[ INLET_BOUNDARY ].triangle[n].pressure_amp,
	       vis->boundary[ INLET_BOUNDARY ].triangle[n].pressure_phs);
    }
  
  fprintf (pars, "%i\n", vis->boundary[ OUTLET_BOUNDARY ].triangles);
  
  for (n = 0; n < vis->boundary[ OUTLET_BOUNDARY ].triangles; n++)
    {
      fprintf (pars, "%f %f %f\n",
	       vis->boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_avg,
	       vis->boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_amp,
	       vis->boundary[ OUTLET_BOUNDARY ].triangle[n].pressure_phs);
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
      fprintf (pars, "%f %f %f\n", nor[0], nor[1], nor[2]);
    }
  fclose (pars);
}


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
