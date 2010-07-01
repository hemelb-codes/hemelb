#include "utilityFunctions.h"
#include "rt.h"
#include "lb.h"
#include "xdrFileWriter.h"
#include <math.h>
#include <string.h>

glyphDrawer *myGlypher;

void visRotate (float sin_1, float cos_1,
		float sin_2, float cos_2,
		float  x1, float  y1, float  z1,
		float *x2, float *y2, float *z2)
{
  //sin_1 = sin(longitude), cos_1 = cos(longitude)
  //sin_2 = sin(latitude) , cos_2 = cos(latitude)
  //rotation of latitude  DEG around the X axis followed by a
  //rotation of longitude DEG around the Y axis
  
  float temp = z1 * cos_2 - y1 * sin_2;
  
  *x2 = temp * sin_1 + x1 * cos_1;
  *y2 =   z1 * sin_2 + y1 * cos_2;
  *z2 = temp * cos_1 - x1 * sin_1;
}


void visProject (float p1[], float p2[])
{
  float x1[3], x2[3];
  float temp;
  
  int l;
  
  
  for (l = 0; l < 3; l++)
    {
      x1[l] = p1[l] - viewpoint.x[l];
    }
  temp = viewpoint.cos_1 * x1[2] + viewpoint.sin_1 * x1[0];
  
  x2[0] = viewpoint.cos_1 * x1[0] - viewpoint.sin_1 * x1[2];
  x2[1] = viewpoint.cos_2 * x1[1] - viewpoint.sin_2 * temp;
  x2[2] = viewpoint.cos_2 * temp + viewpoint.sin_2 * x1[1];
  
  temp = viewpoint.dist / (p2[2] = -x2[2]);
  
  p2[0] = temp * x2[0];
  p2[1] = temp * x2[1];
}


void visProjection (float ortho_x, float ortho_y,
		    int pixels_x, int pixels_y,
		    float ctr_x, float ctr_y, float ctr_z,
		    float rad,
		    float longitude, float latitude,
		    float dist,
		    float zoom)
{
  float temp;
  
  
  screen.max_x = ortho_x / zoom;
  screen.max_y = ortho_y / zoom;
  
  screen.pixels_x = pixels_x;
  screen.pixels_y = pixels_y;
  
  temp = longitude * 0.01745329F;
  
  viewpoint.sin_1 = sinf(temp);
  viewpoint.cos_1 = cosf(temp);
  
  temp = latitude * 0.01745329F;
  
  viewpoint.sin_2 = sinf(temp);
  viewpoint.cos_2 = cosf(temp);
  
  temp = rad * viewpoint.cos_2;
  
  viewpoint.x[0] = temp * viewpoint.sin_1 + ctr_x;
  viewpoint.x[1] = rad  * viewpoint.sin_2 + ctr_y;
  viewpoint.x[2] = temp * viewpoint.cos_1 + ctr_z;
  
  viewpoint.dist = dist;
  
  temp = dist / rad;
  
  ctr_x = viewpoint.x[0] + temp * (ctr_x - viewpoint.x[0]);
  ctr_y = viewpoint.x[1] + temp * (ctr_y - viewpoint.x[1]);
  ctr_z = viewpoint.x[2] + temp * (ctr_z - viewpoint.x[2]);
  
  screen.zoom = zoom;
  
  visRotate (viewpoint.sin_1, viewpoint.cos_1,
	     viewpoint.sin_2, viewpoint.cos_2,
	     screen.max_x, 0.0F, 0.0F,
	     &screen.dir1[0], &screen.dir1[1], &screen.dir1[2]);
  
  visRotate (viewpoint.sin_1, viewpoint.cos_1,
	     viewpoint.sin_2, viewpoint.cos_2,
	     0.0F, screen.max_y, 0.0F,
	     &screen.dir2[0], &screen.dir2[1], &screen.dir2[2]);
  
  screen.scale_x = (float)pixels_x / (2.F * screen.max_x);
  screen.scale_y = (float)pixels_y / (2.F * screen.max_y);
  
  screen.vtx[0] = ctr_x - screen.dir1[0] - screen.dir2[0] - viewpoint.x[0];
  screen.vtx[1] = ctr_y - screen.dir1[1] - screen.dir2[1] - viewpoint.x[1];
  screen.vtx[2] = ctr_z - screen.dir1[2] - screen.dir2[2] - viewpoint.x[2];
  
  screen.dir1[0] *= (2.F / (float)pixels_x);
  screen.dir1[1] *= (2.F / (float)pixels_x);
  screen.dir1[2] *= (2.F / (float)pixels_x);
  
  screen.dir2[0] *= (2.F / (float)pixels_y);
  screen.dir2[1] *= (2.F / (float)pixels_y);
  screen.dir2[2] *= (2.F / (float)pixels_y);
}


void visMergePixels (ColPixel *col_pixel1, ColPixel *col_pixel2)
{
  // merge raytracing data
  if ((col_pixel1->i & RT) && (col_pixel2->i & RT))
    {
      col_pixel2->vel_r += col_pixel1->vel_r;
      col_pixel2->vel_g += col_pixel1->vel_g;
      col_pixel2->vel_b += col_pixel1->vel_b;
      
      if (lbm_stress_type != SHEAR_STRESS)
	{
	  col_pixel2->stress_r += col_pixel1->stress_r;
	  col_pixel2->stress_g += col_pixel1->stress_g;
	  col_pixel2->stress_b += col_pixel1->stress_b;
	}
      col_pixel2->dt += col_pixel1->dt;
      
      if (col_pixel1->t < col_pixel2->t)
	{
	  col_pixel2->t       = col_pixel1->t;
	  col_pixel2->density = col_pixel1->density;
	  col_pixel2->stress  = col_pixel1->stress;
	}
    }
  else if ((col_pixel1->i & RT) && !(col_pixel2->i & RT))
    {
      col_pixel2->vel_r = col_pixel1->vel_r;
      col_pixel2->vel_g = col_pixel1->vel_g;
      col_pixel2->vel_b = col_pixel1->vel_b;
      
      if (lbm_stress_type != SHEAR_STRESS)
	{
	  col_pixel2->stress_r = col_pixel1->stress_r;
	  col_pixel2->stress_g = col_pixel1->stress_g;
	  col_pixel2->stress_b = col_pixel1->stress_b;
	}
      col_pixel2->t       = col_pixel1->t;
      col_pixel2->dt      = col_pixel1->dt;
      col_pixel2->density = col_pixel1->density;
      col_pixel2->stress  = col_pixel1->stress;
      
      col_pixel2->i |= RT;
    }
  if (lbm_stress_type != SHEAR_STRESS && (vis_mode == 0 || vis_mode == 1))
    {
      // merge glyph data
      if (col_pixel1->i & GLYPH)
	{
	  col_pixel2->i |= GLYPH;
	}
    }
  else
    {
#ifndef NO_STREAKLINES
      // merge streakline data
      if ((col_pixel1->i & STREAKLINE) && (col_pixel2->i & STREAKLINE))
	{
	  if (col_pixel1->particle_z < col_pixel2->particle_z)
	    {
	      col_pixel2->particle_z        = col_pixel1->particle_z;
	      col_pixel2->particle_vel      = col_pixel1->particle_vel;
	      col_pixel2->particle_inlet_id = col_pixel1->particle_inlet_id;
	    }
	}
      else if ((col_pixel1->i & STREAKLINE) && !(col_pixel2->i & STREAKLINE))
	{
	  col_pixel2->particle_z        = col_pixel1->particle_z;
	  col_pixel2->particle_vel      = col_pixel1->particle_vel;
	  col_pixel2->particle_inlet_id = col_pixel1->particle_inlet_id;
	  
	  col_pixel2->i |= STREAKLINE;
	}
#endif
    }
}


void visWritePixel (ColPixel *col_pixel_p)
{
  int *col_pixel_id_p, i, j;
  
  
  i = PixelI(col_pixel_p->i);
  j = PixelJ(col_pixel_p->i);
  
  if (*(col_pixel_id_p = &col_pixel_id[ i*screen.pixels_y+j ]) != -1)
    {
      visMergePixels (col_pixel_p, &col_pixel_send[ *col_pixel_id_p ]);
    }
  else
    {
      if (col_pixels >= COLOURED_PIXELS_MAX)
	{
          return;
	}
      *col_pixel_id_p = col_pixels;
      
      memcpy (&col_pixel_send[ col_pixels ], col_pixel_p, sizeof(ColPixel));
      ++col_pixels;
    }
}


void rawWritePixel (ColPixel *col_pixel_p, unsigned int *pixel_index,
		    unsigned char rgb_data[],
		    void (*ColourPalette) (float value, float col[]))
{
  float density_col[3], stress_col[3], particle_col[3];
  
  int bits_per_char = sizeof(char) * 8;
  int pixel_i, pixel_j;
  
  unsigned char r1, g1, b1;
  unsigned char r2, g2, b2;
  unsigned char r3, g3, b3;
  unsigned char r4, g4, b4;
  
  
  // store pixel id
  pixel_i = PixelI(col_pixel_p->i);
  pixel_j = PixelJ(col_pixel_p->i);
  
  *pixel_index = (pixel_i << (2*bits_per_char)) + pixel_j;
  
  r1 = g1 = b1 = 255;
  r2 = g2 = b2 = 255;
  
  if (col_pixel_p->i & RT)
    {
      // store velocity volume rendering colour
      makePixelColour(r1, g1, b1, 
        (int)(col_pixel_p->vel_r / col_pixel_p->dt),
        (int)(col_pixel_p->vel_g / col_pixel_p->dt),
        (int)(col_pixel_p->vel_b / col_pixel_p->dt));
      
      if (lbm_stress_type != SHEAR_STRESS)
	{
	  // store von Mises stress volume rendering colour
          makePixelColour(r2, g2, b2,
	    (int)(col_pixel_p->stress_r / col_pixel_p->dt),
	    (int)(col_pixel_p->stress_g / col_pixel_p->dt),
	    (int)(col_pixel_p->stress_b / col_pixel_p->dt));
	}
      else if (col_pixel_p->stress < 1.0e+30F)
	{
	  ColourPalette (col_pixel_p->stress, stress_col);

	  // store wall shear stress colour
	  makePixelColour(r2, g2, b2,
            (int)(255.0F * stress_col[0]),
	    (int)(255.0F * stress_col[1]),
	    (int)(255.0F * stress_col[2]));
	}
      else
	{
	  r2 = g2 = b2 = 0;
	}
    }
  if (lbm_stress_type != SHEAR_STRESS && vis_mode == 0)
    {
      ColourPalette (col_pixel_p->density, density_col);
      ColourPalette (col_pixel_p->stress, stress_col);
      
      // store wall pressure colour
      makePixelColour(r3, g3, b3,
        (int)(255.0F * density_col[0]),
        (int)(255.0F * density_col[1]),
        (int)(255.0F * density_col[2]));
      
      // store von Mises stress colour
      makePixelColour(r4, g4, b4,
        (int)(255.0F * stress_col[0]),
        (int)(255.0F * stress_col[1]),
        (int)(255.0F * stress_col[2]));
    }
  else if (lbm_stress_type != SHEAR_STRESS && vis_mode == 1)
    {
      ColourPalette (col_pixel_p->density, density_col);
      ColourPalette (col_pixel_p->stress, stress_col);
      
      if (col_pixel_p->i & RT)
	{
	  if (!(col_pixel_p->i & GLYPH))
	    {
	      density_col[0] += 1.0F;
	      density_col[1] += 1.0F;
	      density_col[2] += 1.0F;
	      
	      stress_col[0] += 1.0F;
	      stress_col[1] += 1.0F;
	      stress_col[2] += 1.0F;
	    }
	  // store wall pressure (+glyph) colour
          makePixelColour(r3, g3, b3,
            (int)(127.5F * density_col[0]),
	    (int)(127.5F * density_col[1]),
	    (int)(127.5F * density_col[2]));
	  
	  // store von Mises stress (+glyph) colour 
          makePixelColour(r4, g4, b4,
            (int)(127.5F * stress_col[0]),
	    (int)(127.5F * stress_col[1]),
	    (int)(127.5F * stress_col[2]));
	}
      else
	{
	  r3 = g3 = b3 = 0;
	  r4 = g4 = b4 = 0;
	}
    }
  else
    {
      if (col_pixel_p->i & STREAKLINE)
	{
	  float scaled_vel = col_pixel_p->particle_vel * vis_velocity_threshold_max_inv;
	  
	  ColourPalette (scaled_vel, particle_col);
	  
	  // store particle colour
          makePixelColour(r3, g3, b3,
            (int)(255.0F * particle_col[0]),
	    (int)(255.0F * particle_col[1]),
	    (int)(255.0F * particle_col[2]));

          r4 = r3;
          g4 = g3;
          b4 = b3;
	}
      else
	{
	  // store pressure colour
	  r3 = g3 = b3 = (unsigned char)UtilityFunctions::enforceBounds((int)(127.5F * col_pixel_p->density), 0, 127);
	  
	  // store shear stress or von Mises stress
	  if (col_pixel_p->stress < 1.0e+30F)
	    {
	      r4 = g4 = b4 = (unsigned char)UtilityFunctions::enforceBounds((int)(127.5F * col_pixel_p->stress), 0, 127);
	    }
	  else
	    {
	      r4 = g4 = b4 = 0;
	    }
	} 
    }
  rgb_data[0] = r1; rgb_data[1] = g1; rgb_data[2] = b1;
  rgb_data[3] = r2; rgb_data[4] = g2; rgb_data[5] = b2;
  rgb_data[6] = r3; rgb_data[7] = g3; rgb_data[8] = b3;
  rgb_data[9] = r4; rgb_data[10] = g4; rgb_data[11] = b4;
  
  //pix_data[1] = (r1<<(3*bits_per_char)) + (g1<<(2*bits_per_char)) + (b1<<bits_per_char) + r2;
  //pix_data[2] = (g2<<(3*bits_per_char)) + (b2<<(2*bits_per_char)) + (r3<<bits_per_char) + g3;
  //pix_data[3] = (b3<<(3*bits_per_char)) + (r4<<(2*bits_per_char)) + (g4<<bits_per_char) + b4;
}

void makePixelColour(unsigned char& red, unsigned char& green, unsigned char& blue,
  int rawRed, int rawGreen, int rawBlue)
{
  red   = (unsigned char)UtilityFunctions::enforceBounds(rawRed, 0, 255);
  green = (unsigned char)UtilityFunctions::enforceBounds(rawGreen, 0, 255);
  blue  = (unsigned char)UtilityFunctions::enforceBounds(rawBlue, 0, 255);
}

void visRenderLine (float p1[], float p2[])
{
  int pixels_x, pixels_y;
  int x, y;
  int x1, y1;
  int x2, y2;
  int dx, dy;
  int incE, incNE;
  int d;
  int m;
  
  ColPixel col_pixel;
  
  
  pixels_x = screen.pixels_x;
  pixels_y = screen.pixels_y;
  
  x1 = (int)p1[0];
  y1 = (int)p1[1];
  
  x2 = (int)p2[0];
  y2 = (int)p2[1];
  
  if (x2 < x1)
    {
      x = x1;
      y = y1;
      x1 = x2;
      y1 = y2;
      x2 = x;
      y2 = y;
    }
  x = x1;
  y = y1;
  
  if (y1 < y2)
    {
      dy = y2 - y1;
      m = 1;
    }
  else
    {
      dy = y1 - y2;
      m = -1;
    }
  dx = x2 - x1;
    
  if (dx > dy)
    {
      incE = dy;
      d = dy - dx;
      incNE = d;
      
      while (x <= x2)
	{
	  if (!(x < 0 || x >= pixels_x ||
		y < 0 || y >= pixels_y))
	    {
	      col_pixel.i = PixelId(x, y) | GLYPH;
	      
	      visWritePixel (&col_pixel);
	    }
	  if (d < 0)
	    {
	      d += incE;
	    }
	  else
	    {
	      d += incNE;
	      y += m;
	    }
	  ++x;
	}
    }
  else if (y1 < y2)
    {
      incE = dx;
      d = dx - dy;
      incNE = d;
      
      while (y <= y2)
	{
	  if (!(x < 0 || x >= pixels_x ||
		y < 0 || y >= pixels_y))
	    {
	      col_pixel.i = PixelId(x, y) | GLYPH;
	      
	      visWritePixel (&col_pixel);
	    }
	  if (d < 0)
	    {
	      d += incE;
	    }
	  else
	    {
	      d += incNE;
	      ++x;
	    }
	  ++y;
	}
    }
  else
    {
      incE = dx;
      d = dx - dy;
      incNE = d;
      
      while (y >= y2)
	{
	  if (!(x < 0 || x >= pixels_x ||
		y < 0 || y >= pixels_y))
	    {
	      col_pixel.i = PixelId(x, y) | GLYPH;
	      
	      visWritePixel (&col_pixel);
	    }
	  if (d < 0)
	    {
	      d += incE;
	    }
	  else
	    {
	      d += incNE;
	      ++x;
	    }
	  --y;
	}
    }
}


void visReadParameters (char *parameters_file_name, LBM *lbm, Net *net, Vis *vis)
{
  FILE *parameters_file;
  
  float par_to_send[9];
  float ctr_x, ctr_y, ctr_z;
  float longitude, latitude;
  float zoom;
  float density_min, density_max, velocity_max, stress_max;
  float physical_velocity_max, physical_stress_max;
  
  int i;
  
  
  if (net->id == 0)
    {
      fprintf(stderr, "opening ray tracing configuration file %s\n", parameters_file_name);

      parameters_file = fopen (parameters_file_name, "r");

      if( parameters_file == NULL ) {
        fprintf(stderr, "unable to open file %s, exiting\n", parameters_file_name);
        fflush(0x0);
        exit(0x0);
      } else {
        fprintf(stderr, "done\n");
      }

      fflush(NULL);
      
      fscanf (parameters_file, "%e \n", &ctr_x);
      fscanf (parameters_file, "%e \n", &ctr_y);
      fscanf (parameters_file, "%e \n", &ctr_z);
      fscanf (parameters_file, "%e \n", &longitude);
      fscanf (parameters_file, "%e \n", &latitude);
      fscanf (parameters_file, "%e \n", &zoom);
      fscanf (parameters_file, "%e \n", &vis_brightness);
      fscanf (parameters_file, "%e \n", &physical_velocity_max);
      fscanf (parameters_file, "%e \n", &physical_stress_max);
      fclose (parameters_file);
      
      velocity_max = lbm->lbmConvertVelocityToLatticeUnits (physical_velocity_max);
      stress_max   = lbm->lbmConvertStressToLatticeUnits (physical_stress_max);
      
      par_to_send[ 0 ] = ctr_x;
      par_to_send[ 1 ] = ctr_y;
      par_to_send[ 2 ] = ctr_z;
      par_to_send[ 3 ] = longitude;
      par_to_send[ 4 ] = latitude;
      par_to_send[ 5 ] = zoom;
      par_to_send[ 6 ] = vis_brightness;
      par_to_send[ 7 ] = velocity_max;
      par_to_send[ 8 ] = stress_max;
    }
#ifndef NOMPI
  net->err = MPI_Bcast (par_to_send, 9, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  
  ctr_x          = par_to_send[ 0 ];
  ctr_y          = par_to_send[ 1 ];
  ctr_z          = par_to_send[ 2 ];
  longitude      = par_to_send[ 3 ];
  latitude       = par_to_send[ 4 ];
  zoom           = par_to_send[ 5 ];
  vis_brightness = par_to_send[ 6 ];
  velocity_max   = par_to_send[ 7 ];
  stress_max     = par_to_send[ 8 ];
  
  visProjection (0.5F * vis->system_size, 0.5F * vis->system_size,
		 512, 512,
		 ctr_x, ctr_y, ctr_z,
		 5.F * vis->system_size,
		 longitude, latitude,
		 0.5F * (5.F * vis->system_size),
		 zoom);
  
  density_min = +1.0e+30F;
  density_max = -1.0e+30F;
  
  for (i = 0; i < lbm->inlets; i++)
    {
      density_min = fminf(density_min, inlet_density_avg[ i ] - inlet_density_amp[ i ]);
      density_max = fmaxf(density_max, inlet_density_avg[ i ] + inlet_density_amp[ i ]);
    }
  for (i = 0; i < lbm->outlets; i++)
    {
      density_min = fminf(density_min, outlet_density_avg[ i ] - outlet_density_amp[ i ]);
      density_max = fmaxf(density_max, outlet_density_avg[ i ] + outlet_density_amp[ i ]);
    }
  vis_density_threshold_min = density_min;
  
  if (!is_bench)
    {
      vis_density_threshold_minmax_inv = 1.0F / (density_max - density_min);
      vis_velocity_threshold_max_inv   = 1.0F / velocity_max;
      vis_stress_threshold_max_inv     = 1.0F / stress_max;
    }
  else
    {
      vis_density_threshold_minmax_inv = 1.0F;
      vis_velocity_threshold_max_inv   = 1.0F;
      vis_stress_threshold_max_inv     = 1.0F;
    }
}
 

void visInit (Net *net, Vis *vis, streaklineDrawer *sl)
{
  blocks_yz = blocks_y * blocks_z;
  
  vis->half_dim[0] = 0.5F * (float)sites_x;
  vis->half_dim[1] = 0.5F * (float)sites_y;
  vis->half_dim[2] = 0.5F * (float)sites_z;
  
  vis->system_size = 2.F * fmaxf(vis->half_dim[0], fmaxf(vis->half_dim[1], vis->half_dim[2]));
  
  block_size_f = (float)block_size;
  block_size2 = block_size * block_size;
  block_size3 = block_size * block_size2;
  block_size_1 = block_size - 1;
  
  block_size_inv = 1.F / (float)block_size;
  
#ifndef NOMPI
  int col_pixel_count = 15;
  int col_pixel_blocklengths[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  MPI_Datatype col_pixel_types[15] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				      MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
				      MPI_FLOAT, MPI_FLOAT,
				      MPI_FLOAT,
				      MPI_FLOAT,
				      MPI_FLOAT,
				      MPI_FLOAT,
				      MPI_INT,
				      MPI_INT,
				      MPI_UB};
  
  MPI_Aint col_pixel_disps[15];
#endif
  
  col_pixels_max = COLOURED_PIXELS_MAX;
  
  vis_pixels_max = COLOURED_PIXELS_MAX;
  col_pixel_id = (int *)malloc(sizeof(int) * vis_pixels_max);
  
  for (int i = 0; i < COLOURED_PIXELS_MAX; i++)
    {
      col_pixel_id[i] = -1;
    }
#ifndef NOMPI
  // create the derived datatype for the MPI communications
  
  col_pixel_disps[0] = 0;
  
  for (int i = 1; i < col_pixel_count; i++)
    {
      if (col_pixel_types[i - 1] == MPI_FLOAT)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(float) * col_pixel_blocklengths[i - 1]);
	}
      else if (col_pixel_types[i - 1] == MPI_INT)
	{
	  col_pixel_disps[i] = col_pixel_disps[i - 1] + (sizeof(int) * col_pixel_blocklengths[i - 1]);
	}
    }
  MPI_Type_struct (col_pixel_count, col_pixel_blocklengths, col_pixel_disps, col_pixel_types, &MPI_col_pixel_type);
  MPI_Type_commit (&MPI_col_pixel_type);
#endif
  
  rtInit (net);
  
  if (!is_bench)
    {
      myGlypher = new glyphDrawer(net);
    }
#ifndef NO_STREAKLINES
  sl->slInit (net);
#endif
  vis_ctr_x -= vis->half_dim[0];
  vis_ctr_y -= vis->half_dim[1];
  vis_ctr_z -= vis->half_dim[2];
}


void visUpdateImageSize (int pixels_x, int pixels_y)
{
  if (pixels_x * pixels_y > screen.pixels_x * screen.pixels_y)
    {
      vis_pixels_max = pixels_x * pixels_y;
      col_pixel_id = (int *)realloc(col_pixel_id, sizeof(int) * vis_pixels_max);
    }
  for (int i = 0; i < pixels_x * pixels_y; i++)
    {
      col_pixel_id[ i ] = -1;
    }
}


void visCompositeImage (int recv_buffer_id, Net *net)
{
  // here, the communications needed to composite the image are
  // handled through a binary tree pattern and parallel pairwise
  // blocking communications.
  
  int *col_pixel_id_p;
  int col_pixels_temp;
  int comm_inc, send_id, recv_id;
  int i, j;
  int m, n;
  
  ColPixel *col_pixel1, *col_pixel2;
  
  
#ifndef NEW_COMPOSITING
  memcpy (col_pixel_recv[recv_buffer_id], col_pixel_send, col_pixels * sizeof(ColPixel));
#else
  if (net->id != 0)
    {
      memcpy (col_pixel_recv[recv_buffer_id], col_pixel_send, col_pixels * sizeof(ColPixel));
    }
#endif
  comm_inc = 1;
  m = 1;
  
  while (m < net->procs)
    {
      m <<= 1;
#ifndef NEW_COMPOSITING
      int start_id = 0;
#else
      int start_id = 1;
#endif
      for (recv_id = start_id; recv_id < net->procs;)
	{
	  send_id = recv_id + comm_inc;
	  
	  if (net->id != recv_id && net->id != send_id)
	    {
	      recv_id += comm_inc << 1;
	      continue;
	    }
	  if (send_id >= net->procs || recv_id == send_id)
	    {
	      recv_id += comm_inc << 1;
	      continue;
	    }
	  if (net->id == send_id)
	    {
#ifndef NOMPI
	      MPI_Send (&col_pixels, 1, MPI_INT, recv_id, 20, MPI_COMM_WORLD);
#endif
	      if (col_pixels > 0)
		{
#ifndef NOMPI
		  MPI_Send (col_pixel_send, col_pixels, MPI_col_pixel_type,
			    recv_id, 20, MPI_COMM_WORLD);
#endif
		}
	    }
	  else
	    {
#ifndef NOMPI
	      MPI_Recv (&col_pixels_temp, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD, net->status);
	      
	      if (col_pixels_temp > 0)
		{
		  MPI_Recv (col_pixel_send, col_pixels_temp, MPI_col_pixel_type,
			    send_id, 20, MPI_COMM_WORLD, net->status);
		}
#else
	      col_pixels_temp = 0;
#endif
	      for (n = 0; n < col_pixels_temp; n++)
		{
		  col_pixel1 = &col_pixel_send[ n ];
		  i = PixelI(col_pixel1->i);
		  j = PixelJ(col_pixel1->i);
		  
		  if (*(col_pixel_id_p = &col_pixel_id[ i * screen.pixels_y + j ]) == -1)
		    {
		      col_pixel2 = &col_pixel_recv[recv_buffer_id][ *col_pixel_id_p = col_pixels ];
		      
		      memcpy (col_pixel2, col_pixel1, sizeof(ColPixel));
		      ++col_pixels;
		    }
		  else
		    {
		      col_pixel2 = &col_pixel_recv[recv_buffer_id][ *col_pixel_id_p ];
		      
		      visMergePixels (col_pixel1, col_pixel2);
		    }
		}
	    }
	  if (m < net->procs && net->id == recv_id)
	    {
	      memcpy (col_pixel_send, col_pixel_recv[recv_buffer_id], col_pixels * sizeof(ColPixel));
	    }
	  recv_id += comm_inc << 1;
	}
      comm_inc <<= 1;
    }
#ifdef NEW_COMPOSITING
  if (net->id == 1)
    {
#ifndef NOMPI
      MPI_Send (&col_pixels, 1, MPI_INT, 0, 20, MPI_COMM_WORLD);
#endif
      if (col_pixels > 0)
	{
#ifndef NOMPI
	  MPI_Send (col_pixel_recv[recv_buffer_id], col_pixels, MPI_col_pixel_type,
		    0, 20, MPI_COMM_WORLD);
#endif
	}
    }
  else if (net->id == 0)
    {
      MPI_Recv (&col_pixels, 1, MPI_INT, 1, 20, MPI_COMM_WORLD, net->status);
      
      if (col_pixels > 0)
	{
	  MPI_Recv (col_pixel_recv[recv_buffer_id], col_pixels, MPI_col_pixel_type,
		    1, 20, MPI_COMM_WORLD, net->status);
	}
    }
#endif // NEW_COMPOSITING
}


void visRender (int recv_buffer_id, void (*ColourPalette) (float value, float col[]), Net *net, streaklineDrawer *sl)
{
  if (screen.pixels_x * screen.pixels_y > vis_pixels_max)
    {
      vis_pixels_max = UtilityFunctions::max(2 * vis_pixels_max, screen.pixels_x * screen.pixels_y);
      
      col_pixel_id = (int *)realloc(col_pixel_id, sizeof(int) * vis_pixels_max);
    }
  col_pixels = 0;
  
  rtRayTracing (ColourPalette);
  
  if (!is_bench && vis_mode == 1)
    {
      myGlypher->render();
    }
#ifndef NO_STREAKLINES
  if (vis_streaklines &&
      (lbm_stress_type == SHEAR_STRESS || vis_mode == 2))
    {
      sl->render();
    }
#endif
  visCompositeImage (recv_buffer_id, net);
  
  col_pixels_recv[recv_buffer_id] = col_pixels;
#ifndef NEW_COMPOSITING
  if (net->id == 0)
    {
      return;
    }
#endif
  for (int m = 0; m < col_pixels_recv[recv_buffer_id]; m++)
    {
      col_pixel_id[ (PixelI(col_pixel_send[m].i)) * screen.pixels_y + (PixelJ(col_pixel_send[m].i)) ] = -1;
    }
}


void visWriteImage (int recv_buffer_id, char *image_file_name,
		    void (*ColourPalette) (float value, float col[]))
{
  XdrFileWriter myWriter = XdrFileWriter(image_file_name);
  
  myWriter.writeInt(&vis_mode);

  myWriter.writeFloat(&vis_physical_pressure_threshold_min);
  myWriter.writeFloat(&vis_physical_pressure_threshold_max);
  myWriter.writeFloat(&vis_physical_velocity_threshold_max);
  myWriter.writeFloat(&vis_physical_stress_threshold_max);
  
  myWriter.writeInt(&screen.pixels_x);
  myWriter.writeInt(&screen.pixels_y);
  myWriter.writeInt(&col_pixels_recv[recv_buffer_id]);
  
  for (int n = 0; n < col_pixels_recv[recv_buffer_id]; n++)
  {
    myWriter.writePixel (&col_pixel_recv[recv_buffer_id][n], ColourPalette);
  }
}


void visCalculateMouseFlowField (ColPixel *col_pixel_p, LBM *lbm)
{
  double density = vis_density_threshold_min + col_pixel_p->density / vis_density_threshold_minmax_inv;
  double stress = col_pixel_p->stress / vis_stress_threshold_max_inv;
  
  vis_mouse_pressure = lbm->lbmConvertPressureToPhysicalUnits (density * Cs2);
  vis_mouse_stress = lbm->lbmConvertStressToPhysicalUnits (stress);
}


void visEnd (streaklineDrawer*sl)
{
#ifndef NO_STREAKLINES
  sl->slEnd ();
#endif
  if (!is_bench)
    {
      delete myGlypher;
    }
  rtEnd ();

  free(col_pixel_id);
}

