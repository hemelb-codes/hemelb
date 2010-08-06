#include <math.h>
#include <string.h>

#include "utilityFunctions.h"
#include "rt.h"
#include "lb.h"
#include "io/xdrFileWriter.h"

namespace vis {

  glyphDrawer *myGlypher;
  streaklineDrawer *myStreaker;

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
    if (lbm_stress_type != SHEAR_STRESS && (mode == 0 || mode == 1))
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



  void makePixelColour(unsigned char& red, unsigned char& green, unsigned char& blue,
		       int rawRed, int rawGreen, int rawBlue)
  {
    red   = (unsigned char)util::enforceBounds(rawRed, 0, 255);
    green = (unsigned char)util::enforceBounds(rawGreen, 0, 255);
    blue  = (unsigned char)util::enforceBounds(rawBlue, 0, 255);
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
	fscanf (parameters_file, "%e \n", &brightness);
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
	par_to_send[ 6 ] = brightness;
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
    brightness = par_to_send[ 6 ];
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
	density_min = fminf(density_min, lbm->inlet_density_avg[ i ] - lbm->inlet_density_amp[ i ]);
	density_max = fmaxf(density_max, lbm->inlet_density_avg[ i ] + lbm->inlet_density_amp[ i ]);
      }
    for (i = 0; i < lbm->outlets; i++)
      {
	density_min = fminf(density_min, lbm->outlet_density_avg[ i ] - lbm->outlet_density_amp[ i ]);
	density_max = fmaxf(density_max, lbm->outlet_density_avg[ i ] + lbm->outlet_density_amp[ i ]);
      }
    density_threshold_min = density_min;
  
    if (!is_bench)
      {
	density_threshold_minmax_inv = 1.0F / (density_max - density_min);
	velocity_threshold_max_inv   = 1.0F / velocity_max;
	stress_threshold_max_inv     = 1.0F / stress_max;
      }
    else
      {
	density_threshold_minmax_inv = 1.0F;
	velocity_threshold_max_inv   = 1.0F;
	stress_threshold_max_inv     = 1.0F;
      }
  }
 

  void visInit (Net *net, Vis *vis)
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
  
    pixels_max = COLOURED_PIXELS_MAX;
    col_pixel_id = (int *)malloc(sizeof(int) * pixels_max);
  
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
    myStreaker = new streaklineDrawer (net);
#endif
    ctr_x -= vis->half_dim[0];
    ctr_y -= vis->half_dim[1];
    ctr_z -= vis->half_dim[2];
  }


  void visUpdateImageSize (int pixels_x, int pixels_y)
  {
    if (pixels_x * pixels_y > screen.pixels_x * screen.pixels_y)
      {
	pixels_max = pixels_x * pixels_y;
	col_pixel_id = (int *)realloc(col_pixel_id, sizeof(int) * pixels_max);
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


  void visRender (int recv_buffer_id, void (*ColourPalette) (float value, float col[]), Net *net)
  {
    if (screen.pixels_x * screen.pixels_y > pixels_max)
      {
	pixels_max = util::max(2 * pixels_max, screen.pixels_x * screen.pixels_y);
      
	col_pixel_id = (int *)realloc(col_pixel_id, sizeof(int) * pixels_max);
      }
    col_pixels = 0;
  
    rtRayTracing (ColourPalette);
  
    if (!is_bench && mode == 1)
      {
	myGlypher->render();
      }
#ifndef NO_STREAKLINES
    if (streaklines &&
	(lbm_stress_type == SHEAR_STRESS || mode == 2))
      {
	myStreaker->render();
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
    io::XdrFileWriter writer = io::XdrFileWriter(image_file_name);
  
    writer << mode;
  
    writer << physical_pressure_threshold_min 
	   << physical_pressure_threshold_max
	   << physical_velocity_threshold_max
	   << physical_stress_threshold_max;
  
    writer << screen.pixels_x
	   << screen.pixels_y
	   << col_pixels_recv[recv_buffer_id];
  
    for (int n = 0; n < col_pixels_recv[recv_buffer_id]; n++)
      {
	writer.writePixel(&col_pixel_recv[recv_buffer_id][n],
			  ColourPalette);
      }
  }


  void visCalculateMouseFlowField (ColPixel *col_pixel_p, LBM *lbm)
  {
    double density = density_threshold_min + col_pixel_p->density / density_threshold_minmax_inv;
    double stress = col_pixel_p->stress / stress_threshold_max_inv;
  
    mouse_pressure = lbm->lbmConvertPressureToPhysicalUnits (density * Cs2);
    mouse_stress = lbm->lbmConvertStressToPhysicalUnits (stress);
  }

  void visStreaklines(int time_step, int period, Net *net)
  {
    myStreaker->streakLines (time_step, period, net);
  }

  void visRestart()
  {
    myStreaker->restart();
  }

  void visEnd ()
  {
#ifndef NO_STREAKLINES
    delete myStreaker;
#endif
    if (!is_bench)
      {
	delete myGlypher;
      }
    rtEnd ();

    free(col_pixel_id);
  }

}
