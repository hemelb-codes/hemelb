#include <vector>
#include <math.h>

#include "utilityFunctions.h"
#include "vis/Control.h"
#include "vis/Layer.h"
#include "vis/RayTracer.h"
#include "vis/GlyphDrawer.h"

#include "io/XdrFileWriter.h"

/* Hywel's bit
// Constructor for the controller.
Control::Control()
{
myLayers = std::vector<Layer*>();
}

// Destructor, which includes deleting all added layers.
// This is done in reverse order of addition.
Control::~Control()
{
  for(int ii = myLayers.size() - 1; ii >= 0; ii--)
    {
      delete myLayers[ii];
    }
}

// Adds a layer to the visualisation. Note that rendering is done in
// the order in which layers are added.
void Control::addLayer(Layer *newLayer)
{
  myLayers.push_back(newLayer);
}

// Combines the output of all visualisation methods in the order added.
void Control::render()
{
  for(int i = 0; i < myLayers.size(); i++) {
    myLayers[i]->render();
  }
}
*/

namespace hemelb
{
  namespace vis
  {
    // make a global controller
    Control *controller;
    
    Control::Control()
    {
      
      this->vis = new Vis;
      
      //sites_x etc are globals declared in net.h
      vis->half_dim[0] = 0.5F * float(sites_x);
      vis->half_dim[1] = 0.5F * float(sites_y);
      vis->half_dim[2] = 0.5F * float(sites_z);
      
      vis->system_size = 2.F * fmaxf(vis->half_dim[0],
				     fmaxf(vis->half_dim[1], vis->half_dim[2]));
      col_pixels_max = COLOURED_PIXELS_MAX;
      
      col_pixel_recv[0] = new ColPixel[col_pixels_max];
      col_pixel_recv[1] = new ColPixel[col_pixels_max];
      
      pixels_max = COLOURED_PIXELS_MAX;
      col_pixel_id = (int *)malloc(sizeof(int) * pixels_max);
  
      for (int i = 0; i < COLOURED_PIXELS_MAX; i++)
	{
	  col_pixel_id[i] = -1;
	}

    }
  
    
    void Control::initLayers(Net *net)
    {
      rtInit (net);
  
      if (!is_bench)
	{
	  myGlypher = new GlyphDrawer(net);
	}
#ifndef NO_STREAKLINES
      myStreaker = new StreaklineDrawer (net);
#endif
      // Note that rtInit does stuff to this->ctr_x (because this has
      // to be global)
      this->ctr_x -= vis->half_dim[0];
      this->ctr_y -= vis->half_dim[1];
      this->ctr_z -= vis->half_dim[2];
    }
    
    
    void Control::rotate (float sin_1, float cos_1,
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


    void Control::project (float p1[], float p2[])
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


    void Control::setProjection (int pixels_x, int pixels_y,
				 float local_ctr_x, float local_ctr_y, float local_ctr_z,
				 float longitude, float latitude,
				 float zoom)
    {
      float temp;
      float ortho_x = 0.5 * vis->system_size;
      float ortho_y = 0.5 * vis->system_size;
      float rad = 5.F * vis->system_size;
      float dist = 0.5 * rad;
      
      screen.max_x = ortho_x / zoom;
      screen.max_y = ortho_y / zoom;
      
      screen.pixels_x = pixels_x;
      screen.pixels_y = pixels_y;
  
      // Convert to radians
      temp = longitude * 0.01745329F;
      
      viewpoint.sin_1 = sinf(temp);
      viewpoint.cos_1 = cosf(temp);
      
      // Convert to radians
      temp = latitude * 0.01745329F;
  
      viewpoint.sin_2 = sinf(temp);
      viewpoint.cos_2 = cosf(temp);
  
      temp = rad * viewpoint.cos_2;
  
      viewpoint.x[0] = temp * viewpoint.sin_1 + local_ctr_x;
      viewpoint.x[1] = rad  * viewpoint.sin_2 + local_ctr_y;
      viewpoint.x[2] = temp * viewpoint.cos_1 + local_ctr_z;
  
      viewpoint.dist = dist;
  
      temp = dist / rad;
  
      local_ctr_x = viewpoint.x[0] + temp * (local_ctr_x - viewpoint.x[0]);
      local_ctr_y = viewpoint.x[1] + temp * (local_ctr_y - viewpoint.x[1]);
      local_ctr_z = viewpoint.x[2] + temp * (local_ctr_z - viewpoint.x[2]);
  
      screen.zoom = zoom;
  
      rotate (viewpoint.sin_1, viewpoint.cos_1,
	      viewpoint.sin_2, viewpoint.cos_2,
	      screen.max_x, 0.0F, 0.0F,
	      &screen.dir1[0], &screen.dir1[1], &screen.dir1[2]);
  
      rotate (viewpoint.sin_1, viewpoint.cos_1,
	      viewpoint.sin_2, viewpoint.cos_2,
	      0.0F, screen.max_y, 0.0F,
	      &screen.dir2[0], &screen.dir2[1], &screen.dir2[2]);
  
      screen.scale_x = (float)pixels_x / (2.F * screen.max_x);
      screen.scale_y = (float)pixels_y / (2.F * screen.max_y);
  
      screen.vtx[0] = local_ctr_x - screen.dir1[0] - screen.dir2[0] - viewpoint.x[0];
      screen.vtx[1] = local_ctr_y - screen.dir1[1] - screen.dir2[1] - viewpoint.x[1];
      screen.vtx[2] = local_ctr_z - screen.dir1[2] - screen.dir2[2] - viewpoint.x[2];
  
      screen.dir1[0] *= (2.F / (float)pixels_x);
      screen.dir1[1] *= (2.F / (float)pixels_x);
      screen.dir1[2] *= (2.F / (float)pixels_x);
  
      screen.dir2[0] *= (2.F / (float)pixels_y);
      screen.dir2[1] *= (2.F / (float)pixels_y);
      screen.dir2[2] *= (2.F / (float)pixels_y);
    }


    void Control::mergePixels (ColPixel *col_pixel1, ColPixel *col_pixel2)
    {
      // Merge raytracing data
    
      if (col_pixel1->i.isRt && col_pixel2->i.isRt) {
	// Both are ray-tracing
	col_pixel2->vel_r += col_pixel1->vel_r;
	col_pixel2->vel_g += col_pixel1->vel_g;
	col_pixel2->vel_b += col_pixel1->vel_b;
      
	if (lbm_stress_type != SHEAR_STRESS) {
	  col_pixel2->stress_r += col_pixel1->stress_r;
	  col_pixel2->stress_g += col_pixel1->stress_g;
	  col_pixel2->stress_b += col_pixel1->stress_b;
	}
      
	col_pixel2->dt += col_pixel1->dt;
      
	if (col_pixel1->t < col_pixel2->t) {
	  col_pixel2->t       = col_pixel1->t;
	  col_pixel2->density = col_pixel1->density;
	  col_pixel2->stress  = col_pixel1->stress;
	}
      
      } else if (col_pixel1->i.isRt && !col_pixel2->i.isRt) {
	// Only 1 is ray-tracing
	col_pixel2->vel_r = col_pixel1->vel_r;
	col_pixel2->vel_g = col_pixel1->vel_g;
	col_pixel2->vel_b = col_pixel1->vel_b;
      
	if (lbm_stress_type != SHEAR_STRESS) {
	  col_pixel2->stress_r = col_pixel1->stress_r;
	  col_pixel2->stress_g = col_pixel1->stress_g;
	  col_pixel2->stress_b = col_pixel1->stress_b;
	}
      
	col_pixel2->t       = col_pixel1->t;
	col_pixel2->dt      = col_pixel1->dt;
	col_pixel2->density = col_pixel1->density;
	col_pixel2->stress  = col_pixel1->stress;
      
	col_pixel2->i.isRt = true;
      }
      // Done merging ray-tracing
    
      // Now merge glyph data
      if (lbm_stress_type != SHEAR_STRESS && (mode == 0 || mode == 1)) {
      
	if (col_pixel1->i.isGlyph) {
	  col_pixel2->i.isGlyph = true;
	}
      
      } else {
      
#ifndef NO_STREAKLINES
	// merge streakline data
	if (col_pixel1->i.isStreakline &&
	    col_pixel2->i.isStreakline) {
	
	  if (col_pixel1->particle_z < col_pixel2->particle_z) {
	    col_pixel2->particle_z        = col_pixel1->particle_z;
	    col_pixel2->particle_vel      = col_pixel1->particle_vel;
	    col_pixel2->particle_inlet_id = col_pixel1->particle_inlet_id;
	  }
	
	}	else if (col_pixel1->i.isStreakline &&
			 !col_pixel2->i.isStreakline) {
	  col_pixel2->particle_z        = col_pixel1->particle_z;
	  col_pixel2->particle_vel      = col_pixel1->particle_vel;
	  col_pixel2->particle_inlet_id = col_pixel1->particle_inlet_id;
	
	  col_pixel2->i.isStreakline = true;
	}
#endif
      }
    
    }


    void Control::writePixel (ColPixel *col_pixel_p) {
      int *col_pixel_id_p, i, j;
    
    
      i = col_pixel_p->i.i;
      j = col_pixel_p->i.j;
    
      col_pixel_id_p = &col_pixel_id[ i*screen.pixels_y+j ];
    
      if (*col_pixel_id_p != -1) {
	mergePixels(col_pixel_p,
		    &col_pixel_send[ *col_pixel_id_p ]);
      
      } else { // col_pixel_id_p == -1
      
	if (col_pixels >= COLOURED_PIXELS_MAX) {
	  return;
	}
      
	*col_pixel_id_p = col_pixels;
      
	memcpy (&col_pixel_send[ col_pixels ],
		col_pixel_p,
		sizeof(ColPixel));
	++col_pixels;
      }
    
    }



    

    void Control::renderLine (float p1[], float p2[]) {
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
  
      x1 = int(p1[0]);
      y1 = int(p1[1]);
  
      x2 = int(p2[0]);
      y2 = int(p2[1]);
  
      if (x2 < x1) {
	x = x1;
	y = y1;
	x1 = x2;
	y1 = y2;
	x2 = x;
	y2 = y;
      }
    
      x = x1;
      y = y1;
    
      if (y1 < y2) {
	dy = y2 - y1;
	m = 1;
      
      } else {
	dy = y1 - y2;
	m = -1;
      }
    
      dx = x2 - x1;
    
      if (dx > dy) {
	incE = dy;
	d = dy - dx;
	incNE = d;
      
	while (x <= x2) {
	  if (!(x < 0 || x >= pixels_x ||
		y < 0 || y >= pixels_y)) {
	  
	    col_pixel.i = PixelId(x, y);
	    col_pixel.i.isGlyph = true;
	  
	    writePixel (&col_pixel);
	  }
	
	  if (d < 0) {
	    d += incE;
	  } else {
	    d += incNE;
	    y += m;
	  }
	  ++x;
	
	} // end while
      
      } else if (y1 < y2) {
	incE = dx;
	d = dx - dy;
	incNE = d;
      
	while (y <= y2) {
	  if (!(x < 0 || x >= pixels_x ||
		y < 0 || y >= pixels_y)) {
	  
	    col_pixel.i = PixelId(x, y);
	    col_pixel.i.isGlyph = true;
	  
	    writePixel (&col_pixel);
	  }
	
	  if (d < 0) {
	    d += incE;
	  } else {
	    d += incNE;
	    ++x;
	  }
	  ++y;
	
	} // while
      
      } else {
	incE = dx;
	d = dx - dy;
	incNE = d;
      
	while (y >= y2) {
	  if (!(x < 0 || x >= pixels_x ||
		y < 0 || y >= pixels_y)) {
	  
	    col_pixel.i = PixelId(x, y);
	    col_pixel.i.isGlyph = true;
	  
	    writePixel (&col_pixel);
	  }
	
	  if (d < 0) {
	    d += incE;
	  } else {
	    d += incNE;
	    ++x;
	  }
	  --y;
	
	} // while
      }
    }


    void Control::readParameters (char *parameters_file_name, LBM *lbm, Net *net)
    {
      FILE *parameters_file;
  
      float par_to_send[9];
      float local_ctr_x, local_ctr_y, local_ctr_z;
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
      
	  fscanf (parameters_file, "%e \n", &local_ctr_x);
	  fscanf (parameters_file, "%e \n", &local_ctr_y);
	  fscanf (parameters_file, "%e \n", &local_ctr_z);
	  fscanf (parameters_file, "%e \n", &longitude);
	  fscanf (parameters_file, "%e \n", &latitude);
	  fscanf (parameters_file, "%e \n", &zoom);
	  fscanf (parameters_file, "%e \n", &brightness);
	  fscanf (parameters_file, "%e \n", &physical_velocity_max);
	  fscanf (parameters_file, "%e \n", &physical_stress_max);
	  fclose (parameters_file);
      
	  velocity_max = lbm->lbmConvertVelocityToLatticeUnits (physical_velocity_max);
	  stress_max   = lbm->lbmConvertStressToLatticeUnits (physical_stress_max);
      
	  par_to_send[ 0 ] = local_ctr_x;
	  par_to_send[ 1 ] = local_ctr_y;
	  par_to_send[ 2 ] = local_ctr_z;
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
  
      local_ctr_x    = par_to_send[ 0 ];
      local_ctr_y    = par_to_send[ 1 ];
      local_ctr_z    = par_to_send[ 2 ];
      longitude      = par_to_send[ 3 ];
      latitude       = par_to_send[ 4 ];
      zoom           = par_to_send[ 5 ];
      brightness     = par_to_send[ 6 ];
      velocity_max   = par_to_send[ 7 ];
      stress_max     = par_to_send[ 8 ];
  
      setProjection(512, 512,
		    local_ctr_x, local_ctr_y, local_ctr_z,
		    longitude, latitude,
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
    

    void Control::updateImageSize (int pixels_x, int pixels_y)
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


    void Control::compositeImage (int recv_buffer_id, Net *net)
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
      memcpy (col_pixel_recv[recv_buffer_id],
	      col_pixel_send,
	      col_pixels * sizeof(ColPixel));
#else
      if (net->id != 0) {
	memcpy (col_pixel_recv[recv_buffer_id],
		col_pixel_send,
		col_pixels * sizeof(ColPixel));
      }
#endif
    
      comm_inc = 1;
      m = 1;
    
      while (m < net->procs) {
	m <<= 1;
#ifndef NEW_COMPOSITING
	int start_id = 0;
#else
	int start_id = 1;
#endif
	for (recv_id = start_id; recv_id < net->procs;) {
	  send_id = recv_id + comm_inc;
	
	  if (net->id != recv_id && net->id != send_id) {
	    recv_id += comm_inc << 1;
	    continue;
	  }
	
	  if (send_id >= net->procs || recv_id == send_id) {
	    recv_id += comm_inc << 1;
	    continue;
	  }
	
	  if (net->id == send_id) {
#ifndef NOMPI
	    MPI_Send(&col_pixels, 1,
		     MPI_INT, recv_id, 20, MPI_COMM_WORLD);
#endif
	    if (col_pixels > 0) {
#ifndef NOMPI
	      MPI_Send (col_pixel_send, col_pixels, ColPixel::getMpiType(),
			recv_id, 20, MPI_COMM_WORLD);
#endif
	    }
	  
	  } else {
#ifndef NOMPI
	    MPI_Recv (&col_pixels_temp, 1, MPI_INT, send_id, 20, MPI_COMM_WORLD, net->status);
	  
	    if (col_pixels_temp > 0) {
	      MPI_Recv (col_pixel_send, col_pixels_temp, ColPixel::getMpiType(),
			send_id, 20, MPI_COMM_WORLD, net->status);
	    }
#else
	    col_pixels_temp = 0;
#endif
	    for (n = 0; n < col_pixels_temp; n++) {
	      col_pixel1 = &col_pixel_send[ n ];
	      i = col_pixel1->i.i;
	      j = col_pixel1->i.j;
	    
	      if (*(col_pixel_id_p = &col_pixel_id[ i * screen.pixels_y + j ]) == -1) {
		col_pixel2 = &col_pixel_recv[recv_buffer_id][ *col_pixel_id_p = col_pixels ];
	      
		memcpy (col_pixel2, col_pixel1, sizeof(ColPixel));
		++col_pixels;
	      
	      } else {
		col_pixel2 = &col_pixel_recv[recv_buffer_id][ *col_pixel_id_p ];
	      
		mergePixels (col_pixel1, col_pixel2);
	      }
	    
	    }
	  }
	  if (m < net->procs && net->id == recv_id) {
	    memcpy(col_pixel_send,
		   col_pixel_recv[recv_buffer_id],
		   col_pixels * sizeof(ColPixel));
	  }
	
	  recv_id += comm_inc << 1;
	}
	comm_inc <<= 1;
      }
#ifdef NEW_COMPOSITING
      if (net->id == 1) {
#ifndef NOMPI
	MPI_Send (&col_pixels, 1, MPI_INT, 0, 20, MPI_COMM_WORLD);
#endif
	if (col_pixels > 0) {
#ifndef NOMPI
	  MPI_Send (col_pixel_recv[recv_buffer_id], col_pixels, MPI_col_pixel_type,
		    0, 20, MPI_COMM_WORLD);
#endif
	}
      
      } else if (net->id == 0) {
	MPI_Recv (&col_pixels, 1, MPI_INT, 1, 20, MPI_COMM_WORLD, net->status);
      
	if (col_pixels > 0) {
	  MPI_Recv (col_pixel_recv[recv_buffer_id], col_pixels, MPI_col_pixel_type,
		    1, 20, MPI_COMM_WORLD, net->status);
	}
      
      }
#endif // NEW_COMPOSITING
    }
  

    void Control::render (int recv_buffer_id, void (*ColourPalette) (float value, float col[]), Net *net)
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
      if (shouldDrawStreaklines &&
	  (lbm_stress_type == SHEAR_STRESS || mode == 2))
	{
	  myStreaker->render();
	}
#endif
      compositeImage (recv_buffer_id, net);
  
      col_pixels_recv[recv_buffer_id] = col_pixels;
#ifndef NEW_COMPOSITING
      if (net->id == 0)
	{
	  return;
	}
#endif
      for (int m = 0; m < col_pixels_recv[recv_buffer_id]; m++)
	{
	  col_pixel_id[ col_pixel_send[m].i.i * screen.pixels_y + col_pixel_send[m].i.j ] = -1;
	}
    }


    void Control::writeImage (int recv_buffer_id, char *image_file_name,
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


    void Control::calculateMouseFlowField (ColPixel *col_pixel_p, LBM *lbm)
    {
      double density = density_threshold_min + col_pixel_p->density / density_threshold_minmax_inv;
      double stress = col_pixel_p->stress / stress_threshold_max_inv;
  
      mouse_pressure = lbm->lbmConvertPressureToPhysicalUnits (density * Cs2);
      mouse_stress = lbm->lbmConvertStressToPhysicalUnits (stress);
    }

    void Control::streaklines(int time_step, int period, Net *net)
    {
      myStreaker->streakLines (time_step, period, net);
    }

    void Control::restart()
    {
      myStreaker->restart();
    }

    Control::~Control()
    {
#ifndef NO_STREAKLINES
      delete myStreaker;
#endif
      if (!is_bench)
	{
	  delete myGlypher;
	}
      rtEnd ();

      delete col_pixel_recv[0];
      delete col_pixel_recv[1];

      free(col_pixel_id);
    }

  } // namespace vis
} // namespace hemelb
