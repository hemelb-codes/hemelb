#include "glyphDrawer.h"

double glyphDrawer::vis_glyph_length = -1.F;

glyphDrawer::glyphDrawer(Net *net)
{
  int glyphs_max;
  int site_i, site_j, site_k;
  int i, j, k;
  int m, n;
  int c1, c2;
  
  DataBlock *map_block_p;
  
  glyphs_max = 0;
  
  for (n = 0; n < blocks; n++)
    {
      if (net->map_block[ n ].site_data != NULL)
	{
	  ++glyphs_max;
	}
    }
  glyph = new Glyph[glyphs_max];
  
  glyphs = 0;
  
  if (check_conv)
    {
      c1 = 30;
      c2 = 15;
    }
  else
    {
      c1 = 15;
      c2 = 0;
    }
  n = -1;
  
  for (i = 0; i < sites_x; i += block_size)
    {
      for (j = 0; j < sites_y; j += block_size)
	{
	  for (k = 0; k < sites_z; k += block_size)
	    {
	      map_block_p = &net->map_block[ ++n ];
	      
	      if (map_block_p->site_data == NULL) continue;
	      
	      site_i = (block_size >> 1);
	      site_j = (block_size >> 1);
	      site_k = (block_size >> 1);
	      
	      m = (((site_i << shift) + site_j) << shift) + site_k;
	      
	      if (map_block_p->site_data[ m ] & (1U << 31U)) continue;
	      
	      glyph[ glyphs ].x = (float)(i + site_i) - 0.5F * (float)sites_x;
	      glyph[ glyphs ].y = (float)(j + site_j) - 0.5F * (float)sites_y;
	      glyph[ glyphs ].z = (float)(k + site_k) - 0.5F * (float)sites_z;
	      
	      glyph[ glyphs ].f = &f_old[ map_block_p->site_data[m]*c1+c2 ];
	      ++glyphs;
	    }
	}
    }
}

glyphDrawer::~glyphDrawer()
{
  delete[] glyph;
}

void glyphDrawer::render()
{
  double density;
  double vx, vy, vz;
  double temp;
  
  float screen_max[4];
  float scale[4];
  float p1[3], p2[3], p3[3], p4[3];
  
  int n;
  
  
  screen_max[0] = screen.max_x;
  screen_max[1] = screen.max_x;
  screen_max[2] = screen.max_y;
  screen_max[3] = screen.max_y;
  
  scale[0] = scale[1] = screen.scale_x;
  scale[2] = scale[3] = screen.scale_y;
  
  for (n = 0; n < glyphs; n++)
    {
      lbmDensityAndVelocity (glyph[ n ].f, &density, &vx, &vy, &vz);
      
      temp = vis_glyph_length * block_size * vis_velocity_threshold_max_inv / density;
      
      vx *= temp;
      vy *= temp;
      vz *= temp;
      
      p1[0] = glyph[n].x;
      p1[1] = glyph[n].y;
      p1[2] = glyph[n].z;
      
      p2[0] = glyph[n].x + vx;
      p2[1] = glyph[n].y + vy;
      p2[2] = glyph[n].z + vz;
      
      visProject (p1, p3);
      visProject (p2, p4);
      
      p3[0] = scale[0] * (p3[0] + screen_max[0]);
      p3[1] = scale[1] * (p3[1] + screen_max[1]);
      p4[0] = scale[2] * (p4[0] + screen_max[2]);
      p4[1] = scale[3] * (p4[1] + screen_max[3]);
      
      visRenderLine (p3, p4);
    }
}     