#include "vis/colPixel.h"
#include "vis/rt.h"
#include "utilityFunctions.h"

namespace vis {

#ifndef NOMPI
  MPI_Datatype MPI_col_pixel_type;
#endif

  void ColPixel::rawWritePixel (unsigned int *pixel_index,
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
    pixel_i = PixelI(this->i);
    pixel_j = PixelJ(this->i);
  
    *pixel_index = (pixel_i << (2*bits_per_char)) + pixel_j;
  
    r1 = g1 = b1 = 255;
    r2 = g2 = b2 = 255;
  
    if (this->i & RT)
      {
	// store velocity volume rendering colour
	makePixelColour(r1, g1, b1, 
			(int)(this->vel_r / this->dt),
			(int)(this->vel_g / this->dt),
			(int)(this->vel_b / this->dt));
      
	if (lbm_stress_type != SHEAR_STRESS)
	  {
	    // store von Mises stress volume rendering colour
	    makePixelColour(r2, g2, b2,
			    (int)(this->stress_r / this->dt),
			    (int)(this->stress_g / this->dt),
			    (int)(this->stress_b / this->dt));
	  }
	else if (this->stress < 1.0e+30F)
	  {
	    ColourPalette (this->stress, stress_col);

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
    if (lbm_stress_type != SHEAR_STRESS && mode == 0)
      {
	ColourPalette (this->density, density_col);
	ColourPalette (this->stress, stress_col);
      
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
    else if (lbm_stress_type != SHEAR_STRESS && mode == 1)
      {
	ColourPalette (this->density, density_col);
	ColourPalette (this->stress, stress_col);
      
	if (this->i & RT)
	  {
	    if (!(this->i & GLYPH))
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
	if (this->i & STREAKLINE)
	  {
	    float scaled_vel = this->particle_vel * velocity_threshold_max_inv;
	  
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
	    r3 = g3 = b3 = (unsigned char)util::enforceBounds((int)(127.5F * this->density), 0, 127);
	  
	    // store shear stress or von Mises stress
	    if (this->stress < 1.0e+30F)
	      {
		r4 = g4 = b4 = (unsigned char)util::enforceBounds((int)(127.5F * this->stress), 0, 127);
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

}
