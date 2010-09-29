#include "vis/ColPixel.h"

#include "lb.h"
#include "utilityFunctions.h"
#include "vis/Control.h"

namespace hemelb
{
  namespace vis
  {

#ifndef NOMPI
    MPI_Datatype MPI_col_pixel_type;
#endif
    PixelId::PixelId(int i_, int j_) :isRt(false),
				      isGlyph(false),
				      isStreakline(false),
				      i(i_),
				      j(j_) {}
    
    PixelId::PixelId() :  isRt(false),
			  isGlyph(false),
			  isStreakline(false),
			  i(0),
			  j(0) { }
    
#ifndef NOMPI
    MPI_Datatype ColPixel::mpiType = MPI_DATATYPE_NULL;
#endif// NOMPI
    
    // create the derived datatype for the MPI communications
    void ColPixel::registerMpiType() {
#ifndef NOMPI
      int col_pixel_count = 15;
      int col_pixel_blocklengths[15] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
, 1};
  
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
      MPI_Type_struct (col_pixel_count, col_pixel_blocklengths,
                       col_pixel_disps, col_pixel_types,
                       &mpiType);
      MPI_Type_commit (&mpiType);
#endif
    }
    
    const MPI_Datatype& ColPixel::getMpiType()
    {
      if (mpiType == MPI_DATATYPE_NULL) {
        registerMpiType();
      }
      return mpiType;
    }


    void ColPixel::makePixelColour(unsigned char& red, unsigned char& green, unsigned char& blue,
				   int rawRed, int rawGreen, int rawBlue)
    {
      red   = (unsigned char)util::enforceBounds(rawRed, 0, 255);
      green = (unsigned char)util::enforceBounds(rawGreen, 0, 255);
      blue  = (unsigned char)util::enforceBounds(rawBlue, 0, 255);
    }
    
    void ColPixel::rawWritePixel (unsigned int *pixel_index,
				  unsigned char rgb_data[],
				  ColourPaletteFunction *colourPalette)
    {
      float density_col[3], stress_col[3], particle_col[3];
  
      int bits_per_char = sizeof(char) * 8;
      int pixel_i, pixel_j;
  
      unsigned char r1, g1, b1;
      unsigned char r2, g2, b2;
      unsigned char r3, g3, b3;
      unsigned char r4, g4, b4;
  
  
      // store pixel id
      pixel_i = i.i;
      pixel_j = i.j;
  
      *pixel_index = (pixel_i << (2*bits_per_char)) + pixel_j;
  
      r1 = g1 = b1 = 255;
      r2 = g2 = b2 = 255;
  
      if (i.isRt) {
	// store velocity volume rendering colour
	makePixelColour(r1, g1, b1, 
			int(vel_r / dt),
			int(vel_g / dt),
			int(vel_b / dt));
      
	if (lbm_stress_type != SHEAR_STRESS) {
	  // store von Mises stress volume rendering colour
	  makePixelColour(r2, g2, b2,
			  int(stress_r / dt),
			  int(stress_g / dt),
			  int(stress_b / dt));
	
	} else if (stress < 1.0e+30F) {
	  colourPalette(stress, stress_col);

	  // store wall shear stress colour
	  makePixelColour(r2, g2, b2,
			  int(255.0F * stress_col[0]),
			  int(255.0F * stress_col[1]),
			  int(255.0F * stress_col[2]));
	
	} else {
	  r2 = g2 = b2 = 0;
	}
      
      } // if (isRt)
    
      if (lbm_stress_type != SHEAR_STRESS && controller->mode == 0) {
	colourPalette(density, density_col);
	colourPalette(stress, stress_col);
      
	// store wall pressure colour
	makePixelColour(r3, g3, b3,
			int(255.0F * density_col[0]),
			int(255.0F * density_col[1]),
			int(255.0F * density_col[2])
			);
      
	// store von Mises stress colour
	makePixelColour(r4, g4, b4,
			int(255.0F * stress_col[0]),
			int(255.0F * stress_col[1]),
			int(255.0F * stress_col[2])
			);
      
      } else if (lbm_stress_type != SHEAR_STRESS && controller->mode == 1) {
	colourPalette(density, density_col);
	colourPalette(stress, stress_col);
      
	if (i.isRt) {
	  if (!i.isGlyph) {
	    density_col[0] += 1.0F;
	    density_col[1] += 1.0F;
	    density_col[2] += 1.0F;
	  
	    stress_col[0] += 1.0F;
	    stress_col[1] += 1.0F;
	    stress_col[2] += 1.0F;
	  }
	
	  // store wall pressure (+glyph) colour
	  makePixelColour(r3, g3, b3,
			  int(127.5F * density_col[0]),
			  int(127.5F * density_col[1]),
			  int(127.5F * density_col[2]));
	
	  // store von Mises stress (+glyph) colour 
	  makePixelColour(r4, g4, b4,
			  int(127.5F * stress_col[0]),
			  int(127.5F * stress_col[1]),
			  int(127.5F * stress_col[2]));
	} else {
	  r3 = g3 = b3 = 0;
	  r4 = g4 = b4 = 0;
	}
      
      } else {
      
	if (i.isStreakline) {
	  float scaled_vel = particle_vel * controller->velocity_threshold_max_inv;
	
	  colourPalette(scaled_vel, particle_col);
	
	  // store particle colour
	  makePixelColour(r3, g3, b3,
			  int(255.0F * particle_col[0]),
			  int(255.0F * particle_col[1]),
			  int(255.0F * particle_col[2]));
	
	  r4 = r3;
	  g4 = g3;
	  b4 = b3;
	
	} else {
	  // store pressure colour
	  r3 = g3 = b3 = (unsigned char)util::enforceBounds(int(127.5F * density), 0, 127);
	
	  // store shear stress or von Mises stress
	  if (stress < 1.0e+30F) {
	    r4 = g4 = b4 = (unsigned char)util::enforceBounds(int(127.5F * stress), 0, 127);
	  
	  } else {
	    r4 = g4 = b4 = 0;
	  }
	}
      
      }
    
      rgb_data[0] = r1; rgb_data[1] = g1; rgb_data[2] = b1;
      rgb_data[3] = r2; rgb_data[4] = g2; rgb_data[5] = b2;
      rgb_data[6] = r3; rgb_data[7] = g3; rgb_data[8] = b3;
      rgb_data[9] = r4; rgb_data[10] = g4; rgb_data[11] = b4;
    
    }
    
    
  }
}
