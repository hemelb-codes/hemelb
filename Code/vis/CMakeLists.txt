# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.


file(GLOB_RECURSE HEMELB_CU_FILES ${CMAKE_SOURCE_DIR}/*.cc)
set_source_files_properties( ${HEMELB_CU_FILES} PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ )
cuda_add_library(hemelb_vis
  GlyphDrawer.cc Control.cc Screen.cc Viewpoint.cc BasicPixel.cc Rendering.cc ResultPixel.cc
  rayTracer/ClusterNormal.cc rayTracer/ClusterWithWallNormals.cc rayTracer/HSLToRGBConverter.cc
  rayTracer/RayDataEnhanced.cc rayTracer/RayDataNormal.cc
  streaklineDrawer/NeighbouringProcessor.cc streaklineDrawer/Particle.cc streaklineDrawer/ParticleManager.cc
  streaklineDrawer/StreaklineDrawer.cc streaklineDrawer/VelocityField.cc streaklineDrawer/StreakPixel.cc
  )
