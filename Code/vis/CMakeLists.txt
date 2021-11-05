# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_vis
  GlyphDrawer.cc Control.cc Screen.cc Viewpoint.cc BasicPixel.cc Rendering.cc ResultPixel.cc
  rayTracer/ClusterNormal.cc rayTracer/ClusterWithWallNormals.cc rayTracer/HSLToRGBConverter.cc
  rayTracer/RayDataEnhanced.cc rayTracer/RayDataNormal.cc
  streaklineDrawer/NeighbouringProcessor.cc streaklineDrawer/Particle.cc streaklineDrawer/ParticleManager.cc
  streaklineDrawer/StreaklineDrawer.cc streaklineDrawer/VelocityField.cc streaklineDrawer/StreakPixel.cc
  )
target_link_libraries(hemelb_vis PRIVATE hemelb_configuration)
