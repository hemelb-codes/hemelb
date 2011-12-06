add_library(hemelb_vis
	GlyphDrawer.cc Control.cc Screen.cc Viewpoint.cc BasicPixel.cc Rendering.cc ResultPixel.cc
	rayTracer/ClusterNormal.cc rayTracer/ClusterWithWallNormals.cc rayTracer/HSLToRGBConverter.cc
	rayTracer/RayDataEnhanced.cc rayTracer/RayDataNormal.cc
	streaklineDrawer/NeighbouringProcessor.cc streaklineDrawer/Particle.cc streaklineDrawer/ParticleManager.cc
	streaklineDrawer/StreaklineDrawer.cc streaklineDrawer/VelocityField.cc streaklineDrawer/StreakPixel.cc
	)