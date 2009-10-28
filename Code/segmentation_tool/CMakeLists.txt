PROJECT(segtool C CXX)
#SET(ITK_DIR "/home/steven/ITK/install-06-02-2009/lib/InsightToolkit")
#SET(CMAKE_CXX_COMPILER "/usr/bin/g++")
cmake_minimum_required(VERSION 2.6)
INCLUDE(FindGLUT)
INCLUDE(FindOpenGL)
ADD_DEFINITIONS(-DMESH)
#FIND_PACKAGE(ITK)
FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(GLUT REQUIRED)
#INCLUDE(${ITK_USE_FILE}) 
include_directories(${OPENGL_INCLUDE_DIR})
ADD_EXECUTABLE(segtool config.cc  editing.cc  io.cc  main.cc  math.cc  menu.cc
rt.cc segmentation.cc timing.cc usage.cc vis.cc io_dicom.cc)
#TARGET_LINK_LIBRARIES(segtool ITKAlgorithms ITKBasicFilters ITKCommon ITKFEM
#ITKIO ITKNrrdIO ITKNumerics ITKSpatialObject ITKStatistics ITKniftiio ITKznz)
target_link_libraries(segtool ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES})

