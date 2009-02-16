#include "config.h"
#include "vis.h"
#include "editing.h"


int main (int argc, char *argv[])
{
  visInit (argc, argv, &vis);
  
  
  glutReshapeFunc (Reshape);
  glutIdleFunc (Visualise);
  glutDisplayFunc (Visualise);
  glutMouseFunc (MouseFunction);
  glutMotionFunc (MotionFunction);
  glutPassiveMotionFunc (PassiveMotionFunction);
  glutMainLoop ();

  return(0);
}
