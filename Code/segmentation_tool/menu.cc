#include "menu.h"


void menuProcessMenuEvents (int option)
{
  block_id b_id;
  triangle_id t_id;
  
  if (option & (SEGMENT_1X|SEGMENT_2X|SEGMENT_3X|SEGMENT_4X|SEGMENT_5X|SEGMENT_6X))
    {
      if      (option == SEGMENT_1X) { vis.res_factor = 1; }
      else if (option == SEGMENT_2X) { vis.res_factor = 2; }
      else if (option == SEGMENT_3X) { vis.res_factor = 3; }
      else if (option == SEGMENT_4X) { vis.res_factor = 4; }
      else if (option == SEGMENT_5X) { vis.res_factor = 5; }
      else if (option == SEGMENT_6X) { vis.res_factor = 6; }

      editRescaleGrid (&vis);
      
      if (segSegmentation (&vis) == SUCCESS)
	{
	  visDeleteCubeDisplayList ();
	  visCreateCubeDisplayList (&vis);
	}
      else
	{
	  editRescaleGrid (&vis);
	  segSegmentation (&vis);
	}
    }
  else if (option & (ZOOM_SCENE|ROTATE_SCENE))
    {
      vis.menu.option = option;
      vis.mouse.b_id = -1;
    }
  else if (option & (CREATE_INLET|CREATE_OUTLET|CREATE_WALL))
    {
      if      (option == CREATE_INLET)  { b_id = INLET_BOUNDARY; }
      else if (option == CREATE_OUTLET) { b_id = OUTLET_BOUNDARY; }
      else if (option == CREATE_WALL)   { b_id = WALL_BOUNDARY; }
      t_id = segCreateOptimisedTriangle (b_id, &vis);
    }
  else if (option & (SCALE_BOUNDARY|ROTATE_BOUNDARY))
    {
      vis.menu.option = option;
      
      if (vis.mouse.b_id < 0) return;
    }
  else if (option & (REVERSE_INLET_NORMAL|DELETE_BOUNDARY))
    {
      if (vis.mouse.b_id < 0) return;
      
      if (option == REVERSE_INLET_NORMAL && ((unsigned int)vis.mouse.b_id) == INLET_BOUNDARY)
	{
	  editInvertTriangleNormal (vis.mouse.b_id, vis.mouse.t_id, &vis);
	}
      else if (option == DELETE_BOUNDARY)
	{
	  editDeleteTriangle (vis.mouse.b_id, vis.mouse.t_id, &vis);
	  
	  vis.mouse.b_id = -1;
	}
    }
  else if (option & (CHANGE_PRESSURE_AMPLITUDE|CHANGE_MEAN_PRESSURE|CHANGE_PRESSURE_PHASE))
    {
      vis.menu.option = option;
    }
  else if (option == SAVE_DATA)
    {
      printf("Opening ppm file: ./snapshot.ppm\n");
      std::string lSavePath = "./snapshot.ppm";
      ioSaveWindowImage (lSavePath.c_str());
      
      segSetBoundaryConfigurations (&vis);
      
      printf("Opening output pars file: %s\n", vis.output_pars);
      ioWritePars (&vis);
      
      printf("Opening output config file: %s\n", vis.output_config);
      ioWriteConfig (&vis);
      
      printf("Opening output coordinates file: %s\n", vis.output_coords);
      ioWriteCoords (&vis);
      
      printf("Opening checkpoint file: %s\n", vis.checkpoint);
      ioWriteCheckpoint (&vis);
    }
  else if (option == QUIT)
    {
      visEndBoundaries (&vis);
      rtEndRayTracing (&vis.mesh);
      visEnd (&vis);
      
      exit(0);
    }
}


void menuCreateMenu (Vis *vis)
{
  vis->menu.option = NULL_MENU_OPTION;
  vis->mouse.state = !ACTIVE;
  vis->mouse.b_id = -1;
  
  vis->menu.id = glutCreateMenu (menuProcessMenuEvents);

    {
      glutAddMenuEntry ("Segmentation 1X", SEGMENT_1X);
      glutAddMenuEntry ("Segmentation 2X", SEGMENT_2X);
      glutAddMenuEntry ("Segmentation 3X", SEGMENT_3X);
      glutAddMenuEntry ("Segmentation 4X", SEGMENT_4X);
      glutAddMenuEntry ("Segmentation 5X", SEGMENT_5X);
      glutAddMenuEntry ("Segmentation 6X", SEGMENT_6X);
      glutAddMenuEntry ("Zoom scene", ZOOM_SCENE);
      glutAddMenuEntry ("Rotate scene", ROTATE_SCENE);
      glutAddMenuEntry ("Create inlet", CREATE_INLET);
      glutAddMenuEntry ("Create outlet", CREATE_OUTLET);
      glutAddMenuEntry ("Create wall", CREATE_WALL);
      glutAddMenuEntry ("Scale boundary", SCALE_BOUNDARY);
      glutAddMenuEntry ("Rotate boundary", ROTATE_BOUNDARY);
      glutAddMenuEntry ("Reverse inlet normal", REVERSE_INLET_NORMAL);
      glutAddMenuEntry ("Change pressure amplitude", CHANGE_PRESSURE_AMPLITUDE);
      glutAddMenuEntry ("Change mean pressure", CHANGE_MEAN_PRESSURE);
      glutAddMenuEntry ("Change pressure phase", CHANGE_PRESSURE_PHASE);
      glutAddMenuEntry ("Delete boundary", DELETE_BOUNDARY);
      glutAddMenuEntry ("Save data", SAVE_DATA);
      glutAddMenuEntry ("Quit", QUIT);
    }
  glutAttachMenu (GLUT_RIGHT_BUTTON);
}
