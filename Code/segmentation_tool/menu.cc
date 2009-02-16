#include "menu.h"


void menuProcessMenuEvents (int option)
{
  int b_id, t_id;
  
  ScreenVoxel *voxel_p;
  
  
  if (option == CHANGE_SLICE)
    {
      vis.menu.option = CHANGE_SLICE;
    }
  else if (option & (SEGMENT_1X|SEGMENT_2X|SEGMENT_3X|SEGMENT_4X|SEGMENT_5X|SEGMENT_6X))
    {
      float res_factor_temp = vis.res_factor;
      
      if      (option == SEGMENT_1X) { vis.res_factor = 1; }
      else if (option == SEGMENT_2X) { vis.res_factor = 2; }
      else if (option == SEGMENT_3X) { vis.res_factor = 3; }
      else if (option == SEGMENT_4X) { vis.res_factor = 4; }
      else if (option == SEGMENT_5X) { vis.res_factor = 5; }
      else if (option == SEGMENT_6X) { vis.res_factor = 6; }
      
      if (vis.mode == 0)
	{
	  vis.selected_voxel[0] = (vis.mouse.x[0] * vis.input_voxels[0]) / vis.viewport_pixels[0];
	  vis.selected_voxel[1] = (vis.mouse.x[1] * vis.input_voxels[1]) / vis.viewport_pixels[1];
	  
	  vis.selected_voxel[0] = max(0, min(vis.input_voxels[0]-1, vis.selected_voxel[0]));
	  vis.selected_voxel[1] = max(0, min(vis.input_voxels[1]-1, vis.selected_voxel[1]));
	  
	  vis.selected_grey = vis.voxel[ VoxelId(vis.selected_voxel,vis.input_voxels) ];
	}
      editRescaleSystem (&vis);
      editRescaleTriangles (vis.res_factor / res_factor_temp, &vis);
      
      if (segSegmentation (&vis) == SUCCESS)
	{
	  visCalculateSceneCenter (&vis);
	  visRescaleViewpoint (vis.res_factor / res_factor_temp, &vis);
	}
      else
	{
	  float temp = vis.res_factor;
	  vis.res_factor = res_factor_temp;
	  res_factor_temp = temp;
	  
	  editRescaleSystem (&vis);
	  editRescaleTriangles (vis.res_factor / res_factor_temp, &vis);
	  segSegmentation (&vis);
	}
    }
  else if (option == CHANGE_THRESHOLD)
    {
      vis.menu.option = CHANGE_THRESHOLD;
    }
  else if (option & (ZOOM_SCENE|ROTATE_SCENE))
    {
      vis.menu.option = option;
      vis.mouse.b_id = -1;
    }
  else if (option & (CREATE_INLET|CREATE_OUTLET|CREATE_WALL))
    {
      voxel_p = visScreenVoxelPtr (vis.mouse.x, &vis);
      
      if (voxel_p->site[0] < 0) return;
      
      if      (option == CREATE_INLET)  { b_id = INLET_BOUNDARY; }
      else if (option == CREATE_OUTLET) { b_id = OUTLET_BOUNDARY; }
      else if (option == CREATE_WALL)   { b_id = WALL_BOUNDARY; }
      
      t_id = segCreateOptimisedTriangle (b_id, voxel_p->site, &vis);
      
      if (t_id == -1) return;
      
      vis.mouse.b_id = b_id;
      vis.mouse.t_id = t_id;
    }
  else if (option & (ZOOM_BOUNDARY|ROTATE_BOUNDARY))
    {
      vis.menu.option = option;
      
      if (vis.mouse.b_id < 0) return;
    }
  else if (option & (REVERSE_INLET_NORMAL|DELETE_BOUNDARY))
    {
      if (vis.mouse.b_id < 0) return;
      
      if (option == REVERSE_INLET_NORMAL && vis.mouse.b_id == INLET_BOUNDARY)
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
  else if (option == CHANGE_VIS_MODE)
    {
      vis.menu.option = NULL_MENU_OPTION;
      
      if (vis.mode == 0)
	{
	  vis.mode = 1;
	  glutChangeToMenuEntry (8,  "Zoom scene", ZOOM_SCENE);
	  glutChangeToMenuEntry (9,  "Rotate scene", ROTATE_SCENE);
	  glutChangeToMenuEntry (10, "Create inlet", CREATE_INLET);
	  glutChangeToMenuEntry (11, "Create outlet", CREATE_OUTLET);
	  glutAddMenuEntry ("Create wall", CREATE_WALL);
	  glutAddMenuEntry ("Scale boundary", ZOOM_BOUNDARY);
	  glutAddMenuEntry ("Rotate boundary", ROTATE_BOUNDARY);
	  glutAddMenuEntry ("Reverse inlet normal", REVERSE_INLET_NORMAL);
	  glutAddMenuEntry ("Change pressure amplitude", CHANGE_PRESSURE_AMPLITUDE);
	  glutAddMenuEntry ("Change mean pressure", CHANGE_MEAN_PRESSURE);
	  glutAddMenuEntry ("Change pressure phase", CHANGE_PRESSURE_PHASE);
	  glutAddMenuEntry ("Delete boundary", DELETE_BOUNDARY);
	  glutAddMenuEntry ("2D rendering", CHANGE_VIS_MODE);
	  glutAddMenuEntry ("Save data", SAVE_DATA);
	  glutAddMenuEntry ("Quit", QUIT);
	}
      else
	{
	  vis.mode = 0;
	  vis.mouse.state = !ACTIVE;
	  vis.mouse.b_id = -1;
	  
	  glutChangeToMenuEntry (8,  "Change slice", CHANGE_SLICE);
	  glutChangeToMenuEntry (9,  "3D rendering", CHANGE_VIS_MODE);
	  glutChangeToMenuEntry (10, "Save data", SAVE_DATA);
	  glutChangeToMenuEntry (11, "Quit", QUIT);
	  
	  for (int i = glutGet(GLUT_MENU_NUM_ITEMS); i >= 12; i--)
	    {
	      glutRemoveMenuItem (i);
	    }
	}
    }
  else if (option == SAVE_DATA)
    {
      printf("Opening ppm file: ./snapshot.ppm\n");
      ioSaveWindowImage ("./snapshot.ppm");
      
      segSetBoundaryConfigurations (&vis);
      
      printf("Opening output pars file: %s\n", vis.output_pars);
      ioWritePars (&vis);
      
      printf("Opening output config file: %s\n", vis.output_config);
      ioWriteConfig (&vis);
      
      printf("Opening checkpoint file: %s\n", vis.checkpoint);
      ioWriteCheckpoint (&vis);
    }
  else if (option == QUIT)
    {
      visEndBoundaries (&vis);
      
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
  
  if (vis->mode == 0)
    {
      glutAddMenuEntry ("Segmentation 1X", SEGMENT_1X);
      glutAddMenuEntry ("Segmentation 2X", SEGMENT_2X);
      glutAddMenuEntry ("Segmentation 3X", SEGMENT_3X);
      glutAddMenuEntry ("Segmentation 4X", SEGMENT_4X);
      glutAddMenuEntry ("Segmentation 5X", SEGMENT_5X);
      glutAddMenuEntry ("Segmentation 6X", SEGMENT_6X);
      glutAddMenuEntry ("Change threshold", CHANGE_THRESHOLD);
      glutAddMenuEntry ("Change slice", CHANGE_SLICE);
      glutAddMenuEntry ("3D rendering", CHANGE_VIS_MODE);
      glutAddMenuEntry ("Save data", SAVE_DATA);
      glutAddMenuEntry ("Quit", QUIT);
    }
  else
    {
      glutAddMenuEntry ("Segmentation 1X", SEGMENT_1X);
      glutAddMenuEntry ("Segmentation 2X", SEGMENT_2X);
      glutAddMenuEntry ("Segmentation 3X", SEGMENT_3X);
      glutAddMenuEntry ("Segmentation 4X", SEGMENT_4X);
      glutAddMenuEntry ("Segmentation 5X", SEGMENT_5X);
      glutAddMenuEntry ("Segmentation 6X", SEGMENT_6X);
      glutAddMenuEntry ("Change threshold", CHANGE_THRESHOLD);
      glutAddMenuEntry ("Zoom scene", ZOOM_SCENE);
      glutAddMenuEntry ("Rotate scene", ROTATE_SCENE);
      glutAddMenuEntry ("Create inlet", CREATE_INLET);
      glutAddMenuEntry ("Create outlet", CREATE_OUTLET);
      glutAddMenuEntry ("Create wall", CREATE_WALL);
      glutAddMenuEntry ("Scale boundary", ZOOM_BOUNDARY);
      glutAddMenuEntry ("Rotate boundary", ROTATE_BOUNDARY);
      glutAddMenuEntry ("Reverse inlet normal", REVERSE_INLET_NORMAL);
      glutAddMenuEntry ("Change pressure amplitude", CHANGE_PRESSURE_AMPLITUDE);
      glutAddMenuEntry ("Change mean pressure", CHANGE_MEAN_PRESSURE);
      glutAddMenuEntry ("Change pressure phase", CHANGE_PRESSURE_PHASE);
      glutAddMenuEntry ("Delete boundary", DELETE_BOUNDARY);
      glutAddMenuEntry ("2D rendering", CHANGE_VIS_MODE);
      glutAddMenuEntry ("Save data", SAVE_DATA);
      glutAddMenuEntry ("Quit", QUIT);
    }
  glutAttachMenu (GLUT_RIGHT_BUTTON);
}
