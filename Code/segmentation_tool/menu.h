#ifndef MENU
#define MENU

#include "rt.h"
#include "vis.h"


#define NULL_MENU_OPTION            0
#define SEGMENT_1X                  (1<<0)
#define SEGMENT_2X                  (1<<1)
#define SEGMENT_3X                  (1<<2)
#define SEGMENT_4X                  (1<<3)
#define SEGMENT_5X                  (1<<4)
#define SEGMENT_6X                  (1<<5)
#define ZOOM_SCENE                  (1<<8)
#define ROTATE_SCENE                (1<<9)
#define CREATE_INLET                (1<<10)
#define CREATE_OUTLET               (1<<11)
#define CREATE_WALL                 (1<<12)
#define SCALE_BOUNDARY              (1<<13)
#define ROTATE_BOUNDARY             (1<<14)
#define REVERSE_INLET_NORMAL        (1<<15)
#define CHANGE_PRESSURE_AMPLITUDE   (1<<16)
#define CHANGE_MEAN_PRESSURE        (1<<17)
#define CHANGE_PRESSURE_PHASE       (1<<18)
#define DELETE_BOUNDARY             (1<<19)
#define SAVE_DATA                   (1<<21)
#define QUIT                        (1<<22)
#define ACTIVE                      1


void menuProcessMenuEvents (int option);
void menuCreateMenu (Vis *vis);

#endif // MENU
