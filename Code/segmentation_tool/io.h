#ifndef IO
#define IO


#include <rpc/types.h>
#include <rpc/xdr.h>
#include <dirent.h>
#include <errno.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "config.h"
#include "editing.h"
#include "segmentation.h"

using namespace std;


#ifndef MESH
// int ioGetDir (string dir, vector<string> &files);
// void ioGetFileNames (Vis *vis);
// void ioReadFirstSlice (Vis *vis);
// void ioReadSlice (int slice_id, Vis *vis);
void ioReadConfig (Vis *vis);
#else
void ioReadConfig (Vis *vis);
#endif
void ioWriteCheckpoint (Vis *vis);
void ioReadCheckpoint (Vis *vis);
void ioWriteConfig (Vis *vis);
void ioWriteCoords (Vis *vis);
void ioWritePars (Vis *vis);
void ioSaveWindowImage (const char *file_name);


#endif // IO
