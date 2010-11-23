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

void ioReadConfig (Vis *vis);
void ioWriteCheckpoint(Vis *vis);
void ioReadCheckpoint(Vis *vis);
void ioWriteConfig(Vis *vis);
void ioWriteCoords(Vis *vis);
void ioWritePars(Vis *vis);
void ioSaveWindowImage(const char *file_name);

#endif // IO
