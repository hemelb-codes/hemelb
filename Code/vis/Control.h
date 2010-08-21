#ifndef HEME_VIS_CONTROL_H
#define HEME_VIS_CONTROL_H

#include <vector>
#include "vis/Layer.h"

namespace heme {
  namespace vis {
  
    // Class to control and use the effects of different visualisation methods.
    class Control
    {
    public:
      Control();
      ~Control();
    
      // Adds a layer to the visualisation. Note that rendering is done
      // in the order in which layers are added.
      void addLayer(Layer *newLayer);
    
      // Combines the output of all visualisation methods in the order
      // added.
      void render();
    
    private:
      int nLayers;
      std::vector<Layer *> myLayers;
    
    };

  }
}

#endif // HEME_VIS_CONTROL_H
