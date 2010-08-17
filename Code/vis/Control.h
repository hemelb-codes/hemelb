#ifndef __vis_Control_h_
#define __vis_Control_h_

#include <vector>
#include "vis/Layer.h"

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
#endif //__vis_Control_h_
