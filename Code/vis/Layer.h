#ifndef __vis_Layer_h_
#define __vis_Layer_h_

namespace vis {
  
  // Base class for any method of adding to the visualisation,
  // e.g. Glyphs, Streamlines, Ray-tracing etc.
  class Layer
  {
  public:
    // Method to render the layer's output to a visualisation.
    virtual void render() = 0;
  };
  
}

#endif //__vis_Layer_h_
