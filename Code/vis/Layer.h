#ifndef HEMELB_VIS_LAYER_H
#define HEMELB_VIS_LAYER_H

namespace hemelb
{
  namespace vis
  {
  
    // Base class for any method of adding to the visualisation,
    // e.g. Glyphs, Streamlines, Ray-tracing etc.
    class Layer
    {
    public:
      // Method to render the layer's output to a visualisation.
      virtual void render() = 0;
      virtual ~Layer() = 0;
    protected:
      Layer();
    };
  
  }
}

#endif // HEMELB_VIS_LAYER_H
