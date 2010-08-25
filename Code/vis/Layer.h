#ifndef HEME_VIS_LAYER_H
#define HEME_VIS_LAYER_H

namespace heme
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

#endif // HEME_VIS_LAYER_H
