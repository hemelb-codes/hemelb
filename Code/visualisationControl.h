#ifndef __visualisationControl_h_
#define __visualisationControl_h_

#include <vector>

// Base class for any method of adding to the visualisation, e.g. Glyphs, Streamlines, Ray-tracing etc.
class visualisationLayer
{
  public:
    // Method to render the layer's output to a visualisation.
    virtual void render() = 0;
};

// Class to control and use the effects of different visualisation methods.
class visualisationControl
{
  private:
    int numberOfLayers;
    std::vector<visualisationLayer*> myLayers;
  
  public:
    visualisationControl();
    ~visualisationControl();

    // Adds a layer to the visualisation. Note that rendering is done in the order in which
    // layers are added.
    void addLayer(visualisationLayer *newLayer);
    // Combines the output of all visualisation methods in the order added.
    void render();      
};

#endif //__visualisationControl_h_