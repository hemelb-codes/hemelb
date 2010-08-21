#include <vector>

#include "vis/Control.h"
#include "vis/Layer.h"

using namespace heme::vis;

// Constructor for the controller.
Control::Control()
{
  myLayers = std::vector<Layer*>();
}

// Destructor, which includes deleting all added layers.
// This is done in reverse order of addition.
Control::~Control()
{
  for(int ii = myLayers.size() - 1; ii >= 0; ii--)
    {
      delete myLayers[ii];
    }
}

// Adds a layer to the visualisation. Note that rendering is done in
// the order in which layers are added.
void Control::addLayer(Layer *newLayer)
{
  myLayers.push_back(newLayer);
}

// Combines the output of all visualisation methods in the order added.
void Control::render()
{
  for(int i = 0; i < myLayers.size(); i++) {
    myLayers[i]->render();
  }
}
  
