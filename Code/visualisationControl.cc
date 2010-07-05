#include "visualisationControl.h"

// Constructor for the controller.
visualisationControl::visualisationControl()
{
  myLayers = std::vector<visualisationLayer*>();
}

// Destructor, which includes deleting all added layers.
// This is done in reverse order of addition.
visualisationControl::~visualisationControl()
{
  for(int ii = myLayers.size() - 1; ii >= 0; ii--)
  {
    delete myLayers[ii];
  }
}

// Adds a layer to the visualisation. Note that rendering is done in the order in which
// layers are added.
void visualisationControl::addLayer(visualisationLayer *newLayer)
{
   myLayers.push_back(newLayer);
}

// Combines the output of all visualisation methods in the order added.
void visualisationControl::render()
{
  for(int ii = 0; ii < myLayers.size(); ii++)
  {
    myLayers[ii]->render();
  }
}