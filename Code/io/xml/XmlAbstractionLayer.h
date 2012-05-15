#ifndef HEMELB_CONFIGURATION_XMLABSTRACTIONLAYER_H
#define HEMELB_CONFIGURATION_XMLABSTRACTIONLAYER_H

#include <cstdlib>
#include <stack>
#include "tinyxml.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace configuration
  {
    class XmlAbstractionLayer
    {
      public:
        XmlAbstractionLayer(std::string path);
        
        ~XmlAbstractionLayer();

        void ResetToTopLevel();
        bool MoveToChild();
        bool MoveToChild(std::string name);
        bool NextSibling();
        bool NextSibling(std::string name);
        bool MoveToParent();

        bool GetUnsignedLongValue(std::string name, unsigned long& value);
        bool GetDoubleValue(std::string name, double& value);
        bool GetDoubleVector(std::string name, util::Vector3D<double>& vector);

      private:
        TiXmlDocument *xmlDocument;
        TiXmlElement  *currentNode;
        std::stack<TiXmlElement*> parentNodes;
    };
  }
}
#endif /* HEMELB_CONFIGURATION_XMLABSTRACTIONLAYER_H */
