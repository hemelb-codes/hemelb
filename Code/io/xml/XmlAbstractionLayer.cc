#include "io/xml/XmlAbstractionLayer.h"
#include "util/fileutils.h"

namespace hemelb
{
  namespace configuration
  {
    XmlAbstractionLayer::XmlAbstractionLayer(std::string path) :
      xmlDocument(new TiXmlDocument())
    {
      util::check_file(path.c_str());
      xmlDocument->LoadFile(path);
      ResetToTopLevel();
    }

    XmlAbstractionLayer::~XmlAbstractionLayer()
    {
      delete xmlDocument;
    }

    void XmlAbstractionLayer::ResetToTopLevel()
    {
      while (!parentNodes.empty())
        parentNodes.pop();
      currentNode = xmlDocument->RootElement();
    }

    bool XmlAbstractionLayer::MoveToChild()
    {
      TiXmlElement* nextNode = currentNode->FirstChildElement();
      if (nextNode == NULL)
        return false;

      parentNodes.push(currentNode);
      currentNode = nextNode;
      return true;
    }

    bool XmlAbstractionLayer::MoveToChild(std::string name)
    {
      TiXmlElement* nextNode = currentNode->FirstChildElement(name);
      if (nextNode == NULL)
        return false;

      parentNodes.push(currentNode);
      currentNode = nextNode;
      return true;
    }

    bool XmlAbstractionLayer::NextSibling()
    {
      currentNode = currentNode->NextSiblingElement();
      return (currentNode != NULL);
    }

    bool XmlAbstractionLayer::NextSibling(std::string name)
    {
      currentNode = currentNode->NextSiblingElement(name);
      return (currentNode != NULL);
    }

    bool XmlAbstractionLayer::MoveToParent()
    {
      if (parentNodes.empty())
        return false;
      currentNode = parentNodes.top();
      parentNodes.pop();
      return true;
    }

    bool XmlAbstractionLayer::GetUnsignedLongValue(std::string name, unsigned long& value)
    {
      const std::string* const stringValue = currentNode->Attribute(name);
      if (stringValue == NULL)
      {
        value = 0;
        return false;
      }
      char *dummy;
      // Read in, in base 10.
      value = std::strtoul(stringValue->c_str(), &dummy, 10);
      return true;
    }

    bool XmlAbstractionLayer::GetDoubleValue(std::string name, double& value)
    {
      const std::string* const stringValue = currentNode->Attribute(name);
      if (stringValue == NULL)
      {
        value = 0.0;
        return false;
      }
      char *dummy;
      value = std::strtod(stringValue->c_str(), &dummy);
      return true;
    }

    bool XmlAbstractionLayer::GetDoubleVector(std::string name, util::Vector3D<double>& vector)
    {
      bool ok = true;

      ok &= MoveToChild(name);
      if (!ok)
        return false;

      double x,y,z;
      ok &= GetDoubleValue("xValue", x);
      ok &= GetDoubleValue("yValue", y);
      ok &= GetDoubleValue("zValue", z);
      if (ok)
      {
        vector.x = x;
        vector.y = y;
        vector.z = z;
      }

      MoveToParent();
      return ok;
    }
  }
}
