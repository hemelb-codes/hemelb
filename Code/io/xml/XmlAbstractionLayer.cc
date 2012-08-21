// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "io/xml/XmlAbstractionLayer.h"
#include "util/fileutils.h"

namespace hemelb
{
  namespace io
  {
    namespace xml
    {
      // TODO: use the converter object to convert from physical units to lattice units
      XmlAbstractionLayer::XmlAbstractionLayer(
                             const std::string path,
                             const util::UnitConverter& converter) :
        xmlDocument(new TiXmlDocument()), converter(converter)
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

      bool XmlAbstractionLayer::MoveToChild(const std::string name)
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

      bool XmlAbstractionLayer::NextSibling(const std::string name)
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

      bool XmlAbstractionLayer::GetUnsignedLongValue(const std::string name, unsigned long& value)
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

      bool XmlAbstractionLayer::GetDoubleValue(const std::string name, double& value)
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

      bool XmlAbstractionLayer::GetDoubleValueAndConvert(const std::string name, double& value)
      {
        bool ok = GetDoubleValue(name, value);
        if (!ok)
          return false;
        
        const std::string* const stringUnits = currentNode->Attribute((std::string)"units");
        if (stringUnits == NULL)
        {
          value = 0.0;
          return false;
        }
        ok &= converter.Convert(*stringUnits, value);
        return ok;
      }

      bool XmlAbstractionLayer::GetDoubleVector(const std::string name, util::Vector3D<double>& vector)
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

      bool XmlAbstractionLayer::GetDoubleVectorAndConvert(const std::string name, util::Vector3D<double>& vector)
      {
        bool ok = true;

        ok &= MoveToChild(name);
        if (!ok)
          return false;

        double x,y,z;
        ok &= GetDoubleValueAndConvert("xValue", x);
        ok &= GetDoubleValueAndConvert("yValue", y);
        ok &= GetDoubleValueAndConvert("zValue", z);
        if (ok)
        {
          vector.x = x;
          vector.y = y;
          vector.z = z;
        }

        MoveToParent();
        return ok;

      }

      bool XmlAbstractionLayer::GetLatticePosition(const std::string name, LatticePosition& vector)
      {
        bool ok = GetDoubleVectorAndConvert(name, vector);
        if (!ok) return false;

        LatticePosition physicalOrigin = converter.GetPhysicalOrigin();
        vector.x += physicalOrigin.x;
        vector.y += physicalOrigin.y;
        vector.z += physicalOrigin.z;

        return true;
      }

      bool XmlAbstractionLayer::GetString(const std::string name, std::string& value)
      {
        const std::string* const stringValue = currentNode->Attribute(name);
        if (stringValue == NULL)
        {
          value = std::string();
          return false;
        }
        value = *stringValue;
        return true;
      }

    }
  }
}
