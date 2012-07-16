// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_IO_XML_XMLABSTRACTIONLAYER_H
#define HEMELB_IO_XML_XMLABSTRACTIONLAYER_H

#include <cstdlib>
#include <stack>
#include "tinyxml.h"
#include "util/Vector3D.h"

namespace hemelb
{
  namespace io
  {
    namespace xml
    {
      /** an abstraction for an XML document
       *
       * this class localises the dependency on an external XML library
       */
      class XmlAbstractionLayer
      {
        public:
          /**
           * constructor
           *
           * @param $path
           *   the path to the XML file to be read by this object
           */
          XmlAbstractionLayer(std::string path);
        
          /** destructor */
          ~XmlAbstractionLayer();

          /** resets the internal position pointer to the top-level or root node */
          void ResetToTopLevel();

          /**
           * moves to the first child element irrespective of name
           *
           * moves the internal position pointer
           * to the first child element of the current element, if any
           * if no suitable child element exists then
           * the internal position pointer is not moved
           *
           * @return
           *   returns true if a suitable child element was found
           *   or false if the internal position pointer was not moved
           */
          bool MoveToChild();

          /**
           * moves to the first child element with the specified name
           *
           * moves the internal position pointer
           * to the first child element with the specified name, if any
           * if no suitable child element exists then
           * the internal position pointer is not moved
           *
           * @param $name
           *   the name of the child element to move to
           *
           * @return
           *   returns true if a suitable child element was found
           *   or false if the internal position pointer was not moved
           */
          bool MoveToChild(std::string name);

          /**
           * moves to the next sibling element irrespective of name
           *
           * moves the internal position pointer
           * to the next sibling element of the current element, if any
           * if no suitable sibling element exists then
           * the internal position pointer is not moved
           *
           * @return
           *   returns true if a suitable child element was found
           *   or false if the internal position pointer was not moved
           */
          bool NextSibling();

          /**
           * moves to the next sibling element with the specified name
           *
           * moves the internal position pointer
           * to the next sibling element with the specified name, if any
           * if no suitable sibling element exists then
           * the internal position pointer is not moved
           *
           * @param $name
           *   the name of the sibling element to move to
           *
           * @return
           *   returns true if a suitable child element was found
           *   or false if the internal position pointer was not moved
           */
          bool NextSibling(std::string name);

          /**
           * moves to the parent element irrespective of name
           *
           * moves the internal position pointer
           * to the parent element of the current element, if any
           * if no parent element exists (i.e. the current node is the root) then
           * the internal position pointer is not moved
           *
           * @return
           *   returns true if a parent element was found
           *   or false if the internal position pointer was not moved
           */
          bool MoveToParent();

          /**
           * reads an integer valued attribute from the current element
           *
           * reads the value of the attribute with the specified name
           * and converts it to an unsigned long before returning it
           *
           * @param $name
           *   the name of the attribute within the current element
           * @param $value
           *   the converted value of the attribute (output)
           *
           * @return
           *   returns true if the attribute was found and converted successfully
           *   or false if the attribute was not found or could not be converted
           */
          bool GetUnsignedLongValue(std::string name, unsigned long& value);

          /**
           * reads a floating-point valued attribute from the current element
           *
           * reads the value of the attribute with the specified name
           * and converts it to a double before returning it
           *
           * @param $name
           *   the name of the attribute within the current element
           * @param $value
           *   the converted value of the attribute (output)
           *
           * @return
           *   returns true if the attribute was found and converted successfully
           *   or false if the attribute was not found or could not be converted
           */
          bool GetDoubleValue(std::string name, double& value);

          /**
           * reads a 3D vector of double values from the current element
           *
           * reads the value of the vector with the specified name
           * the name specifies the name of a child element
           * that contains attributes for each co-ordinate value
           * the attributes are named: xValue, yValue, and zValue
           * each co-ordinate attribute value is converted to a double
           * before setting the x, y, and z properties of the vector parameter
           * if this function returns false then no change was made to vector
           *
           * @param $name
           *   the name of the child vector element within the current element
           * @param $value
           *   the converted value of the child vector element (output)
           *
           * @return
           *   returns true if the vector was found and converted successfully
           *   or false if the vector was not found or could not be converted
           */
          bool GetDoubleVector(std::string name, util::Vector3D<double>& vector);

        private:
          /** a pointer to the xml file being abstracted by this object */
          TiXmlDocument *xmlDocument;

          /** a pointer to the current node - the current position pointer */
          TiXmlElement  *currentNode;

          /** a stack containing pointers to all parents of the current node */
          std::stack<TiXmlElement*> parentNodes;
      };
    }
  }
}
#endif /* HEMELB_CONFIGURATION_XMLABSTRACTIONLAYER_H */
