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
#include <sstream>
#include "Exception.h"

// Forward declare the TinyXML types needed.
class TiXmlDocument;
class TiXmlElement;

namespace hemelb
{
  namespace io
  {
    namespace xml
    {

      class Element
      {
        public:
          static Element Missing();

          Element(TiXmlElement* el);
          ~Element();

          /**
           * Get the name of the element
           * @return
           */
          const std::string& GetName() const;
          void GetName(const std::string& str);

          /**
           * Gets the first child element with the specified name
           *
           * @param $name
           *   the name of the child element to return
           *
           * @return
           *   returns the child element if it was found or
           *   Element::Missing() if not
           */
          Element GetChildOrNull(const std::string& name);
          /**
           * Gets the first child element with the specified name or throw
           * ChildError if it does not exist.
           *
           * @param $name
           *   the name of the child element to return
           *
           * @return
           *   returns the child element
           */
          Element GetChildOrThrow(const std::string& name);

          /**
           * Return the next sibling element with the specified name, if any
           * If no suitable sibling element exists then return Element::Missing()
           *
           * @param $name
           *   the name of the sibling element to get
           *
           * @return
           *   Next sibling or Element::Missing()
           */
          Element NextSiblingOrNull(const std::string name);

          /**
           * Return the next sibling element with the specified name, if any
           * If no suitable sibling element exists then throw SiblingError
           *
           * @param $name
           *   the name of the sibling element to get
           *
           * @return
           *   Next sibling.
           */
          Element NextSiblingOrThrow(const std::string name);

          /**
           * Return the parent element unless this is the root element.
           * If so then return Element::Missing()
           *
           * @return
           *   Parent or Element::Missing()
           */
          Element GetParentOrNull();
          /**
           * Return the parent element unless this is the root element.
           * If so then throws ParentError
           *
           * @return
           *   Parent
           */
          Element GetParentOrThrow();

          /**
           * Get the value (as a string) contained in the specified attribute.
           * If it does not exist, return NULL
           * @param $name
           *   The name of the attribute to get
           * @return
           *   A pointer to a string containing the attribute value (or NULL
           *   on failure)
           */
          const std::string* GetAttributeOrNull(const std::string& name);
          /**
           * Get the value (as a string) contained in the specified attribute.
           * If it does not exist, throws AttributeError
           * @param $name
           *   The name of the attribute to get
           * @return
           *   A reference to a string containing the attribute value
           */
          const std::string& GetAttributeOrThrow(const std::string& name);

          /**
           * Get the value (as a string) contained in the specified attribute.
           * If it does not exist, return NULL
           * @param $name
           *   The name of the attribute to get
           * @return
           *   A pointer to a string containing the attribute value
           */

          /**
           * Get the value contained in the specified attribute. This function
           * will attempt to convert the string to the type of the second
           * argument using std::istream::operator>> (i.e. the standard library's
           * formatted input operator) so this can work for arbitrary types,
           * as long as you supply an implementation for:
           *   std::istream& operator>>(std::istream&, T&);
           *
           * This also returns the attribute as a pointer to the string value
           * of the attribute, which will be Element::Missing() if it does not
           * exist.
           *
           * If the attribute exists but does not convert correctly to type T,
           * then this will throw std::stringstream::failure.
           *
           * If the attribute exists, and converts to type T but does not use
           * the whole attribute string, then this will throw ParseError.
           *
           * @param $name
           *   Attribute to read and convert
           * @param $out
           *   Variable in which to store the converted attribute.
           * @return
           *   A pointer to a string containing the attribute value (or NULL
           *   if it does not exist)
           */
          template<class T>
          const std::string* GetAttributeOrNull(const std::string& name, T& out);

          /**
           * Get the value contained in the specified attribute. This function
           * will attempt to convert the string to the type of the second
           * argument using std::istream::operator>> (i.e. the standard library's
           * formatted input operator) so this can work for arbitrary types,
           * as long as you supply an implementation for:
           *   std::istream& operator>>(std::istream&, T&);
           *
           * This also returns the attribute as a reference to the string value
           * of the attribute.
           *
           * If the attribute does not exists, throw AttributeError.
           *
           * If the attribute exists but does not convert correctly to type T,
           * then this will throw std::stringstream::failure.
           *
           * If the attribute exists, and converts to type T but does not use
           * the whole attribute string, then this will throw ParseError.
           *
           * @param $name
           *   Attribute to read and convert
           * @param $out
           *   Variable in which to store the converted attribute.
           * @return
           *   A reference to a string containing the attribute value
           */
          template<class T>
          const std::string& GetAttributeOrThrow(const std::string& name, T& out);

          /**
           * Return a string giving a full path to the element.
           * @return
           */
          std::string GetPath() const;

        private:
          TiXmlElement* el;

          /**
           * Recursive function used by GetPath
           * @param el
           * @param ans
           */
          static void GetPathWorker(const TiXmlElement* el, std::string& ans);

          /**
           * Equality and inequality operators for
           * @param left
           * @param right
           * @return
           */
          friend bool operator==(const Element& left, const Element& right);
          friend bool operator!=(const Element& left, const Element& right);
      };
      bool operator==(const Element& left, const Element& right);
      bool operator!=(const Element& left, const Element& right);
      /** an abstraction for an XML document
       *
       * this class localises the dependency on an external XML library
       */
      class Document
      {
        public:
          /**
           * constructor
           *
           * @param $path
           *   the path to the XML file to be read by this object
           */
          Document(const std::string path);

          /** destructor */
          ~Document();
          Element GetRoot();
        private:
          TiXmlDocument* xmlDoc;
      };

      /**
       * Base class for XML errors. Should not be instantiated.
       */
      class XmlError : public hemelb::Exception
      {
        protected:
          //
          /**
           * Should not be instantiated, only its subclasses.
           *
           * Requires the Element where the error occured.
           *
           * @param element
           */
          XmlError(const Element& element);

        public:
          /**
           * Access the Element where the error occurred.
           * @return
           */
          inline const Element& GetNode() const
          {
            return elem;
          }
          /**
           * D'tor because the default one does not have the correct exception
           * specification.
           */
          virtual ~XmlError() throw ()
          {
          }
          /**
           * Supply an human readable error message.
           * @return
           */
          //virtual const char* what() const throw ();

        protected:
          const Element elem;
          // Full path to the element
          const std::string elemPath;
      };

      /**
       * Indicate that an element does not have a requested attribute
       */
      class AttributeError : public XmlError
      {
        public:
          AttributeError(const Element& n, const std::string& attr_);

          virtual ~AttributeError() throw ()
          {
          }

        private:
          const std::string attr;
      };

      class ParseError : public XmlError
      {
        public:
          ParseError(const Element& el, const std::string& attrName, const std::string& attrVal);
          virtual ~ParseError() throw ()
          {
          }

        private:
          const std::string name;
          const std::string val;
      };

      /**
       * Indicate that a requested element does not exist. Should not be used
       * directly; use ChildError, ParentError, or SiblingError.
       */
      class ElementError : public XmlError
      {
        protected:
          ElementError(const Element& n, const std::string& elemName);
        public:
          virtual ~ElementError() throw ()
          {
          }
        protected:
          const std::string elemName;
      };
      /**
       * Indicate that an element lacks the requested child
       */
      class ChildError : public ElementError
      {
        public:
          ChildError(const Element& elem, const std::string& subElemName);
      };

      /**
       * Indicate that an element lacks a parent.
       */
      class ParentError : public ElementError
      {
        public:
          ParentError(const Element& elem);
      };

      /**
       * Indicate that an element lacks the requested sibling.
       */
      class SiblingError : public ElementError
      {
        public:
          SiblingError(const Element& elem, const std::string& subElemName);
      };


      // Implement the template member functions declared above, now that the
      // declarations of the exceptions are available.
      template<class T>
      const std::string* Element::GetAttributeOrNull(const std::string& name, T& out)
      {
        const std::string* attrString = GetAttributeOrNull(name);
        if (attrString != NULL)
        {
          std::stringstream attrStream(*attrString, std::ios_base::in);

          attrStream >> out;
          if (attrStream.fail())
          {
            throw ParseError(*this, name, *attrString) << " error in extraction operator";
          }

          size_t pos = attrStream.tellg();
          if (pos != attrString->size())
          {
            throw ParseError(*this, name, *attrString) << " not all characters consumed";
          }
        }
        return attrString;
      }
      template<class T>
      const std::string& Element::GetAttributeOrThrow(const std::string& name, T& out)
      {
        const std::string* ans = GetAttributeOrNull(name, out);
        if (ans == NULL)
          throw AttributeError(*this, name);
        return *ans;
      }
    }
  }
}
#endif /* HEMELB_CONFIGURATION_XMLABSTRACTIONLAYER_H */
