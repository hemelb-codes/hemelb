
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_XML_XMLABSTRACTIONLAYER_H
#define HEMELB_IO_XML_XMLABSTRACTIONLAYER_H

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <limits>
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
      // Forward declare
      class ChildIterator;

      class Element
      {
        public:
          static const Element Missing();

          Element(TiXmlElement* el);
          ~Element();

          /**
           * Get the name of the element
           * @return
           */
          const std::string& GetName() const;
          void GetName(const std::string& str);
          /**
           * Get the line number
           * @return
           */
          int GetLine() const;
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
          const Element GetChildOrNull(const std::string& name) const;
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
          const Element GetChildOrThrow(const std::string& name) const;

          /**
           * Return iterator over children with the specified name.
           * @param name
           * @return
           */
          ChildIterator IterChildren(const std::string& name) const;

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
          const std::string* GetAttributeOrNull(const std::string& name) const;
          /**
           * Get the value (as a string) contained in the specified attribute.
           * If it does not exist, throws AttributeError
           * @param $name
           *   The name of the attribute to get
           * @return
           *   A reference to a string containing the attribute value
           */
          const std::string& GetAttributeOrThrow(const std::string& name) const;

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
          const std::string* GetAttributeOrNull(const std::string& name, T& out) const;

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
          const std::string& GetAttributeOrThrow(const std::string& name, T& out) const;

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
          static void GetPathWorker(const TiXmlElement* el, std::ostringstream& ans);

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

      /**
       * Want to model the ForwardIterator concept
       */
      class ChildIterator : public std::iterator<std::forward_iterator_tag, Element>
      {
        public:
          /**
           * Default constructor
           */
          ChildIterator();

          /**
           * Constructor that will iterate over subelements with the given name.
           * @param elem
           * @param subElemName
           */
          ChildIterator(const Element& elem, const std::string& subElemName);

          /**
           * Copy constructor
           * @param other
           */
          ChildIterator(const ChildIterator& other);

          // Default is fine
          // /**
          //  * Destructor.
          //  */
          // ~ChildIterator();

          /**
           * Copy assignment
           * @param other
           * @return
           */
          ChildIterator& operator=(const ChildIterator& other);


          /**
           * Dereference
           * @return
           */
          reference operator*();

          /**
           * Dereference
           * @return
           */
          pointer operator->();

          /**
           * Prefix increment
           * @return
           */
          ChildIterator& operator++();
          /**
           * Postfix increment
           * @param
           * @return
           */
          ChildIterator operator++(int);

          bool AtEnd() const;

        private:
          Element parent;
          Element current;
          std::string name;
          friend bool operator==(const ChildIterator& a, const ChildIterator& b);
      };
      /**
       * Equality comparable
       * @param
       * @return
       */
      bool operator==(const ChildIterator& a, const ChildIterator& b);

      /**
       * Inequality
       * @param
       * @return
       */
      bool operator!=(const ChildIterator& a, const ChildIterator& b);

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
          template<typename T>
          XmlError& operator<<(const T& t)
          {
            return static_cast<XmlError&> (Exception::operator<<(t));
          }

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
          template<typename T>
          AttributeError& operator<<(const T& t)
          {
            return static_cast<AttributeError&> (XmlError::operator<<(t));
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

          template<typename T>
          ParseError& operator<<(const T& t)
          {
            return static_cast<ParseError&> (XmlError::operator<<(t));
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
          template<typename T>
          ElementError& operator<<(const T& t)
          {
            return static_cast<ElementError&> (XmlError::operator<<(t));
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
          template<typename T>
          ChildError& operator<<(const T& t)
          {
            return static_cast<ChildError&> (ElementError::operator<<(t));
          }

      };

      /**
       * Indicate that an element lacks a parent.
       */
      class ParentError : public ElementError
      {
        public:
          ParentError(const Element& elem);
          template<typename T>
          ParentError& operator<<(const T& t)
          {
            return static_cast<ParentError&> (ElementError::operator<<(t));
          }

      };

      /**
       * Indicate that an element lacks the requested sibling.
       */
      class SiblingError : public ElementError
      {
        public:
          SiblingError(const Element& elem, const std::string& subElemName);
          template<typename T>
          SiblingError& operator<<(const T& t)
          {
            return static_cast<SiblingError&> (ElementError::operator<<(t));
          }

      };

      // Implement the template member functions declared above, now that the
      // declarations of the exceptions are available.
      template<class T>
      const std::string* Element::GetAttributeOrNull(const std::string& name, T& out) const
      {
        const std::string* attrString = GetAttributeOrNull(name);
        if (attrString != NULL)
        {
          /*
           * So, basically parsing of unsigned values varies across platforms
           * such that on OS X 10.8
           *   unsigned ans
           *   stream = istringstream("-1");
           *   stream >> ans;
           * will NOT set ans but will set stream's failbit.
           *
           * On various Linux boxes, this will set ans == 2**32 - 1 and NOT set
           * the fail bit.
           *
           * Here, for types numeric_limits knows are unsigned we explicitly
           * check for a "-" and throw a suitable error.
           */
          if (std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_signed)
          {
            if (attrString->at(0) == '-')
            {
              throw ParseError(*this, name, *attrString)
                  << " attempt to convert negative number to unsigned type";
            }
          }
          std::stringstream attrStream(*attrString, std::ios_base::in);

          // Don't skip whitespace as that could indicate a malformed value
          attrStream >> std::noskipws;

          attrStream >> out;
          if (attrStream.fail())
          {
            throw ParseError(*this, name, *attrString) << " error in extraction operator";
          }
          bool eof = attrStream.eof();
          std::stringstream::pos_type pos = attrStream.tellg();
          if (!eof && pos != int(attrString->size()))
          {
            throw ParseError(*this, name, *attrString) << " not all characters consumed";
          }
        }
        return attrString;
      }
      template<class T>
      const std::string& Element::GetAttributeOrThrow(const std::string& name, T& out) const
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
