
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <stdexcept>
#include "io/xml/XmlAbstractionLayer.h"
#include "util/fileutils.h"

#include "tinyxml.h"

namespace hemelb
{
  namespace io
  {
    namespace xml
    {
      Document::Document(const std::string path)
      {
        util::check_file(path.c_str());
        xmlDoc = new ::TiXmlDocument();
        xmlDoc->LoadFile(path);
      }
      Document::~Document()
      {
        delete xmlDoc;
        xmlDoc = NULL;
      }

      Element Document::GetRoot()
      {
        return Element(xmlDoc->RootElement());
      }

      Element::Element(TiXmlElement* el_) :
        el(el_)
      {
      }

      Element::~Element()
      {
      }
      const Element Element::Missing()
      {
        return Element(NULL);
      }
      const std::string& Element::GetName() const
      {
        return el->ValueStr();
      }

      int Element::GetLine() const
      {
        return el->Row();
      }

      Element Element::GetChildOrNull(const std::string& name)
      {
        TiXmlElement* ans = el->FirstChildElement(name);
        return Element(ans);
      }
      const Element Element::GetChildOrNull(const std::string& name) const
      {
        TiXmlElement* ans = el->FirstChildElement(name);
        return Element(ans);
      }

      Element Element::GetChildOrThrow(const std::string& name)
      {
        TiXmlElement* ans = el->FirstChildElement(name);
        if (ans == NULL)
          throw ChildError(*this, name);

        return Element(ans);
      }
      const Element Element::GetChildOrThrow(const std::string& name) const
      {
        TiXmlElement* ans = el->FirstChildElement(name);
        if (ans == NULL)
          throw ChildError(*this, name);

        return Element(ans);
      }

      ChildIterator Element::IterChildren(const std::string& name) const
      {
        return ChildIterator(*this, name);
      }

      Element Element::NextSiblingOrNull(const std::string name)
      {
        TiXmlElement* ans = el->NextSiblingElement(name);
        if (ans == NULL)
          return NULL;

        return Element(ans);
      }

      Element Element::NextSiblingOrThrow(const std::string name)
      {
        TiXmlElement* ans = el->NextSiblingElement(name);
        if (ans == NULL)
          throw SiblingError(*this, name);

        return Element(ans);
      }
      const std::string* Element::GetAttributeOrNull(const std::string& name) const
      {
        return el->Attribute(name);
      }
      const std::string& Element::GetAttributeOrThrow(const std::string& name) const
      {
        const std::string* ans = el->Attribute(name);
        if (ans == NULL)
          throw AttributeError(*this, name);

        return *ans;
      }

      Element Element::GetParentOrNull()
      {
        return Element(el->Parent()->ToElement());
      }
      Element Element::GetParentOrThrow()
      {
        TiXmlElement* parent = el->Parent()->ToElement();
        if (parent == NULL)
          throw ParentError(*this);
        return Element(parent);
      }

      bool operator==(const Element& left, const Element& right)
      {
        return (left.el == right.el);
      }
      bool operator!=(const Element& left, const Element& right)
      {
        return (left.el != right.el);
      }

      std::string Element::GetPath() const
      {
        std::ostringstream ans;
        GetPathWorker(el, ans);
        return ans.str();
      }
      void Element::GetPathWorker(const TiXmlElement* el, std::ostringstream& ans)
      {
        const TiXmlNode* parent = el->Parent();
        const TiXmlElement* parentEl = parent->ToElement();
        if (parentEl != NULL)
        {
          GetPathWorker(parentEl, ans);
        }
        else
        {
          const TiXmlDocument* doc = parent->ToDocument();
          if (doc != NULL)
          {
            ans << doc->Value() << ":";
          }
          else
          {
            ans << "?:";
          }
        }

        ans << "/" << el->Value() << "(" << el->Row() << ")";
      }

      /**
       * Default constructor
       */
      ChildIterator::ChildIterator() :
        parent(Element::Missing()), current(Element::Missing()), name()
      {
      }

      /**
       * Constructor that will iterate over subelements with the given name.
       * @param elem
       * @param subElemName
       */
      ChildIterator::ChildIterator(const Element& elem, const std::string& subElemName) :
        parent(elem), current(elem.GetChildOrNull(subElemName)), name(subElemName)
      {
      }

      /**
       * Copy constructor
       * @param other
       */
      ChildIterator::ChildIterator(const ChildIterator& other) :
        parent(other.parent), current(other.current), name(other.name)
      {
      }

      /**
       * Copy assignment
       * @param other
       * @return
       */
      ChildIterator& ChildIterator::operator=(const ChildIterator& other)
      {
        parent = other.parent;
        current = other.current;
        name = other.name;
        return *this;
      }

      /**
       * Equality comparable
       * @param
       * @return
       */
      bool operator==(const ChildIterator& a, const ChildIterator& b)
      {
        return (a.parent == b.parent) && (a.name == b.name) && (a.current == b.current);
      }

      /**
       * Inequality
       * @param
       * @return
       */
      bool operator!=(const ChildIterator& a, const ChildIterator& b)
      {
        return ! (a == b);
      }

      /**
       * Dereference
       * @return
       */
      ChildIterator::reference ChildIterator::operator*()
      {
        return current;
      }

      /**
       * Dereference
       * @return
       */
      ChildIterator::pointer ChildIterator::operator->()
      {
        return &current;
      }

      /**
       * Prefix increment
       * @return
       */
      ChildIterator& ChildIterator::operator++()
      {
        // increment and return the updated version
        current = current.NextSiblingOrNull(name);
        return *this;
      }

      /**
       * Postfix increment
       * @param
       * @return
       */
      ChildIterator ChildIterator::operator++(int n)
      {
        // increment and return the value pre-increment
        ChildIterator ans = *this;
        ++(*this);
        return ans;
      }

      bool ChildIterator::AtEnd() const
      {
        return current == Element::Missing();
      }

      // XML exception base class
      XmlError::XmlError(const Element& el) :
        elem(el), elemPath(el.GetPath())
      {
        *this << "xml::";
      }

      // Missing attribute
      AttributeError::AttributeError(const Element& n, const std::string& attr_) :
        XmlError(n), attr(attr_)
      {
        *this << "AttributeError: '" << elemPath << "' has no attribute '" << attr << "'";
      }

      // Attribute parsing error
      ParseError::ParseError(const Element& el, const std::string& attrName,
                             const std::string& attrVal) :
        XmlError(el), name(attrName), val(attrVal)
      {
        *this << "ParseError: '" << elemPath << "' Cannot convert attribute '" << name << "=\""
            << val << "\"'";
      }

      ElementError::ElementError(const Element& el, const std::string& elName) :
        XmlError(el), elemName(elName)
      {
      }
      ChildError::ChildError(const Element& elem, const std::string& subElemName) :
        ElementError(elem, subElemName)
      {
        *this << "ChildError: '" << elemPath << "' has no child '" << elemName << "'";
      }

      ParentError::ParentError(const Element& elem) :
        ElementError(elem, "")
      {
        *this << "ParentError: '" << elemPath << "' is root element";
      }

      SiblingError::SiblingError(const Element& elem, const std::string& subElemName) :
        ElementError(elem, subElemName)
      {
        *this << "SiblingError: '" << elemPath << "' has no later sibling '" << elemName << "'";
      }

    }
  }
}

