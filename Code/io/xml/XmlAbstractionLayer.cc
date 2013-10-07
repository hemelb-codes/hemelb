// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 
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
      Element Element::Missing()
      {
        return Element(NULL);
      }
      const std::string& Element::GetName() const
      {
        return el->ValueStr();
      }

      Element Element::GetChildOrNull(const std::string& name)
      {
        TiXmlElement* ans = el->FirstChildElement(name);
        if (ans == NULL)
          return NULL;
        return Element(ans);
      }

      Element Element::GetChildOrThrow(const std::string& name)
      {
        TiXmlElement* ans = el->FirstChildElement(name);
        if (ans == NULL)
          throw ChildError(*this, name);

        return Element(ans);
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
      const std::string* Element::GetAttributeOrNull(const std::string& name)
      {
        return el->Attribute(name);
      }
      const std::string& Element::GetAttributeOrThrow(const std::string& name)
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
        std::string ans;
        GetPathWorker(el, ans);
        return ans;
      }
      void Element::GetPathWorker(const TiXmlElement* el, std::string& ans)
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
            ans.append(doc->Value());
            ans.append(":");
          }
          else
          {
            ans.append("?");
          }
        }

        ans.append("/");
        ans.append(el->Value());
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

