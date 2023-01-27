// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include <stdexcept>

#include "tinyxml.h"

#include "io/xml.h"

namespace hemelb::io::xml
{
      Document::Document()
      {
        xmlDoc = std::make_unique<::TiXmlDocument>();
      }

      Document::Document(const std::filesystem::path& path) : Document()
      {
          LoadFile(path);
      }

      Document::~Document() {
      }

      Element Document::GetRoot()
      {
        return {xmlDoc->RootElement()};
      }

      void Document::LoadFile(const std::filesystem::path& path) {
          if (!xmlDoc->LoadFile(path.c_str())) {
              throw ParseError(xmlDoc.get());
          }
      }

    void Document::LoadString(std::string_view data) {
        if (xmlDoc->Parse(data.data()) == nullptr) {
            throw ParseError(xmlDoc.get());
        }
    }

      Element::Element(TiXmlElement const* el_) :
          el(el_)
      {
      }

      Element::~Element() {
      }

      const Element Element::Missing()
      {
        return Element(nullptr);
      }
      const std::string& Element::GetName() const
      {
        return el->ValueStr();
      }

      int Element::GetLine() const
      {
        return el->Row();
      }

      Element Element::GetChildOrNull(std::string_view name) const
      {
        return Element(el->FirstChildElement(name.data()));
      }

      Element Element::GetChildOrThrow(std::string_view name) const
      {
        auto ans = el->FirstChildElement(name.data());
        if (ans == nullptr)
          throw ChildError(*this, name);

        return Element(ans);
      }

      ChildIterator Element::IterChildren(std::string_view name) const
      {
        return ChildIterator(*this, std::string{name});
      }

      Element Element::NextSiblingOrNull(std::string_view name) const
      {
        auto ans = el->NextSiblingElement(name.data());
        if (ans == nullptr)
          return nullptr;

        return Element(ans);
      }

      Element Element::NextSiblingOrThrow(std::string_view name) const
      {
        auto ans = el->NextSiblingElement(name.data());
        if (ans == nullptr)
          throw SiblingError(*this, name);

        return Element(ans);
      }
      char const* Element::GetAttrOrNull(std::string_view attr) const {
          return el->Attribute(attr.data());
      }

    std::optional<std::string_view> Element::GetAttributeMaybe(std::string_view name) const
    {
        auto attr_p = el->Attribute(name.data());
        if (attr_p == nullptr) {
            return std::nullopt;
        } else {
            return std::make_optional<std::string_view>(attr_p);
        }
    }

    std::string_view Element::GetAttributeOrThrow(std::string_view name) const
    {
        const char* ans = el->Attribute(name.data());
        if (ans == nullptr)
            throw AttributeError(*this, name);

        return {ans};
    }

      Element Element::GetParentOrNull() const
      {
        return Element(el->Parent()->ToElement());
      }
      Element Element::GetParentOrThrow() const
      {
        auto parent = el->Parent()->ToElement();
        if (parent == nullptr)
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
        if (parentEl != nullptr)
        {
          GetPathWorker(parentEl, ans);
        }
        else
        {
          const TiXmlDocument* doc = parent->ToDocument();
          if (doc != nullptr)
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

      Element::operator bool() const {
	return el != nullptr;
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
      ChildIterator::ChildIterator(const Element& elem, std::string subElemName) :
          parent(elem), current(elem.GetChildOrNull(subElemName)), name(std::move(subElemName))
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
        ++ (*this);
        return ans;
      }

      bool ChildIterator::AtEnd() const
      {
        return current == Element::Missing();
      }

      // XML exception base class
      XmlError::XmlError()
      {
        *this << "xml::";
      }

      ParseError::ParseError(const TiXmlDocument* node) : XmlError() {
	*this << "Error parsing XML. TiXml says: " << node->ErrorDesc();
      }

      // Missing attribute
      AttributeError::AttributeError(const Element& n, std::string_view attr) :
        XmlError()
      {
        *this << "AttributeError: '" << n.GetPath() << "' has no attribute '" << attr << "'";
      }

      // Attribute parsing error
      DeserialisationError::DeserialisationError(const Element& el, std::string_view attrName,
                             std::string_view attrVal) :
        XmlError()
      {
        *this << "ParseError: '" << el.GetPath() << "' Cannot convert attribute '" << attrName << "=\""
            << attrVal << "\"'";
      }

      ElementError::ElementError(const Element& el) :
        XmlError(), elemPath(el.GetPath())
      {
      }
      ChildError::ChildError(const Element& elem, std::string_view subElemName) :
        ElementError(elem)
      {
        *this << "ChildError: '" << elemPath << "' has no child '" << subElemName << "'";
      }

      ParentError::ParentError(const Element& elem) :
        ElementError(elem)
      {
        *this << "ParentError: '" << elemPath << "' is root element";
      }

      SiblingError::SiblingError(const Element& elem, std::string_view sibElemName) :
        ElementError(elem)
      {
        *this << "SiblingError: '" << elemPath << "' has no later sibling '" << sibElemName << "'";
      }

}

