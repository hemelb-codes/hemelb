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

    Document::~Document() = default;

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

    Element::~Element() = default;

    Element Element::Missing()
    {
        return {};
    }
    const std::string& Element::GetName() const
    {
        return el->ValueStr();
    }

    int Element::GetLine() const
    {
        return el->Row();
    }

    Element Element::GetChildOrNull() const
    {
        return Element(el->FirstChildElement());
    }

    Element Element::GetChildOrNull(std::string_view name) const
    {
        return Element(el->FirstChildElement(name.data()));
    }

    Element Element::GetChildOrThrow() const
    {
        auto ans = el->FirstChildElement();
        if (ans == nullptr)
            throw ChildError(*this, "of any name");

        return Element(ans);
    }
    Element Element::GetChildOrThrow(std::string_view name) const
    {
        auto ans = el->FirstChildElement(name.data());
        if (ans == nullptr)
            throw ChildError(*this, name);

        return Element(ans);
    }

    NamedChildIterator Element::IterChildren(std::string_view name) const
    {
        return {*this, std::string{name}};
    }

    NamedIterationRange Element::Children(std::string_view name) const {
        return {*this, std::string{name}};
    }

    UnnamedChildIterator Element::IterChildren() const {
        return {*this};
    }

    UnnamedIterationRange Element::Children() const {
        return {*this};
    }

    Element Element::NextSiblingOrNull() const
    {
        auto ans = el->NextSiblingElement();
        if (ans == nullptr)
            return nullptr;

        return Element(ans);
    }

    Element Element::NextSiblingOrNull(std::string_view name) const
    {
        auto ans = el->NextSiblingElement(name.data());
        if (ans == nullptr)
            return nullptr;

        return Element(ans);
    }

    Element Element::NextSiblingOrThrow() const
    {
        auto ans = el->NextSiblingElement();
        if (ans == nullptr)
            throw SiblingError(*this, "of any name");

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
        next();
        return *this;
    }

    UnnamedChildIterator::UnnamedChildIterator(const Element &elem) {
        parent = elem;
        start();
    }

    UnnamedChildIterator::UnnamedChildIterator(const Element &elem, const Element &pos) {
        parent = elem;
        current = pos;
    }

    void UnnamedChildIterator::start() {
        current = parent.GetChildOrNull();
    }
    void UnnamedChildIterator::next() {
        current = current.NextSiblingOrNull();
    }

    /**
     * Constructor that will iterate over subelements with the given name.
     * @param elem
     * @param subElemName
     */
    NamedChildIterator::NamedChildIterator(const Element& elem, std::string subElemName) :
            name(std::move(subElemName))
    {
        parent = elem;
        start();
    }
    NamedChildIterator::NamedChildIterator(const Element &elem, std::string subElemName, const Element &pos) :
            name(std::move(subElemName))
    {
        parent = elem;
        current = pos;
    }

    void NamedChildIterator::start() {
        current = parent.GetChildOrNull(name);
    }
    void NamedChildIterator::next() {
        current = current.NextSiblingOrNull(name);
    }

    /**
     * Equality comparable
     * @param
     * @return
     */
    bool operator==(const NamedChildIterator& a, const NamedChildIterator& b)
    {
        return (a.parent == b.parent) && (a.name == b.name) && (a.current == b.current);
    }

    /**
     * Inequality
     * @param
     * @return
     */
    bool operator!=(const NamedChildIterator& a, const NamedChildIterator& b)
    {
        return ! (a == b);
    }

    bool operator==(const UnnamedChildIterator& a, const UnnamedChildIterator& b)
    {
        return (a.parent == b.parent) && (a.current == b.current);
    }

    bool operator!=(const UnnamedChildIterator& a, const UnnamedChildIterator& b) {
        return !(a == b);
    }

    NamedChildIterator NamedIterationRange::begin() const {
        return {parent, name};
    }

    NamedChildIterator NamedIterationRange::end() const {
        return {parent, name, Element::Missing()};
    }

    UnnamedChildIterator UnnamedIterationRange::begin() const {
        return {parent};
    }

    UnnamedChildIterator UnnamedIterationRange::end() const {
        return {parent, Element::Missing()};
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

