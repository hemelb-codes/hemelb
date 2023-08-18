// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <tinyxml2.h>

#include "io/xml.h"

namespace hemelb::io::xml
{
    // Concept checks for the iterators
    static_assert(std::forward_iterator<UnnamedChildIterator>);
    static_assert(std::forward_iterator<NamedChildIterator>);
    static_assert(std::sentinel_for<ChildIteratorSentinel, NamedChildIterator>);
    static_assert(std::sentinel_for<ChildIteratorSentinel, UnnamedChildIterator>);

    using namespace tinyxml2;

    Document::Document()
    {
        xmlDoc = std::make_unique<XMLDocument>();
        // Associate our Document object with TinyXML's to allow
        // look up of the filename in error reporting.
        xmlDoc->SetUserData(this);
    }

    Document::Document(const std::filesystem::path& path) : Document()
    {
        LoadFile(path);
    }

    // Supply default implementations of destructor/move
    Document::~Document() = default;
    Document::Document(Document &&) = default;
    Document &Document::operator=(Document &&) = default;

    Element Document::GetRoot() const
    {
        return {xmlDoc->RootElement()};
    }

    Element Document::AddChild(char const* name) {
        auto el = xmlDoc->NewElement(name);
        xmlDoc->InsertEndChild(el);
        return el;
    }

    void Document::LoadFile(const std::filesystem::path& path) {
        filename = path;
        if (xmlDoc->LoadFile(path.c_str()) != tinyxml2::XML_SUCCESS) {
            throw ParseError(xmlDoc.get());
        }
    }

    void Document::LoadString(char const* data) {
        filename = "";
        if (xmlDoc->Parse(data) != XML_SUCCESS) {
            throw ParseError(xmlDoc.get());
        }
    }

    void Document::LoadString(const std::string &data) {
        LoadString(data.c_str());
    }

    void Document::SaveFile(const std::filesystem::path &path) const {
        if (xmlDoc->SaveFile(path.c_str()) != XML_SUCCESS) {
            throw (Exception() << "XML save error: " << xmlDoc->ErrorStr());
        }
    }

    Element::Element(XMLElement* el_) :
            el(el_)
    {
    }

    Element::~Element() = default;

    Element Element::Missing()
    {
        return {};
    }
    std::string_view Element::GetName() const
    {
        return el->Name();
    }

    int Element::GetLine() const
    {
        return el->GetLineNum();
    }

    Element Element::GetChildOrNull() const
    {
        return Element(el->FirstChildElement());
    }

    Element Element::GetChildOrNull(char const* name) const
    {
        return Element(el->FirstChildElement(name));
    }

    Element Element::GetChildOrThrow() const
    {
        auto ans = el->FirstChildElement();
        if (ans == nullptr)
            throw ChildError(*this, "of any name");

        return Element(ans);
    }
    Element Element::GetChildOrThrow(char const* name) const
    {
        auto ans = el->FirstChildElement(name);
        if (ans == nullptr)
            throw ChildError(*this, name);

        return Element(ans);
    }

    NamedIterationRange Element::Children(char const* name) const {
        return {*this, std::string{name}};
    }

    NamedIterationRange Element::Children(std::string name) const {
        return {*this, std::move(name)};
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

    Element Element::NextSiblingOrNull(char const* name) const
    {
        auto ans = el->NextSiblingElement(name);
        if (ans == nullptr)
            return nullptr;

        return Element(ans);
    }

    char const* Element::GetAttrOrNull(char const* attr) const {
        return el->Attribute(attr);
    }

    void Element::SetAttribute(const char *name, const char *value) {
        el->SetAttribute(name, value);
    }

    std::optional<std::string_view> Element::GetAttributeMaybe(char const* name) const
    {
        auto attr_p = el->Attribute(name);
        if (attr_p == nullptr) {
            return std::nullopt;
        } else {
            return std::make_optional<std::string_view>(attr_p);
        }
    }

    std::string_view Element::GetAttributeOrThrow(char const* name) const
    {
        const char* ans = el->Attribute(name);
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
    void Element::GetPathWorker(const XMLElement* el, std::ostringstream& ans)
    {
        const XMLNode* parent = el->Parent();
        const XMLElement* parentEl = parent->ToElement();
        if (parentEl != nullptr)
        {
            GetPathWorker(parentEl, ans);
        }
        else
        {
            // Documents owned by our Document class have a user data
            // pointer to the instance.
            const XMLDocument* doc = parent->ToDocument();
            auto heme_doc = static_cast<Document const*>(doc->GetUserData());
            auto fn = heme_doc->filename;
            if (fn != "")
            {
                ans << fn << ":";
            }
            else
            {
                ans << "?:";
            }
        }

        ans << "/" << el->Value() << "(" << el->GetLineNum() << ")";
    }

    Element::operator bool() const {
        return el != nullptr;
    }

    Element Element::AddChild(const char *name) {
        return {el->InsertNewChildElement(name)};
    }

    Element Element::CopyAsChild(const Element & source) {
        auto new_el = source.el->DeepClone(el->GetDocument());
        return {el->InsertEndChild(new_el)->ToElement()};
    }

    void Element::DeleteChild(const Element &elem_to_del) {
        el->DeleteChild(elem_to_del.el);
    }

    void Element::Delete() {
        GetParentOrThrow().DeleteChild(*this);
    }

    /**
     * Dereference
     * @return
     */
    ChildIterator::reference ChildIterator::operator*() const
    {
        return current;
    }

    /**
     * Dereference
     * @return
     */
    Element const* ChildIterator::operator->() const
    {
        return &current;
    }

    UnnamedChildIterator::UnnamedChildIterator(const Element &elem) {
        parent = elem;
        current = parent.GetChildOrNull();
    }

    UnnamedChildIterator::UnnamedChildIterator(const Element &elem, const Element &pos) {
        parent = elem;
        current = pos;
    }

    // Prefix increment
    UnnamedChildIterator& UnnamedChildIterator::operator++()
    {
        // increment and return the updated version
        current = current.NextSiblingOrNull();
        return *this;
    }

    // Postfix increment
    UnnamedChildIterator UnnamedChildIterator::operator++(int)
    {
        auto old = *this;
        // increment this
        ++*this;
        // return old
        return old;
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
        current = parent.GetChildOrNull(name.c_str());
    }
    NamedChildIterator::NamedChildIterator(const Element &elem, std::string subElemName, const Element &pos) :
            name(std::move(subElemName))
    {
        parent = elem;
        current = pos;
    }

    // Prefix increment
    NamedChildIterator& NamedChildIterator::operator++()
    {
        // increment and return the updated version
        current = current.NextSiblingOrNull(name.c_str());
        return *this;
    }

    // Postfix increment
    NamedChildIterator NamedChildIterator::operator++(int)
    {
        auto old = *this;
        // increment this
        ++*this;
        // return old
        return old;
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

    ChildIteratorSentinel NamedIterationRange::end() const {
        return {};
    }

    UnnamedChildIterator UnnamedIterationRange::begin() const {
        return {parent};
    }

    ChildIteratorSentinel UnnamedIterationRange::end() const {
        return {};
    }

    // XML exception base class
    XmlError::XmlError()
    {
        *this << "xml::";
    }

    ParseError::ParseError(const XMLDocument* node) : XmlError() {
        *this << "Error parsing XML. TinXML2 says: " << node->ErrorStr();
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

