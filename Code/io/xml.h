// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_IO_XML_H
#define HEMELB_IO_XML_H

#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <sstream>
#include <limits>
#include <memory>
#include <optional>
#include <filesystem>

#include "Exception.h"
#include "util/traits.h"

// Forward declare the TinyXML types needed.
class TiXmlDocument;
class TiXmlElement;

namespace hemelb::io::xml
{
    // Forward declare
    class NamedChildIterator;
    class UnnamedChildIterator;
    struct NamedIterationRange;
    struct UnnamedIterationRange;

    // Represent an element in an XML file, or no element (if equal to
    // the result of Missing());
    //
    // Exposes a basic monadic interface.
    class Element
    {
        TiXmlElement const* el = nullptr;

    public:
        static Element Missing();

        Element() = default;
        Element(TiXmlElement const* el);
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
         * Gets the first child element, whatever its name.
         * @return element or missing if none
         */
        Element GetChildOrNull() const;

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
        Element GetChildOrNull(std::string_view name) const;

        /**
         * Gets the first child element whatever its name or throw
         * ChildError if it does not exist.
         *
         * @param $name
         *   the name of the child element to return
         *
         * @return
         *   returns the child element
         */
        Element GetChildOrThrow() const;

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
        Element GetChildOrThrow(std::string_view name) const;

        /**
         * Return iterator over children with the specified name.
         * @param name
         * @return
         */
        NamedChildIterator IterChildren(std::string_view name) const;
        NamedIterationRange Children(std::string_view name) const;

        /**
         * Return iterator over all children.
         * @param name
         * @return
         */
        UnnamedChildIterator IterChildren() const;
        UnnamedIterationRange Children() const;

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
        Element NextSiblingOrNull() const;

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
        Element NextSiblingOrNull(std::string_view name) const;

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
        Element NextSiblingOrThrow() const;

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
        Element NextSiblingOrThrow(std::string_view name) const;

        /**
         * Return the parent element unless this is the root element.
         * If so then return Element::Missing()
         *
         * @return
         *   Parent or Element::Missing()
         */
        Element GetParentOrNull() const;
        /**
         * Return the parent element unless this is the root element.
         * If so then throws ParentError
         *
         * @return
         *   Parent
         */
        Element GetParentOrThrow() const;

        /**
         * Get the value (as a string) contained in the specified attribute.
         * If it does not exist, return std::nullopt
         * @param $name
         *   The name of the attribute to get
         * @return
         *   An optional string containing the attribute value (or nullopt
         *   on failure)
         */
        std::optional<std::string_view> GetAttributeMaybe(std::string_view name) const;
        /**
         * Get the value (as a string) contained in the specified attribute.
         * If it does not exist, throws AttributeError
         * @param $name
         *   The name of the attribute to get
         * @return
         *   A reference to a string containing the attribute value
         */
        std::string_view GetAttributeOrThrow(std::string_view name) const;

        /**
         * Get the value contained in the specified attribute. This function
         * will attempt to convert the string to the type of templace
         * argument using std::istream::operator>> (i.e. the standard library's
         * formatted input operator) so this can work for arbitrary types,
         * as long as you supply an implementation for:
         *   std::istream& operator>>(std::istream&, T&);
         *
         * If the attribute is missing, you will get std::nullopt.
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
         *   bool indicating if the attribute was present.
         */
        template<class T>
        std::optional<T> GetAttributeMaybe(std::string_view name) const;

        template<class T>
        T GetAttributeOrThrow(std::string_view name) const;

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
         * the whole attribute string, then this will throw DeserialisationError.
         *
         * @param $name
         *   Attribute to read and convert
         * @param $out
         *   Variable in which to store the converted attribute.
         * @return
         *   A reference to a string containing the attribute value
         */
        template<class T>
        void GetAttributeOrThrow(std::string_view name, T& out) const;

        /**
         * Return a string giving a full path to the element.
         * @return
         */
        std::string GetPath() const;

        operator bool() const;


        // Returns the result of invocation of f on an Element that does
        // contain a value. Otherwise, returns std::nullopt. f must return
        // Element or a specialisation of std::optional.
        template <typename F>
        auto and_then(F&& f) const {
            using R = std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<F, const Element&>>>;
            static_assert(std::disjunction_v<
                    std::is_same<R, Element>,
                    util::is_optional<R>
            >);
            if (*this) {
                return std::invoke(std::forward<F>(f), *this);
            } else {
                if constexpr (std::is_same_v<R, Element>) {
                    return Missing();
                } else {
                    return R{};
                }
            }
        }

        // Returns an std::optional containing the result of invoking f on
        // an Element that contains a value. Otherwise returns
        // std::nullopt.
        template <typename F>
        auto transform(F&& f) const {
            using R = std::remove_cv_t<std::remove_reference_t<std::invoke_result_t<F, const Element&>>>;
            if (*this) {
                return std::optional<R>{std::invoke(std::forward<F>(f), *this)};
            } else {
                return std::optional<R>{};
            }
        }


    private:

        // For use in templates to hide TinyXML API
        char const* GetAttrOrNull(std::string_view attr) const;

        // Convert a string to a value
        template <typename T>
        void StringToVal(std::string_view name, char const* s, T& val) const;

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

    /**
     * Want to model the ForwardIterator concept
     */
    class ChildIterator
    {
    public:
        using value_type = Element;
        using reference = Element&;
        using pointer = Element*;
        using iterator_category = std::forward_iterator_tag;

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

    protected:
        virtual void start() = 0;
        virtual void next() = 0;

        Element parent = Element::Missing();
        Element current = Element::Missing();
    };

    class UnnamedChildIterator final : public ChildIterator {
        void start() override;
        void next() override;
    public:
        UnnamedChildIterator() = default;
        UnnamedChildIterator(const Element& elem);
        UnnamedChildIterator(const Element& elem, const Element& pos);

        friend bool operator==(const UnnamedChildIterator& a, const UnnamedChildIterator& b);
    };

    class NamedChildIterator final : public ChildIterator {
        std::string name;
        void start() override;
        void next() override;
    public:
        NamedChildIterator() = default;
        NamedChildIterator(const Element& elem, std::string subElemName);
        NamedChildIterator(const Element& elem, std::string subElemName, const Element& pos);
        friend bool operator==(const NamedChildIterator& a, const NamedChildIterator& b);
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

    struct NamedIterationRange {
        Element parent;
        std::string name;

        NamedChildIterator begin() const;
        NamedChildIterator end() const;
    };
    struct UnnamedIterationRange {
        Element parent;
        std::string name;

        UnnamedChildIterator begin() const;
        UnnamedChildIterator end() const;
    };

    /** an abstraction for an XML document
     *
     * this class localises the dependency on an external XML library
     */
    class Document
    {
    public:
        // Create an empty document
        Document();
        /**
         * constructor
         *
         * @param $path
         *   the path to the XML file to be read by this object
         */
        Document(const std::filesystem::path& path);

        /** destructor - needed to avoid polluting all users of this with TinyXML */
        ~Document();
        Element GetRoot();

        // Load some XML from a file
        void LoadFile(const std::filesystem::path& path);

        // Load some XML from a string
        void LoadString(std::string_view path);

    private:
        std::unique_ptr<TiXmlDocument> xmlDoc;
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
        XmlError();

    public:
        /**
         * Supply an human readable error message.
         * @return
         */
        template<typename T>
        XmlError& operator<<(const T& t)
        {
            return static_cast<XmlError&>(Exception::operator<<(t));
        }

    };

    // Indicate that some XML failed to parse
    class ParseError : public XmlError
    {
    public:
        ParseError(const TiXmlDocument*);
    };

    /**
     * Indicate that an element does not have a requested attribute
     */
    class AttributeError : public XmlError
    {
    public:
        AttributeError(const Element& n, std::string_view attr_);

        template<typename T>
        AttributeError& operator<<(const T& t)
        {
            return static_cast<AttributeError&>(XmlError::operator<<(t));
        }
    };

    class DeserialisationError : public XmlError
    {
    public:
        DeserialisationError(const Element& el, std::string_view attrName, std::string_view attrVal);

        template<typename T>
        DeserialisationError& operator<<(const T& t)
        {
            return static_cast<DeserialisationError&> (XmlError::operator<<(t));
        }
    };

    /**
     * Indicate that a requested element does not exist. Should not be used
     * directly; use ChildError, ParentError, or SiblingError.
     */
    class ElementError : public XmlError
    {
    protected:
        ElementError(const Element& el);
    public:
        template<typename T>
        ElementError& operator<<(const T& t)
        {
            return static_cast<ElementError&>(XmlError::operator<<(t));
        }

    protected:
        std::string elemPath;
    };
    /**
     * Indicate that an element lacks the requested child
     */
    class ChildError : public ElementError
    {
    public:
        ChildError(const Element& elem, std::string_view subElemName);
        template<typename T>
        ChildError& operator<<(const T& t)
        {
            return static_cast<ChildError&>(ElementError::operator<<(t));
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
            return static_cast<ParentError&>(ElementError::operator<<(t));
        }

    };

    /**
     * Indicate that an element lacks the requested sibling.
     */
    class SiblingError : public ElementError
    {
    public:
        SiblingError(const Element& elem, std::string_view subElemName);
        template<typename T>
        SiblingError& operator<<(const T& t)
        {
            return static_cast<SiblingError&>(ElementError::operator<<(t));
        }

    };

    template<class T>
    void Element::StringToVal(std::string_view name, char const* s, T& out) const
    {
        const auto N = std::strlen(s);
        if (N == 0)
            throw DeserialisationError(*this, name, s) << "Zero length string for value";

        // So, basically parsing of unsigned values varies across platforms
        // such that on OS X 10.8
        //   unsigned ans
        //   stream = istringstream("-1");
        //   stream >> ans;
        // will NOT set ans but will set stream's failbit.
        //
        // On various Linux boxes, this will set ans == 2**32 - 1 and NOT set
        // the fail bit.
        //
        // Here, for types numeric_limits knows are unsigned we explicitly
        // check for a "-" and throw a suitable error.
        if constexpr (std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_signed)
        {
            if (s[0] == '-')
            {
                throw DeserialisationError(*this, name, s)
                        << " attempt to convert negative number to unsigned type";
            }
        }
        std::stringstream attrStream(s, std::ios_base::in);

        // Don't skip whitespace as that could indicate a malformed value
        attrStream >> std::noskipws;

        attrStream >> out;
        if (attrStream.fail())
        {
            throw DeserialisationError(*this, name, s) << " error in extraction operator";
        }
        bool eof = attrStream.eof();
        if (!eof && attrStream.tellg() != std::streampos(N))
        {
            throw DeserialisationError(*this, name, s) << " not all characters consumed";
        }
    }

    template<class T>
    std::optional<T> Element::GetAttributeMaybe(std::string_view name) const
    {
        const char* attrString = GetAttrOrNull(name);
        if (attrString == nullptr) {
            return std::nullopt;
        } else {
            // Requires T is  default constructible
            std::optional<T> ans{std::in_place};
            StringToVal(name, attrString, *ans);
            return ans;
        }
    }

    template<class T>
    void Element::GetAttributeOrThrow(std::string_view name, T& out) const
    {
        const char* attrString = GetAttrOrNull(name);
        if (attrString == nullptr)
            throw AttributeError(*this, name);
        StringToVal(name, attrString, out);
    }

    template<class T>
    T Element::GetAttributeOrThrow(std::string_view name) const {
        T out;
        GetAttributeOrThrow(name, out);
        return out;
    }

}
#endif  // HEMELB_IO_XML_H
