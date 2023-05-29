// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_CLONE_PTR_H
#define HEMELB_UTIL_CLONE_PTR_H

#include <memory>
#include <type_traits>
#include "util/concepts.h"

namespace hemelb::util {

    template <typename T>
    concept cloneable = requires(T const& t) {
        // Objects must have a clone member function, that returns a
        // pointer-to-base-class (which could be T)
        { t.clone() } -> util::pointer_base_of<T*>;
    };

    // A unique_ptr that, when copied, will return a new instance owning
    // the result of calling `clone` on the current object.
    //
    // Worth stressing that each instance uniquely owns the pointer-to
    // object.
    //
    // This is helpful for polymorphic class hierarchies.
    //
    // If clone() is virtual, then so should be the destructor.
    template <cloneable T>
    class clone_ptr {
        // Get the type of calling clone() on a const T
        using clone_return_t = decltype(std::declval<std::add_const_t<T>>().clone());
        using base_type = std::remove_pointer_t<clone_return_t>;

        static T* cast(base_type* b) {
            if constexpr (std::is_same_v<T, base_type>) {
                return b;
            } else {
                return dynamic_cast<T*>(b);
            }
        }

        // Any other instantiation is a friend.
        template <cloneable U>
        friend class clone_ptr;
    public:
        using unique_ptr = std::unique_ptr<T>;

        using pointer = typename unique_ptr::pointer;
        using element_type = typename unique_ptr::element_type;

        // Repeat a minimal API of unique_ptr

        // Default c'tor - null
        clone_ptr() noexcept: m_ptr{} {
        }
        // From a pointer to T
        explicit clone_ptr(pointer p) noexcept : m_ptr{p} {
        }
        // Move from this instantiation.
        clone_ptr(clone_ptr&& u ) noexcept : m_ptr{std::move(u.m_ptr)} {
        }
        // Move from any type instantation that unique_ptr allows.
        template <cloneable U>
        clone_ptr(clone_ptr<U>&& u) noexcept : m_ptr{std::move(u.m_ptr)} {
        }

        // Copy construct by calling clone on object
        //
        // Note that for a class hierarchy, if we are holding a derived
        // type, then T::clone should probably return pointer-to-base.
        // i.e.
        //
        clone_ptr(clone_ptr const& other) : m_ptr{cast(other->clone())} {
        }

        // And for other types with have the same return type of clone().
        template <cloneable U>
        requires std::same_as<base_type, typename clone_ptr<U>::base_type>
        clone_ptr(clone_ptr<U> const& other) : m_ptr{cast(other->clone())} {
        }

        // Move assign
        clone_ptr& operator=(clone_ptr&& other) noexcept {
            m_ptr = std::move(other.m_ptr);
            return *this;
        }

        // Move assign from any allowed by unique_ptr
        template <cloneable U>
        clone_ptr& operator=(clone_ptr<U>&& u) noexcept {
            m_ptr = std::move(u.m_ptr);
            return *this;
        }

        // Assign nullptr
        clone_ptr& operator=(std::nullptr_t) noexcept {
            m_ptr = nullptr;
            return *this;
        }

        // Copy assign by cloning other object
        clone_ptr& operator=(clone_ptr const& other) {
            m_ptr.reset(cast(other->clone()));
            return *this;
        }

        // And for compatible types
        template <cloneable U>
        requires std::same_as<base_type, typename clone_ptr<U>::base_type>
        clone_ptr& operator=(clone_ptr<U> const& other) {
            m_ptr.reset(cast(other->clone()));
            return *this;
        }

        // Satisfy the rule of 5
        ~clone_ptr() = default;

        // Modifiers
        pointer release() noexcept {
            return m_ptr.release();
        }
        void reset(pointer ptr ) noexcept {
            m_ptr.reset(ptr);
        }
        void swap( clone_ptr& other ) noexcept {
            swap(other.m_ptr);
        }

        // Observers
        pointer get() const noexcept {
            return m_ptr.get();
        }
        explicit operator bool() const noexcept {
            return (bool)m_ptr;
        }

        // Pointer access
        typename std::add_lvalue_reference<T>::type
        operator*() const {
            return *m_ptr;
        }
        pointer operator->() const noexcept {
            return m_ptr.get();
        }

        // Return a clone_ptr owning the result of calling clone on our
        // object.
        clone_ptr clone() const {
            return clone_ptr{m_ptr->clone()};
        }

    private:
        template <cloneable L, cloneable R>
        friend bool operator==(clone_ptr<L> const& lhs, clone_ptr<R> const& rhs);
        template <cloneable L, cloneable R>
        friend bool operator!=(clone_ptr<L> const& lhs, clone_ptr<R> const& rhs);
        unique_ptr m_ptr;
    };

    template <cloneable L, cloneable R>
    bool operator==(clone_ptr<L> const& lhs, clone_ptr<R> const& rhs) {
        return lhs.m_ptr == rhs.m_ptr;
    }
    template <cloneable L, cloneable R>
    bool operator!=(clone_ptr<L> const& lhs, clone_ptr<R> const& rhs) {
        return !(lhs == rhs);
    }

    template <cloneable T>
    bool operator!=(clone_ptr<T> const& v, std::nullptr_t) {
        return bool(v);
    }
    template <cloneable T>
    bool operator!=(std::nullptr_t, clone_ptr<T> const& v) {
        return bool(v);
    }
    template <cloneable T>
    bool operator==(clone_ptr<T> const& v, std::nullptr_t) {
        return !v;
    }
    template <cloneable T>
    bool operator==(std::nullptr_t, clone_ptr<T> const& v) {
        return !v;
    }

    // Helper like std::make_unique
    template <cloneable T, typename... Args>
    clone_ptr<T> make_clone_ptr(Args&&... args) {
        return clone_ptr<T>{new T{std::forward<Args>(args)...}};
    }

    // dynamic_cast (move)
    template <cloneable T, cloneable U>
    clone_ptr<T> clone_dynamic_cast(clone_ptr<U>&& p) {
        // Object still owned by p
        T* derived = dynamic_cast<T*>(p.get());
        if (derived) {
            // We own the object now
            p.release();
            // Pass to return val
            return clone_ptr<T>{derived};
        } else {
            return {};
        }
    }

    // dynamic_cast (clone)
    template <cloneable T, cloneable U>
    clone_ptr<T> clone_dynamic_cast(clone_ptr<U> const& p) {
        // Object still owned by p
        T* derived = dynamic_cast<T*>(p.get());
        if (derived) {
            return clone_ptr<T>{static_cast<T*>(derived->clone())};
        } else {
            return {};
        }
    }
}

#endif
