// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_UTIL_CLONE_PTR_H
#define HEMELB_UTIL_CLONE_PTR_H

#include <memory>
#include <type_traits>

namespace hemelb::util {

  // A unique_ptr that, when copied, will return a new instance owning
  // the result of calling `clone` on the current object.
  //
  // Worth stressing that each instance uniquely owns the pointer-to
  // object.
  //
  // This is helpful for polymorphic class hierarchies.
  //
  // If clone() is virtual, then so should be the destructor.
  template <typename T>
  class clone_ptr {
    // Get the type of calling clone() on a const T
    using clone_return_t = decltype(std::declval<std::add_const_t<T>>().clone());
    // A T* must be implicitly convertible to that type.
    static_assert(std::is_convertible_v<T*, clone_return_t>,
		  "T* must be compatible with the return type of T::clone");
    // Let's be more specific for our use case: we want clone_return_t
    // to be a pointer to (say) B and for B to either be T or to be a
    // base of T.
    using B = std::remove_pointer_t<clone_return_t>;
    static_assert(std::is_base_of_v<B, T>,
		  "T::clone must return pointer-to-T or pointer-to-base-of-T");

    static T* cast(B* b) {
      if constexpr (std::is_same_v<T, B>) {
	return b;
      } else {
	return dynamic_cast<T*>(b);
      }
    }

    // Any other instantiation is a friend.
    template <typename U>
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
    template <typename U>
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
    template <
      typename U,
      typename=std::enable_if<std::is_same_v<B, typename clone_ptr<U>::B>>
    >
    clone_ptr(clone_ptr<U> const& other) : m_ptr{cast(other->clone())} {
    }

    // Move assign
    clone_ptr& operator=(clone_ptr&& other) noexcept {
      m_ptr = std::move(other.m_ptr);
      return *this;
    }
    // Move assign from any allowed by unique_ptr
    template <typename U>
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
    template <
      typename U,
      typename=std::enable_if<std::is_same_v<B, typename clone_ptr<U>::B>>
      >
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
    template <typename L, typename R>
    friend bool operator==(clone_ptr<L> const& lhs, clone_ptr<R> const& rhs);
    template <typename L, typename R>
    friend bool operator!=(clone_ptr<L> const& lhs, clone_ptr<R> const& rhs);
    unique_ptr m_ptr;
  };

  template <typename L, typename R>
  bool operator==(clone_ptr<L> const& lhs, clone_ptr<R> const& rhs) {
    return lhs.m_ptr == rhs.m_ptr;
  }
  template <typename L, typename R>
  bool operator!=(clone_ptr<L> const& lhs, clone_ptr<R> const& rhs) {
    return !(lhs == rhs);
  }

  template <typename T>
  bool operator!=(clone_ptr<T> const& v, std::nullptr_t) {
    return bool(v);
  }
  template <typename T>
  bool operator!=(std::nullptr_t, clone_ptr<T> const& v) {
    return bool(v);
  }
  template <typename T>
  bool operator==(clone_ptr<T> const& v, std::nullptr_t) {
    return !v;
  }
  template <typename T>
  bool operator==(std::nullptr_t, clone_ptr<T> const& v) {
    return !v;
  }

  // Helper like std::make_unique
  template <typename T, typename... Args>
  clone_ptr<T> make_clone_ptr(Args&&... args) {
    return clone_ptr<T>{new T{std::forward<Args>(args)...}};
  }

  // dynamic_cast (move)
  template <typename T, typename U>
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
  template <typename T, typename U>
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
