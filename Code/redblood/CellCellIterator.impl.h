//
// Copyright (C) University College London, 2007-2012, all rights reserved.
//
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
//

// Defines both const and non-const iterators for DivideConquerCells.
// Could be done via templates, but compilation is sufficiently slow as it is.

#ifdef HEMELB_DOING_NONCONST
#  define HEMELB_ITERATOR iterator
#  define HEMELB_CONST
#else
#  define HEMELB_ITERATOR const_iterator
#  define HEMELB_CONST const
#endif

#ifndef NDEBUG
#  define HEMELB_GET(X) .at(X)
#else
#  define HEMELB_GET(X) [X]
#endif


class HEMELB_ITERATOR {
    typedef base_type::HEMELB_ITERATOR wrappee_iterator;
  public:
    typedef LatticePosition value_type;
    typedef value_type HEMELB_CONST & reference;
    typedef value_type HEMELB_CONST * pointer_type;
    typedef value_type const & const_reference;
    typedef value_type const * const_pointer_type;

    HEMELB_ITERATOR(
        DivideConquerCells HEMELB_CONST &_owner,
        wrappee_iterator const &_w
    ) : owner_(_owner), wrappee_(_w) {}
#   ifndef HEMELB_DOING_NONCONST
      HEMELB_ITERATOR(iterator const &_in)
        : owner_(_in.owner_), wrappee_(_in.wrappee_) {}
#   else
      friend class DivideConquerCells::const_iterator;
#   endif

    reference operator*() {
        return (*owner_.cells_) HEMELB_GET(wrappee_->second.cellIndex)
              .GetVertices()
              HEMELB_GET(wrappee_->second.nodeIndex);
    }
    pointer_type operator->() { return &this->operator*(); }

    const_reference operator*() const {
        return (*owner_.cells_) HEMELB_GET(wrappee_->second.cellIndex)
              .GetVertices()
              HEMELB_GET(wrappee_->second.nodeIndex);
    }
    const_pointer_type operator->() const { return &this->operator*(); }
    LatticeVector const & key() const { return wrappee_->first; }
    HEMELB_ITERATOR & operator++() { ++wrappee_; return *this; }
    HEMELB_ITERATOR & operator--() { --wrappee_; return *this; }
    HEMELB_ITERATOR operator++(int) {
      return HEMELB_ITERATOR(owner_, wrappee_++);
    }
    HEMELB_ITERATOR operator--(int) {
      return HEMELB_ITERATOR(owner_, wrappee_--);
    }

    bool operator==(HEMELB_ITERATOR const &_in) const {
      return _in.wrappee_ == wrappee_;
    }
    bool operator!=(HEMELB_ITERATOR const &_in) const {
      return not operator==(_in);
    }

  protected:
    //! Container from which this object was obtained
    DivideConquerCells HEMELB_CONST &owner_;
    //! Iterator over the mapped objects
    wrappee_iterator wrappee_;
};

#undef HEMELB_ITERATOR
#undef HEMELB_CONST
#undef HEMELB_GET
