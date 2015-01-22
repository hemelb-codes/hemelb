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
    typedef value_type HEMELB_CONST * pointer;
    typedef value_type const & const_reference;
    typedef value_type const * const_pointer;
    typedef wrappee_iterator::difference_type difference_type;
    typedef wrappee_iterator::iterator_category iterator_category;

    HEMELB_ITERATOR(
        DivideConquerCells HEMELB_CONST &_owner,
        wrappee_iterator const &_w
    ) : owner_(_owner), wrappee_(_w) {}
    HEMELB_ITERATOR(HEMELB_ITERATOR const &_in)
      : owner_(_in.owner_), wrappee_(_in.wrappee_) {}
#   ifndef HEMELB_DOING_NONCONST
      HEMELB_ITERATOR(iterator const &_in)
        : owner_(_in.owner_), wrappee_(_in.wrappee_) {}
#   else
      friend class DivideConquerCells::const_iterator;
#   endif

    reference operator*() {
        return (*wrappee_->second.cellIterator)->GetVertices()
              HEMELB_GET(wrappee_->second.nodeIndex);
    }
    pointer operator->() { return &this->operator*(); }

    const_reference operator*() const {
        return (*wrappee_->second.cellIterator)->GetVertices()
              HEMELB_GET(wrappee_->second.nodeIndex);
    }
    const_pointer operator->() const { return &this->operator*(); }
    LatticeVector const & GetKey() const { return wrappee_->first; }
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

    CellReference const & GetCellReference() const { return wrappee_->second; }
    CellReference HEMELB_CONST & GetCellReference() {
      return wrappee_->second;
    }

    //! Gets integer coding for whether node is close to boundary
    int GetNearBorder() const { return GetCellReference().isNearBorder; }
    //! True if close to given boundary
    bool IsNearBorder(CellReference::Borders _border) const {
      return GetNearBorder() bitand _border;
    }
    //! True if close to any boundary
    bool IsNearBorder() const { return GetNearBorder() != 0; }

    //! Returns shared pointer to cell pointed to by this object
    CellContainer::const_reference GetCell() const {
      return *GetCellReference().cellIterator;
    }

    operator base_type::HEMELB_ITERATOR() const { return wrappee_; }

    void operator=(iterator const &_c) { wrappee_ = _c.wrappee_; }
#   ifndef HEMELB_DOING_NONCONST
      void operator=(const_iterator const &_c) { wrappee_ = _c.wrappee_; }
#   endif

  protected:
    //! Container from which this object was obtained
    DivideConquerCells HEMELB_CONST &owner_;
    //! Iterator over the mapped objects
    wrappee_iterator wrappee_;
};

#undef HEMELB_ITERATOR
#undef HEMELB_CONST
#undef HEMELB_GET
