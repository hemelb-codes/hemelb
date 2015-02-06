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
#define HEMELB_ITERATOR iterator
#define HEMELB_CONST
#else
#define HEMELB_ITERATOR const_iterator
#define HEMELB_CONST const
#endif

#ifndef NDEBUG
#define HEMELB_GET(X) .at(X)
#else
#define HEMELB_GET(X) [X]
#endif

class HEMELB_ITERATOR
{
    typedef base_type::HEMELB_ITERATOR wrappee_iterator;

  public:
    typedef LatticePosition value_type;
    typedef value_type HEMELB_CONST &reference;
    typedef value_type HEMELB_CONST *pointer;
    typedef value_type const &const_reference;
    typedef value_type const *const_pointer;
    typedef wrappee_iterator::difference_type difference_type;
    typedef wrappee_iterator::iterator_category iterator_category;

    HEMELB_ITERATOR (DivideConquerCells HEMELB_CONST &owner, wrappee_iterator const &w) :
        owner(owner), wrappee(w)
    {
    }
    HEMELB_ITERATOR (HEMELB_ITERATOR const &in) :
        owner(in.owner), wrappee(in.wrappee)
    {
    }
#ifndef HEMELB_DOING_NONCONST
    HEMELB_ITERATOR (iterator const &in) :
        owner(in.owner), wrappee(in.wrappee)
    {
    }
#else
    friend class DivideConquerCells::const_iterator;
#endif

    reference operator*()
    {
      // @formatter:off
      return (*wrappee->second.cellIterator)->GetVertices()
            HEMELB_GET(wrappee->second.nodeIndex);
      // @formatter:on
    }
    pointer operator->()
    {
      return &this->operator*();
    }

    const_reference operator*() const
    {
      // @formatter:off
      return (*wrappee->second.cellIterator)->GetVertices()
            HEMELB_GET(wrappee->second.nodeIndex);
      // @formatter:on
    }
    const_pointer operator->() const
    {
      return &this->operator*();
    }
    LatticeVector const &GetKey() const
    {
      return wrappee->first;
    }
    HEMELB_ITERATOR &operator++()
    {
      ++wrappee;
      return *this;
    }
    HEMELB_ITERATOR &operator--()
    {
      --wrappee;
      return *this;
    }
    HEMELB_ITERATOR operator++(int)
    {
      return HEMELB_ITERATOR(owner, wrappee++);
    }
    HEMELB_ITERATOR operator--(int)
    {
      return HEMELB_ITERATOR(owner, wrappee--);
    }

    bool operator==(HEMELB_ITERATOR const &in) const
    {
      return in.wrappee == wrappee;
    }
    bool operator!=(HEMELB_ITERATOR const &in) const
    {
      return not operator==(in);
    }

    CellReference const &GetCellReference() const
    {
      return wrappee->second;
    }
    CellReference HEMELB_CONST &GetCellReference()
    {
      return wrappee->second;
    }

    //! Gets integer coding for whether node is close to boundary
    int GetNearBorder() const
    {
      return GetCellReference().isNearBorder;
    }
    //! True if close to given boundary
    bool IsNearBorder(CellReference::Borders border) const
    {
      return GetNearBorder() bitand border;
    }
    //! True if close to any boundary
    bool IsNearBorder() const
    {
      return GetNearBorder() != 0;
    }

    //! Returns shared pointer to cell pointed to by this object
    CellContainer::const_reference GetCell() const
    {
      return *GetCellReference().cellIterator;
    }

    operator base_type::HEMELB_ITERATOR() const
    {
      return wrappee;
    }

    void operator=(iterator const &c)
    {
      wrappee = c.wrappee;
    }
#ifndef HEMELB_DOING_NONCONST
    void operator=(const_iterator const &c)
    {
      wrappee = c.wrappee;
    }
#endif

  protected:
    //! Container from which this object was obtained
    DivideConquerCells HEMELB_CONST &owner;
    //! Iterator over the mapped objects
    wrappee_iterator wrappee;
};

#undef HEMELB_ITERATOR
#undef HEMELB_CONST
#undef HEMELB_GET
