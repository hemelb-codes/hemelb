// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REDBLOOD_CELLCELL_H
#define HEMELB_REDBLOOD_CELLCELL_H

#include <vector>
#include <cassert>
#include <memory>
#include <initializer_list>
#include <type_traits>

#include "units.h"
#include "Exception.h"
#include "geometry/LatticeData.h"
#include "redblood/DivideConquer.h"
#include "redblood/Cell.h"
#include "redblood/Node2Node.h"
#include "redblood/stencil.h"
#include "redblood/Interpolation.h"
#include "redblood/Borders.h"
#include "redblood/parallel/CellParallelization.h"

namespace hemelb
{
  namespace redblood
  {
    // References a node of a mesh in the divide-and-conquer box
    class CellReference
    {
      public:
        //! Index of cell in input container
        CellContainer::const_iterator cellIterator;
        //! Index of node in mesh
        site_t nodeIndex;
        //! Id of the nearest borders
        size_t nearBorder;
    };

    class DivideConquerCells;
    namespace detail {
#     ifndef NDEBUG
#       define HEMELB_GET(X) .at(X)
#     else
#       define HEMELB_GET(X) [X]
#     endif

      // Use const base_type for const_iterator, (mutable) base_type for (mutable) iterator.
      template <typename BASE>
      class iterator_base
      {
	// This is the unqualified type
	using base_type = std::remove_const_t<BASE>;
	static constexpr bool am_const = std::is_const<BASE>::value;
	// This is possibly const qualified type we are instantiated on
	using owner_type = BASE;
	// This is the other one
	using other_type = std::conditional_t<am_const,
					      base_type,
					      base_type const>;

	// This is the iterator type we are wrapping
	using wrappee_iterator = std::conditional_t<am_const,
						    typename base_type::const_iterator,
						    typename base_type::iterator>;

      public:
	using value_type = LatticePosition;
	using reference = std::conditional_t<am_const,
					     value_type const &,
					     value_type&>;
	using pointer = std::conditional_t<am_const,
					   value_type const*,
					   value_type*>;
	using const_reference = value_type const &;
	using const_pointer = value_type const *;
	using difference_type = typename wrappee_iterator::difference_type;
	using iterator_category = typename wrappee_iterator::iterator_category;

	// Iterators really should be default constructible
	iterator_base() = default;

	iterator_base(owner_type const &owner, wrappee_iterator const &w) :
	  owner(&owner), wrappee(w)
	{
	}

	// Copy(ish) constructors
	//
	// Want: const{const}, const{mutable}, mutable{mutable}
	// Don't want: mutable{const}
	//
	// Allow self{self}
	iterator_base(iterator_base const& in) = default;
	// Allow self{other} iff self == const
	iterator_base(iterator_base<other_type> const& in) :
	  owner(in.owner), wrappee(in.wrappee)
	{
	  static_assert(am_const, "Cannot convert const->mutable");
	}

	// Want mutable to declare that const is it's friend.
	friend class iterator_base<const base_type>;

	reference operator*()
	{
	  // @formatter:off
	  return (*wrappee->second.cellIterator)->GetVertices()HEMELB_GET(wrappee->second.nodeIndex);
	  // @formatter:on
	}
	pointer operator->()
	{
	  return &this->operator*();
	}

	const_reference operator*() const
	{
	  // @formatter:off
	  return (*wrappee->second.cellIterator)->GetVertices()HEMELB_GET(wrappee->second.nodeIndex);
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
	iterator_base &operator++()
	{
	  ++wrappee;
	  return *this;
	}
	iterator_base &operator--()
	{
	  --wrappee;
	  return *this;
	}
	iterator_base operator++(int)
	{
	  return iterator_base(owner, wrappee++);
	}
	iterator_base operator--(int)
	{
	  return iterator_base(owner, wrappee--);
	}

	friend bool operator==(iterator_base const& a, iterator_base const& b)
	{
	  return (a.owner == b.owner) && (a.wrappee == b.wrappee);
	}
	friend bool operator!=(iterator_base const& a, iterator_base const& b)
	{
	  return !(a == b);
	}

	CellReference const &GetCellReference() const
	{
	  return wrappee->second;
	}

	std::conditional_t<am_const, CellReference const&, CellReference&> GetCellReference()
	{
	  return wrappee->second;
	}

	//! Gets integer coding for whether node is close to boundary
	size_t GetNearBorder() const
	{
	  return GetCellReference().nearBorder;
	}
	//! True if close to given boundary
	bool IsNearBorder(Borders border) const
	{
	  return GetNearBorder() bitand size_t(border);
	}
	//! True if close to any boundary
	bool IsNearBorder() const
	{
	  return GetNearBorder() != static_cast<size_t>(Borders::CENTER);
	}

	//! Returns shared pointer to cell pointed to by this object
	CellContainer::const_reference GetCell() const
	{
	  return *GetCellReference().cellIterator;
	}

	operator wrappee_iterator() const
	{
	  return wrappee;
	}

	// Want: C = C, C = M, M = M
	// Don't want: M = C
	// So always allow self = self
	iterator_base& operator=(iterator_base const &c)
	{
	  owner = c.owner;
	  wrappee = c.wrappee;
	  return *this;
	}

      protected:
	//! Container from which this object was obtained
	owner_type const* owner = nullptr;
	//! Iterator over the mapped objects
	wrappee_iterator wrappee;
      };
#     undef HEMELB_GET
    } // namespace detail

    //! Organizes nodes in cells in boxes
    //! The object is to easily check nodes that are within interaction distance
    class DivideConquerCells : public DivideConquer<CellReference>
    {
        //! Type of the base class
        typedef DivideConquer<CellReference> base_type;

      public:
        //! Iterates over vertices
        //! Wraps a multimap iterator. As such it provides the exact same
        //! garantees.
        using iterator = detail::iterator_base<base_type>;
        //! Iterates over vertices
        //! Wraps a multimap iterator. As such it provides the exact same
        //! garantees.
        using const_iterator = detail::iterator_base<base_type const>;

        typedef std::reverse_iterator<iterator> reverse_iterator;
        typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
        //! Iterator pair over single box
        typedef std::pair<const_iterator, const_iterator> const_range;
        //! Iterates over pairs of nodes that are close to one another
        //! Each pair is visited only once.
        class pair_range;

        //! Constructor
        DivideConquerCells(CellContainer const &cells, LatticeDistance boxsize,
                           LatticeDistance halosize);

        //! Gets all nodes in a box
        const_range operator()(LatticeVector const &pos) const;
        //! Gets all nodes in a box
        const_range operator()(LatticePosition const &pos) const
        {
          return this->operator()(base_type::DowngradeKey(pos));
        }
        template<class T>
        typename std::enable_if<std::is_integral<T>::value, const_range>::type operator()(
            std::initializer_list<T> pos) const
        {
          assert(pos.size() == 3);
          return operator()(LatticeVector(*pos.begin(),
                                          *std::next(pos.begin()),
                                          *std::next(std::next(pos.begin()))));
        }
        template<class T>
        typename std::enable_if<std::is_floating_point<T>::value, const_range>::type operator()(
            std::initializer_list<T> pos) const
        {
          assert(pos.size() == 3);
          return operator()(LatticePosition(*pos.begin(),
                                            *std::next(pos.begin()),
                                            *std::next(std::next(pos.begin()))));
        }

        iterator begin()
        {
          return iterator(*this, base_type::begin());
        }
        iterator end()
        {
          return iterator(*this, base_type::end());
        }
        const_iterator begin() const
        {
          return const_iterator(*this, base_type::begin());
        }
        const_iterator end() const
        {
          return const_iterator(*this, base_type::end());
        }
        reverse_iterator rbegin()
        {
          return reverse_iterator(end());
        }
        reverse_iterator rend()
        {
          return reverse_iterator(begin());
        }
        const_reverse_iterator rbegin() const
        {
          return const_reverse_iterator(end());
        }
        const_reverse_iterator rend() const
        {
          return const_reverse_iterator(begin());
        }
        base_type::size_type size() const
        {
          return base_type::size();
        }

        //! Loops over pair of vertices closer than input distance
        pair_range pair_begin(LatticeDistance maxdist) const;

        //! Distance from border below which an object is in the halo
        LatticeDistance GetHaloLength() const
        {
          return haloLength;
        }
        //! Size of each box
        LatticeDistance GetBoxSize() const
        {
          return base_type::GetBoxSize();
        }
        //! Cells that are present in this object
        CellContainer const &GetCells() const
        {
          return cells;
        }

        //! After vertices have moved, update mapping and whether it is near
        //! boundary
        void update();
        //! recomputes using current cells
        void SetBoxSizeAndHalo(LatticeDistance boxSize, LatticeDistance halo);

        //! In a parallel simulation, some cells will leave/enter the domain.
        //! \param[in] 3-tuple with the newly owned cells, the disowned cells, and the lent cells
        void update(parallel::ExchangeCells::ChangedCells const& changedCells);

        //! Insert a new cell
        //! Returns true if the cell was inserted, false if it already existed.
        bool insert(CellContainer::value_type cell);
        //! Removes a cell
        //! Returns true if the cell did exist.
        bool remove(CellContainer::value_type cell);

      protected:
        //! Distance from border below which an object is in the halo
        LatticeDistance haloLength;
        //! Container of cells
        CellContainer cells;
        //! Keeps track of which cells in the box are lent (currentlyLentCells must always be a subset of DivideConquerCells::cells)
        CellContainer currentlyLentCells;
        //! Helper method that combines all the CellContainers in a LentCells object into a single CellContainer
        CellContainer LentCellsToSingleContainer(parallel::LentCells const &lentCells) const;
    };

    class DivideConquerCells::pair_range
    {
        //! Parent iterator
        typedef DivideConquerCells::const_iterator iterator;

      public:
        typedef std::pair<iterator, iterator> value_type;

        //! Constructs a pair range iterator
        pair_range(DivideConquerCells const &owner, iterator const &begin, iterator const &end,
                   LatticeDistance maxdist);

        //! Whether current iteration is valid
        bool is_valid() const
        {
          return currents.first != ends.first;
        }
        //! Whether current iteration is valid
        operator bool() const
        {
          return is_valid();
        }

        value_type const &operator*() const
        {
#ifndef NDEBUG

          if (not is_valid())
          {
            throw Exception() << "Iterator is invalid\n";
          }

#endif
          return currents;
        }
        value_type const *operator->() const
        {
          return & (this->operator*());
        }

        //! Goes to next value and returns true if valid
        bool operator++();

      protected:
        //! Maximum distance for which to report pair
        LatticeDistance maxdist;
        //! Iterates over boxes we want to see
        BorderBoxIterator box_iterator;
        //! Iterator for main item
        value_type currents;
        //! range for iteration over second item
        value_type ends;
        //! Container over which this one iterates
        DivideConquerCells const &owner;

        bool doBox();
        bool nextDist();
    };

    namespace
    {
      //! Spread force from given vertices to lattice-sites
      template<class STENCIL>
      void spreadForce(LatticePosition const &vertex, geometry::LatticeData &latticeData,
                       LatticeForceVector const &force)
      {
        proc_t procid;
        site_t siteid;
        InterpolationIterator<STENCIL> spreader = interpolationIterator<STENCIL>(vertex);

        for (; spreader; ++spreader)
        {
          if (latticeData.GetContiguousSiteId(*spreader, procid, siteid))
          {
            latticeData.GetSite(siteid).AddToForce(force * spreader.weight());
          }
        }
      }
    }
    //! \brief Computes cell <-> cell interactions and spread to grid
    //! \details Given a partition of the cells' nodes and node <-> node interaction
    //! functional, computes the short-range that can occur between cells that are
    //! too close to one another. The interaction forces are computed and spread to
    //! the lattice.
    template<class STENCIL>
    void addCell2CellInteractions(DivideConquerCells const &dnc, Node2NodeForce const &functional,
                                  geometry::LatticeData &latticeData)
    {
      auto range = dnc.pair_begin(functional.cutoff);

      for (; range.is_valid(); ++range)
      {
        LatticeForceVector const force(functional(*range->first, *range->second));
        // spread to the grid from from one node and from the other
        spreadForce<STENCIL>(*range->first, latticeData, force);
        spreadForce<STENCIL>(*range->second, latticeData, -force);
      }
    }

  }
} // hemelb::redblood

#endif
