// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_LOOKUPTREE_H
#define HEMELB_GEOMETRY_LOOKUPTREE_H

#include <algorithm>
#include <cstdint>
#include <compare>
#include <memory>
#include <vector>

#include "Exception.h"
#include "units.h"
#include "util/Vector3D.h"
#include "net/MpiCommunicator.h"
#include "net/MpiWindow.h"

namespace hemelb::geometry::octree
{

    //  Children are stored in the standard C++ order with the
    //  associated level-local index being their 3D local coordinate
    //  concatenated and treated as a binary number.
    //
    //  I.e.:
    //
    // (0,0,0) - 000 - 0
    // (0,0,1) - 001 - 1
    // (0,1,0) - 010 - 2
    // (0,1,1) - 011 - 3
    // (1,0,0) - 100 - 4
    // (1,0,1) - 101 - 5
    // (1,1,0) - 110 - 6
    // (1,1,1) - 111 - 7

    // This is repeated down each level.
    // E.g. a point at (2, 0, 7)  = (010, 000, 111)
    // would go root node -> child(0,0,1) -> child (1,0,1) -> child(0,0,1)
    // effectively we are interleaving the bits.

    // Restrict each dimension to a 16 bit int and the oct index to 64 bits.
    using U16 = std::uint16_t;
    using U64 = std::uint64_t;
    using Vec16 = util::Vector3D<U16>;

    // We can squeeze an extra one in, because the root node in level zero
    // always has index zero.
    constexpr U16 MAX_LEVELS = 17;

    // These constants used for interleaving the bits of the dimension coords
    constexpr std::array<U64, 4> shifts = {2, 4, 8, 16};
    constexpr std::array<U64, 5> masks = {
            0b0000000000000000001001001001001001001001001001001001001001001001,
            0b0000000000000000000011000011000011000011000011000011000011000011,
            0b0000000000000000000000001111000000001111000000001111000000001111,
            0b0000000000000000000000000000000011111111000000000000000011111111,
            0b0000000000000000000000000000000000000000000000001111111111111111
    };

    // Spread the first 16 bits of the input out. I.e.:
    // 0b0000000000000000000000000000000000000000000000001111111111111111
    // goes to:
    // 0b0000000000000000001001001001001001001001001001001001001001001001
    constexpr U64 bitspread_one(U64 x) {
        x = (x | (x << shifts[3])) & masks[3];
        x = (x | (x << shifts[2])) & masks[2];
        x = (x | (x << shifts[1])) & masks[1];
        x = (x | (x << shifts[0])) & masks[0];
        return x;
    }

    // Condense the bits of a spread number back to the lowest 16 bits
    constexpr U16 bitcondense_one(U64 x) {
        // mask out irrelevant bits
        x = x & masks[0];
        x = (x | x >> shifts[0]) & masks[1];
        x = (x | x >> shifts[1]) & masks[2];
        x = (x | x >> shifts[2]) & masks[3];
        x = (x | x >> shifts[3]) & masks[4];
        return std::uint16_t(x);
    }

    // Interleave the bits of the input
    constexpr U64 ijk_to_oct(Vec16 const & ijk) {
        auto x = bitspread_one(ijk.x());
        auto y = bitspread_one(ijk.y());
        auto z = bitspread_one(ijk.z());
        return (x << 2U) ^ (y << 1U) ^ z;
    }

    // Separate the bits of the input and return a Vector
    constexpr Vec16 oct_to_ijk(U64 oct) {
        auto x = bitcondense_one(oct >> 2U);
        auto y = bitcondense_one(oct >> 1U);
        auto z = bitcondense_one(oct);
        return {x, y, z};
    }

    class LookupTree;

    // Single level of a condensed octree
    // Have a structure of arrays view on the data - vectors must have the same length
    class Level {
    public:
        // Constant for no child
        static constexpr std::size_t NC = ~0U;
        static constexpr std::array<std::size_t, 8> NOCHILDREN = {NC, NC, NC, NC, NC, NC, NC, NC};

        // The octree ID of the point - unique at the level and can be converted to 3D coordinates
        std::vector<U64> node_ids;
        // The number of fluid sites under the node (is the sum of child sites_per_nodes)
        std::vector<U64> sites_per_node;
        // The indexes of child nodes in the next Level's arrays
        std::vector<std::array<std::size_t, 8>> child_indices;
        // This level's ID (root node has level == 0)
        U16 level;
    };

    class LeafRef {
        Level const* level;
        std::size_t idx;
        friend class NodeRef;
        friend class LookupTree;
        friend class LeafIterator;

        inline LeafRef(Level const& l, std::size_t i) : level(&l), idx(i) {}
    public:
        inline explicit operator bool() const{
            return idx != Level::NC;
        }

        [[nodiscard]] inline std::size_t index() const {
            return idx;
        }
        [[nodiscard]] inline U64 const& id() const {
            return level->node_ids[idx];
        }
        [[nodiscard]] inline U64 const& sites() const {
            return level->sites_per_node[idx];
        }
        [[nodiscard]] inline Vec16 coords() const {
            return oct_to_ijk(id());
        }

        inline friend std::partial_ordering operator<=>(const LeafRef& a, const LeafRef& b) {
            if (a.level != b.level)
                return std::partial_ordering::unordered;
            return a.idx <=> b.idx;
        }
        inline friend bool operator==(const LeafRef& a, const LeafRef& b) {
            return (a <=> b) == std::partial_ordering::equivalent;
        }

        };

    // Refer to a particular node in the tree
    class NodeRef {
    public:
        std::array<std::size_t, MAX_LEVELS> path;
        LookupTree const* tree;

        inline explicit NodeRef(LookupTree const* t) : tree{t} {
            std::fill(path.begin(), path.end(), Level::NC);
        }

        [[nodiscard]] std::size_t leaf() const;
    };

    class LeafIterator;
    using LeafRange = std::pair<LeafIterator, LeafIterator>;

    // A tree with one level holds only the root node
    class LookupTree {
    public:
        std::vector<Level> levels;
        U16 n_levels;

        explicit LookupTree(U16 N);

        // Get from the lowest level in the tree
        [[nodiscard]] NodeRef GetPath(Vec16 ijk) const;

        [[nodiscard]] LeafRef GetLeaf(Vec16 ijk) const;

        [[nodiscard]] LeafRange IterLeaves() const;

        [[nodiscard]] Vec16 GetLeafCoords(std::size_t) const;
    };

    // Iterate over tree leaf nodes in storage order
    class LeafIterator {
        LeafRef r;
    public:
        inline explicit LeafIterator(LeafRef lr) : r(std::move(lr)) {}

        inline LeafIterator& operator++() {
            r.idx++;
            return *this;
        }
        inline LeafRef const& operator*() const {
            return r;
        }
        friend auto operator<=>(const LeafIterator&, const LeafIterator&) = default;
    };

    inline LeafIterator begin(LeafRange const& r) {
        return r.first;
    }
    inline LeafIterator end(LeafRange const& r) {
        return r.second;
    }
    // Iterate over all the octree leaf coordinates that are within a simple cuboid
    // bounding box starting at origin up to (not including) the bounds.
    //
    // Used to iterate a GMY blocks in octree order
    class WithinBoundsIterator {
        Vec16 mPos;
        Vec16 mBounds;
        U64 mEnd;

        friend struct IterBounds;
        explicit WithinBoundsIterator(Vec16 bounds);
    public:

        [[nodiscard]] bool within_dims(Vec16 const& coord) const;

        [[nodiscard]] Vec16 advance(U64 prev) const;

        WithinBoundsIterator& operator++();

        Vec16 const& operator*() const;

        friend bool operator==(WithinBoundsIterator const& lhs, WithinBoundsIterator const& rhs);
        friend bool operator!=(WithinBoundsIterator const& lhs, WithinBoundsIterator const& rhs);
    };

    struct IterBounds {
        Vec16 bounds;

        [[nodiscard]] WithinBoundsIterator begin() const;
        [[nodiscard]] WithinBoundsIterator end() const;
    };

    LookupTree build_block_tree(const Vec16& dimensionsInBlocks, std::vector<site_t> const& fluidSitesPerBlock);



    // Store something (in this case the rank that owns a given site)
    // distributed across MPI processes/

    class DistributedStore {
        // Needed for offset calculations
        MPI_Aint sites_per_block;
        // Octree describing layout of blocks
        LookupTree block_tree;
        // Which MPI rank holds data for which block, **in octree order**.
        std::vector<int> storage_rank;

        net::MpiCommunicator comm;
        // MPI_Win and allocated data - holds which MPI rank owns a site
        using WinData = net::WinData<int>;
        WinData rank_that_owns_site_win;

        // Compute the index within a partition's array where a block lives (block given by its flat index).
        MPI_Aint ComputeBlockStart(std::size_t block_idx) const;
    public:
        // Construct - collective on the communicator
        DistributedStore(site_t sites_per_block, LookupTree tree, std::vector<int> ranks, net::MpiCommunicator c);

        [[nodiscard]] inline LookupTree const& GetTree() const {
            return block_tree;
        }

        [[nodiscard]] inline auto GetBlockCount() const {
            return storage_rank.size();
        }

        // Allow many writes to the same block, wherever it is.
        // Constructed by a WriteSession below
        struct BlockWriteChunk {
            // where the block lives
            int rank;
            // offset of the block in the owner's window
            MPI_Aint block_start_idx;
            // window wrapper
            WinData* window;

            // Return a write-object for the given site
            WinData::reference operator()(site_t site_id);
        };

        // Hide some access details behind in a more natural access style
        // Recall that the destructor has a fence so is collective.
        struct WriteSession : public WinData::WriteSession {
            DistributedStore* store;
            explicit WriteSession(DistributedStore* s);
            // Do the block level calculations once
            BlockWriteChunk operator()(Vec16 const& block_coord);
            BlockWriteChunk operator()(std::size_t block_idx);
            // Do everything
            WinData::reference operator()(Vec16 const& block_coord, site_t site_id);
        };

        // Start a load of writes
        WriteSession begin_writes();

        [[nodiscard]] int GetSiteRank(Vec16 const& blockIjk, site_t siteIdx) const;
        [[nodiscard]] int GetSiteRank(std::size_t blockIdx, site_t siteIdx) const;
    };

}


#endif
