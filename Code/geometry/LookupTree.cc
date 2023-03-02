// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/LookupTree.h"

#include <cassert>

namespace hemelb::geometry::octree {

    LookupTree::LookupTree(U16 N) : levels(N+1), n_levels(N) {
        for (U16 i = 0; i <= N; ++i) {
            levels[i].level = i;
        }
    }

    // Get from the lowest level in the tree
    NodeRef LookupTree::operator()(Vec16 ijk) const {
        auto ans = NodeRef{this};
        auto leaf_oct = ijk_to_oct(ijk);

        std::size_t current_idx = 0U;
        U16 current_level = 0U;
        U64 current_oct = 0U;
        ans.path[current_level] = current_idx;
        for (; current_level < n_levels; ++current_level) {
            U64 next_oct = leaf_oct >> (3U * (n_levels - current_level - 1));
            auto local = next_oct & 7U;

            auto& L = levels[current_level];
            if (L.node_ids[current_idx] != current_oct)
                throw Exception() << "Tree fail!";

            current_idx = ans.path[current_level + 1] = L.child_indices[current_idx][local];
            if (current_idx == Level::NC) {
                // The leaf doesn't exist, return the path to the last extant node
                return ans;
            }
            current_oct = next_oct;
        }
        return ans;
    }

    constexpr U64 bounds_to_end(Vec16 bounds) {
        return ijk_to_oct(bounds -= Vec16::Ones()) + 1;
    }

    WithinBoundsIterator IterBounds::begin() const {
        return WithinBoundsIterator{bounds};
    }
    WithinBoundsIterator IterBounds::end() const {
        WithinBoundsIterator ans{bounds};
        ans.mPos = oct_to_ijk(ans.mEnd);
        return ans;
    }

    WithinBoundsIterator::WithinBoundsIterator(Vec16 bounds) :
            mPos(Vec16::Zero()), mBounds(bounds), mEnd(bounds_to_end(mBounds)) {
    }

    bool WithinBoundsIterator::within_dims(Vec16 const& coord) const {
        return (coord[0] < mBounds[0]) && (coord[1] < mBounds[1]) && (coord[2] < mBounds[2]);
    }

    [[nodiscard]] Vec16 WithinBoundsIterator::advance(U64 prev) const {
        U64 next = prev + 1U;
        if (next >= mEnd)
            return oct_to_ijk(mEnd);

        auto next_pos = oct_to_ijk(next);
        if (within_dims(next_pos))
            return next_pos;

        unsigned N = 0;
        auto tmp = next;
        // Pop off the trios of zero, keeping count.
        // Effectively, we are moving back up levels of the tree
        for (; (tmp & 7U) == 0U; tmp >>= 3U)
            ++N;

        // Fill the popped zero trios with 7 = 0b111
        // This moves us to the last leaf node in the subtree
        for (unsigned i = 0; i < N; ++i)
            tmp = (tmp << 3U) ^ 7U;

        // Now advance to the next leaf point
        return advance(tmp);
    }

    WithinBoundsIterator& WithinBoundsIterator::operator++() {
        U64 prev = ijk_to_oct(mPos);
        mPos = advance(prev);
        return *this;
    }

    Vec16 const& WithinBoundsIterator::operator*() const {
        return mPos;
    }

    bool operator==(WithinBoundsIterator const& lhs, WithinBoundsIterator const& rhs) {
        return lhs.mPos == rhs.mPos;
    }
    bool operator!=(WithinBoundsIterator const& lhs, WithinBoundsIterator const& rhs) {
        return !(lhs.mPos == rhs.mPos);
    }

    LookupTree build_block_tree(const Vec16& dimensionsInBlocks, std::vector<site_t> const& fluidSitesPerBlock) {
        auto biggest_dim = *std::max_element(dimensionsInBlocks.begin(), dimensionsInBlocks.end());
        // What power of two is greater than or equal to the biggest dimension of the domain?
        U16 N = 1;
        U16 cube_size = 2U;
        while (cube_size < biggest_dim) {
            ++N;
            cube_size *= 2;
        }
//        U64 max_blocks = cube_size;
//        max_blocks *= cube_size;
//        max_blocks *= cube_size;

        LookupTree ans(N);

        // Geometry file format decrees this layout of blocks
        auto const block_strides = util::Vector3D<std::size_t>(
                dimensionsInBlocks.y() * dimensionsInBlocks.z(),
                dimensionsInBlocks.z(),
                1
        );

        // Now iterate over the blocks, IN OCTREE ORDER
        for (auto block_ijk: IterBounds{dimensionsInBlocks}) {
            auto block_i = Dot(block_ijk, block_strides);
            if (auto nsites = fluidSitesPerBlock[block_i]) {
                // Add to the tree
                auto oct = ijk_to_oct(block_ijk);
                // ll = leaf level - start here
                auto ll = N;
                ans.levels[ll].node_ids.push_back(oct);
                ans.levels[ll].sites_per_node.push_back(nsites);

                // pl = parent level
                U16 pl = ll;
                while (pl != 0) {
                    // Walk up the tree - doing the "increment" first to have an easy to express condition
                    --pl;
                    // The index of the child node relative to its parent
                    auto local = oct & 7U;
                    // Level index of the parent node by popping off the local part
                    oct >>= 3U;

                    auto& lvl = ans.levels[pl];
                    if (lvl.node_ids.empty() || lvl.node_ids.back() != oct) {
                        // Have either:
                        // - never reached this level of the tree, or
                        // - need a new node (this is true because iterating over the block in octree order)
                        lvl.node_ids.push_back(oct);
                        lvl.sites_per_node.push_back(nsites);
                        lvl.child_indices.push_back(Level::NOCHILDREN);
                    } else {
                        // Update existing node
                        lvl.sites_per_node.back() += nsites;
                    }
                    // Add index of child to parent
                    if (lvl.child_indices.back()[local] == Level::NC) {
                        lvl.child_indices.back()[local] = ans.levels[pl + 1].node_ids.size() - 1;
                    } else {
                        assert(lvl.child_indices.back()[local] == ans.levels[pl + 1].node_ids.size() - 1);
                    }
                }
            }
        }
        return ans;
    }
}