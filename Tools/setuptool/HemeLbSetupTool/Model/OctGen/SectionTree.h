// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_SECTIONTREE_H
#define HEMELBSETUPTOOL_SECTIONTREE_H

#include <vector>
#include <list>

#include "Oct.h"
#include "MaskTree.h"
#include "FluidSiteTree.h"

// This is a flattened Octree that at the lowest level stores zero or
// more Sections, in the meaning of PETSc. This is a pair of integers:
// an offset and a count. These will be used to index into other
// arrays. The section are stored by index but have an attached name
// (which can be an empty string.

// The tree data is stored in nLevel vectors

// Each branch node in the tree is 8 indices giving the offsets of the
// children in the next level down's vector.

// Leaf nodes are the same but stores the offset into the section

class SectionTree {
public:
  typedef std::shared_ptr<SectionTree> Ptr;

  typedef MaskTree::Int Int;
  typedef uint64_t IndT;
  
  typedef std::vector<IndT> TreeLevel;
  typedef std::vector<TreeLevel> Tree;
  
  typedef std::pair<IndT, IndT> RangeT;
  typedef std::vector<RangeT> SectionT;

  static constexpr IndT NA() {return ~0;};
  
  static inline Int LocalOffset(Int i, Int j, Int k, Int lvl) {
    Int xbit = (i >> lvl) & 1;
    Int ybit = (j >> lvl) & 1;
    Int zbit = (k >> lvl) & 1;
    return (xbit << 2) | (ybit << 1) | zbit;
  }
  
  IndT FindIndex(Int i, Int j, Int k) const;
  
  void AddSection(std::string& name);

  const Tree& GetTree() const;
  // SectionT& GetSection(unsigned i);
  // SectionT& GetSection(const std::string& name);
  // const SectionT& GetSection(unsigned i) const;
  // const SectionT& GetSection(const std::string& name) const;
  
  void Write(const std::string& fn) const;
  
  
private:
  
  friend class SectionTreeBuilder;
  
  SectionTree(size_t nl);
  
  Int nLevels;
  Tree indices;
  IndT section_size;
  Tree counts;
  
  std::list<SectionT> secs;
  std::list<std::string> sec_names;
  
};


#endif
