// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_SECTIONTREEBUILDER_H
#define HEMELBSETUPTOOL_SECTIONTREEBUILDER_H

#include "MaskTree.h"
#include "SectionTree.h"

class SectionTreeBuilder : public MaskTree::ConstVisitor {
public:
  using Int = MaskTree::Int;
  
  static inline Int LocalOffset(const MaskTree::Node& n) {
    return SectionTree::LocalOffset(n.X(), n.Y(), n.Z(), n.Level());
  }
  
  SectionTreeBuilder(const MaskTree& mask, const FluidTree& edges);
  
  SectionTree::Ptr operator()();
  
  virtual void Arrive(MaskTree::ConstNodePtr np);
  virtual void Depart(MaskTree::ConstNodePtr n);
  
  inline SectionTree::IndT GetSectionSize() const {
    return offsets[0];
  }

private:
  const MaskTree& maskTree;
  const FluidTree& edgeTree;
  
  const Int nLevels;
  const unsigned nEdgeSites;
  // Current insertion index for each level
  // Has size == nLevels +1
  std::vector<SectionTree::IndT> offsets;
  std::vector<FluidTree::ConstNodePtr> edge_ptrs;
  SectionTree::Ptr output;
};

#endif
