// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_SECTIONTREEBUILDER_H
#define HEMELBSETUPTOOL_SECTIONTREEBUILDER_H

#include "MaskTree.h"
#include "SectionTree.h"

class SectionTreeBuilder : public MaskTree::ConstVisitor {
public:
  using Int = MaskTree::Int;
  
  static inline Int LocalOffset(const MaskTree::Node& n) {
    Int lvl = n.Level();
    Int xbit = (n.X() >> lvl) & 1;
    Int ybit = (n.Y() >> lvl) & 1;
    Int zbit = (n.Z() >> lvl) & 1;
    return (xbit << 2) | (ybit << 1) | zbit;
  }
  
  SectionTreeBuilder(const MaskTree& mask);
  
  SectionTree::Ptr operator()();
  
  virtual void Arrive(MaskTree::ConstNodePtr np);
  virtual void Depart(MaskTree::ConstNodePtr n);
  
  inline SectionTree::IndT GetSectionSize() const {
    return offsets[0];
  }

private:
  const MaskTree& maskTree;
  
  Int nLevels;
  
  // Current insertion index for each level
  // Has size == nLevels +1
  std::vector<SectionTree::IndT> offsets;
  
  SectionTree::Ptr output;
};

#endif
