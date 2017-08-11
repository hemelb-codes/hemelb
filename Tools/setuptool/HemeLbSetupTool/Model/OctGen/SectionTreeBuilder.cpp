#include "SectionTreeBuilder.h"

SectionTreeBuilder::SectionTreeBuilder(const MaskTree& mask) :
  maskTree(mask), nLevels(mask.Level()), offsets(mask.Level()+1)
{
}
  
SectionTree::Ptr SectionTreeBuilder::operator()() {
  output = SectionTree::Ptr(new SectionTree(nLevels));
  maskTree.Root()->Accept(*this);
  return output;
}
  
void SectionTreeBuilder::Arrive(MaskTree::ConstNodePtr np) {
  const MaskTree::Node& n = *np;
  Int lvl = n.Level();
    
  if (lvl < nLevels) {
    // I'm not the root node
      
    auto parent_current_i = offsets[lvl+1];
    auto loffset = LocalOffset(n);
      
    // Write my level's offset value into parent's current data oct
    // We will then increment it on Departure.
    output->indices[lvl][parent_current_i + loffset] = offsets[lvl];

    // Write the current total sites to output
    // Will store the difference on departure
    output->counts[lvl][parent_current_i + loffset] = GetSectionSize();
  }
    
  if (lvl) {
    // We are not a leaf node
    // increase storage by 8 and fill with NA's
    auto new_size = output->indices[lvl-1].size() + 8;
    output->indices[lvl-1].resize(new_size, SectionTree::NA());
    // Or zeros in case of counts
    output->counts[lvl-1].resize(new_size, 0);
  } else {
    // We are at a leaf node

    // Add storage for one more RangeT
    // This is implicitly done by departing the leaf
  }
}

void SectionTreeBuilder::Depart(MaskTree::ConstNodePtr n) {
  Int lvl = n->Level();
  // increment is 2^3 for branch and 1 for leaf nodes
  SectionTree::IndT inc = lvl ? 8 : 1;
  // Increment my level's current offset value
  offsets[lvl] += inc;

  if (lvl < nLevels) {
    // I'm not the root node
    auto parent_current_i = offsets[lvl+1];
    auto loffset = LocalOffset(*n);
    auto arrival_size = output->counts[lvl][parent_current_i + loffset];
    auto departure_size = GetSectionSize();
    output->counts[lvl][parent_current_i + loffset] = departure_size - arrival_size;
  } else {
    output->section_size = GetSectionSize();
  }    
}
 
