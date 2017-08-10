// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_SECTIONTREEBUILDER_HPP
#define HEMELBSETUPTOOL_SECTIONTREEBUILDER_HPP

template <class... OtherTreeTs>
class SectionTreeBuilder : public MaskTree::ConstVisitor {
public:
  using Int = MaskTree::Int;
  
  static Int LocalOffset(const MaskTree::Node& n) {
    Int lvl = n.Level();
    Int xbit = (n.X() >> lvl) & 1;
    Int ybit = (n.Y() >> lvl) & 1;
    Int zbit = (n.Z() >> lvl) & 1;
    return (xbit << 2) | (ybit << 1) | zbit;
  }
  
  SectionTreeBuilder(const MaskTree& mask, const OtherTreeTs&... otherTrees) :
    maskTree(mask), otherInputTrees(otherTrees...), nLevels(mask.Level()), counters(mask.Level()+1)
  {
  }
  
  SectionTree::Ptr operator()() {
    output = SectionTree::Ptr(new SectionTree(nLevels));
    maskTree.Root()->Accept(*this);
    output->section_size = GetSectionSize();
    return output;
  }
  
  virtual void Arrive(MaskTree::ConstNodePtr np) {
    const MaskTree::Node& n = *np;
    Int lvl = n.Level();
    
    if (lvl < nLevels) {
      // I'm not the root node
      
      // Write my level's counter value into parent's current data oct
      // We will then increment it on Departure.
      auto& pdata = output->indices[lvl];
      auto pdata_current_i = counters[lvl+1];
      auto loffset = LocalOffset(n);
      pdata[pdata_current_i + loffset] = counters[lvl];
    }
    
    if (lvl) {
      // We are not a leaf node
      // increase storage by 8 and fill with NA's
      auto& my_data = output->indices[lvl-1];
      my_data.resize(my_data.size() + 8, SectionTree::NA());
    
    } else {
      // We are at a leaf node

      // Add storage for one more RangeT
      // This is implicitly done by departing the leaf
    }
  }

  virtual void Depart(MaskTree::ConstNodePtr n) {
    Int lvl = n->Level();
    // Increment my level's counter value
    counters[lvl] += lvl ? 8:1;
  }
  SectionTree::IndT GetSectionSize() const {
    return counters[0];
  }

private:
  const MaskTree& maskTree;

  typedef std::tuple<const OtherTreeTs&...> OTT;
  OTT otherInputTrees;
  
  Int nLevels;
  std::vector<SectionTree::IndT> counters;
  SectionTree::Ptr output;
};

#endif
