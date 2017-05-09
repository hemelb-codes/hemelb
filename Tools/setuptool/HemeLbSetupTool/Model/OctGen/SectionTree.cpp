#include <boost/filesystem.hpp>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "SectionTree.h"
#include "enumerate.hpp"

SectionTreeBuilder::SectionTreeBuilder(const MaskTree& mask) : maskTree(mask), nLevels(mask.Level()), counters(mask.Level()+1) {
  
}

auto SectionTreeBuilder::LocalOffset(const MaskTree::Node& n) -> Int{
  Int lvl = n.Level();
  Int xbit = (n.X() >> lvl) & 1;
  Int ybit = (n.Y() >> lvl) & 1;
  Int zbit = (n.Z() >> lvl) & 1;
  return (xbit << 2) | (ybit << 1) | zbit;
}
  
void SectionTreeBuilder::Arrive(MaskTree::ConstNodePtr np) {
  const MaskTree::Node& n = *np;
  Int lvl = n.Level();
  //std::cout << lvl << " [" << n.X() << ", " << n.Y() << ", " << n.Z() << "]\n";
  
  if (lvl < nLevels) {
    // I'm not the root node

    // Write my level's counter value into parent's current data oct
    // We will then increment it on Departure.
    auto& pdata = output->indices[lvl];
    auto pdata_current_i = counters[lvl+1];
    auto loffset = LocalOffset(n);
    //std::cout << "output.indices["<<lvl <<"]["<< pdata_current_i << " + " << loffset << "] = counters[" <<lvl <<"] = " << counters[lvl] << std::endl;
    pdata[pdata_current_i + loffset] = counters[lvl];
  }
    
  if (lvl) {
    // We are not a leaf node
    // increase storage by 8 and fill with NA's
    //    std::cout << "Alloc my node storage" << std::endl;
    auto& my_data = output->indices[lvl-1];
    my_data.resize(my_data.size() + 8, SectionTree::NA());
    
  } else {
    // We are at a leaf node

    // Add storage for one more RangeT
    // This is implicitly done by departing the leaf
  }
}

void SectionTreeBuilder::Depart(MaskTree::ConstNodePtr n) {
  Int lvl = n->Level();
  // Increment my level's counter value
  counters[lvl] += lvl ? 8:1;
}
  
SectionTree::IndT SectionTreeBuilder::GetSectionSize() const {
  return counters[0];
}

SectionTree::Ptr SectionTreeBuilder::operator()() {
  output = SectionTree::Ptr(new SectionTree(nLevels));
  maskTree.Root()->Accept(*this);
  output->section_size = GetSectionSize();
  return output;
}

SectionTree::SectionTree(size_t nl) : nLevels(nl), indices(nl) {
}

const SectionTree::Tree& SectionTree::GetTree() const {
  return indices;
}

void SectionTree::Write(const std::string& fn) const {
  namespace fs = boost::filesystem;
  fs::path dirname(fn);
  
  fs::create_directory(dirname);
  for (auto x: const_enumerate(indices)) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(4) << x.first << ".txt";
    auto ifilename = dirname / ss.str();
    std::ofstream outfile(ifilename.native());
    for (auto el: x.second) {
      outfile << el << std::endl;
    }
  }
}
