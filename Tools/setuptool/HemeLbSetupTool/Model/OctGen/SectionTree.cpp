#include <boost/filesystem.hpp>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "SectionTree.h"
#include "enumerate.hpp"


SectionTree::SectionTree(size_t nl) : nLevels(nl), indices(nl+1), counts(nl+1) {
}

const SectionTree::Tree& SectionTree::GetTree() const {
  return indices;
}

auto SectionTree::FindIndex(Int i, Int j, Int k) const -> IndT {
  Int cur_level = nLevels;
  // Root level (==nLevel) has implicit offset of zero
  IndT cur_offset = 0;
  while (cur_level) {
    // Get the local index
    auto lIndex = LocalOffset(i,j,k, cur_level-1);
    cur_offset = indices[cur_level][cur_offset + lIndex];
    if (cur_offset == NA()) {
      break;
    }
    --cur_level;
  }
  return cur_offset;
}

void SectionTree::Write(const std::string& fn) const {
  // Simple text file in directory format for testing
  namespace fs = boost::filesystem;
  
  fs::path dirname(fn);
  fs::create_directory(dirname);
  
  for (auto x: const_enumerate(indices)) {
    // Each level is a single file 000n.txt
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(4) << x.first << ".txt";
    auto ifilename = dirname / ss.str();
    std::ofstream outfile(ifilename.native());
    // Each file just contains the indices one per line
    for (auto el: x.second) {
      outfile << el << std::endl;
    }

    links.write((dirname / "links").native());
    wall_normals.write((dirname / "wall_normals").native());
  }

  
}
