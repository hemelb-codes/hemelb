#include <boost/filesystem.hpp>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "SectionTree.h"
#include "enumerate.hpp"


SectionTree::SectionTree(size_t nl) : nLevels(nl), indices(nl) {
}

const SectionTree::Tree& SectionTree::GetTree() const {
  return indices;
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
  }
}
