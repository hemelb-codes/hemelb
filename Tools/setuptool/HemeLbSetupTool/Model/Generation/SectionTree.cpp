#include <boost/filesystem.hpp>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "SectionTree.h"
#include "enumerate.hpp"
#include "H5.h"

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

namespace H5 {
  template<>
  DataTypeSharedPtr DataTypeTraits<SVector>::GetType()
  {
    return DataType::Array<float>({3});
  }
  
  template<>
  DataTypeSharedPtr DataTypeTraits<std::array<Link,26>>::GetType()
  {
    return DataType::Array<Link>({26});
  }

  template<>
  DataTypeSharedPtr DataTypeTraits<Link>::GetType()
  {
    hid_t link_id = H5Tcreate(H5T_COMPOUND, sizeof(Link));
    H5Tinsert (link_id, "type", HOFFSET(Link, type),
	       H5T_NATIVE_INT);
    H5Tinsert (link_id, "dist", HOFFSET(Link, dist),
	       H5T_NATIVE_FLOAT);
    H5Tinsert (link_id, "id", HOFFSET(Link, id),
	       H5T_NATIVE_INT);
    
    return DataTypeSharedPtr(new DataType(link_id));
  }
}
template <class T>
void Section<T>::write(H5::GroupPtr grp) const {
  grp->CreateWriteDataSet("offsets", offsets);
  grp->CreateWriteDataSet("counts", counts);
  grp->CreateWriteDataSet("data", data);
}

void SectionTree::Write(const std::string& fn) const {
  auto outfile =  H5::File::Create(fn, H5F_ACC_TRUNC);
  auto root = outfile->CreateGroup("hemelb")->CreateGroup("geometry");
  root->SetAttribute("version", 1);
  root->SetAttribute("levels", int(nLevels));
  for (auto x: const_enumerate(indices)) {
    // Each level is a single dataset
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(4) << x.first;
    root->CreateWriteDataSet(ss.str(), x.second);
  }

  auto wng = root->CreateGroup("wall_normals");
  wall_normals.write(wng);
  links.write(root->CreateGroup("links"));
  
}
