// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_OCT_HPP
#define HEMELBSETUPTOOL_OCT_HPP

#include "Oct.h"
#include <stdexcept>

template<class T>
Octree<T>::Node::Node(Int i, Int j, Int k, Int l) : x(i), y(j), z(k), level(l) {
}
template<class T>
T& Octree<T>::Node::Data() {
  return value;
}
 
template<class T>
auto Octree<T>::Node::X() -> Int{
  return x;
}
template<class T>
auto Octree<T>::Node::Y() -> Int{
  return y;
}
template<class T>
auto Octree<T>::Node::Z() -> Int{
  return z;
}
template<class T>
auto Octree<T>::Node::Level() -> Int{
  return level;
}

// Get a node without creating - returns null pointer if doensn't exist
template<class T>
auto Octree<T>::Branch::Get(Int i, Int j, Int k, Int l) -> NodePtr {
  if (l >= this->level) {
    // a parent - error
    throw std::out_of_range("trying to get a parent node");
  }
      
  // the least significant this->level bits represent the relative index
  // the other bits must match my bits
  // tl = this level
  Int tl_mask = (~0) << this->level;
  // logical not(0) == all the ones
      
  if ((tl_mask & i) != (tl_mask & this->x) ||
      (tl_mask & j) != (tl_mask & this->y) ||
      (tl_mask & k) != (tl_mask & this->z)) {
    throw std::out_of_range("requested node not in my range");
  }
      
  // OK
  return get_nocreate_internal(i, j, k, l);
}

// Get a node without creating - returns null pointer if doensn't exist
template<class T>
auto Octree<T>::Branch::GetCreate(Int i, Int j, Int k, Int l) -> NodePtr {
  if (l >= this->level) {
    // a parent - error
    throw std::out_of_range("trying to get a parent node");
  }
      
  // the least significant this->level bits represent the relative index
  // the other bits must match my bits
  // tl = this level
  Int tl_mask = (~0) << this->level;
  // logical not(0) == all the ones
      
  if ((tl_mask & i) != (tl_mask & this->x) ||
      (tl_mask & j) != (tl_mask & this->y) ||
      (tl_mask & k) != (tl_mask & this->z)) {
    throw std::out_of_range("requested node not in my range");
  }
      
  // OK
  return get_create_internal(i, j, k, l);
}

template<class T>
void Octree<T>::Branch::Accept(Visitor& v) {
  if (v.ShouldDescend(*this))
    for (auto i: {0, 1})
      for (auto j: {0, 1})
	for (auto k: {0, 1}) {
	  auto child = children[i][j][k];
	  if (child)
	    child->Accept(v);
	}
  
  v.Do(*this);
}

template<class T>
auto Octree<T>::Branch::get_create_internal(Int i, Int j, Int k, Int l) -> NodePtr {
  // Get the local index
  Int child_level = this->level - 1;
  Int li = (i >> child_level) & 1;
  Int lj = (j >> child_level) & 1;
  Int lk = (k >> child_level) & 1;
  auto child = children[li][lj][lk];
      
  if (!child) {
    // Create the child node
    Int ci = this->x + (li << child_level);
    Int cj = this->y + (lj << child_level);
    Int ck = this->z + (lk << child_level);
	
    if (child_level) {
      // child is branch
      child = children[li][lj][lk] = NodePtr(new Branch(ci, cj, ck, child_level));
    } else {
      // child is leaf
      child = children[li][lj][lk] = NodePtr(new Leaf(ci, cj, ck, child_level));
    }
  }
      
  if (child_level == l) {
    return child;
  } else {
    auto branch = dynamic_cast<Branch*>(child.get());
    return branch->get_create_internal(i, j, k, l);
  }
}
template<class T>
auto Octree<T>::Branch::get_nocreate_internal(Int i, Int j, Int k, Int l) -> NodePtr {
  // Get the local index
  Int child_level = this->level - 1;
  Int li = (i >> child_level) & 1;
  Int lj = (j >> child_level) & 1;
  Int lk = (k >> child_level) & 1;
  auto child = children[li][lj][lk];
      
  if (!child) {
    return child;
  }
      
  if (child_level == l) {
    return child;
  } else {
    auto branch = dynamic_cast<Branch*>(child.get());
    return branch->get_nocreate_internal(i, j, k, l);
  }
}
  
// Get a node without creating
template<class T>
auto Octree<T>::Leaf::Get(Int i, Int j, Int k, Int l) -> NodePtr {
  throw std::out_of_range("trying to get a child of a leaf node");
}
template<class T>
auto Octree<T>::Leaf::GetCreate(Int i, Int j, Int k, Int l) -> NodePtr {
  throw std::out_of_range("trying to get a child of a leaf node");
}
template<class T>
void Octree<T>::Leaf::Accept(Visitor& v) {
  v.Do(*this);
}

template<class T>
Octree<T>::Octree(Int level_) : level(level_) {
  if (level) {
    root = NodePtr(new Branch(0,0,0, level));
  } else {
    root = NodePtr(new Leaf(0,0,0,0));
  }
}
template<class T>
auto Octree<T>::Root() -> NodePtr {
  return root;
}
template<class T>
auto Octree<T>::Get(Int i, Int j, Int k, Int l) -> NodePtr {
  return root->Get(i,j,k,l);
}
template<class T>
auto Octree<T>::GetCreate(Int i, Int j, Int k, Int l) -> NodePtr {
  return root->GetCreate(i,j,k,l);
}

#endif
