// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_OCT_HPP
#define HEMELBSETUPTOOL_OCT_HPP

#include "Oct.h"
#include <stdexcept>
#include <string>
#include <ostream>

template<class T>
Octree<T>::Node::Node(Int i, Int j, Int k, Int l) : x(i), y(j), z(k), level(l) {
}
template<class T>
T& Octree<T>::Node::Data() {
  return value;
}
template<class T>
const T& Octree<T>::Node::Data() const {
  return value;
}

template<class T>
auto Octree<T>::Node::X() const -> Int{
  return x;
}
template<class T>
auto Octree<T>::Node::Y() const -> Int{
  return y;
}
template<class T>
auto Octree<T>::Node::Z() const -> Int{
  return z;
}
template<class T>
auto Octree<T>::Node::Level() const -> Int{
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
auto Octree<T>::Branch::Get(Int i, Int j, Int k, Int l) const -> ConstNodePtr {
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
void Octree<T>::Branch::Set(Int i, Int j, Int k, Int l, NodePtr n) {
	Int pl = l + 1;
	Int mask = ~(1 << l);
	Int pi = i & mask;
	Int pj = j & mask;
	Int pk = k & mask;

	auto parent = std::dynamic_pointer_cast<Branch>(GetCreate(pi, pj, pk, pl));
	if (!parent)
		throw std::out_of_range("parent of request index not a branch node");

	Int li = (i >> l) & 1;
	Int lj = (j >> l) & 1;
	Int lk = (k >> l) & 1;
	parent->children[li][lj][lk] = n;
}

template<class T>
void Octree<T>::Branch::Accept(Visitor& v) {
  v.Arrive(this->shared_from_this());
  if (v.ShouldDescend(this->shared_from_this()))
    for (auto i: {0, 1})
      for (auto j: {0, 1})
	for (auto k: {0, 1}) {
	  auto child = children[i][j][k];
	  if (child)
	    child->Accept(v);
	}
  
  v.Depart(this->shared_from_this());
}
template<class T>
void Octree<T>::Branch::Accept(ConstVisitor& v) const {
  v.Arrive(this->shared_from_this());
  if (v.ShouldDescend(this->shared_from_this()))
    for (auto i: {0, 1})
      for (auto j: {0, 1})
	for (auto k: {0, 1}) {
	  auto child = children[i][j][k];
	  if (child)
	    child->Accept(v);
	}

  v.Depart(this->shared_from_this());
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
template<class T>
auto Octree<T>::Branch::get_nocreate_internal(Int i, Int j, Int k, Int l) const -> ConstNodePtr {
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
auto Octree<T>::Leaf::Get(Int i, Int j, Int k, Int l) const -> ConstNodePtr {
  throw std::out_of_range("trying to get a child of a leaf node");
}
template<class T>
auto Octree<T>::Leaf::GetCreate(Int i, Int j, Int k, Int l) -> NodePtr {
  throw std::out_of_range("trying to get a child of a leaf node");
}
template<class T>
void Octree<T>::Leaf::Set(Int i, Int j, Int k, Int l, NodePtr n) {
	throw std::out_of_range("trying to set a child of a leaf node");
}

template<class T>
void Octree<T>::Leaf::Accept(Visitor& v) {
  v.Arrive(this->shared_from_this());
  v.Depart(this->shared_from_this());
}

template<class T>
void Octree<T>::Leaf::Accept(ConstVisitor& v) const {
  v.Arrive(this->shared_from_this());
  v.Depart(this->shared_from_this());
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
auto Octree<T>::Root() const -> ConstNodePtr {
	return root;
}

template<class T>
auto Octree<T>::Get(Int i, Int j, Int k, Int l) -> NodePtr {
  return root->Get(i,j,k,l);
}
template<class T>
auto Octree<T>::Get(Int i, Int j, Int k, Int l) const -> ConstNodePtr {
  return root->Get(i,j,k,l);
}

template<class T>
auto Octree<T>::GetCreate(Int i, Int j, Int k, Int l) -> NodePtr {
  return root->GetCreate(i,j,k,l);
}
template<class T>
void Octree<T>::Set(Int i, Int j, Int k, Int l, NodePtr n) {
  root->Set(i,j,k,l, n);
}

template<class T>
void WriteNodeData(std::ostream& os, const T& data);

template<class T>
class Writer : public Octree<T>::Visitor {
public:
  typedef Octree<T> Tree;
  typedef typename Tree::Node Node;
  
  Writer(std::ostream& stream) : os(stream), padding(0) {}
  
  virtual void Arrive(Node& node) {
    os << std::string(padding, ' ') << node.X() << ", " << node.Y() << ", " << node.Z()
       << " = ";
    WriteNodeData(os, node.Data());
    os << std::endl;
    padding++;
  }
  virtual void Depart(Node& node) {
    padding--;
  }
  
private:
  std::ostream& os;
  int padding;
};

template<class T>
std::ostream& operator<<(std::ostream& os, Octree<T>& tree) {
  os << "Octree levels = " << tree.Level() << std::endl;
  Writer<T> w(os);
  tree.Root()->Accept(w);
  return os;
}
template<class T>
std::ostream& operator<<(std::ostream& os, typename Octree<T>::Node& node) {
  Writer<T> w(os);
  node.Accept(w);
  return os;
}

template<class T>
void WriteNodeData(std::ostream& os, const T& data) {
  os << data;
}

#endif
