// -*- mode: c++; -*-
#ifndef HEMELBSETUPTOOL_OCT_H
#define HEMELBSETUPTOOL_OCT_H

#include <memory>
 
template<class T>
class Octree {
public:
  // int type for 3D indices into our box
  typedef uint16_t Int;
    
  class Node;
  typedef std::shared_ptr<Node> NodePtr;
  class Branch;
  class Leaf;
  
  class Visitor;
  
  class Node {
  public:
    Node(Int i, Int j, Int k, Int l);
    T& Data();
    
    Int X();
    Int Y();
    Int Z();
    Int Level();
    
    // Get a node without creating
    virtual NodePtr Get(Int i, Int j, Int k, Int l) = 0;
    // Get a node with creating
    virtual NodePtr GetCreate(Int i, Int j, Int k, Int l) = 0;
    virtual void Accept(Visitor& v) = 0;
    
    template<class FuncT>
    void IterDepthFirst(Int bot, Int top, FuncT f) {
      LevelVisitor<FuncT> v(bot, top, f);
      Accept(v);
    }
    template<class FuncT>
    void IterDepthFirst(FuncT f) {
      IterDepthFirst(0, level, f);
    }
    template<class FuncT>
    void IterDepthFirst(Int bot, FuncT f) {
      IterDepthFirst(bot, level, f);
    }
    
  protected:
    Int x, y, z, level;
    T value;
  };
  
  class Branch : public Node {
  public:
    using Node::Node;
    
    // Get a node without creating - returns null pointer if doensn't exist
    virtual NodePtr Get(Int i, Int j, Int k, Int l);
      
    // Get a node without creating - returns null pointer if doensn't exist
    virtual NodePtr GetCreate(Int i, Int j, Int k, Int l);
    
    virtual void Accept(Visitor& v);
  private:
    NodePtr get_create_internal(Int i, Int j, Int k, Int l);
    NodePtr get_nocreate_internal(Int i, Int j, Int k, Int l);
    
    NodePtr children[2][2][2];
  };
  
  class Leaf : public Node {
  public:
    using Node::Node;
    // Get a node without creating
    virtual NodePtr Get(Int i, Int j, Int k, Int l);
    virtual NodePtr GetCreate(Int i, Int j, Int k, Int l);
    virtual void Accept(Visitor& v);
  };

  class Visitor {
  public:
    virtual void Do(Node& n) = 0;
    virtual bool ShouldDescend(Node& n) {
      return true;
    }
  };
  
  template <class FuncT>
  class LevelVisitor : public Visitor {
    FuncT f;
    Int lowest;
    Int highest;
  public:
    LevelVisitor(Int bot, Int top, FuncT func) : f(func), lowest(bot), highest(top) {}
    virtual void Do(Node& node) {
      f(node);
    }
    virtual bool ShouldDescend(Node& n) {
      return lowest < n.Level();
    }
    
  };
  
  template<class FuncT>
  void IterDepthFirst(Int bot, Int top, FuncT f) {
    root->IterDepthFirst(bot, top, f);
  }
  template<class FuncT>
  void IterDepthFirst(FuncT f) {
    IterDepthFirst(0, level, f);
  }
  template<class FuncT>
  void IterDepthFirst(Int bot, FuncT f) {
    IterDepthFirst(bot, level, f);
  }
  
  Octree(Int level_);
  //T& GetValue(Int level, Int i, Int j, Int k);
  NodePtr Root();
  // Get a node without creating
  NodePtr Get(Int i, Int j, Int k, Int l);
  // Get a node with creating
  NodePtr GetCreate(Int i, Int j, Int k, Int l);

  Int Level() const {
    return level;
  }
private:
  NodePtr root;
  Int level;
};

#include "Oct.hpp"
#endif
