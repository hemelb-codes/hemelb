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
  typedef std::shared_ptr<const Node> ConstNodePtr;

  class Branch;
  class Leaf;
  
  class Visitor;
  class ConstVisitor;
  
  class Node : public std::enable_shared_from_this<Node> {
  public:
    Node(Int i, Int j, Int k, Int l);
    T& Data();
    const T& Data() const;
    
    Int X() const;
    Int Y() const;
    Int Z() const;
    Int Level() const;
    
    // Get a node without creating
    virtual NodePtr Get(Int i, Int j, Int k, Int l) = 0;
    virtual ConstNodePtr Get(Int i, Int j, Int k, Int l) const = 0;

    // Get a node with creating
    virtual NodePtr GetCreate(Int i, Int j, Int k, Int l) = 0;
    virtual void Accept(Visitor& v) = 0;
    virtual void Accept(ConstVisitor& v) const = 0;
    
    virtual void Set(Int i, Int j, Int k, Int l, NodePtr n) = 0;

    template<class FuncT>
    void IterDepthFirst(Int bot, Int top, FuncT f) const {
      ConstLevelVisitor<FuncT> v(bot, top, f);
      Accept(v);
    }
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
    virtual ConstNodePtr Get(Int i, Int j, Int k, Int l) const;

    // Get a node without creating - returns null pointer if doensn't exist
    virtual NodePtr GetCreate(Int i, Int j, Int k, Int l);
    
    virtual void Accept(Visitor& v);
    virtual void Accept(ConstVisitor& v) const;
    virtual void Set(Int i, Int j, Int k, Int l, NodePtr n);

  private:
    NodePtr get_create_internal(Int i, Int j, Int k, Int l);
    NodePtr get_nocreate_internal(Int i, Int j, Int k, Int l);
    ConstNodePtr get_nocreate_internal(Int i, Int j, Int k, Int l) const;
    
    NodePtr children[2][2][2];
  };
  
  class Leaf : public Node {
  public:
    using Node::Node;
    // Get a node without creating
    virtual NodePtr Get(Int i, Int j, Int k, Int l);
    virtual ConstNodePtr Get(Int i, Int j, Int k, Int l) const;

    virtual NodePtr GetCreate(Int i, Int j, Int k, Int l);
    virtual void Accept(Visitor& v);
    virtual void Accept(ConstVisitor& v) const;
    virtual void Set(Int i, Int j, Int k, Int l, NodePtr n);

  };

  class Visitor {
  public:
    virtual void Arrive(NodePtr n) = 0;
    virtual void Depart(NodePtr n) = 0;
    virtual bool ShouldDescend(NodePtr n) {
      return true;
    }
  };
  
  class ConstVisitor {
  public:
    virtual void Arrive(ConstNodePtr n) = 0;
    virtual void Depart(ConstNodePtr n) = 0;
    virtual bool ShouldDescend(ConstNodePtr n) {
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
    virtual void Arrive(NodePtr node) {}
    virtual void Depart(NodePtr node) {
      if (node->Level() <= highest)
	f(node);
    }
    virtual bool ShouldDescend(NodePtr n) {
      return lowest < n->Level();
    }
    
  };
  template <class FuncT>
  class ConstLevelVisitor : public ConstVisitor {
    FuncT f;
    Int lowest;
    Int highest;
  public:
    ConstLevelVisitor(Int bot, Int top, FuncT func) : f(func), lowest(bot), highest(top) {}
    virtual void Arrive(ConstNodePtr node) {}
    virtual void Depart(ConstNodePtr node) {
      if (node->Level() <= highest)
	f(node);
    }
    virtual bool ShouldDescend(ConstNodePtr n) {
      return lowest < n->Level();
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
  template<class FuncT>
  void IterDepthFirst(Int bot, Int top, FuncT f) const {
    root->IterDepthFirst(bot, top, f);
  }
  
  Octree(Int level_);
  //T& GetValue(Int level, Int i, Int j, Int k);
  NodePtr Root();
  ConstNodePtr Root() const;
  // Get a node without creating
  NodePtr Get(Int i, Int j, Int k, Int l);
  ConstNodePtr Get(Int i, Int j, Int k, Int l) const;

  // Get a node with creating
  NodePtr GetCreate(Int i, Int j, Int k, Int l);

  void Set(Int i, Int j, Int k, Int l, NodePtr n);
  Int Level() const {
    return level;
  }
private:
  NodePtr root;
  Int level;
};

template<class T>
std::ostream& operator<<(std::ostream&, const Octree<T>& obj);
template<class T>
std::ostream& operator<<(std::ostream&, const typename Octree<T>::Node& obj);


#include "Oct.hpp"
#endif
