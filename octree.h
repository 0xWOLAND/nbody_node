#ifndef OCTREE_H
#define OCTREE_H
#include "imports.h"
#include "point.h"

struct Octree {
  bool isleaf;
  virtual ~Octree() { }
//  virtual bool isLeaf() const = 0;
  bool isLeaf() const{ return isleaf; }
};

struct OctreeInternal : Octree {
  Octree* child[8];
  Point pos;
  double mass;
  OctreeInternal(Point _pos) : pos(_pos), mass(0.0) {
	isleaf = false;
    bzero(child, sizeof(*child) * 8);
  }
  OctreeInternal(){
  }
  void assign(Point _pos){
	pos = _pos;
	mass = 0;
	isleaf = false;
    bzero(child, sizeof(*child) * 8);
  }
  bool isLeaf() const {
    return false;
  }
  virtual ~OctreeInternal() {
    for (int i = 0; i < 8; i++) {
      if (child[i] != NULL && !child[i]->isLeaf()) {
        delete child[i];
      }
    }
  }
};

struct Body : Octree {
  Point pos;
  Point vel;
  Point acc;
  double mass;
  Body() { isleaf = true; }
  bool isLeaf() const { return true; }
  ~Body(){}
};

#endif