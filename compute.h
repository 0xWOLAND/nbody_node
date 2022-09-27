#ifndef COMPUTE_H
#define COMPUTE_H
#include "imports.h"
#include "point.h"
#include "octree.h"
#include "box.h"
#include "compute.h"

typedef std::vector<Body> Bodies;

std::ostream& operator<<(std::ostream& os, const Point& p) {
  os << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
  return os;
}
std::ostream& operator<<(std::ostream& os, const Body& b) {
  os << "(pos:" << b.pos
     << " vel:" << b.vel
     << " acc:" << b.acc
     << " mass:" << b.mass << ")";
  return os;
}

std::ostream& operator<<(std::ostream& os, const BoundingBox& b) {
  os << "(min:" << b.min << " max:" << b.max << ")";
  return os;
}
struct Config {
  const double dtime; // length of one time step
  const double eps; // potential softening parameter
  const double tol; // tolerance for stopping recursion, <0.57 to bound error
  const double dthf, epssq, itolsq;
  Config() :
    dtime(0.5),
    eps(0.05),
    tol(0.025),
    dthf(dtime * 0.5),
    epssq(eps * eps),
    itolsq(1.0 / (tol * tol))  { }
};

Config config;

inline int getIndex(const Point& a, const Point& b) {
  int index = 0;
  if (a.x < b.x)
    index += 1;
  if (a.y < b.y)
    index += 2;
  if (a.z < b.z)
    index += 4;
  return index;
}

inline void updateCenter(Point& p, int index, double radius) {
  for (int i = 0; i < 3; i++) {
    double v = (index & (1 << i)) > 0 ? radius : -radius;
    p[i] += v;
  }
}

template<typename T, int N>
class Stack{
	T elem[N+1];
	int sz;
public:
	typedef T* iterator;
public:
	Stack() : sz(0){
	}
	iterator begin(){
		return elem;
	}
	iterator end(){
		return elem+sz;
	}
	bool empty()const{
		return sz == 0;
	}
	int size()const{
		return sz;
	}
	T& operator[](int i){
		return elem[i];
	}
	void push_back(const T& e){
		if(sz >= N) throw std::runtime_error("too many elements");
		elem[sz++] = e;
	}
	void pop_back(){
		sz--;
	}
	const T& back()const{
		return elem[sz-1];
	}
};

struct BuildOctree {
  typedef int tt_does_not_need_stats;

  OctreeInternal* root;
  double root_radius;

  BuildOctree(OctreeInternal* _root, double radius) :
    root(_root),
    root_radius(radius) { }

  template<typename Context>
  void operator()(Body* b, Context&) {
    insert(b, root, root_radius);
  }

  void insert(Body* b, OctreeInternal* node, double radius) {
    int index = getIndex(node->pos, b->pos);
    assert(!node->isLeaf());

    Octree *child = node->child[index];
    
    if (child == NULL) {
      node->child[index] = b;
      return;
    }
    
    radius *= 0.5;
    if (child->isLeaf()) {
      Body* n = static_cast<Body*>(child);
      Point new_pos(node->pos);
      updateCenter(new_pos, index, radius);
      OctreeInternal* new_node = new OctreeInternal(new_pos);

      assert(n->pos != b->pos);
      
      insert(b, new_node, radius);
      insert(n, new_node, radius);
      node->child[index] = new_node;
    } else {
      OctreeInternal* n = static_cast<OctreeInternal*>(child);
      insert(b, n, radius);
    }
  }
};

struct ComputeCenterOfMass {
  typedef int tt_does_not_need_stats;
  OctreeInternal* root;
  int cnt;

  ComputeCenterOfMass(OctreeInternal* _root) : root(_root) { }

  void operator()() {
    cnt = 0;
    root->mass = recurse(root);
  }

private:
  double recurse(OctreeInternal* node) {
    double mass = 0.0;
    int index = 0;
    Point accum;
    
    for (int i = 0; i < 8; i++) {
      Octree* child = node->child[i];
      if (child == NULL)
        continue;

      // Reorganize leaves to be denser up front 
      if (index != i) {
        node->child[index] = child;
        node->child[i] = NULL;
      }
      index++;
      
      double m;
      cnt++;
      const Point* p;
      if (child->isLeaf()) {
        Body* n = static_cast<Body*>(child);
        m = n->mass;
        p = &n->pos;
      } else {
        OctreeInternal* n = static_cast<OctreeInternal*>(child);
        m = recurse(n);
        p = &n->pos;
      }

      mass += m;
      for (int j = 0; j < 3; j++){
        accum[j] += (*p)[j] * m;
      }
    }

    node->mass = mass;
    
    if (mass > 0.0) {
      double inv_mass = 1.0 / mass;
      for (int j = 0; j < 3; j++){
        node->pos[j] = accum[j] * inv_mass;
      }
    }

    return mass;
  }
};

struct ComputeForces {
  typedef int tt_does_not_need_context;

  OctreeInternal* top;
  double diameter;
  double root_dsq;

  ComputeForces(OctreeInternal* _top, double _diameter) :
    top(_top),
    diameter(_diameter) {
    root_dsq = diameter * diameter * config.itolsq;
  }
  
  template<typename Context>
  void operator()(Body* bb, Context&) {
    Body& b = *bb;
    Point p = b.acc;
    for (int i = 0; i < 3; i++){
      b.acc[i] = 0;
    }

    //recurse(b, top, root_dsq);
    iterate(b, root_dsq);
    for (int i = 0; i < 3; i++){
      b.vel[i] += (b.acc[i] - p[i]) * config.dthf;
    }
  }

  void recurse(Body& b, Body* node, double dsq) {
    Point p;
    for (int i = 0; i < 3; i++){
      p[i] = node->pos[i] - b.pos[i];
    }

    double psq = p.x * p.x + p.y * p.y + p.z * p.z;
    psq += config.epssq;
    double idr = 1 / sqrt(psq);
    // b.mass is fine because every body has the same mass
    double nphi = b.mass * idr;
    double scale = nphi * idr * idr;
    for (int i = 0; i < 3; i++){
      b.acc[i] += p[i] * scale;
    }
  }

  struct Frame {
    double dsq;
    OctreeInternal* node;
	Frame(){}
    Frame(OctreeInternal* _node, double _dsq) : dsq(_dsq), node(_node) { }
  };

  void iterate(Body& b, double root_dsq) {
//    std::vector<Frame> stack;
    Stack<Frame, 100> stack;
    stack.push_back(Frame(top, root_dsq));

    Point p;
    while (!stack.empty()) {
      Frame f = stack.back();
      stack.pop_back();

      for (int i = 0; i < 3; i++){
        p[i] = f.node->pos[i] - b.pos[i];
      }

      double psq = p.x * p.x + p.y * p.y + p.z * p.z;
      if (psq >= f.dsq) {
        // Node is far enough away, summarize contribution
        psq += config.epssq;
        double idr = 1 / sqrt(psq);
        double nphi = f.node->mass * idr;
        double scale = nphi * idr * idr;
        for (int i = 0; i < 3; i++){
#pragma tls forkpoint id 1 cannot
          b.acc[i] += p[i] * scale;
#pragma tls joinpoint id 1
        }
#pragma tls barrierpoint id 1
        continue;
      }

      double dsq = f.dsq * 0.25;

      for (int i = 0; i < 8; i++) {
        Octree *next = f.node->child[i];
        if (next == NULL)
          break;
        if (next->isLeaf()) {
          // Check if it is me
          if (&b != next) {
            recurse(b, static_cast<Body*>(next), dsq);
          }
        } else {
          stack.push_back(Frame(static_cast<OctreeInternal*>(next), dsq));
        }
      }
    }
  }
};


struct ReduceBoxes {
  // NB: only correct when run sequentially or tree-like reduction
  typedef int tt_does_not_need_stats;
  BoundingBox& initial;

  ReduceBoxes(BoundingBox& _initial): initial(_initial) { }

  template<typename Context>
  void operator()(Body* b, Context&) {
    initial.merge(b->pos);
  }
};

double nextDouble() {
  return rand() / (double) RAND_MAX;
}

void generateInput(Bodies& bodies, int nbodies, int seed) {
  double v, sq, scale;
  Point p;
  double PI = acos(-1.0);

  srand(seed);

  double rsc = (3 * PI) / 16;
  double vsc = sqrt(1.0 / rsc);

  for (int body = 0; body < nbodies; body++) {
    double r = 1.0 / sqrt(pow(nextDouble() * 0.999, -2.0 / 3.0) - 1);
    do {
      for (int i = 0; i < 3; i++)
        p[i] = nextDouble() * 2.0 - 1.0;
      sq = p.x * p.x + p.y * p.y + p.z * p.z;
    } while (sq > 1.0);
    scale = rsc * r / sqrt(sq);

    Body b;
    b.mass = 1.0 / nbodies;
    for (int i = 0; i < 3; i++)
      b.pos[i] = p[i] * scale;

    do {
      p.x = nextDouble();
      p.y = nextDouble() * 0.1;
    } while (p.y > p.x * p.x * pow(1 - p.x * p.x, 3.5));
    v = p.x * sqrt(2.0 / sqrt(1 + r * r));
    do {
      for (int i = 0; i < 3; i++)
        p[i] = nextDouble() * 2.0 - 1.0;
      sq = p.x * p.x + p.y * p.y + p.z * p.z;
    } while (sq > 1.0);
    scale = vsc * v / sqrt(sq);
    for (int i = 0; i < 3; i++)
      b.vel[i] = p[i] * scale;

    bodies.push_back(b);
  }
}

struct timeval tm0, tm1, tm2;

double gettime(struct timeval t1, struct timeval t2){
	return (t2.tv_sec-t1.tv_sec)*1000.0+(t2.tv_usec-t1.tv_usec)/1000.0;
}


void run(int nbodies, int ntimesteps, int seed) {
  static Bodies bodies;
  generateInput(bodies, nbodies, seed);

  for (int step = 0; step < ntimesteps; step++) {
    // Do tree building sequentially

    size_t i, k, n=bodies.size(), block=64;

    gettimeofday(&tm1, NULL);
    BoundingBox box;
    ReduceBoxes reduceBoxes(box);
	for(Bodies::iterator it = bodies.begin(); it != bodies.end(); ++it){
        (ReduceBoxes(box))(&*it,bodies);
	}

	OctreeInternal* top = new OctreeInternal(box.center());
	for(k=0;k<block;k++){
	    for(i=n*k/block; i<n*(k+1)/block; i++){
		BuildOctree(top, box.radius())(&bodies[i],bodies);
	    }
	}

    ComputeCenterOfMass computeCenterOfMass(top);

    computeCenterOfMass();

    for(k=0;k<block;k++){
        for(i=n*k/block; i<n*(k+1)/block; i++){
   	    ComputeForces(top, box.diameter())(&bodies[i],bodies);
	}
    }

    std::cout 
      << "Timestep " << step
      << " Center of Mass = " << top->pos << "\n";
    delete top;
    cout << "CNT: " << top->cnt << "\n";
  }
}
#endif