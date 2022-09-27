#include "imports.h"
#include "point.h"
#include "octree.h"
#include "box.h"
#include "compute.h"

int nbodies=1200;
int ntimesteps=1;
int seed=0;
int numThreads=1;


int main(int argc, char** argv) {
  std::cout.setf(std::ios::right|std::ios::scientific|std::ios::showpoint);

  std::cerr << "configuration: "
            << nbodies << " bodies, "
            << ntimesteps << " time steps" << std::endl << std::endl;

  run(nbodies, ntimesteps, seed);
}
