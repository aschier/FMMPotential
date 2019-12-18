#ifndef buildtree_h
#define buildtree_h
#include "exafmm.h"
#include <list>

namespace exafmm {
  //! Get bounding box of bodies
  void getBounds(Bodies & bodies, real_t & R0, vec3 & X0) {
    vec3 Xmin = bodies[0].X;
    vec3 Xmax = bodies[0].X;
    for (size_t b=0; b<bodies.size(); b++) {
      Xmin = min(bodies[b].X, Xmin);
      Xmax = max(bodies[b].X, Xmax);
    }
    X0 = (Xmax + Xmin) / 2;
    R0 = fmax(max(X0-Xmin), max(Xmax-X0));
    R0 *= 1.00001;
  }

  struct BuildQueueItem {
	  int begin;
	  int end;
	  int cell_index;
	  const vec3 X;
	  int level;
	  bool direction;
  };

  //! Build cells of tree adaptively using a top-down approach based on tail-recursion implemented using a loop
  void buildCells(Bodies &bodies_in, Bodies &buffer_in, Cells & cells, const vec3 & X0, real_t R0) {
	  std::list<BuildQueueItem> queue;
	  queue.push_front({
		  0,
		  bodies_in.size(),
		  0,
		  X0,
		  0,
		  false
	  });
	  while(!queue.empty()) {
		  BuildQueueItem item = queue.front();
		  queue.pop_front();
		  Bodies * bodies = &bodies_in;
		  Bodies * buffer = &buffer_in;
		  if(item.direction) {
			  std::swap(bodies, buffer);
		  }
		  //! Create a tree cell
		  cells[item.cell_index].body = &bodies->at(item.begin);
		  if(item.direction) {
			  cells[item.cell_index].body = &buffer->at(item.begin);
		  }
		  cells[item.cell_index].numBodies = item.end - item.begin;
		  cells[item.cell_index].numChilds = 0;
		  cells[item.cell_index].X = item.X;
		  cells[item.cell_index].R = R0 / (1 << item.level);
		  //! If cell is a leaf
		  if(item.end - item.begin <= NCRIT) {
			  if(item.direction) {
				  for(int i = item.begin; i < item.end; i++) {
					  buffer->at(i).X = bodies->at(i).X;
					  buffer->at(i).q = bodies->at(i).q;
#ifdef EXAFMM_INDEXED_BODIES
					  buffer->at(i).index = bodies->at(i).index;
#endif
				  }
			  }
			  continue;
		  }
		  //! Count number of bodies in each octant
		  int size[8] = {0,0,0,0,0,0,0,0};
		  vec3 x;
		  for(int i = item.begin; i < item.end; i++) {
			  x = bodies->at(i).X;
			  int octant = (x[0] > item.X[0]) + ((x[1] > item.X[1]) << 1) + ((x[2] > item.X[2]) << 2);
			  size[octant]++;
		  }
		  //! Exclusive scan to get offsets
		  int offset = item.begin;
		  int offsets[8], counter[8];
		  for(int i = 0; i < 8; i++) {
			  offsets[i] = offset;
			  offset += size[i];
			  if(size[i] > 0) {
				  cells[item.cell_index].numChilds++;
			  }
		  }
		  //! Sort bodies by octant
		  for(int i = 0; i < 8; i++) {
			  counter[i] = offsets[i];
		  }
		  for(int i = item.begin; i < item.end; i++) {
			  x = bodies->at(i).X;
			  int octant = (x[0] > item.X[0]) + ((x[1] > item.X[1]) << 1) + ((x[2] > item.X[2]) << 2);
			  buffer->at(counter[octant]).X = bodies->at(i).X;
			  buffer->at(counter[octant]).q = bodies->at(i).q;
#ifdef EXAFMM_INDEXED_BODIES
			  buffer->at(counter[octant]).index = bodies->at(i).index;
#endif
			  counter[octant]++;
		  }
		  //! Loop over children
		  vec3 Xchild;
		  int first_child_index = cells.size();
		  cells.resize(cells.size() + cells[item.cell_index].numChilds);
		  cells[item.cell_index].child = &cells[first_child_index];
		  int c = 0;
		  for(int i = 0; i < 8; i++) {
			  Xchild = item.X;
			  real_t r = R0 / (1 << (item.level + 1));
			  for(int d = 0; d < 3; d++) {
				  Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);
			  }
			  if(size[i] > 0) {
				  queue.push_back({
					  offsets[i],
					  offsets[i] + size[i],
					  first_child_index + c,
					  Xchild,
					  item.level+1,
					  !item.direction
				  });
				  c++;
			  }
		  }
	  }
  }

  Cells buildTree(Bodies & bodies) {
    real_t R0;
    vec3 X0;
    getBounds(bodies, R0, X0);
    Bodies buffer = bodies;
    Cells cells(1);
    cells.reserve(bodies.size());
    buildCells(bodies, buffer, cells, X0, R0);
    return cells;
  }
}
#endif
