#ifndef LMP_DECOMP_GRID_H
#define LMP_DECOMP_GRID_H

#include "subgrid.h"

#include <vector>

#ifdef NUFEB_DEBUG_COMM
#include <fstream>
#include <sstream>
#endif

namespace LAMMPS_NS {
template <class Derived>
class DecompGrid {
 public:
  void setup_exchange(const Grid<double, 3> &grid,
		      const Box<int, 3> &box,
		      const std::array<bool, 3> &periodic) {
    Derived *derived = static_cast<Derived *>(this);

#ifdef NUFEB_DEBUG_COMM
    std::stringstream ss;
    ss << "debug_comm_" << derived->comm->me << ".txt";
    debug.open(ss.str().c_str(), std::fstream::app);
    debug << ">>> Entering setup_exchange" << std::endl;
    debug << "Grid:    [origin](" << grid.get_origin()[0] << ", " << grid.get_origin()[1] << ", " << grid.get_origin()[2]
	  << ") [dimensions](" << grid.get_dimensions()[0] << ", " << grid.get_dimensions()[1] << ", " << grid.get_dimensions()[2] << ")" << std::endl;
    debug << "Subgrid: [lower] (" << box.lower[0] << ", " << box.lower[1] << ", " << box.lower[2]
	  << ") [upper](" << box.upper[0] << ", " << box.upper[1] << ", " << box.upper[2] << ")" << std::endl;
#endif

    // communicate grid extent
    std::vector<int> boxlo(3 * derived->comm->nprocs); 
    MPI_Allgather(const_cast<int *>(&box.lower[0]), 3, MPI_INT, boxlo.data(), 3, MPI_INT, derived->world);
    std::vector<int> boxhi(3 * derived->comm->nprocs);
    MPI_Allgather(const_cast<int *>(&box.upper[0]), 3, MPI_INT, boxhi.data(), 3, MPI_INT, derived->world);

    clear();
    recv_begin.resize(derived->comm->nprocs);
    send_begin.resize(derived->comm->nprocs);
    recv_end.resize(derived->comm->nprocs);
    send_end.resize(derived->comm->nprocs);
    requests.resize(2 * derived->comm->nprocs);
    // look for intersetions
    Box<int, 3> ext_box = extend(box);
    Subgrid<double, 3> subgrid(grid, ext_box);
    int epc = derived->get_elem_per_cell();
    for (int p = 0; p < derived->comm->nprocs; p++) {
      recv_begin[p] = recv_cells.size() * epc;
      send_begin[p] = send_cells.size() * epc;
      Box<int, 3> other(&boxlo[3 * p], &boxhi[3 * p]);
      if (p != derived->comm->me) {
#ifdef NUFEB_DEBUG_COMM
      debug << "Checking for intersections with proc " << p
	    << ", box: [lower](" << box.lower[0] << ", " << box.lower[1] << ", " << box.lower[2]
	    << ") [upper](" << box.upper[0] << ", " << box.upper[1] << ", " << box.upper[2] << ")" << std::endl;
#endif
	// identify which cell we are going to send and receive
	setup_comm_cells(subgrid, box, other);
	// check for periodic boundary conditions
	if (ext_box.lower[0] < 0 && periodic[0]) {
	  setup_comm_cells(subgrid, box, translate(other, {-grid.get_dimensions()[0], 0, 0}));
	}
	if (ext_box.upper[0] > grid.get_dimensions()[0] && periodic[0]) {      
	  setup_comm_cells(subgrid, box, translate(other, {grid.get_dimensions()[0], 0, 0}));
	}
	if (ext_box.lower[1] < 0 && periodic[1]) {
	  setup_comm_cells(subgrid, box, translate(other, {0, -grid.get_dimensions()[1], 0}));
	}
	if (ext_box.upper[1] > grid.get_dimensions()[1] && periodic[1]) {      
	  setup_comm_cells(subgrid, box, translate(other, {0, grid.get_dimensions()[1], 0}));
	}
	if (ext_box.lower[2] < 0 && periodic[2]) {
	  setup_comm_cells(subgrid, box, translate(other, {0, 0, -grid.get_dimensions()[2]}));
	}
	if (ext_box.upper[2] > grid.get_dimensions()[2] && periodic[2]) {      
	  setup_comm_cells(subgrid, box, translate(other, {0, 0, grid.get_dimensions()[2]}));
	}
      }
      recv_end[p] = recv_cells.size() * epc;
      send_end[p] = send_cells.size() * epc;
    }
    recv_buff.resize(recv_end.back());
    send_buff.resize(send_end.back());

#ifdef NUFEB_DEBUG_COMM
    debug << "<<< Leaving setup_exchange" << std::endl;
    debug.close();
#endif
  }

  void exchange() {
    Derived *derived = static_cast<Derived *>(this);
    int epc = derived->get_elem_per_cell();

    // pack data to send buffer
    derived->pack_cells(send_cells.begin(), send_cells.end(), send_buff.begin());

#ifdef NUFEB_DEBUG_COMM
    std::stringstream ss;
    ss << "debug_comm_" << derived->comm->me << ".txt";
    debug.open(ss.str().c_str(), std::fstream::app);
    debug << ">>> Entering exchange" << std::endl;
    debug << "Packing cells: ";
    {
      auto data = send_buff.begin();
      for (auto cell = send_cells.begin(); cell != send_cells.end(); ++cell) {
        if (cell != send_cells.begin())
          debug << ", ";
        debug << "[" << *cell << "](";
        for (int i = 0; i < epc; i++, ++data) {
	  if (i != 0)
	    debug << ", ";
	  debug << *data;
	}
	debug << ")";
      }
      debug << std::endl;
    }
#endif

    // send and recv grid data
    int nrequests = 0;
    for (int p = 0; p < derived->comm->nprocs; p++) {
      if (p == derived->comm->me)
	continue;
      if (recv_begin[p] < recv_end[p]) {
	MPI_Irecv(&recv_buff[recv_begin[p]], recv_end[p] - recv_begin[p], MPI_DOUBLE, p, 0, derived->world, &requests[nrequests++]);
      }
      if (send_begin[p] < send_end[p]) {
	MPI_Isend(&send_buff[send_begin[p]], send_end[p] - send_begin[p], MPI_DOUBLE, p, 0, derived->world, &requests[nrequests++]);
      }
    }
    // wait for all MPI requests
    MPI_Waitall(nrequests, requests.data(), MPI_STATUS_IGNORE);
    // unpack data from recv buffer
    derived->unpack_cells(recv_cells.begin(), recv_cells.end(), recv_buff.begin());

#ifdef NUFEB_DEBUG_COMM
    debug << "Unpacking cells: ";
    {
      auto data = recv_buff.begin();
      for (auto cell = recv_cells.begin(); cell != recv_cells.end(); ++cell) {
        if (cell != recv_cells.begin())
          debug << ", ";
        debug << "[" << *cell << "](";
        for (int i = 0; i < epc; i++, ++data) {
	  if (i != 0)
	    debug << ", ";
	  debug << *data;
	}
	debug << ")";
      }
      debug << std::endl;
    }
    debug << "<<< Leaving exchange" << std::endl;
    debug.close();
#endif
  }


  void migrate(const Grid<double, 3> &grid, const Box<int, 3> &from, const Box<int, 3> &to) {
    migrate(grid, from, to, from, to);
  }

  void migrate(const Grid<double, 3> &grid, const Box<int, 3> &from, const Box<int, 3> &to, const Box<int, 3> &from_base, const Box<int, 3> &to_base) {
    // TODO: check if from_base contains from and to_base contains to
    Derived *derived = static_cast<Derived *>(this);

    std::vector<Box<int, 3>> old_boxes(derived->comm->nprocs);
    std::vector<Box<int, 3>> new_boxes(derived->comm->nprocs);

    std::vector<int> old_boxlo(3 * derived->comm->nprocs); 
    MPI_Allgather(const_cast<int *>(&from.lower[0]), 3, MPI_INT, old_boxlo.data(), 3, MPI_INT, derived->world);
    std::vector<int> old_boxhi(3 * derived->comm->nprocs);
    MPI_Allgather(const_cast<int *>(&from.upper[0]), 3, MPI_INT, old_boxhi.data(), 3, MPI_INT, derived->world);

    std::vector<int> new_boxlo(3 * derived->comm->nprocs); 
    MPI_Allgather(const_cast<int *>(&to.lower[0]), 3, MPI_INT, new_boxlo.data(), 3, MPI_INT, derived->world);
    std::vector<int> new_boxhi(3 * derived->comm->nprocs);
    MPI_Allgather(const_cast<int *>(&to.upper[0]), 3, MPI_INT, new_boxhi.data(), 3, MPI_INT, derived->world);

    for (int p = 0; p < derived->comm->nprocs; p++) {
      old_boxes[p] = Box<int, 3>(&old_boxlo[3 * p], &old_boxhi[3 * p]);
      new_boxes[p] = Box<int, 3>(&new_boxlo[3 * p], &new_boxhi[3 * p]);
    }
    
    std::vector<int> mig_recv_cells;
    std::vector<int> mig_send_cells;
    std::vector<int> mig_recv_begin(derived->comm->nprocs);
    std::vector<int> mig_send_begin(derived->comm->nprocs);
    std::vector<int> mig_recv_end(derived->comm->nprocs);
    std::vector<int> mig_send_end(derived->comm->nprocs);
    std::vector<MPI_Request> mig_requests(2 * derived->comm->nprocs);

    int me = derived->comm->me;
    int epc = derived->get_elem_per_cell();
    Subgrid<double, 3> from_subgrid(grid, from_base);
    Subgrid<double, 3> to_subgrid(grid, to_base);
    for (int p = 0; p < derived->comm->nprocs; p++) {
      mig_recv_begin[p] = mig_recv_cells.size() * epc;
      Box<int, 3> recv_box = intersect(new_boxes[me], old_boxes[p]);
      if (!is_empty(recv_box)) {
	add_cells(to_subgrid, recv_box, mig_recv_cells);
      }
      mig_recv_end[p] = mig_recv_cells.size() * epc;
      mig_send_begin[p] = mig_send_cells.size() * epc;
      Box<int, 3> send_box = intersect(old_boxes[me], new_boxes[p]);
      if (!is_empty(send_box)) {
	add_cells(from_subgrid, send_box, mig_send_cells);
      }
      mig_send_end[p] = mig_send_cells.size() * epc;
    }

    std::vector<double> mig_send_buff(mig_send_end.back());
    derived->pack_cells(mig_send_cells.begin(), mig_send_cells.end(), mig_send_buff.begin());

    int nrequests = 0;
    std::vector<double> mig_recv_buff(mig_recv_end.back());
    for (int p = 0; p < derived->comm->nprocs; p++) {
      if (p == derived->comm->me) {
	if (mig_recv_begin[p] != mig_recv_end[p]) {
	  std::copy(mig_send_buff.begin() + mig_send_begin[p],
		    mig_send_buff.begin() + mig_send_end[p],
		    mig_recv_buff.begin() + mig_recv_begin[p]);
	}
      }
      else {
	if (mig_recv_begin[p] != mig_recv_end[p]) {
	  MPI_Irecv(&mig_recv_buff[mig_recv_begin[p]], mig_recv_end[p] - mig_recv_begin[p], MPI_DOUBLE, p, 0, derived->world, &mig_requests[nrequests++]);
	}
	if (mig_send_begin[p] != mig_send_end[p]) {
	  MPI_Isend(&mig_send_buff[mig_send_begin[p]], mig_send_end[p] - mig_send_begin[p], MPI_DOUBLE, p, 0, derived->world, &mig_requests[nrequests++]);
	}
      }
    }
    MPI_Waitall(nrequests, mig_requests.data(), MPI_STATUS_IGNORE);
    derived->resize(to_subgrid);
    derived->unpack_cells(mig_recv_cells.begin(), mig_recv_cells.end(), mig_recv_buff.begin());
  }

 private:
  void add_cells(const Subgrid<double, 3> &subgrid, const Box<int, 3> &box, std::vector<int> &cells)
  {
    for (int k = box.lower[2]; k < box.upper[2]; k++) {
      for (int j = box.lower[1]; j < box.upper[1]; j++) {
	for (int i = box.lower[0]; i < box.upper[0]; i++) {
	  cells.push_back(subgrid.get_linear_index({i, j, k}));
#ifdef NUFEB_DEBUG_COMM
	  if (i != box.lower[0] || j != box.lower[1] || k != box.lower[2])
	    debug << ", ";
	  debug << "[" << subgrid.get_linear_index({i, j, k}) 
		<< "](" << i << ", " << j << ", " << k << ")";
#endif
	}
      }
    }
  }

  bool check_intersection(const Box<int, 3> &g)
  {
    // check if the intersection is empty
    if (is_empty(g))
      return false;
    // check if the intersection is a corner
    int n[3];
    for (int i = 0; i < 3; i++)
      n[i] = g.upper[i] - g.lower[i];
    if ((n[0] > 1 && n[1] > 1) || (n[0] > 1 && n[2] > 1) || (n[1] > 1 && n[2] > 1))
      return true;
    return false;
  }

  int cell_count(const Box<int, 3> &box) {
    std::array<int, 3> s = size(box);
    return s[0] * s[1] * s[2];
  }

  void setup_comm_cells(const Subgrid<double, 3> &subgrid, const Box<int, 3> &box, const Box<int, 3> &other) {
    // identify which cells we need to recv
    Box<int, 3> recvbox = intersect(extend(box), other);
    int n = cell_count(recvbox);
    if (check_intersection(recvbox)) {
#ifdef NUFEB_DEBUG_COMM
      debug << "Receiving cells from proc: ";
#endif
      add_cells(subgrid, recvbox, recv_cells);
#ifdef NUFEB_DEBUG_COMM
    debug << std::endl;
#endif	  
    }
    // identify which cells we need to send
    Box<int, 3> sendbox = intersect(extend(other), box);
    n = cell_count(sendbox);
    if (check_intersection(sendbox)) {
#ifdef NUFEB_DEBUG_COMM
      debug << "Sending cells from proc: ";
#endif
      add_cells(subgrid, sendbox, send_cells);
#ifdef NUFEB_DEBUG_COMM
    debug << std::endl;
#endif	  
    }
  }

  void clear() {
    recv_cells.clear();
    send_cells.clear();
    recv_begin.clear();
    recv_end.clear();
    send_begin.clear();
    send_end.clear();
    recv_buff.clear();
    send_buff.clear();
    requests.clear();
  }

  std::vector<int> recv_cells;
  std::vector<int> send_cells;
  std::vector<int> recv_begin;
  std::vector<int> recv_end;
  std::vector<int> send_begin;
  std::vector<int> send_end;
  std::vector<double> recv_buff;
  std::vector<double> send_buff;
  std::vector<MPI_Request> requests;

#ifdef NUFEB_DEBUG_COMM
  std::ofstream debug;
#endif
};
}

#endif // LMP_DECOMP_GRID_H
