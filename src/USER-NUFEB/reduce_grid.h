#ifndef LMP_REDUCE_GRID_H
#define LMP_REDUCE_GRID_H

#include "subgrid.h"

#include <vector>

namespace LAMMPS_NS {
template <class Derived>
class ReduceGrid {
public:
  void setup() {
    Derived *derived = static_cast<Derived *>(this);
    // communicate grid extent
    Subgrid<double, 2> subgrid = derived->get_subgrid();
    Box<int, 2> box = subgrid.get_box();
    std::vector<int> boxlo(2 * derived->comm->nprocs);
    MPI_Allgather(&box.lower[0], 2, MPI_INT, boxlo.data(), 2, MPI_INT, derived->world);
    std::vector<int> boxhi(2 * derived->comm->nprocs);
    MPI_Allgather(&box.upper[0], 2, MPI_INT, boxhi.data(), 2, MPI_INT, derived->world);
    std::vector<double> sublo_z(derived->comm->nprocs);
    MPI_Allgather(&derived->domain->sublo[2], 1, MPI_DOUBLE, sublo_z.data(), 1, MPI_DOUBLE, derived->world);

    clear();
    recv_begin.resize(derived->comm->nprocs);
    send_begin.resize(derived->comm->nprocs);
    recv_end.resize(derived->comm->nprocs);
    send_end.resize(derived->comm->nprocs);
    requests.resize(derived->comm->nprocs);
    // look for intersections
    int nrecv = 0;
    int nsend = 0;
    int me = derived->comm->me;
    for (int p = 0; p < derived->comm->nprocs; p++) {
      recv_begin[p] = nrecv;
      send_begin[p] = nsend;
      Box<int, 2> other(&boxlo[2 * p], &boxhi[2 * p]);
      if (p != me) {
        // identify which cells we need to recv if we are the bottom most proc
        Box<int, 2> intersection = intersect(box, other);
        int n = cell_count(intersection);
        if (n > 0 && sublo_z[me] == derived->domain->boxlo[2] && sublo_z[p] != derived->domain->boxlo[2]) {
          add_cells(subgrid, intersection, recv_cells);
          nrecv += n;
        }
        // identify which cells we need to send if we are not the bottom most proc
        if (n > 0 && sublo_z[me] != derived->domain->boxlo[2] && sublo_z[p] == derived->domain->boxlo[2]) {
          add_cells(subgrid, intersection, send_cells);
          nsend += n;
        }
      }
      recv_end[p] = nrecv;
      send_end[p] = nsend;
    }
    recv_buff.resize(derived->get_cell_data_size(nrecv));
    send_buff.resize(derived->get_cell_data_size(nsend));
  }

  void exchange() {
    Derived *derived = static_cast<Derived *>(this);
    int nrequests = 0;
    int recv_offset = 0;
    int send_offset = 0;
    for (int p = 0; p < derived->comm->nprocs; p++) {
      if (p == derived->comm->me)
        continue;
      int nrecv = recv_end[p] - recv_begin[p];
      if (nrecv > 0) {
        int count = derived->get_cell_data_size(nrecv);
        MPI_Irecv(&recv_buff[recv_offset], count, MPI_DOUBLE, p, 0, derived->world, &requests[nrequests++]);
        recv_offset += count;
      }
    }
    for (int i = 0; i < nrequests; i++) {
      MPI_Wait(&requests[i], MPI_STATUS_IGNORE);
      derived->unpack_cells_reduce(recv_cells.begin(), recv_cells.end(), recv_buff.begin(), [](double v0, double v1) {return MAX(v0, v1);});
    }
    // pack data to send buffer
    derived->pack_cells(send_cells.begin(), send_cells.end(), send_buff.begin());
    for (int p = 0; p < derived->comm->nprocs; p++) {
      if (p == derived->comm->me)
        continue;
      int nsend = send_end[p] - send_begin[p];
      if (nsend > 0) {
        int count = derived->get_cell_data_size(nsend);
        MPI_Send(&send_buff[send_offset], count, MPI_DOUBLE, p, 0, derived->world);
        send_offset += count;
      }
    }
  }

private:
  int cell_count(const Box<int, 2> &box) {
    std::array<int, 2> s = size(box);
    return s[0] * s[1];
  }

  void add_cells(const Subgrid<double, 2> &subgrid, const Box<int, 2> &box, std::vector<int> &cells) {
    for (int j = box.lower[1]; j < box.upper[1]; j++) {
      for (int i = box.lower[0]; i < box.upper[0]; i++) {
        cells.push_back(subgrid.get_linear_index( { i, j }));
      }
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

  Grid<double, 3> grid;
  Grid<double, 3> subgrid;
  std::vector<int> recv_cells;
  std::vector<int> send_cells;
  std::vector<int> recv_begin;
  std::vector<int> recv_end;
  std::vector<int> send_begin;
  std::vector<int> send_end;
  std::vector<double> recv_buff;
  std::vector<double> send_buff;
  std::vector<MPI_Request> requests;
};
}

#endif // LMP_REDUCE_GRID_H
