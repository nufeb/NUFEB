/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ENABLE_DUMP_BIO_HDF5
#include "dump_bio_hdf5.h"

#include "atom_vec_bio.h"
#include "bio.h"
#include "fix_bio_kinetics.h"
#include "fix_bio_kinetics_energy.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "modify.h"
#include "update.h"

#include <regex>
#include <sstream>

using namespace LAMMPS_NS;

DumpBioHDF5::DumpBioHDF5(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg), bio(nullptr), kinetics(nullptr), energy(nullptr) {
  AtomVecBio *avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Dump bio/hdf5 requires atom style bio");

  bio = avec->bio;

  parse_fields(narg, arg);
}

DumpBioHDF5::~DumpBioHDF5() {}

void DumpBioHDF5::init_style() {
  for (int j = 0; j < modify->nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    }
    if (strcmp(modify->fix[j]->style, "kinetics/growth/energy") == 0) {
      energy = static_cast<FixKineticsEnergy *>(lmp->modify->fix[j]);
    }
  }
}

hid_t DumpBioHDF5::create_filespace_atom(bool oneperproc) {
  hsize_t dims = atom->nlocal;
  if (!oneperproc)
    dims = atom->natoms;
  return H5Screate_simple(1, &dims, NULL); 
}

hid_t DumpBioHDF5::create_filespace_grid(bool oneperproc) {
  hsize_t dims[3];
  if (oneperproc) {
    for (int i = 0; i < 3; i++)
      dims[i] = kinetics->subgrid.get_dimensions()[i];
  } else {
    for (int i = 0; i < 3; i++)
      dims[i] = kinetics->grid.get_dimensions()[i];
  }
  return H5Screate_simple(3, dims, NULL); 
}

void DumpBioHDF5::write() {
  std::string str(filename);
  auto perc = str.find('%');
  bool oneperproc = false; // one file per proc?
  if (perc != std::string::npos) {
    oneperproc = true;
  }
  auto star = str.find('*');
  if (star == std::string::npos) {
    error->warning(FLERR, "dump bio/hdf5 output file will be ovewriten each timestep. Consider using '*' special character");
  }

  str = std::regex_replace(str, std::regex("%"), std::to_string(comm->me));
  str = std::regex_replace(str, std::regex("\\*"), std::to_string(update->ntimestep));

  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    proplist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(proplist, MPI_COMM_WORLD, MPI_INFO_NULL);
  }
  hid_t file = H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, proplist);
  H5Pclose(proplist);

  int offset = 0;
  if (!oneperproc)
    MPI_Scan(&atom->nlocal, &offset, 1, MPI_INT, MPI_SUM, world); 

  for (auto it = fields.begin(); it != fields.end(); ++it) {
    if (*it == "id") {
      write_atoms_scalar(file, "id", H5T_NATIVE_INT, atom->tag, oneperproc, offset);
    } else if (*it == "type") {
      write_atoms_scalar(file, "type", H5T_NATIVE_INT, atom->type, oneperproc, offset);
    } else if (*it == "x") {
      write_atoms_comp(file, "x", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 0);
    } else if (*it == "y") {
      write_atoms_comp(file, "y", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 1);
    } else if (*it == "z") {
      write_atoms_comp(file, "z", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 2);
    } else if (*it == "vx") {
      write_atoms_comp(file, "vx", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 0);
    } else if (*it == "vy") {
      write_atoms_comp(file, "vy", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 1);
    } else if (*it == "vz") {
      write_atoms_comp(file, "vz", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 2);
    } else if (*it == "fx") {
      write_atoms_comp(file, "fx", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 0);
    } else if (*it == "fy") {
      write_atoms_comp(file, "fy", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 1);
    } else if (*it == "fz") {
      write_atoms_comp(file, "fz", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 2);
    } else if (*it == "radius") {
      write_atoms_scalar(file, "radius", H5T_NATIVE_DOUBLE, atom->radius, oneperproc, offset);
    } else if (*it == "con") {
      hid_t group = H5Gcreate(file, "concentration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "concentration/" << bio->nuname[i];
	write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->nus[0], oneperproc, bio->nnu + 1, i);
      }
      H5Gclose(group);
    } else if (*it == "upt") {
      hid_t group = H5Gcreate(file, "uptake", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "uptake/" << bio->nuname[i];
	write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->nur[0], oneperproc, bio->nnu + 1, i);
      }
      H5Gclose(group);
    } else if (*it == "act") {
    } else if (*it == "yie") {
      hid_t group = H5Gcreate(file, "yield", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= atom->ntypes; i++) {
	std::ostringstream oss;
	oss << "yield/" << bio->nuname[i];
	write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->grid_yield[0], oneperproc, atom->ntypes + 1, i);
      }
      H5Gclose(group);
    } else if (*it == "cat") {
      hid_t group = H5Gcreate(file, "catabolism", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= atom->ntypes; i++) {
	std::ostringstream oss;
	oss << "catabolism/" << bio->tname[i];
	write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->gibbs_cata[0], oneperproc, atom->ntypes + 1, i);
      }
      H5Gclose(group);
    } else if (*it == "ana") {
      hid_t group = H5Gcreate(file, "anabolism", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= atom->ntypes; i++) {
	std::ostringstream oss;
	oss << "anabolism/" << bio->tname[i];
	write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->gibbs_anab[0], oneperproc, atom->ntypes + 1, i);
      }
      H5Gclose(group);
    } else if (*it == "hyd") {
      write_grid(file, "hydronium", H5T_NATIVE_DOUBLE, kinetics->sh, oneperproc);
    }
  }
  H5Fclose(file);
}

int DumpBioHDF5::parse_fields(int narg, char **arg) {
  int i;
  for (int iarg = 0; iarg < narg; iarg++) {
    i = iarg;
    if (strcmp(arg[iarg],"id") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "type") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "x") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "y") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "z") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "vx") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "vy") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "vz") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "fx") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "fy") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "fz") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "radius") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "con") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "upt") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "act") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "yie") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "cat") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "ana") == 0) {
      fields.push_back(arg[iarg]);
    } else if (strcmp(arg[iarg], "hyd") == 0) {
      fields.push_back(arg[iarg]);
    }
  }
}

template <class T>
herr_t DumpBioHDF5::write_atoms_scalar(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int offset) {
  hsize_t dims = atom->nlocal;
  if (!oneperproc)
    dims = atom->natoms;
  hid_t filespace = H5Screate_simple(1, &dims, NULL);
  hid_t dataset = H5Dcreate(file, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const hsize_t count = atom->nlocal;
  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    const hsize_t start = offset - atom->nlocal; 
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL);
    proplist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(proplist, H5FD_MPIO_INDEPENDENT);
  }
  hid_t memspace = H5Screate_simple(1, &count, NULL);
  herr_t err = H5Dwrite(dataset, type, memspace, filespace, proplist, buf);
  H5Sclose(memspace);
  H5Pclose(proplist);
  H5Dclose(dataset);
  H5Sclose(filespace);
  return err;
}

template <class T>
herr_t DumpBioHDF5::write_atoms_comp(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int offset, int comp) {
  hsize_t dims = atom->nlocal;
  if (!oneperproc)
    dims = atom->natoms;
  hid_t filespace = H5Screate_simple(1, &dims, NULL);
  hid_t dataset = H5Dcreate(file, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    const hsize_t start = offset - atom->nlocal; 
    const hsize_t count = atom->nlocal;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &start, NULL, &count, NULL);
    proplist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(proplist, H5FD_MPIO_INDEPENDENT);
  }
  hsize_t memcount[2];
  memcount[0] = atom->nlocal;
  memcount[1] = 3;
  hid_t memspace = H5Screate_simple(2, memcount, NULL);
  hsize_t memstart[2];
  memstart[0] = 0;
  memstart[1] = comp;
  memcount[0] = atom->nlocal;
  memcount[1] = 1;
  hsize_t stride[2];
  stride[0] = 1;
  stride[1] = 3;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, stride, memcount, NULL);
  herr_t err = H5Dwrite(dataset, type, memspace, filespace, proplist, buf);
  H5Sclose(memspace);
  H5Pclose(proplist);
  H5Dclose(dataset);
  H5Sclose(filespace);
  return err;
}

template <class T>
herr_t DumpBioHDF5::write_grid(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc) {
  hsize_t dims[3];
  if (oneperproc) {
    for (int i = 0; i < 3; i++)
      dims[2 - i] = kinetics->subgrid.get_dimensions()[i];
  } else {
    for (int i = 0; i < 3; i++)
      dims[2 - i] = kinetics->grid.get_dimensions()[i];
  }
  hid_t filespace = H5Screate_simple(3, dims, NULL); 
  hid_t dataset = H5Dcreate(file, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    hsize_t start[3];
    for (int i = 0; i < 3; i++)
      start[2 - i] = kinetics->subgrid.get_origin()[i];
    hsize_t count[3];
    for (int i = 0; i < 3; i++)
      count[2 - i] = kinetics->subgrid.get_dimensions()[i];
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
    proplist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(proplist, H5FD_MPIO_INDEPENDENT);
  }
  hsize_t memcount = kinetics->ngrids;
  hid_t memspace = H5Screate_simple(1, &memcount, NULL);
  herr_t err = H5Dwrite(dataset, type, memspace, filespace, proplist, buf);
  H5Sclose(memspace);
  H5Pclose(proplist);
  H5Dclose(dataset);
  H5Sclose(filespace);
  return err;
}

template <class T>
herr_t DumpBioHDF5::write_grid(hid_t file, const char *name, hid_t type, T *buf, bool oneperproc, int n, int index) {
  hsize_t dims[3];
  if (oneperproc) {
    for (int i = 0; i < 3; i++)
      dims[2 - i] = kinetics->subgrid.get_dimensions()[i];
  } else {
    for (int i = 0; i < 3; i++)
      dims[2 - i] = kinetics->grid.get_dimensions()[i];
  }
  hid_t filespace = H5Screate_simple(3, dims, NULL); 
  hid_t dataset = H5Dcreate(file, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  hid_t proplist = H5P_DEFAULT;
  if (!oneperproc) {
    hsize_t start[3];
    for (int i = 0; i < 3; i++)
      start[2 - i] = kinetics->subgrid.get_origin()[i];
    hsize_t count[3];
    for (int i = 0; i < 3; i++)
      count[2 - i] = kinetics->subgrid.get_dimensions()[i];
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
    proplist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(proplist, H5FD_MPIO_INDEPENDENT);
  }
  hsize_t memcount[2];
  memcount[0] = n;
  memcount[1] = kinetics->ngrids;
  hid_t memspace = H5Screate_simple(2, memcount, NULL);
  hsize_t memstart[2];
  memstart[0] = index;
  memstart[1] = 0;
  memcount[0] = 1;
  memcount[1] = kinetics->ngrids;
  hsize_t stride[2];
  stride[0] = 1;
  stride[1] = 1;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, memstart, stride, memcount, NULL);
  herr_t err = H5Dwrite(dataset, type, memspace, filespace, proplist, buf);
  H5Sclose(memspace);
  H5Pclose(proplist);
  H5Dclose(dataset);
  H5Sclose(filespace);
  return err;
}

#endif // ENABLE_DUMP_BIO_HDF5
