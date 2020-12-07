/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
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

  if (!multifile)
    create_one_file();
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
  hid_t file;
  hid_t proplist = H5P_DEFAULT;
  std::string str(filename);
  int offset = 0;
  bool oneperproc = false;

  auto perc = str.find('%');
  if (perc != std::string::npos) {
    oneperproc = true;
    str = std::regex_replace(str, std::regex("%"), std::to_string(comm->me));
  }

  if (!oneperproc) {
    proplist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(proplist, MPI_COMM_WORLD, MPI_INFO_NULL);
  }

  if (multifile) {
    str = std::regex_replace(str, std::regex("\\*"), std::to_string(update->ntimestep));
    file = H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, proplist);
  } else {
    file = H5Fopen(str.c_str(), H5F_ACC_RDWR, proplist);
  }

  H5Pclose(proplist);

  if (!oneperproc)
    MPI_Scan(&atom->nlocal, &offset, 1, MPI_INT, MPI_SUM, world);

  for (auto it = fields.begin(); it != fields.end(); ++it) {
    if (*it == "id") {
      if (multifile) {
	write_atoms_scalar(file, "id", H5T_NATIVE_INT, atom->tag, oneperproc, offset);
      } else {
	hid_t group = H5Gopen(file, "id", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "id/" << std::to_string(update->ntimestep);
	write_atoms_scalar(file, oss.str().c_str(), H5T_NATIVE_INT, atom->tag, oneperproc, offset);
	H5Gclose(group);
      }
    } else if (*it == "type") {
      if (multifile) {
	write_atoms_scalar(file, "type", H5T_NATIVE_INT, atom->type, oneperproc, offset);
      } else {
	hid_t group = H5Gopen(file, "type", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "type/" << std::to_string(update->ntimestep);
	write_atoms_scalar(file, oss.str().c_str(), H5T_NATIVE_INT, atom->type, oneperproc, offset);
	H5Gclose(group);
      }
    } else if (*it == "x") {
      if (multifile) {
	write_atoms_comp(file, "x", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 0);
      } else {
	hid_t group = H5Gopen(file, "x", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "x/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 0);
	H5Gclose(group);
      }
    } else if (*it == "y") {
      if (multifile) {
	write_atoms_comp(file, "y", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 1);
      } else {
	hid_t group = H5Gopen(file, "y", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "y/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 1);
	H5Gclose(group);
      }
    } else if (*it == "z") {
      if (multifile) {
	write_atoms_comp(file, "z", H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 2);
      } else {
	hid_t group = H5Gopen(file, "z", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "z/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->x[0], oneperproc, offset, 2);
	H5Gclose(group);
      }
    } else if (*it == "vx") {
      if (multifile) {
         write_atoms_comp(file, "vx", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 0);
      } else {
	hid_t group = H5Gopen(file, "vx", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "vx/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 0);
	H5Gclose(group);
      }
    } else if (*it == "vy") {
      if (multifile) {
	 write_atoms_comp(file, "vy", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 1);
      } else {
	hid_t group = H5Gopen(file, "vy", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "vy/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 1);
	H5Gclose(group);
      }
    } else if (*it == "vz") {
      if (multifile) {
	 write_atoms_comp(file, "vz", H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 2);
      } else {
	hid_t group = H5Gopen(file, "vz", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "vz/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->v[0], oneperproc, offset, 2);
	H5Gclose(group);
      }
    } else if (*it == "fx") {
      if (multifile) {
	write_atoms_comp(file, "fx", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 0);
      } else {
	hid_t group = H5Gopen(file, "fx", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "fx/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 0);
	H5Gclose(group);
      }
    } else if (*it == "fy") {
      if (multifile) {
	write_atoms_comp(file, "fy", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 1);
      } else {
	hid_t group = H5Gopen(file, "fy", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "fy/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 1);
	H5Gclose(group);
      }
    } else if (*it == "fz") {
      if (multifile) {
	write_atoms_comp(file, "fz", H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 2);
      } else {
	hid_t group = H5Gopen(file, "fz", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "fz/" << std::to_string(update->ntimestep);
	write_atoms_comp(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->f[0], oneperproc, offset, 2);
	H5Gclose(group);
      }
    } else if (*it == "radius") {
      if (multifile) {
	write_atoms_scalar(file, "radius", H5T_NATIVE_DOUBLE, atom->radius, oneperproc, offset);
      } else {
	hid_t group = H5Gopen(file, "radius", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "radius" << std::to_string(update->ntimestep);
	write_atoms_scalar(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, atom->radius, oneperproc, offset);
	H5Gclose(group);
      }
    } else if (*it == "con") {
      hid_t group;
      if (multifile)
	group = H5Gcreate(file, "concentration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	group = H5Gopen(file, "concentration", H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "concentration/" << bio->nuname[i];
	if (multifile) {
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->nus[0], oneperproc, bio->nnu + 1, i);
	} else {
	  hid_t subgroup = H5Gopen(file, oss.str().c_str(), H5P_DEFAULT);
	  oss <<"/" << std::to_string(update->ntimestep);
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->nus[0], oneperproc, bio->nnu + 1, i);
	  H5Gclose(subgroup);
	}
      }
      H5Gclose(group);
    } else if (*it == "upt") {
      hid_t group;
      if (multifile)
	group = H5Gcreate(file, "uptake", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	group = H5Gopen(file, "uptake", H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "uptake/" << bio->nuname[i];
	if (multifile) {
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->nur[0], oneperproc, bio->nnu + 1, i);
	} else {
	  hid_t subgroup = H5Gopen(file, oss.str().c_str(), H5P_DEFAULT);
	  oss <<"/" << std::to_string(update->ntimestep);
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->nur[0], oneperproc, bio->nnu + 1, i);
	  H5Gclose(subgroup);
	}
      }
      H5Gclose(group);
    } else if (*it == "act") {
    } else if (*it == "yie") {
      hid_t group;
      if (multifile)
	group = H5Gcreate(file, "yield", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	group = H5Gopen(file, "yield", H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "yield/" << bio->nuname[i];
	if (multifile) {
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->grid_yield[0], oneperproc, atom->ntypes + 1, i);
	} else {
	  hid_t subgroup = H5Gopen(file, oss.str().c_str(), H5P_DEFAULT);
	  oss <<"/" << std::to_string(update->ntimestep);
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->grid_yield[0], oneperproc, atom->ntypes + 1, i);
	  H5Gclose(subgroup);
	}
      }
      H5Gclose(group);
    } else if (*it == "cat") {
      hid_t group;
      if (multifile)
	group = H5Gcreate(file, "catabolism", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	group = H5Gopen(file, "catabolism", H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "catabolism/" << bio->nuname[i];
	if (multifile) {
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->gibbs_cata[0], oneperproc, atom->ntypes + 1, i);
	} else {
	  hid_t subgroup = H5Gopen(file, oss.str().c_str(), H5P_DEFAULT);
	  oss <<"/" << std::to_string(update->ntimestep);
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->gibbs_cata[0], oneperproc, atom->ntypes + 1, i);
	  H5Gclose(subgroup);
	}
      }
      H5Gclose(group);
    } else if (*it == "ana") {
      hid_t group;
      if (multifile)
	group = H5Gcreate(file, "anabolism", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      else
	group = H5Gopen(file, "anabolism", H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "anabolism/" << bio->nuname[i];
	if (multifile) {
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->gibbs_anab[0], oneperproc, atom->ntypes + 1, i);
	} else {
	  hid_t subgroup = H5Gopen(file, oss.str().c_str(), H5P_DEFAULT);
	  oss <<"/" << std::to_string(update->ntimestep);
	  write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->gibbs_anab[0], oneperproc, atom->ntypes + 1, i);
	  H5Gclose(subgroup);
	}
      }
      H5Gclose(group);
    } else if (*it == "hyd") {
      if (multifile) {
	write_grid(file, "hydronium", H5T_NATIVE_DOUBLE, kinetics->sh, oneperproc);
      } else {
	hid_t group = H5Gopen(file, "hydronium", H5P_DEFAULT);
	std::ostringstream oss;
	oss << "hydronium/" << std::to_string(update->ntimestep);
	write_grid(file, oss.str().c_str(), H5T_NATIVE_DOUBLE, kinetics->sh, oneperproc);
	H5Gclose(group);
      }
    }
  }
  H5Fclose(file);
}

void DumpBioHDF5::create_one_file() {
  std::string str(filename);
  hid_t proplist = H5P_DEFAULT;
  bool oneperproc = false;

  auto perc = str.find('%');
  if (perc != std::string::npos) {
    str = std::regex_replace(str, std::regex("%"), std::to_string(comm->me));
    oneperproc = true;
  }

  if (!oneperproc) {
    proplist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(proplist, MPI_COMM_WORLD, MPI_INFO_NULL);
  }
  hid_t file = H5Fcreate(str.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, proplist);
  H5Pclose(proplist);

  int offset = 0;
  MPI_Scan(&atom->nlocal, &offset, 1, MPI_INT, MPI_SUM, world);
  // create data structure
  for (auto it = fields.begin(); it != fields.end(); ++it) {
    if (*it == "id") {
      hid_t group = H5Gcreate(file, "id", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "type") {
      hid_t group = H5Gcreate(file, "type", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "x") {
      hid_t group = H5Gcreate(file, "x", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "y") {
      hid_t group = H5Gcreate(file, "y", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "z") {
      hid_t group = H5Gcreate(file, "z", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "vx") {
      hid_t group = H5Gcreate(file, "vx", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "vy") {
      hid_t group = H5Gcreate(file, "vy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "vz") {
      hid_t group = H5Gcreate(file, "vz", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "fx") {
      hid_t group = H5Gcreate(file, "fx", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "fy") {
      hid_t group = H5Gcreate(file, "fy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "fz") {
      hid_t group = H5Gcreate(file, "fz", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "radius") {
      hid_t group = H5Gcreate(file, "radius", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
    } else if (*it == "con") {
      hid_t group = H5Gcreate(file, "concentration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "concentration/" << bio->nuname[i];
	hid_t group = H5Gcreate(file, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group);
      }
      H5Gclose(group);
    } else if (*it == "upt") {
      hid_t group = H5Gcreate(file, "uptake", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= bio->nnu; i++) {
	std::ostringstream oss;
	oss << "uptake/" << bio->nuname[i];
	hid_t group = H5Gcreate(file, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group);
      }
      H5Gclose(group);
    } else if (*it == "act") {
    } else if (*it == "yie") {
      hid_t group = H5Gcreate(file, "yield", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= atom->ntypes; i++) {
	std::ostringstream oss;
	oss << "yield/" << bio->nuname[i];
	hid_t group = H5Gcreate(file, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group);
      }
      H5Gclose(group);
    } else if (*it == "cat") {
      hid_t group = H5Gcreate(file, "catabolism", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= atom->ntypes; i++) {
	std::ostringstream oss;
	oss << "catabolism/" << bio->tname[i];
	hid_t group = H5Gcreate(file, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group);
      }
      H5Gclose(group);
    } else if (*it == "ana") {
      hid_t group = H5Gcreate(file, "anabolism", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      for (int i = 1; i <= atom->ntypes; i++) {
	std::ostringstream oss;
	oss << "anabolism/" << bio->tname[i];
	hid_t group = H5Gcreate(file, oss.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group);
      }
      H5Gclose(group);
    } else if (*it == "ph") {
      hid_t group = H5Gcreate(file, "ph", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Gclose(group);
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
  return i;
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
