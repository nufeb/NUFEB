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

#ifdef ENABLE_DUMP_GRID
#include "dump_grid.h"

#include "atom_vec_bio.h"
#include "bio.h"
#include "fix_bio_kinetics.h"
#include "fix_bio_kinetics_energy.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "modify.h"
#include "update.h"

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkXMLImageDataWriter.h>

#include <regex>
#include <sstream>

using namespace LAMMPS_NS;

DumpGrid::DumpGrid(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg), bio(nullptr), kinetics(nullptr), energy(nullptr) {
  AtomVecBio *avec = (AtomVecBio *) atom->style_match("bio");
  if (!avec)
    error->all(FLERR, "Dump grid requires atom style bio");

  bio = avec->bio;

  parse_fields(narg, arg);
}

DumpGrid::~DumpGrid() {}

void DumpGrid::init_style() {
  using std::placeholders::_1;
  for (int j = 0; j < modify->nfix; j++) {
    if (strcmp(modify->fix[j]->style, "kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
    }
    if (strcmp(modify->fix[j]->style, "kinetics/growth/energy") == 0) {
      energy = static_cast<FixKineticsEnergy *>(lmp->modify->fix[j]);
    }
  }

  for (auto it = fields.begin(); it != fields.end(); ++it) {
    if (*it == "con") {
      packs.push_back(std::bind(&DumpGrid::pack_concentration, this, _1));
    }
    else if (*it == "upt") {
      packs.push_back(std::bind(&DumpGrid::pack_uptake, this, _1));
    }
    else if (*it == "act") {
      if (energy)
	packs.push_back(std::bind(&DumpGrid::pack_activity, this, _1));
      else
	error->warning(FLERR, "dump grid 'act' argument only available when using kinetics/energy");
    }
    else if (*it == "yie") {
      if (energy)
	packs.push_back(std::bind(&DumpGrid::pack_yield, this, _1));
      else
	error->warning(FLERR, "dump grid 'yie' argument only available when using kinetics/energy");
    }
    else if (*it == "cat") {
      if (energy)
	packs.push_back(std::bind(&DumpGrid::pack_catabolism, this, _1));
      else
	error->warning(FLERR, "dump grid 'cat' argument only available when using kinetics/energy");
    }
    else if (*it == "ana") {
      if (energy)
	packs.push_back(std::bind(&DumpGrid::pack_anabolism, this, _1));
      else
	error->warning(FLERR, "dump grid 'ana' argument only available when using kinetics/energy");
    }
    else if (*it == "hyd") {
      if (energy)
	packs.push_back(std::bind(&DumpGrid::pack_hydronium, this, _1));
      else
	error->warning(FLERR, "dump grid 'hyd' argument only available when using kinetics/energy");
    }
  }
}

void DumpGrid::write() {
  vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
  std::array<int, 3> dim = kinetics->subgrid.get_dimensions();
  image->SetDimensions(dim[0] + 1, dim[1] + 1, dim[2] + 1);
  std::array<double, 3> spacing = kinetics->grid.get_cell_size();
  image->SetSpacing(spacing[0], spacing[1], spacing[2]);
  std::array<int, 3> origin = kinetics->subgrid.get_origin();
  image->SetOrigin(origin[0] * spacing[0], origin[1] * spacing[1], origin[2] * spacing[2]);

  for (auto it = packs.begin(); it != packs.end(); ++it) {
    (*it)(image);
  }

  // filename must contain '%' and '*'
  std::string str(filename);
  auto perc = str.find('%');
  if (perc == std::string::npos) {
    error->all(FLERR, "dump grid filename must contain '%' special character");
  }
  auto star = str.find('*');
  if (star == std::string::npos) {
    error->all(FLERR, "dump grid filename must contain '*' special character");
  }

  str = std::regex_replace(str, std::regex("%"), std::to_string(comm->me));
  str = std::regex_replace(str, std::regex("\\*"), std::to_string(update->ntimestep));

  vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
  writer->SetFileName(str.c_str());
  writer->SetInputData(image);
  writer->Write();
}

int DumpGrid::parse_fields(int narg, char **arg) {
  int i;
  for (int iarg = 0; iarg < narg; iarg++) {
    i = iarg;
    if (strcmp(arg[iarg],"con") == 0) {
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

void DumpGrid::pack_concentration(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "concentration", kinetics->nus, bio->nuname, bio->nnu);
}

void DumpGrid::pack_uptake(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "uptake", kinetics->nur, bio->nuname, bio->nnu);
}

void DumpGrid::pack_activity(vtkSmartPointer<vtkImageData> image) {
  pack_tuple5(image, "activity", kinetics->activity, bio->nuname, bio->nnu);
}

void DumpGrid::pack_yield(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "yield", kinetics->grid_yield, bio->tname, atom->ntypes);
}

void DumpGrid::pack_catabolism(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "catabolism", kinetics->gibbs_cata, bio->tname, atom->ntypes);
}

void DumpGrid::pack_anabolism(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "anabolism", kinetics->gibbs_anab, bio->tname, atom->ntypes);
}

void DumpGrid::pack_hydronium(vtkSmartPointer<vtkImageData> image) {
  pack_tuple1(image, "hydronium", kinetics->sh);
}

void DumpGrid::pack_tuple1(vtkSmartPointer<vtkImageData> image, const char *name, double *data) {
  vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
  array->SetName(name);
  array->SetNumberOfComponents(1);
  for (int i = 0; i < kinetics->subgrid.cell_count(); i++) {
    array->InsertNextTuple1(data[i]);
  }
  image->GetCellData()->AddArray(array);
}

void DumpGrid::pack_tuple1(vtkSmartPointer<vtkImageData> image, const char *name, double **data, char **names, int count) {
  for (int n = 1; n <= count; n++) {
    vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
    std::ostringstream oss;
    oss << names[n] << " " << name;
    array->SetName(oss.str().c_str());
    array->SetNumberOfComponents(1);
    for (int i = 0; i < kinetics->subgrid.cell_count(); i++) {
      array->InsertNextTuple1(data[n][i]);
    }
    image->GetCellData()->AddArray(array);
  }
}

void DumpGrid::pack_tuple5(vtkSmartPointer<vtkImageData> image, const char *name, double ***data, char **names, int count) {
  for (int n = 1; n <= count; n++) {
    vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
    std::ostringstream oss;
    oss << names[n] << " " << name;
    array->SetName(oss.str().c_str());
    array->SetNumberOfComponents(5);
    for (int i = 0; i < kinetics->subgrid.cell_count(); i++) {
      double tuple[5];
      for (int j = 0; j < 5; j++)
	tuple[j] = data[n][j][i];
      array->InsertNextTuple(tuple);
    }
    image->GetCellData()->AddArray(array);
  }
}
#endif // ENABLE_DUMP_GRID
