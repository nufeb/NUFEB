/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "dump_bio.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <dirent.h>

#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "fix.h"
#include "modify.h"
#include "atom.h"
#include "bio.h"
#include "fix_bio_kinetics.h"

#include "compute_bio_diameter.h"
#include "compute_bio_dimension.h"
#include "compute_bio_diversity.h"
#include "compute_bio_height.h"
#include "compute_bio_rough.h"
#include "compute_bio_segregate.h"
#include "compute_bio_ntypes.h"
#include "compute_bio_biomass.h"
#include "compute_bio_avgcon.h"
#include "compute_bio_avgph.h"
#include "compute_bio_gas.h"

struct stat st = {0};

using namespace LAMMPS_NS;

// customize by adding keyword

/* ---------------------------------------------------------------------- */

DumpBio::DumpBio(LAMMPS *lmp, int narg, char **arg) :
  Dump(lmp, narg, arg)
{
  if (narg <= 4) error->all(FLERR,"No dump bio arguments specified");

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery <= 0) error->all(FLERR,"Illegal dump bio command");

  nkeywords = narg - 4;

  nfix = 0;
  id_fix = NULL;
  fix = NULL;
  keywords = NULL;
  fp = NULL;

  anab_flag = 0;
  nus_flag = 0;
  avgnus_flag = 0;
  cata_flag = 0;
  ph_flag = 0;
  mass_flag = 0;
  mass_header = 0;
  div_header = 0;
  type_header = 0;
  bulk_header = 0;
  avgnus_header = 0;
  gas_flag = 0;
  yield_flag = 0;
  bulk_flag = 0;
  avgph_flag = 0;
  avgph_header = 0;

  dia_flag = 0;
  dim_flag = 0;
  div_flag = 0;
  height_flag = 0;
  rough_flag = 0;
  seg_flag = 0;

  // customize for new sections
  keywords = (char **) memory->srealloc(keywords, nkeywords*sizeof(char *), "keywords");

  for (int i = 0; i < nkeywords; i++) {
    int n = strlen(arg[i+4]) + 2;
    keywords[i] = new char[n];
    strcpy(keywords[i], arg[i+4]);
  }
}

/* ---------------------------------------------------------------------- */

DumpBio::~DumpBio()
{
  for (int i = 0; i < nkeywords; i++) {
    delete [] keywords[i];
  }

  memory->sfree(keywords);
}


/* ---------------------------------------------------------------------- */
void DumpBio::init_style()
{
  // register fix kinetics with this class
  kinetics = NULL;
  nfix = modify->nfix;
  ncompute = modify->ncompute;

  for (int j = 0; j < nfix; j++) {
    if (strcmp(modify->fix[j]->style,"kinetics") == 0) {
      kinetics = static_cast<FixKinetics *>(lmp->modify->fix[j]);
      break;
    }
  }

  if (kinetics == NULL)
    lmp->error->all(FLERR,"The fix kinetics command is required");

  bio = kinetics->bio;

  nx = kinetics->nx;
  ny = kinetics->ny;
  nz = kinetics->nz;

  //Get computational domain size
  if (domain->triclinic == 0) {
    xlo = domain->boxlo[0];
    xhi = domain->boxhi[0];
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
  }
  else {
    xlo = domain->boxlo_bound[0];
    xhi = domain->boxhi_bound[0];
    ylo = domain->boxlo_bound[1];
    yhi = domain->boxhi_bound[1];
    zlo = domain->boxlo_bound[2];
    zhi = domain->boxhi_bound[2];
  }

  stepx = (xhi - xlo) / nx;
  stepy = (yhi - ylo) / ny;
  stepz = (zhi - zlo) / nz;

  int i = 0;
  int nnus = kinetics->bio->nnu;
  int ntypes = atom->ntypes;

  //create directory
  if (stat("./Results", &st) == -1) {
      mkdir("./Results", 0700);
  }

  while (i < nkeywords) {
    if (strcmp(keywords[i],"con") == 0) {
      nus_flag = 1;
      if (stat("./Results/S", &st) == -1) {
          mkdir("./Results/S", 0700);
      }
      for (int j = 1; j < nnus + 1; j++) {
        if (bio->nustate[j] == 0 && strcmp(bio->nuname[j], "h") != 0 && strcmp(bio->nuname[j], "h2o") != 0) {
          char *name = bio->nuname[j];
          int len = 13;
          len += strlen(name);
          char path[len];
          strcpy(path, "./Results/S/");
          strcat(path, name);

          if (stat(path, &st) == -1) {
              mkdir(path, 0700);
          }
        }
      }
    } else if (strcmp(keywords[i],"DGRAn") == 0) {
      anab_flag = 1;
      if (stat("./Results/DGRAn", &st) == -1) {
          mkdir("./Results/DGRAn", 0700);
      }
      for (int j = 1; j < ntypes + 1; j++) {
        char *name = bio->tname[j];
        int len = 17;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/DGRAn/");
        strcat(path, name);

        if (stat(path, &st) == -1) {
            mkdir(path, 0700);
        }
      }
    } else if (strcmp(keywords[i],"DGRCat") == 0) {
      cata_flag = 1;
      if (stat("./Results/DGRCat", &st) == -1) {
          mkdir("./Results/DGRCat", 0700);
      }
      for (int j = 1; j < ntypes + 1; j++) {
        char *name = bio->tname[j];
        int len = 18;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/DGRCat/");
        strcat(path, name);

        if (stat(path, &st) == -1) {
            mkdir(path, 0700);
        }
      }
    } else if (strcmp(keywords[i],"yield") == 0) {
      yield_flag = 1;
      if (stat("./Results/yield", &st) == -1) {
          mkdir("./Results/yield", 0700);
      }
      for (int j = 1; j < ntypes + 1; j++) {
        char *name = bio->tname[j];
        int len = 17;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/yield/");
        strcat(path, name);

        if (stat(path, &st) == -1) {
            mkdir(path, 0700);
        }
      }
    } else if (strcmp(keywords[i],"ph") == 0) {
      ph_flag = 1;
      if (stat("./Results/pH", &st) == -1) {
          mkdir("./Results/pH", 0700);
      }
    } else if (strcmp(keywords[i],"biomass") == 0) {
      mass_flag = 1;
    } else if (strcmp(keywords[i],"bulk") == 0) {
      bulk_flag = 1;
    }  else if (strcmp(keywords[i],"diameter") == 0) {
      dia_flag = 1;
    } else if (strcmp(keywords[i],"dimension") == 0) {
      dim_flag = 1;
    } else if (strcmp(keywords[i],"diversity") == 0) {
      div_flag = 1;
    } else if (strcmp(keywords[i],"ave_height") == 0) {
      height_flag = 1;
    } else if (strcmp(keywords[i],"roughness") == 0) {
      rough_flag = 1;
    } else if (strcmp(keywords[i],"segregation") == 0) {
      seg_flag = 1;
    }  else if (strcmp(keywords[i],"ntypes") == 0) {
      ntypes_flag = 1;
    } else if (strcmp(keywords[i],"avg_con") == 0) {
      avgnus_flag = 1;
    } else if (strcmp(keywords[i],"avg_ph") == 0) {
      avgph_flag = 1;
    } else if (strcmp(keywords[i],"gas") == 0) {
      gas_flag = 1;
    }
    else lmp->error->all(FLERR,"Undefined dump_bio keyword");

    i++;
  }

  for (int j = 0; j < ncompute; j++) {
    if (strcmp(modify->compute[j]->style,"diameter") == 0) {
      cdia = static_cast<ComputeNufebDiameter *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"dimension") == 0) {
      cdim = static_cast<ComputeNufebDimension *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"diversity") == 0) {
      cdiv = static_cast<ComputeNufebDiversity *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"ave_height") == 0) {
      cheight = static_cast<ComputeNufebHeight *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"roughness") == 0) {
      crough = static_cast<ComputeNufebRough *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"segregation") == 0) {
      cseg = static_cast<ComputeNufebSegregate *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"ntypes") == 0) {
      ctype = static_cast<ComputeNufebNtypes *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"biomass") == 0) {
      cmass = static_cast<ComputeNufebBiomass *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"avg_con") == 0) {
      cavgs = static_cast<ComputeNufebAvgcon *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"avg_ph") == 0) {
      cavgph = static_cast<ComputeNufebAvgph *>(lmp->modify->compute[j]);
      continue;
    } else if (strcmp(modify->compute[j]->style,"gas") == 0) {
      cgas = static_cast<ComputeNufebGas *>(lmp->modify->compute[j]);
      continue;
    }
  }

}

/* ---------------------------------------------------------------------- */

void DumpBio::write()
{
  if (update-> ntimestep == 0) return;

  if (ntypes_flag == 1) ctype->compute_vector();
  if (mass_flag == 1) cmass->compute_vector();
  if (dia_flag == 1)  cdia->compute_scalar();
  if (dim_flag == 1)   cdim->compute_scalar();
  if (div_flag == 1)  cdiv->compute_scalar();
  if (height_flag == 1)  cheight->compute_scalar();
  if (rough_flag == 1) crough->compute_scalar();
  if (seg_flag == 1) cseg->compute_scalar();
  if (avgnus_flag == 1) cavgs->compute_vector();
  if (avgph_flag == 1) cavgph->compute_scalar();
  if (gas_flag == 1) cgas->compute_scalar();

  int nnus = kinetics->bio->nnu;
  int ntypes = atom->ntypes;

  if (mass_flag == 1 && comm->me == 0) {
    int len = 35;
    char path[len];
    strcpy(path, "./Results/biomass.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_biomass_data();
    fclose(fp);
  }

  if (dia_flag == 1 && comm->me == 0) {
    int len = 36;
    char path[len];
    strcpy(path, "./Results/diameter.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_diameter_data();
    fclose(fp);
  }

  if (dim_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/dimension.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_dimension_data();
    fclose(fp);
  }

  if (div_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/diversity.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_diversity_data();
    fclose(fp);
  }

  if (height_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/ave_height.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_height_data();
    fclose(fp);
  }

  if (rough_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/roughness.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_rough_data();
    fclose(fp);
  }

  if (seg_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/segregation.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_segregate_data();
    fclose(fp);
  }

  if (ntypes_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/ntypes.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_ntype_data();
    fclose(fp);
  }

  if (bulk_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/bulk.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_bulk_data();
    fclose(fp);
  }

  if (avgnus_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/ave_concentration.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_avgcon_data();
    fclose(fp);
  }

  if (avgph_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/ave_ph.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_avgph_data();
    fclose(fp);
  }

  if (gas_flag == 1 && comm->me == 0) {
    int len = 38;
    char path[len];
    strcpy(path, "./Results/gas.csv");

    filename = path;
    fp = fopen(filename,"a");
    write_gas_data();
    fclose(fp);
  }

  if (nus_flag == 1) {
    for (int i = 1; i < nnus+1; i++) {
      if (bio->nustate[i] == 0 && strcmp(bio->nuname[i], "h") != 0 && strcmp(bio->nuname[i], "h2o") != 0) {
        char *name = bio->nuname[i];
        int len = 30;
        len += strlen(name);
        char path[len];
        strcpy(path, "./Results/S/");
        strcat(path, name);
        strcat(path, "/r*.csv");

        filename = path;
        openfile();
        write_nus_data(i);
        fclose(fp);
      }
    }
}
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_header(bigint n)
{
}

void DumpBio::write_data(int n, double *mybuf)
{
}
/* ---------------------------------------------------------------------- */
void DumpBio::openfile()
{
  //replace '*' with current timestep
  char *filecurrent = filename;
  //printf("%s \n", filename);
  char *filestar = filecurrent;
  filecurrent = new char[strlen(filestar) + 16];
  char *ptr = strchr(filestar,'*');
  *ptr = '\0';
  if (nprocs > 1)
  {
    sprintf(filecurrent,"%s_%d_" BIGINT_FORMAT "%s",
            filestar,me,update->ntimestep,ptr+1);
  }
  else
  {
    sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
            filestar,update->ntimestep,ptr+1);
  }

  *ptr = '*';
  fp = fopen(filecurrent,"a");
  delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpBio::pack(tagint *ids)
{
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_nus_data(int nuID)
{
 // fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    fprintf(fp, "%i,%f,%f,%f,%e\n",i, x, y, z, kinetics->nus[nuID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_gibbs_cata_data(int typeID)
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,%f,%f,%f,%e\n",i, x, y, z, kinetics->gibbs_cata[typeID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_yield_data(int typeID)
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,%f,%f,%f,%e\n",i, x, y, z, kinetics->grid_yield[typeID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_gibbs_anab_data(int typeID)
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,%f,%f,%f,%e\n",i, x, y, z, kinetics->gibbs_anab[typeID][i]);
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_ph_data()
{
  fprintf(fp, ",x,y,z,scalar,1,1,1,0.5\n");

  for(int i = 0; i < kinetics->ngrids; i++){
    int zpos = i/(kinetics->subn[0] * kinetics->subn[1]) + 1;
    int ypos = (i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1])) / kinetics->subn[0] + 1;
    int xpos = i - (zpos - 1) * (kinetics->subn[0] * kinetics->subn[1]) - (ypos - 1) * kinetics->subn[0] + 1;

    double x = kinetics->sublo[0] + xpos * stepx - stepx/2;
    double y = kinetics->sublo[1] + ypos * stepy - stepy/2;
    double z = zpos * stepz - stepz/2;

    //average += kinetics->DRGCat[2][i];

    fprintf(fp, "%i,%f,%f,%f,%e\n",i, x, y, z, -log10(kinetics->sh[i]));
  }
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_avgcon_data()
{
  if (!avgnus_header) {
    for(int i = 1; i < bio->nnu+1; i++){
      fprintf(fp, "%s,", kinetics->bio->nuname[i]);
    }
    avgnus_header = 1;
    fprintf(fp, "\n");
  }

  fprintf(fp, "%i,", update->ntimestep);

  for(int i = 1; i < bio->nnu+1; i++){
    fprintf(fp, "%e,", cavgs->vector[i]);
  }
  fprintf(fp, "\n");
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_avgph_data()
{
  fprintf(fp, "%i, %e \n", update->ntimestep, cavgph->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_gas_data()
{
  fprintf(fp, "%i, %e \n", update->ntimestep, cgas->scalar);
}


/* ---------------------------------------------------------------------- */

void DumpBio::write_bulk_data()
{
  if (!bulk_header) {
    for(int i = 1; i < bio->nnu+1; i++){
      fprintf(fp, "%s,", kinetics->bio->nuname[i]);
    }
    bulk_header = 1;
    fprintf(fp, "\n");
  }

  fprintf(fp, "%i,", update->ntimestep);

  for(int nu = 1; nu < bio->nnu+1; nu++){
    if (bio->nustate[nu] != 0) continue;
    fprintf(fp, "%e,",  kinetics->nubs[nu]);
  }
  fprintf(fp, "\n");
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_biomass_data()
{
  if (!mass_header) {
    for(int i = 0; i < atom->ntypes+1; i++){
      if(i == 0) fprintf(fp, "Total,");
      else fprintf(fp, "%s,", kinetics->bio->tname[i]);
    }
    mass_header = 1;
    fprintf(fp, "\n");
  }

  fprintf(fp, "%i,", update->ntimestep);

  for(int i = 0; i < atom->ntypes+1; i++){
    fprintf(fp, "%e,", cmass->vector[i]);
  }
  fprintf(fp, "\n");
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_diameter_data()
{
  fprintf(fp, "%i, %e \n", update->ntimestep, cdia->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_dimension_data()
{
  fprintf(fp, "%i, %e \n", update->ntimestep, cdim->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_diversity_data()
{
  fprintf(fp, "%i, %e \n", update->ntimestep, cdiv->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_height_data()
{
  fprintf(fp, "%i, %e \n", update->ntimestep, cheight->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_rough_data()
{
  fprintf(fp, "%i, %e \n", update->ntimestep, crough->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_segregate_data()
{
  fprintf(fp, "%i, %e \n", update->ntimestep, cseg->scalar);
}

/* ---------------------------------------------------------------------- */

void DumpBio::write_ntype_data()
{
  if (!type_header) {
    for(int i = 1; i < atom->ntypes+1; i++){
      fprintf(fp, "%s,", kinetics->bio->tname[i]);
    }
    type_header = 1;
    fprintf(fp, "\n");
  }

  fprintf(fp, "%i,", update->ntimestep);

  for(int i = 1; i < atom->ntypes+1; i++){
    fprintf(fp, "%i,", (int)ctype->vector[i]);
  }
  fprintf(fp, "\n");
}


bigint DumpBio::memory_usage() {
  return 0;
}

