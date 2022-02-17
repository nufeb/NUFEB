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

#include "write_data.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "universe.h"
#include "comm.h"
#include "output.h"
#include "thermo.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{II,IJ};

/* ---------------------------------------------------------------------- */

WriteData::WriteData(LAMMPS *lmp) : Pointers(lmp)
{
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
}

/* ----------------------------------------------------------------------
   called as write_data command in input script
------------------------------------------------------------------------- */

void WriteData::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Write_data command before simulation box is defined");

  if (narg < 1) error->all(FLERR,"Illegal write_data command");

  // if filename contains a "*", replace with current timestep

  char *ptr;
  int n = strlen(arg[0]) + 16;
  char *file = new char[n];

  if ((ptr = strchr(arg[0],'*'))) {
    *ptr = '\0';
    sprintf(file,"%s" BIGINT_FORMAT "%s",arg[0],update->ntimestep,ptr+1);
  } else strcpy(file,arg[0]);

  // read optional args
  // noinit is a hidden arg, only used by -r command-line switch

  pairflag = II;
  coeffflag = 1;
  fixflag = 1;
  int noinit = 0;

  int iarg = 1;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal write_data command");
      if (strcmp(arg[iarg+1],"ii") == 0) pairflag = II;
      else if (strcmp(arg[iarg+1],"ij") == 0) pairflag = IJ;
      else error->all(FLERR,"Illegal write_data command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"noinit") == 0) {
      noinit = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"nocoeff") == 0) {
      coeffflag = 0;
      iarg++;
    } else if (strcmp(arg[iarg],"nofix") == 0) {
      fixflag = 0;
      iarg++;
    } else error->all(FLERR,"Illegal write_data command");
  }

  // init entire system since comm->exchange is done
  // comm::init needs neighbor::init needs pair::init needs kspace::init, etc
  // exception is when called by -r command-line switch
  //   then write_data immediately follows reading of restart file
  //   assume that read_restart initialized necessary values
  //   if don't make exception:
  //     pair->init() can fail due to various unset values:
  //     e.g. pair hybrid coeffs, dpd ghost-atom velocity setting

  if (noinit == 0) {
    if (comm->me == 0 && screen)
      fprintf(screen,"System init for write_data ...\n");
    lmp->init();

    // move atoms to new processors before writing file
    // do setup_pre_exchange to force update of per-atom info if needed
    // enforce PBC in case atoms are outside box
    // call borders() to rebuild atom map since exchange() destroys map

    modify->setup_pre_exchange();
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    comm->exchange();
    comm->borders();
    if (domain->triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
  }

  write(file);

  delete [] file;
}

/* ----------------------------------------------------------------------
   called from command()
   might later let it be directly called within run/minimize loop
------------------------------------------------------------------------- */

void WriteData::write(char *file)
{
  // special case where reneighboring is not done in integrator
  //   on timestep data file is written (due to build_once being set)
  // if box is changing, must be reset, else data file will have
  //   wrong box size and atoms will be lost when data file is read
  // other calls to pbc and domain and comm are not made,
  //   b/c they only make sense if reneighboring is actually performed

  //if (neighbor->build_once) domain->reset_box();

  // natoms = sum of nlocal = value to write into data file
  // if unequal and thermo lostflag is "error", don't write data file

  bigint nblocal = atom->nlocal;
  bigint natoms;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (natoms != atom->natoms && output->thermo->lostflag == Thermo::ERROR)
    error->all(FLERR,"Atom count is inconsistent, cannot write data file");

  // sum up bond,angle,dihedral,improper counts
  // may be different than atom->nbonds,nangles, etc. if broken/turned-off

  if (atom->molecular == 1 && (atom->nbonds || atom->nbondtypes)) {
    nbonds_local = atom->avec->pack_bond(NULL);
    MPI_Allreduce(&nbonds_local,&nbonds,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }
  if (atom->molecular == 1 && (atom->nangles || atom->nangletypes)) {
    nangles_local = atom->avec->pack_angle(NULL);
    MPI_Allreduce(&nangles_local,&nangles,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  if (atom->molecular == 1 && (atom->ndihedrals || atom->ndihedraltypes)) {
    ndihedrals_local = atom->avec->pack_dihedral(NULL);
    MPI_Allreduce(&ndihedrals_local,&ndihedrals,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  if (atom->molecular == 1 && (atom->nimpropers || atom->nimpropertypes)) {
    nimpropers_local = atom->avec->pack_improper(NULL);
    MPI_Allreduce(&nimpropers_local,&nimpropers,1,MPI_LMP_BIGINT,MPI_SUM,world);
  }

  // check for bonus data.
  if (me == 0) {
    if ((atom->nellipsoids > 0)
        || (atom->nlines > 0)
        || (atom->ntris > 0)
        || (atom->nbodies > 0))
      error->warning(FLERR,"System has ellipsoids, lines, triangles, or bodies. "
                     "Those are not yet supported by write_data. The data file "
                     "will thus be incomplete.");
  }

  // open data file

  if (me == 0) {
    fp = fopen(file,"w");
    if (fp == NULL) {
      char str[128];
      snprintf(str,128,"Cannot open data file %s",file);
      error->one(FLERR,str);
    }
  }

  // proc 0 writes header, ntype-length arrays, force fields

  if (me == 0) {
    header();
    type_arrays();
    if (coeffflag) force_fields();
  }

  // per atom info
  // do not write molecular topology for atom_style template

  if (natoms) atoms();
  if (natoms) velocities();
  if (atom->molecular == 1) {
    if (atom->nbonds && nbonds) bonds();
    if (atom->nangles && nangles) angles();
    if (atom->ndihedrals) dihedrals();
    if (atom->nimpropers) impropers();
  }

  // extra sections managed by fixes
  if (fixflag)
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->wd_section)
        for (int m = 0; m < modify->fix[i]->wd_section; m++) fix(i,m);

  // close data file

  if (me == 0) fclose(fp);
}

/* ----------------------------------------------------------------------
   proc 0 writes out data file header
------------------------------------------------------------------------- */

void WriteData::header()
{
  fprintf(fp,"LAMMPS data file via write_data, version %s, "
          "timestep = " BIGINT_FORMAT "\n",
          universe->version,update->ntimestep);

  fprintf(fp,"\n");

  fprintf(fp,BIGINT_FORMAT " atoms\n",atom->natoms);
  fprintf(fp,"%d atom types\n",atom->ntypes);

  // do not write molecular topology info for atom_style template

  if (atom->molecular == 1) {
    if (atom->nbonds || atom->nbondtypes) {
      fprintf(fp,BIGINT_FORMAT " bonds\n",nbonds);
      fprintf(fp,"%d bond types\n",atom->nbondtypes);
    }
    if (atom->nangles || atom->nangletypes) {
      fprintf(fp,BIGINT_FORMAT " angles\n",nangles);
      fprintf(fp,"%d angle types\n",atom->nangletypes);
    }
    if (atom->ndihedrals || atom->ndihedraltypes) {
      fprintf(fp,BIGINT_FORMAT " dihedrals\n",ndihedrals);
      fprintf(fp,"%d dihedral types\n",atom->ndihedraltypes);
    }
    if (atom->nimpropers || atom->nimpropertypes) {
      fprintf(fp,BIGINT_FORMAT " impropers\n",nimpropers);
      fprintf(fp,"%d improper types\n",atom->nimpropertypes);
    }
  }

  if (fixflag)
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->wd_header)
        for (int m = 0; m < modify->fix[i]->wd_header; m++)
          modify->fix[i]->write_data_header(fp,m);

  fprintf(fp,"\n");

  fprintf(fp,"%-1.16e %-1.16e xlo xhi\n",domain->boxlo[0],domain->boxhi[0]);
  fprintf(fp,"%-1.16e %-1.16e ylo yhi\n",domain->boxlo[1],domain->boxhi[1]);
  fprintf(fp,"%-1.16e %-1.16e zlo zhi\n",domain->boxlo[2],domain->boxhi[2]);

  if (domain->triclinic)
    fprintf(fp,"%-1.16e %-1.16e %-1.16e xy xz yz\n",
            domain->xy,domain->xz,domain->yz);
}

/* ----------------------------------------------------------------------
   proc 0 writes out any type-based arrays that are defined
------------------------------------------------------------------------- */

void WriteData::type_arrays()
{
  if (atom->mass) {
    double *mass = atom->mass;
    fprintf(fp,"\nMasses\n\n");
    for (int i = 1; i <= atom->ntypes; i++) fprintf(fp,"%d %g\n",i,mass[i]);
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes out force field info
------------------------------------------------------------------------- */

void WriteData::force_fields()
{
  if (force->pair && force->pair->writedata) {
    if (pairflag == II) {
      fprintf(fp,"\nPair Coeffs # %s\n\n", force->pair_style);
      force->pair->write_data(fp);
    } else if (pairflag == IJ) {
      fprintf(fp,"\nPairIJ Coeffs # %s\n\n", force->pair_style);
      force->pair->write_data_all(fp);
    }
  }
  if (force->bond && force->bond->writedata && atom->nbondtypes) {
    fprintf(fp,"\nBond Coeffs # %s\n\n", force->bond_style);
    force->bond->write_data(fp);
  }
  if (force->angle && force->angle->writedata && atom->nangletypes) {
    fprintf(fp,"\nAngle Coeffs # %s\n\n", force->angle_style);
    force->angle->write_data(fp);
  }
  if (force->dihedral && force->dihedral->writedata && atom->ndihedraltypes) {
    fprintf(fp,"\nDihedral Coeffs # %s\n\n", force->dihedral_style);
    force->dihedral->write_data(fp);
  }
  if (force->improper && force->improper->writedata && atom->nimpropertypes) {
    fprintf(fp,"\nImproper Coeffs # %s\n\n", force->improper_style);
    force->improper->write_data(fp);
  }
}

/* ----------------------------------------------------------------------
   write out Atoms section of data file
------------------------------------------------------------------------- */

void WriteData::atoms()
{
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = atom->avec->size_data_atom + 3;

  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my atom data into buf

  atom->avec->pack_data(buf);

  // write one chunk of atoms per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp,"\nAtoms # %s\n\n",atom->atom_style);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_data(fp,recvrow,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Velocities section of data file
------------------------------------------------------------------------- */

void WriteData::velocities()
{
  // communication buffer for all my Atom info
  // max_size = largest buffer needed by any proc

  int ncol = atom->avec->size_velocity + 1;

  int sendrow = atom->nlocal;
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my velocity data into buf

  atom->avec->pack_vel(buf);

  // write one chunk of velocities per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp,"\nVelocities\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_vel(fp,recvrow,buf);
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Bonds section of data file
------------------------------------------------------------------------- */

void WriteData::bonds()
{
  // communication buffer for all my Bond info

  int ncol = 3;
  int sendrow = static_cast<int> (nbonds_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my bond data into buf

  atom->avec->pack_bond(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp,"\nBonds\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_bond(fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Angles section of data file
------------------------------------------------------------------------- */

void WriteData::angles()
{
  // communication buffer for all my Angle info

  int ncol = 4;
  int sendrow = static_cast<int> (nangles_local);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my angle data into buf

  atom->avec->pack_angle(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp,"\nAngles\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_angle(fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Dihedrals section of data file
------------------------------------------------------------------------- */

void WriteData::dihedrals()
{
  // communication buffer for all my Dihedral info
  // max_size = largest buffer needed by any proc

  int ncol = 5;

  tagint *tag = atom->tag;
  int *num_dihedral = atom->num_dihedral;
  tagint **dihedral_atom2 = atom->dihedral_atom2;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i,j;
  int sendrow = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      sendrow += num_dihedral[i];
  } else {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_dihedral[i]; j++)
        if (tag[i] == dihedral_atom2[i][j]) sendrow++;
  }

  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my dihedral data into buf

  atom->avec->pack_dihedral(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp,"\nDihedrals\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_dihedral(fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Impropers section of data file
------------------------------------------------------------------------- */

void WriteData::impropers()
{
  // communication buffer for all my Improper info
  // max_size = largest buffer needed by any proc

  int ncol = 5;

  tagint *tag = atom->tag;
  int *num_improper = atom->num_improper;
  tagint **improper_atom2 = atom->improper_atom2;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int i,j;
  int sendrow = 0;
  if (newton_bond) {
    for (i = 0; i < nlocal; i++)
      sendrow += num_improper[i];
  } else {
    for (i = 0; i < nlocal; i++)
      for (j = 0; j < num_improper[i]; j++)
        if (tag[i] == improper_atom2[i][j]) sendrow++;
  }

  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  tagint **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my improper data into buf

  atom->avec->pack_improper(buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    fprintf(fp,"\nImpropers\n\n");
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_LMP_TAGINT,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_LMP_TAGINT,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      atom->avec->write_improper(fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_LMP_TAGINT,0,0,world);
  }

  memory->destroy(buf);
}

/* ----------------------------------------------------------------------
   write out Mth section of data file owned by Fix ifix
------------------------------------------------------------------------- */

void WriteData::fix(int ifix, int mth)
{
  // communication buffer for Fix info

  int sendrow,ncol;
  modify->fix[ifix]->write_data_section_size(mth,sendrow,ncol);
  int maxrow;
  MPI_Allreduce(&sendrow,&maxrow,1,MPI_INT,MPI_MAX,world);

  double **buf;
  if (me == 0) memory->create(buf,MAX(1,maxrow),ncol,"write_data:buf");
  else memory->create(buf,MAX(1,sendrow),ncol,"write_data:buf");

  // pack my fix data into buf

  modify->fix[ifix]->write_data_section_pack(mth,buf);

  // write one chunk of info per proc to file
  // proc 0 pings each proc, receives its chunk, writes to file
  // all other procs wait for ping, send their chunk to proc 0

  int tmp,recvrow;

  int index = 1;
  if (me == 0) {
    MPI_Status status;
    MPI_Request request;

    modify->fix[ifix]->write_data_section_keyword(mth,fp);
    for (int iproc = 0; iproc < nprocs; iproc++) {
      if (iproc) {
        MPI_Irecv(&buf[0][0],maxrow*ncol,MPI_DOUBLE,iproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,iproc,0,world);
        MPI_Wait(&request,&status);
        MPI_Get_count(&status,MPI_DOUBLE,&recvrow);
        recvrow /= ncol;
      } else recvrow = sendrow;

      modify->fix[ifix]->write_data_section(mth,fp,recvrow,buf,index);
      index += recvrow;
    }

  } else {
    MPI_Recv(&tmp,0,MPI_INT,0,0,world,MPI_STATUS_IGNORE);
    MPI_Rsend(&buf[0][0],sendrow*ncol,MPI_DOUBLE,0,0,world);
  }

  memory->destroy(buf);
}
