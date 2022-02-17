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

#include "dump_atom_gz.h"
#include "domain.h"
#include "error.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;

DumpAtomGZ::DumpAtomGZ(LAMMPS *lmp, int narg, char **arg) :
  DumpAtom(lmp, narg, arg)
{
  gzFp = NULL;

  if (!compressed)
    error->all(FLERR,"Dump atom/gz only writes compressed files");
}

/* ---------------------------------------------------------------------- */

DumpAtomGZ::~DumpAtomGZ()
{
  if (gzFp) gzclose(gzFp);
  gzFp = NULL;
  fp = NULL;
}

/* ----------------------------------------------------------------------
   generic opening of a dump file
   ASCII or binary or gzipped
   some derived classes override this function
------------------------------------------------------------------------- */

void DumpAtomGZ::openfile()
{
  // single file, already opened, so just return

  if (singlefile_opened) return;
  if (multifile == 0) singlefile_opened = 1;

  // if one file per timestep, replace '*' with current timestep

  char *filecurrent = filename;
  if (multiproc) filecurrent = multiname;

  if (multifile) {
    char *filestar = filecurrent;
    filecurrent = new char[strlen(filestar) + 16];
    char *ptr = strchr(filestar,'*');
    *ptr = '\0';
    if (padflag == 0)
      sprintf(filecurrent,"%s" BIGINT_FORMAT "%s",
              filestar,update->ntimestep,ptr+1);
    else {
      char bif[8],pad[16];
      strcpy(bif,BIGINT_FORMAT);
      sprintf(pad,"%%s%%0%d%s%%s",padflag,&bif[1]);
      sprintf(filecurrent,pad,filestar,update->ntimestep,ptr+1);
    }
    *ptr = '*';
    if (maxfiles > 0) {
      if (numfiles < maxfiles) {
        nameslist[numfiles] = new char[strlen(filecurrent)+1];
        strcpy(nameslist[numfiles],filecurrent);
        ++numfiles;
      } else {
        remove(nameslist[fileidx]);
        delete[] nameslist[fileidx];
        nameslist[fileidx] = new char[strlen(filecurrent)+1];
        strcpy(nameslist[fileidx],filecurrent);
        fileidx = (fileidx + 1) % maxfiles;
      }
    }
  }

  // each proc with filewriter = 1 opens a file

  if (filewriter) {
    if (append_flag) {
      gzFp = gzopen(filecurrent,"ab9");
    } else {
      gzFp = gzopen(filecurrent,"wb9");
    }

    if (gzFp == NULL) error->one(FLERR,"Cannot open dump file");
  } else gzFp = NULL;

  // delete string with timestep replaced

  if (multifile) delete [] filecurrent;
}

/* ---------------------------------------------------------------------- */

void DumpAtomGZ::write_header(bigint ndump)
{
  if ((multiproc) || (!multiproc && me == 0)) {
    if (unit_flag && !unit_count) {
      ++unit_count;
      gzprintf(gzFp,"ITEM: UNITS\n%s\n",update->unit_style);
    }
    if (time_flag) gzprintf(gzFp,"ITEM: TIME\n%.16g\n",compute_time());

    gzprintf(gzFp,"ITEM: TIMESTEP\n");
    gzprintf(gzFp,BIGINT_FORMAT "\n",update->ntimestep);
    gzprintf(gzFp,"ITEM: NUMBER OF ATOMS\n");
    gzprintf(gzFp,BIGINT_FORMAT "\n",ndump);
    if (domain->triclinic == 0) {
      gzprintf(gzFp,"ITEM: BOX BOUNDS %s\n",boundstr);
      gzprintf(gzFp,"%g %g\n",boxxlo,boxxhi);
      gzprintf(gzFp,"%g %g\n",boxylo,boxyhi);
      gzprintf(gzFp,"%g %g\n",boxzlo,boxzhi);
    } else {
      gzprintf(gzFp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
      gzprintf(gzFp,"%g %g %g\n",boxxlo,boxxhi,boxxy);
      gzprintf(gzFp,"%g %g %g\n",boxylo,boxyhi,boxxz);
      gzprintf(gzFp,"%g %g %g\n",boxzlo,boxzhi,boxyz);
    }
    gzprintf(gzFp,"ITEM: ATOMS %s\n",columns);
  }
}

/* ---------------------------------------------------------------------- */

void DumpAtomGZ::write_data(int n, double *mybuf)
{
  gzwrite(gzFp,mybuf,sizeof(char)*n);
}

/* ---------------------------------------------------------------------- */

void DumpAtomGZ::write()
{
  DumpAtom::write();
  if (filewriter) {
    if (multifile) {
      gzclose(gzFp);
      gzFp = NULL;
    } else {
      if (flush_flag)
        gzflush(gzFp,Z_SYNC_FLUSH);
    }
  }
}

