/* ----------------------------------------------------------------------
   NUFEB package - A LAMMPS user package for Individual-based Modelling of Microbial Communities
   Contributing authors: Bowen Li & Denis Taniguchi (Newcastle University, UK)
   Email: bowen.li2@newcastle.ac.uk & denis.taniguchi@newcastle.ac.uk

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
------------------------------------------------------------------------- */

#include "read_bio_restart.h"

#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "comm.h"
#include "irregular.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_read_restart.h"
#include "group.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "special.h"
#include "universe.h"
#include "mpiio.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

// same as write_restart.cpp

#define MAGIC_STRING "LammpS RestartT"
#define ENDIAN 0x0001
#define ENDIANSWAP 0x1000
#define VERSION_NUMERIC 0

enum{VERSION,SMALLINT,TAGINT,BIGINT,
     UNITS,NTIMESTEP,DIMENSION,NPROCS,PROCGRID,
     NEWTON_PAIR,NEWTON_BOND,
     XPERIODIC,YPERIODIC,ZPERIODIC,BOUNDARY,
     ATOM_STYLE,NATOMS,NTYPES,
     NBONDS,NBONDTYPES,BOND_PER_ATOM,
     NANGLES,NANGLETYPES,ANGLE_PER_ATOM,
     NDIHEDRALS,NDIHEDRALTYPES,DIHEDRAL_PER_ATOM,
     NIMPROPERS,NIMPROPERTYPES,IMPROPER_PER_ATOM,
     TRICLINIC,BOXLO,BOXHI,XY,XZ,YZ,
     SPECIAL_LJ,SPECIAL_COUL,
     MASS,PAIR,BOND,ANGLE,DIHEDRAL,IMPROPER,
     MULTIPROC,MPIIO,PROCSPERFILE,PERPROC,
     IMAGEINT,BOUNDMIN,TIMESTEP,
     ATOM_ID,ATOM_MAP_STYLE,ATOM_MAP_USER,ATOM_SORTFREQ,ATOM_SORTBIN,
     COMM_MODE,COMM_CUTOFF,COMM_VEL};

#define LB_FACTOR 1.1

/* ---------------------------------------------------------------------- */

ReadBIoRestart::ReadBIoRestart(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void ReadBIoRestart::command(int narg, char **arg)
{
  if (narg != 1 && narg != 2) error->all(FLERR,"Illegal read_restart command");

//  if (domain->box_exist)
//    error->all(FLERR,"Cannot read_restart after simulation box is defined");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  // check for remap option

  int remapflag = 0;
  if (narg == 2) {
    if (strcmp(arg[1],"remap") == 0) remapflag = 1;
    else error->all(FLERR,"Illegal read_restart command");
  }

  // if filename contains "*", search dir for latest restart file

  char *file = new char[strlen(arg[0]) + 16];
  if (strchr(arg[0],'*')) {
    int n;
    if (me == 0) {
      file_search(arg[0],file);
      n = strlen(file) + 1;
    }
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(file,n,MPI_CHAR,0,world);
  } else strcpy(file,arg[0]);

  // check for multiproc files and an MPI-IO filename

  if (strchr(arg[0],'%')) multiproc = 1;
  else multiproc = 0;
  if (strstr(arg[0],".mpiio")) mpiioflag = 1;
  else mpiioflag = 0;

  if (multiproc && mpiioflag)
    error->all(FLERR,
               "Read restart MPI-IO input not allowed with % in filename");

  if (mpiioflag) {
    mpiio = new RestartMPIIO(lmp);
    if (!mpiio->mpiio_exists)
      error->all(FLERR,"Reading from MPI-IO filename when "
                 "MPIIO package is not installed");
  }

  // open single restart file or base file for multiproc case

  if (me == 0) {
    if (screen) fprintf(screen,"Reading restart file ...\n");
    char *hfile;
    if (multiproc) {
      hfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(hfile,"%s%s%s",file,"base",ptr+1);
      *ptr = '%';
    } else hfile = file;
    fp = fopen(hfile,"rb");
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open restart file %s",hfile);
      error->one(FLERR,str);
    }
    if (multiproc) delete [] hfile;
  }

  // read magic string, endian flag, numeric version

  magic_string();
  endian();
  int incompatible = version_numeric();

  // read header info which creates simulation box

  header(incompatible);
  domain->box_exist = 1;

  // problem setup using info from header

  int n;
  if (nprocs == 1) n = static_cast<int> (atom->natoms);
  else n = static_cast<int> (LB_FACTOR * atom->natoms / nprocs);

  atom->allocate_type_arrays();
  atom->deallocate_topology();
  atom->avec->grow(n);
  n = atom->nmax;

  domain->print_box("  ");
  domain->set_initial_box(0);
  domain->set_global_box();
  comm->set_proc_grid();
  domain->set_local_box();

  // read groups, ntype-length arrays, force field, fix info from file
  // nextra = max # of extra quantities stored with each atom

  group->read_restart(fp);
  type_arrays();
  force_fields();

  int nextra = modify->read_restart(fp);
  atom->nextra_store = nextra;
  memory->create(atom->extra,n,nextra,"atom:extra");

  // read file layout info

  file_layout();

  // close header file if in multiproc mode

  if (multiproc && me == 0) {
    fclose(fp);
    fp = NULL;
  }

  // read per-proc info

  AtomVec *avec = atom->avec;

  int maxbuf = 0;
  double *buf = NULL;
  int m,flag;

  // MPI-IO input from single file

  if (mpiioflag) {
    mpiio->openForRead(file);
    memory->create(buf,assignedChunkSize,"read_restart:buf");
    mpiio->read((headerOffset+assignedChunkOffset),assignedChunkSize,buf);
    mpiio->close();

    m = 0;
    while (m < assignedChunkSize) m += avec->unpack_restart(&buf[m]);
  }

  // input of single native file
  // nprocs_file = # of chunks in file
  // proc 0 reads a chunk and bcasts it to other procs
  // each proc unpacks the atoms, saving ones in it's sub-domain
  // if remapflag set, remap the atom to box before checking sub-domain
  // check for atom in sub-domain differs for orthogonal vs triclinic box

  else if (multiproc == 0) {

    int triclinic = domain->triclinic;
    imageint *iptr;
    double *x,lamda[3];
    double *coord,*sublo,*subhi;
    if (triclinic == 0) {
      sublo = domain->sublo;
      subhi = domain->subhi;
    } else {
      sublo = domain->sublo_lamda;
      subhi = domain->subhi_lamda;
    }

    for (int iproc = 0; iproc < nprocs_file; iproc++) {
      if (read_int() != PERPROC)
        error->all(FLERR,"Invalid flag in peratom section of restart file");

      n = read_int();
      if (n > maxbuf) {
        maxbuf = n;
        memory->destroy(buf);
        memory->create(buf,maxbuf,"read_restart:buf");
      }
      read_double_vec(n,buf);

      m = 0;
      while (m < n) {
        x = &buf[m+1];
        if (remapflag) {
          iptr = (imageint *) &buf[m+7];
          domain->remap(x,*iptr);
        }

        if (triclinic) {
          domain->x2lamda(x,lamda);
          coord = lamda;
        } else coord = x;

        if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
            coord[1] >= sublo[1] && coord[1] < subhi[1] &&
            coord[2] >= sublo[2] && coord[2] < subhi[2]) {
          m += avec->unpack_restart(&buf[m]);
        } else m += static_cast<int> (buf[m]);
      }
    }

    if (me == 0) {
      fclose(fp);
      fp = NULL;
    }
  }

  // input of multiple native files with procs <= files
  // # of files = multiproc_file
  // each proc reads a subset of files, striding by nprocs
  // each proc keeps all atoms in all perproc chunks in its files

  else if (nprocs <= multiproc_file) {

    char *procfile = new char[strlen(file) + 16];
    char *ptr = strchr(file,'%');

    for (int iproc = me; iproc < multiproc_file; iproc += nprocs) {
      *ptr = '\0';
      sprintf(procfile,"%s%d%s",file,iproc,ptr+1);
      *ptr = '%';
      fp = fopen(procfile,"rb");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open restart file %s",procfile);
        error->one(FLERR,str);
      }

      fread(&flag,sizeof(int),1,fp);
      if (flag != PROCSPERFILE)
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      int procsperfile;
      fread(&procsperfile,sizeof(int),1,fp);

      for (int i = 0; i < procsperfile; i++) {
        fread(&flag,sizeof(int),1,fp);
        if (flag != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        fread(&n,sizeof(int),1,fp);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        fread(buf,sizeof(double),n,fp);

        m = 0;
        while (m < n) m += avec->unpack_restart(&buf[m]);
      }

      fclose(fp);
      fp = NULL;
    }

    delete [] procfile;
  }

  // input of multiple native files with procs > files
  // # of files = multiproc_file
  // cluster procs based on # of files
  // 1st proc in each cluster reads per-proc chunks from file
  // sends chunks round-robin to other procs in its cluster
  // each proc keeps all atoms in its perproc chunks in file

  else {

    // nclusterprocs = # of procs in my cluster that read from one file
    // filewriter = 1 if this proc reads file, else 0
    // fileproc = ID of proc in my cluster who reads from file
    // clustercomm = MPI communicator within my cluster of procs

    int nfile = multiproc_file;
    int icluster = static_cast<int> ((bigint) me * nfile/nprocs);
    int fileproc = static_cast<int> ((bigint) icluster * nprocs/nfile);
    int fcluster = static_cast<int> ((bigint) fileproc * nfile/nprocs);
    if (fcluster < icluster) fileproc++;
    int fileprocnext =
      static_cast<int> ((bigint) (icluster+1) * nprocs/nfile);
    fcluster = static_cast<int> ((bigint) fileprocnext * nfile/nprocs);
    if (fcluster < icluster+1) fileprocnext++;
    int nclusterprocs = fileprocnext - fileproc;
    int filereader = 0;
    if (me == fileproc) filereader = 1;
    MPI_Comm clustercomm;
    MPI_Comm_split(world,icluster,0,&clustercomm);

    if (filereader) {
      char *procfile = new char[strlen(file) + 16];
      char *ptr = strchr(file,'%');
      *ptr = '\0';
      sprintf(procfile,"%s%d%s",file,icluster,ptr+1);
      *ptr = '%';
      fp = fopen(procfile,"rb");
      if (fp == NULL) {
        char str[128];
        sprintf(str,"Cannot open restart file %s",procfile);
        error->one(FLERR,str);
      }
      delete [] procfile;
    }

    int flag,procsperfile;

    if (filereader) {
      fread(&flag,sizeof(int),1,fp);
      if (flag != PROCSPERFILE)
        error->one(FLERR,"Invalid flag in peratom section of restart file");
      fread(&procsperfile,sizeof(int),1,fp);
    }
    MPI_Bcast(&procsperfile,1,MPI_INT,0,clustercomm);

    int tmp,iproc;
    MPI_Request request;

    for (int i = 0; i < procsperfile; i++) {
      if (filereader) {
        fread(&flag,sizeof(int),1,fp);
        if (flag != PERPROC)
          error->one(FLERR,"Invalid flag in peratom section of restart file");

        fread(&n,sizeof(int),1,fp);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        fread(buf,sizeof(double),n,fp);

        if (i % nclusterprocs) {
          iproc = me + (i % nclusterprocs);
          MPI_Send(&n,1,MPI_INT,iproc,0,world);
          MPI_Recv(&tmp,0,MPI_INT,iproc,0,world,MPI_STATUS_IGNORE);
          MPI_Rsend(buf,n,MPI_DOUBLE,iproc,0,world);
        }

      } else if (i % nclusterprocs == me - fileproc) {
        MPI_Recv(&n,1,MPI_INT,fileproc,0,world,MPI_STATUS_IGNORE);
        if (n > maxbuf) {
          maxbuf = n;
          memory->destroy(buf);
          memory->create(buf,maxbuf,"read_restart:buf");
        }
        MPI_Irecv(buf,n,MPI_DOUBLE,fileproc,0,world,&request);
        MPI_Send(&tmp,0,MPI_INT,fileproc,0,world);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
      }

      if (i % nclusterprocs == me - fileproc) {
        m = 0;
        while (m < n) m += avec->unpack_restart(&buf[m]);
      }
    }

    if (filereader && fp != NULL) {
      fclose(fp);
      fp = NULL;
    }
    MPI_Comm_free(&clustercomm);
  }

  // clean-up memory
  delete [] file;
  memory->destroy(buf);
  modify->nfix_restart_peratom = 0;
  // for multiproc or MPI-IO files:
  // perform irregular comm to migrate atoms to correct procs

  if (multiproc || mpiioflag) {

    // if remapflag set, remap all atoms I read back to box before migrating

    if (remapflag) {
      double **x = atom->x;
      imageint *image = atom->image;
      int nlocal = atom->nlocal;

      for (int i = 0; i < nlocal; i++)
        domain->remap(x[i],image[i]);
    }

    // create a temporary fix to hold and migrate extra atom info
    // necessary b/c irregular will migrate atoms

    if (nextra) {
      char cextra[8],fixextra[8];
      sprintf(cextra,"%d",nextra);
      sprintf(fixextra,"%d",modify->nfix_restart_peratom);
      char **newarg = new char*[5];
      newarg[0] = (char *) "_read_restart";
      newarg[1] = (char *) "all";
      newarg[2] = (char *) "READ_RESTART";
      newarg[3] = cextra;
      newarg[4] = fixextra;
      modify->add_fix(5,newarg);
      delete [] newarg;
    }

    // move atoms to new processors via irregular()
    // turn sorting on in migrate_atoms() to avoid non-reproducible restarts
    // in case read by different proc than wrote restart file
    // first do map_init() since irregular->migrate_atoms() will do map_clear()

    if (atom->map_style) {
      atom->map_init();
      atom->map_set();
    }
    if (domain->triclinic) domain->x2lamda(atom->nlocal);
    Irregular *irregular = new Irregular(lmp);
    irregular->migrate_atoms(1);
    delete irregular;
    if (domain->triclinic) domain->lamda2x(atom->nlocal);

    // put extra atom info held by fix back into atom->extra
    // destroy temporary fix

    if (nextra) {
      memory->destroy(atom->extra);
      memory->create(atom->extra,atom->nmax,nextra,"atom:extra");
      int ifix = modify->find_fix("_read_restart");
      FixReadRestart *fix = (FixReadRestart *) modify->fix[ifix];
      int *count = fix->count;
      double **extra = fix->extra;
      double **atom_extra = atom->extra;
      int nlocal = atom->nlocal;
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < count[i]; j++)
          atom_extra[i][j] = extra[i][j];
      modify->delete_fix("_read_restart");
    }
  }

  // check that all atoms were assigned to procs
  bigint natoms;
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

  if (me == 0) {
    if (screen) fprintf(screen,"  " BIGINT_FORMAT " atoms\n",natoms);
    if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atoms\n",natoms);
  }

  if (natoms != atom->natoms)
    error->all(FLERR,"Did not assign all restart atoms correctly");

  if (me == 0) {
    if (atom->nbonds) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " bonds\n",atom->nbonds);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " bonds\n",atom->nbonds);
    }
    if (atom->nangles) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " angles\n",
                          atom->nangles);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " angles\n",
                           atom->nangles);
    }
    if (atom->ndihedrals) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " dihedrals\n",
                          atom->ndihedrals);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " dihedrals\n",
                           atom->ndihedrals);
    }
    if (atom->nimpropers) {
      if (screen) fprintf(screen,"  " BIGINT_FORMAT " impropers\n",
                          atom->nimpropers);
      if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " impropers\n",
                           atom->nimpropers);
    }
  }

  // check that atom IDs are valid

  atom->tag_check();

  // create global mapping of atoms

  if (atom->map_style) {
    atom->map_init();
    atom->map_set();
  }

  // create special bond lists for molecular systems

  if (atom->molecular == 1) {
    Special special(lmp);
    special.build();
  }
}

/* ----------------------------------------------------------------------
   infile contains a "*"
   search for all files which match the infile pattern
   replace "*" with latest timestep value to create outfile name
   search dir referenced by initial pathname of file
   if infile also contains "%", use "base" when searching directory
   only called by proc 0
------------------------------------------------------------------------- */

void ReadBIoRestart::file_search(char *infile, char *outfile)
{
  char *ptr;

  // separate infile into dir + filename

  char *dirname = new char[strlen(infile) + 1];
  char *filename = new char[strlen(infile) + 1];

  if (strchr(infile,'/')) {
    ptr = strrchr(infile,'/');
    *ptr = '\0';
    strcpy(dirname,infile);
    strcpy(filename,ptr+1);
    *ptr = '/';
  } else {
    strcpy(dirname,"./");
    strcpy(filename,infile);
  }

  // if filename contains "%" replace "%" with "base"

  char *pattern = new char[strlen(filename) + 16];

  if ((ptr = strchr(filename,'%'))) {
    *ptr = '\0';
    sprintf(pattern,"%s%s%s",filename,"base",ptr+1);
    *ptr = '%';
  } else strcpy(pattern,filename);

  // scan all files in directory, searching for files that match pattern
  // maxnum = largest int that matches "*"

  int n = strlen(pattern) + 16;
  char *begin = new char[n];
  char *middle = new char[n];
  char *end = new char[n];

  ptr = strchr(pattern,'*');
  *ptr = '\0';
  strcpy(begin,pattern);
  strcpy(end,ptr+1);
  int nbegin = strlen(begin);
  bigint maxnum = -1;

  struct dirent *ep;
  DIR *dp = opendir(dirname);
  if (dp == NULL)
    error->one(FLERR,"Cannot open dir to search for restart file");
  while ((ep = readdir(dp))) {
    if (strstr(ep->d_name,begin) != ep->d_name) continue;
    if ((ptr = strstr(&ep->d_name[nbegin],end)) == NULL) continue;
    if (strlen(end) == 0) ptr = ep->d_name + strlen(ep->d_name);
    *ptr = '\0';
    if (strlen(&ep->d_name[nbegin]) < n) {
      strcpy(middle,&ep->d_name[nbegin]);
      if (ATOBIGINT(middle) > maxnum) maxnum = ATOBIGINT(middle);
    }
  }
  closedir(dp);
  if (maxnum < 0) error->one(FLERR,"Found no restart file matching pattern");

  // create outfile with maxint substituted for "*"
  // use original infile, not pattern, since need to retain "%" in filename

  ptr = strchr(infile,'*');
  *ptr = '\0';
  sprintf(outfile,"%s" BIGINT_FORMAT "%s",infile,maxnum,ptr+1);
  *ptr = '*';

  // clean up

  delete [] dirname;
  delete [] filename;
  delete [] pattern;
  delete [] begin;
  delete [] middle;
  delete [] end;
}

/* ----------------------------------------------------------------------
   read header of restart file
------------------------------------------------------------------------- */

void ReadBIoRestart::header(int incompatible)
{
  int xperiodic(-1),yperiodic(-1),zperiodic(-1);

  // read flags and fields until flag = -1

  int flag = read_int();
  while (flag >= 0) {

    // check restart file version, warn if different

    if (flag == VERSION) {
      char *version = read_string();
      if (me == 0) {
        if (screen) fprintf(screen,"  restart file = %s, LAMMPS = %s\n",
                            version,universe->version);
      }
      if (incompatible)
        error->all(FLERR,"Restart file incompatible with current version");
      delete [] version;

    // check lmptype.h sizes, error if different

    } else if (flag == SMALLINT) {
      int size = read_int();
      if (size != sizeof(smallint))
        error->all(FLERR,"Smallint setting in lmptype.h is not compatible");
    } else if (flag == IMAGEINT) {
      int size = read_int();
      if (size != sizeof(imageint))
        error->all(FLERR,"Imageint setting in lmptype.h is not compatible");
    } else if (flag == TAGINT) {
      int size = read_int();
      if (size != sizeof(tagint))
        error->all(FLERR,"Tagint setting in lmptype.h is not compatible");
    } else if (flag == BIGINT) {
      int size = read_int();
      if (size != sizeof(bigint))
        error->all(FLERR,"Bigint setting in lmptype.h is not compatible");

    // reset unit_style only if different
    // so that timestep,neighbor-skin are not changed

    } else if (flag == UNITS) {
      char *style = read_string();
      if (strcmp(style,update->unit_style) != 0) update->set_units(style);
      delete [] style;

    } else if (flag == NTIMESTEP) {
      update->ntimestep = read_bigint();

    // set dimension from restart file

    } else if (flag == DIMENSION) {
      int dimension = read_int();
      domain->dimension = dimension;
      if (domain->dimension == 2 && domain->zperiodic == 0)
        error->all(FLERR,
                   "Cannot run 2d simulation with nonperiodic Z dimension");

    // read nprocs from restart file, warn if different

    } else if (flag == NPROCS) {
      nprocs_file = read_int();
      if (nprocs_file != comm->nprocs && me == 0)
        error->warning(FLERR,"Restart file used different # of processors");

    // don't set procgrid, warn if different

    } else if (flag == PROCGRID) {
      int procgrid[3];
      read_int();
      read_int_vec(3,procgrid);
      int flag = 0;
      if (comm->user_procgrid[0] != 0 &&
          procgrid[0] != comm->user_procgrid[0]) flag = 1;
      if (comm->user_procgrid[1] != 0 &&
          procgrid[1] != comm->user_procgrid[1]) flag = 1;
      if (comm->user_procgrid[2] != 0 &&
          procgrid[2] != comm->user_procgrid[2]) flag = 1;
      if (flag && me == 0)
        error->warning(FLERR,"Restart file used different 3d processor grid");

    // don't set newton_pair, leave input script value unchanged
    // set newton_bond from restart file
    // warn if different and input script settings are not default

    } else if (flag == NEWTON_PAIR) {
      int newton_pair_file = read_int();
      if (force->newton_pair != 1) {
        if (newton_pair_file != force->newton_pair && me == 0)
          error->warning(FLERR,
                         "Restart file used different newton pair setting, "
                         "using input script value");
      }
    } else if (flag == NEWTON_BOND) {
      int newton_bond_file = read_int();
      if (force->newton_bond != 1) {
        if (newton_bond_file != force->newton_bond && me == 0)
          error->warning(FLERR,
                         "Restart file used different newton bond setting, "
                         "using restart file value");
      }
      force->newton_bond = newton_bond_file;
      if (force->newton_pair || force->newton_bond) force->newton = 1;
      else force->newton = 0;

    // set boundary settings from restart file
    // warn if different and input script settings are not default

    } else if (flag == XPERIODIC) {
      xperiodic = read_int();
    } else if (flag == YPERIODIC) {
      yperiodic = read_int();
    } else if (flag == ZPERIODIC) {
      zperiodic = read_int();
    } else if (flag == BOUNDARY) {
      int boundary[3][2];
      read_int();
      read_int_vec(6,&boundary[0][0]);

      if (domain->boundary[0][0] || domain->boundary[0][1] ||
          domain->boundary[1][0] || domain->boundary[1][1] ||
          domain->boundary[2][0] || domain->boundary[2][1]) {
        if (boundary[0][0] != domain->boundary[0][0] ||
            boundary[0][1] != domain->boundary[0][1] ||
            boundary[1][0] != domain->boundary[1][0] ||
            boundary[1][1] != domain->boundary[1][1] ||
            boundary[2][0] != domain->boundary[2][0] ||
            boundary[2][1] != domain->boundary[2][1]) {
          if (me == 0)
            error->warning(FLERR,
                           "Restart file used different boundary settings, "
                           "using restart file values");
        }
      }

      domain->boundary[0][0] = boundary[0][0];
      domain->boundary[0][1] = boundary[0][1];
      domain->boundary[1][0] = boundary[1][0];
      domain->boundary[1][1] = boundary[1][1];
      domain->boundary[2][0] = boundary[2][0];
      domain->boundary[2][1] = boundary[2][1];

      if (xperiodic < 0 || yperiodic < 0 || zperiodic < 0)
        error->all(FLERR,"Illegal or unset periodicity in restart");

      domain->periodicity[0] = domain->xperiodic = xperiodic;
      domain->periodicity[1] = domain->yperiodic = yperiodic;
      domain->periodicity[2] = domain->zperiodic = zperiodic;

      domain->nonperiodic = 0;
      if (xperiodic == 0 || yperiodic == 0 || zperiodic == 0) {
        domain->nonperiodic = 1;
        if (boundary[0][0] >= 2 || boundary[0][1] >= 2 ||
            boundary[1][0] >= 2 || boundary[1][1] >= 2 ||
            boundary[2][0] >= 2 || boundary[2][1] >= 2)
          domain->nonperiodic = 2;
      }

    } else if (flag == BOUNDMIN) {
      double minbound[6];
      read_int();
      read_double_vec(6,minbound);
      domain->minxlo = minbound[0]; domain->minxhi = minbound[1];
      domain->minylo = minbound[2]; domain->minyhi = minbound[3];
      domain->minzlo = minbound[4]; domain->minzhi = minbound[5];

    // create new AtomVec class using any stored args

    } else if (flag == ATOM_STYLE) {
      char *style = read_string();
      int nargcopy = read_int();
      char **argcopy = new char*[nargcopy];
      for (int i = 0; i < nargcopy; i++)
        argcopy[i] = read_string();
      //atom->create_avec(style,nargcopy,argcopy,1);
      for (int i = 0; i < nargcopy; i++) delete [] argcopy[i];
      delete [] argcopy;
      delete [] style;

    } else if (flag == NATOMS) {
      atom->natoms = read_bigint();
    } else if (flag == NTYPES) {
      atom->ntypes = read_int();
    } else if (flag == NBONDS) {
      atom->nbonds = read_bigint();
    } else if (flag == NBONDTYPES) {
      atom->nbondtypes = read_int();
    } else if (flag == BOND_PER_ATOM) {
      atom->bond_per_atom = read_int();
    } else if (flag == NANGLES) {
      atom->nangles = read_bigint();
    } else if (flag == NANGLETYPES) {
      atom->nangletypes = read_int();
    } else if (flag == ANGLE_PER_ATOM) {
      atom->angle_per_atom = read_int();
    } else if (flag == NDIHEDRALS) {
      atom->ndihedrals = read_bigint();
    } else if (flag == NDIHEDRALTYPES) {
      atom->ndihedraltypes = read_int();
    } else if (flag == DIHEDRAL_PER_ATOM) {
      atom->dihedral_per_atom = read_int();
    } else if (flag == NIMPROPERS) {
      atom->nimpropers = read_bigint();
    } else if (flag == NIMPROPERTYPES) {
      atom->nimpropertypes = read_int();
    } else if (flag == IMPROPER_PER_ATOM) {
      atom->improper_per_atom = read_int();

    } else if (flag == TRICLINIC) {
      domain->triclinic = read_int();
    } else if (flag == BOXLO) {
      read_int();
      read_double_vec(3,domain->boxlo);
    } else if (flag == BOXHI) {
      read_int();
      read_double_vec(3,domain->boxhi);
    } else if (flag == XY) {
      domain->xy = read_double();
    } else if (flag == XZ) {
      domain->xz = read_double();
    } else if (flag == YZ) {
      domain->yz = read_double();

    } else if (flag == SPECIAL_LJ) {
      read_int();
      read_double_vec(3,&force->special_lj[1]);
    } else if (flag == SPECIAL_COUL) {
      read_int();
      read_double_vec(3,&force->special_coul[1]);

    } else if (flag == TIMESTEP) {
      update->dt = read_double();

    } else if (flag == ATOM_ID) {
      atom->tag_enable = read_int();
    } else if (flag == ATOM_MAP_STYLE) {
      atom->map_style = read_int();
    } else if (flag == ATOM_MAP_USER) {
      atom->map_user  = read_int();
    } else if (flag == ATOM_SORTFREQ) {
      atom->sortfreq = read_int();
    } else if (flag == ATOM_SORTBIN) {
      atom->userbinsize = read_double();

    } else if (flag == COMM_MODE) {
      comm->mode = read_int();
    } else if (flag == COMM_CUTOFF) {
      comm->cutghostuser = read_double();
    } else if (flag == COMM_VEL) {
      comm->ghost_velocity = read_int();

    } else error->all(FLERR,"Invalid flag in header section of restart file");

    flag = read_int();
  }
}

/* ---------------------------------------------------------------------- */

void ReadBIoRestart::type_arrays()
{
  int flag = read_int();
  while (flag >= 0) {

    if (flag == MASS) {
      read_int();
      double *mass = new double[atom->ntypes+1];
      read_double_vec(atom->ntypes,&mass[1]);
      atom->set_mass(mass);
      delete [] mass;

    } else error->all(FLERR,
                      "Invalid flag in type arrays section of restart file");

    flag = read_int();
  }
}

/* ---------------------------------------------------------------------- */

void ReadBIoRestart::force_fields()
{
  char *style;

  int flag = read_int();
  while (flag >= 0) {

    if (flag == PAIR) {
      style = read_string();
      force->create_pair(style,1);
      delete [] style;
      force->pair->read_restart(fp);

    } else if (flag == BOND) {
      style = read_string();
      force->create_bond(style,1);
      delete [] style;
      force->bond->read_restart(fp);

    } else if (flag == ANGLE) {
      style = read_string();
      force->create_angle(style,1);
      delete [] style;
      force->angle->read_restart(fp);

    } else if (flag == DIHEDRAL) {
      style = read_string();
      force->create_dihedral(style,1);
      delete [] style;
      force->dihedral->read_restart(fp);

    } else if (flag == IMPROPER) {
      style = read_string();
      force->create_improper(style,1);
      delete [] style;
      force->improper->read_restart(fp);

    } else error->all(FLERR,
                      "Invalid flag in force field section of restart file");

    flag = read_int();
  }
}

/* ---------------------------------------------------------------------- */

void ReadBIoRestart::file_layout()
{
  int flag = read_int();
  while (flag >= 0) {

    if (flag == MULTIPROC) {
      multiproc_file = read_int();
      if (multiproc == 0 && multiproc_file)
        error->all(FLERR,"Restart file is not a multi-proc file");
      if (multiproc && multiproc_file == 0)
        error->all(FLERR,"Restart file is a multi-proc file");

    } else if (flag == MPIIO) {
      int mpiioflag_file = read_int();
      if (mpiioflag == 0 && mpiioflag_file)
        error->all(FLERR,"Restart file is a MPI-IO file");
      if (mpiioflag && mpiioflag_file == 0)
        error->all(FLERR,"Restart file is not a MPI-IO file");

      if (mpiioflag) {
        bigint *nproc_chunk_offsets;
        memory->create(nproc_chunk_offsets,nprocs,
                       "write_restart:nproc_chunk_offsets");
        bigint *nproc_chunk_sizes;
        memory->create(nproc_chunk_sizes,nprocs,
                       "write_restart:nproc_chunk_sizes");

        // on rank 0 read in the chunk sizes that were written out
        // then consolidate them and compute offsets relative to the
        // end of the header info to fit the current partition size
        // if the number of ranks that did the writing is different

        if (me == 0) {
          int *all_written_send_sizes;
          memory->create(all_written_send_sizes,nprocs_file,
                         "write_restart:all_written_send_sizes");
          int *nproc_chunk_number;
          memory->create(nproc_chunk_number,nprocs,
                         "write_restart:nproc_chunk_number");

          fread(all_written_send_sizes,sizeof(int),nprocs_file,fp);

          int init_chunk_number = nprocs_file/nprocs;
          int num_extra_chunks = nprocs_file - (nprocs*init_chunk_number);

          for (int i = 0; i < nprocs; i++) {
            if (i < num_extra_chunks)
              nproc_chunk_number[i] = init_chunk_number+1;
            else
              nproc_chunk_number[i] = init_chunk_number;
          }

          int all_written_send_sizes_index = 0;
          bigint current_offset = 0;
          for (int i=0;i<nprocs;i++) {
            nproc_chunk_offsets[i] = current_offset;
            nproc_chunk_sizes[i] = 0;
            for (int j=0;j<nproc_chunk_number[i];j++) {
              nproc_chunk_sizes[i] +=
                all_written_send_sizes[all_written_send_sizes_index];
              current_offset +=
                (all_written_send_sizes[all_written_send_sizes_index] *
                 sizeof(double));
              all_written_send_sizes_index++;
            }

          }
          memory->destroy(all_written_send_sizes);
          memory->destroy(nproc_chunk_number);
        }

        // scatter chunk sizes and offsets to all procs

        MPI_Scatter(nproc_chunk_sizes, 1, MPI_LMP_BIGINT,
                    &assignedChunkSize , 1, MPI_LMP_BIGINT, 0,world);
        MPI_Scatter(nproc_chunk_offsets, 1, MPI_LMP_BIGINT,
                    &assignedChunkOffset , 1, MPI_LMP_BIGINT, 0,world);

        memory->destroy(nproc_chunk_sizes);
        memory->destroy(nproc_chunk_offsets);
      }
    }

    flag = read_int();
  }

  // if MPI-IO file, broadcast the end of the header offste
  // this allows all ranks to compute offset to their data

  if (mpiioflag) {
    if (me == 0) headerOffset = ftell(fp);
    MPI_Bcast(&headerOffset,1,MPI_LMP_BIGINT,0,world);
  }
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// low-level fread methods
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void ReadBIoRestart::magic_string()
{
  int n = strlen(MAGIC_STRING) + 1;
  char *str = new char[n];

  int count;
  if (me == 0) count = fread(str,sizeof(char),n,fp);
  MPI_Bcast(&count,1,MPI_INT,0,world);
  if (count < n)
    error->all(FLERR,"Invalid LAMMPS restart file");
  MPI_Bcast(str,n,MPI_CHAR,0,world);
  if (strcmp(str,MAGIC_STRING) != 0)
    error->all(FLERR,"Invalid LAMMPS restart file");
  delete [] str;
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void ReadBIoRestart::endian()
{
  int endian = read_int();
  if (endian == ENDIAN) return;
  if (endian == ENDIANSWAP)
    error->all(FLERR,"Restart file byte ordering is swapped");
  else error->all(FLERR,"Restart file byte ordering is not recognized");
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int ReadBIoRestart::version_numeric()
{
  int vn = read_int();
  if (vn != VERSION_NUMERIC) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   read an int from restart file and bcast it
------------------------------------------------------------------------- */

int ReadBIoRestart::read_int()
{
  int value;
  if ((me == 0) && (fread(&value,sizeof(int),1,fp) < 1))
    value = -1;
  MPI_Bcast(&value,1,MPI_INT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a bigint from restart file and bcast it
------------------------------------------------------------------------- */

bigint ReadBIoRestart::read_bigint()
{
  bigint value;
  if ((me == 0) && (fread(&value,sizeof(bigint),1,fp) < 1))
    value = -1;
  MPI_Bcast(&value,1,MPI_LMP_BIGINT,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a double from restart file and bcast it
------------------------------------------------------------------------- */

double ReadBIoRestart::read_double()
{
  double value;
  if ((me == 0) && (fread(&value,sizeof(double),1,fp) < 1))
    value = 0.0;
  MPI_Bcast(&value,1,MPI_DOUBLE,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read a char string (including NULL) and bcast it
   str is allocated here, ptr is returned, caller must deallocate
------------------------------------------------------------------------- */

char *ReadBIoRestart::read_string()
{
  int n = read_int();
  if (n < 0) error->all(FLERR,"Illegal size string or corrupt restart");
  char *value = new char[n];
  if (me == 0) fread(value,sizeof(char),n,fp);
  MPI_Bcast(value,n,MPI_CHAR,0,world);
  return value;
}

/* ----------------------------------------------------------------------
   read vector of N ints from restart file and bcast them
------------------------------------------------------------------------- */

void ReadBIoRestart::read_int_vec(int n, int *vec)
{
  if (n < 0) error->all(FLERR,"Illegal size integer vector read requested");
  if (me == 0) fread(vec,sizeof(int),n,fp);
  MPI_Bcast(vec,n,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   read vector of N doubles from restart file and bcast them
------------------------------------------------------------------------- */

void ReadBIoRestart::read_double_vec(int n, double *vec)
{
  if (n < 0) error->all(FLERR,"Illegal size double vector read requested");
  if (me == 0) fread(vec,sizeof(double),n,fp);
  MPI_Bcast(vec,n,MPI_DOUBLE,0,world);
}
