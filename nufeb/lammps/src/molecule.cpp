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

#include "molecule.h"
#include <mpi.h>
#include <cctype>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "math_extra.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 256
#define EPSILON 1.0e-7
#define BIG 1.0e20

#define SINERTIA 0.4            // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

Molecule::Molecule(LAMMPS *lmp, int narg, char **arg, int &index) :
  Pointers(lmp), id(NULL), x(NULL), type(NULL), q(NULL), radius(NULL),
  rmass(NULL), num_bond(NULL), bond_type(NULL), bond_atom(NULL),
  num_angle(NULL), angle_type(NULL), angle_atom1(NULL), angle_atom2(NULL),
  angle_atom3(NULL), num_dihedral(NULL), dihedral_type(NULL), dihedral_atom1(NULL),
  dihedral_atom2(NULL), dihedral_atom3(NULL), dihedral_atom4(NULL), num_improper(NULL),
  improper_type(NULL), improper_atom1(NULL), improper_atom2(NULL),
  improper_atom3(NULL), improper_atom4(NULL), nspecial(NULL), special(NULL),
  shake_flag(NULL), shake_atom(NULL), shake_type(NULL), avec_body(NULL), ibodyparams(NULL),
  dbodyparams(NULL), dx(NULL), dxcom(NULL), dxbody(NULL), quat_external(NULL),
  fp(NULL), count(NULL)
{
  me = comm->me;

  if (index >= narg) error->all(FLERR,"Illegal molecule command");

  if (domain->box_exist == 0)
    error->all(FLERR,"Molecule command before simulation box is defined");

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,"Molecule template ID must be "
                 "alphanumeric or underscore characters");

  // parse args until reach unknown arg (next file)

  toffset = 0;
  boffset = aoffset = doffset = ioffset = 0;
  sizescale = 1.0;

  int ifile = index;
  int iarg = ifile+1;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"offset") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal molecule command");
      toffset = force->inumeric(FLERR,arg[iarg+1]);
      boffset = force->inumeric(FLERR,arg[iarg+2]);
      aoffset = force->inumeric(FLERR,arg[iarg+3]);
      doffset = force->inumeric(FLERR,arg[iarg+4]);
      ioffset = force->inumeric(FLERR,arg[iarg+5]);
      if (toffset < 0 || boffset < 0 || aoffset < 0 ||
          doffset < 0 || ioffset < 0)
        error->all(FLERR,"Illegal molecule command");
      iarg += 6;
    } else if (strcmp(arg[iarg],"toff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      toffset = force->inumeric(FLERR,arg[iarg+1]);
      if (toffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"boff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      boffset = force->inumeric(FLERR,arg[iarg+1]);
      if (boffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"aoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      aoffset = force->inumeric(FLERR,arg[iarg+1]);
      if (aoffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"doff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      doffset = force->inumeric(FLERR,arg[iarg+1]);
      if (doffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"ioff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      ioffset = force->inumeric(FLERR,arg[iarg+1]);
      if (ioffset < 0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal molecule command");
      sizescale = force->numeric(FLERR,arg[iarg+1]);
      if (sizescale <= 0.0) error->all(FLERR,"Illegal molecule command");
      iarg += 2;
    } else break;
  }

  index = iarg;

  // last molecule if have scanned all args

  if (iarg == narg) last = 1;
  else last = 0;

  // initialize all fields to empty

  initialize();

  // scan file for sizes of all fields and allocate them

  if (me == 0) open(arg[ifile]);
  read(0);
  if (me == 0) fclose(fp);
  allocate();

  // read file again to populate all fields

  if (me == 0) open(arg[ifile]);
  read(1);
  if (me == 0) fclose(fp);

  // stats

  if (me == 0) {
    if (screen)
      fprintf(screen,"Read molecule %s:\n"
              "  %d atoms with max type %d\n  %d bonds with max type %d\n"
              "  %d angles with max type %d\n  %d dihedrals with max type %d\n"
              "  %d impropers with max type %d\n",
              id,natoms,ntypes,
              nbonds,nbondtypes,nangles,nangletypes,
              ndihedrals,ndihedraltypes,nimpropers,nimpropertypes);
    if (logfile)
      fprintf(logfile,"Read molecule %s:\n"
              "  %d atoms with max type %d\n  %d bonds with max type %d\n"
              "  %d angles with max type %d\n  %d dihedrals with max type %d\n"
              "  %d impropers with max type %d\n",
              id,natoms,ntypes,
              nbonds,nbondtypes,nangles,nangletypes,
              ndihedrals,ndihedraltypes,nimpropers,nimpropertypes);
  }
}

/* ---------------------------------------------------------------------- */

Molecule::~Molecule()
{
  delete [] id;
  deallocate();
}

/* ----------------------------------------------------------------------
   compute center = geometric center of molecule
   also compute:
     dx = displacement of each atom from center
     molradius = radius of molecule from center
       including finite-size particles or body particles
------------------------------------------------------------------------- */

void Molecule::compute_center()
{
  if (centerflag) return;
  centerflag = 1;

  center[0] = center[1] = center[2] = 0.0;
  for (int i = 0; i < natoms; i++) {
    center[0] += x[i][0];
    center[1] += x[i][1];
    center[2] += x[i][2];
  }
  center[0] /= natoms;
  center[1] /= natoms;
  center[2] /= natoms;

  memory->destroy(dx);
  memory->create(dx,natoms,3,"molecule:dx");

  for (int i = 0; i < natoms; i++) {
    dx[i][0] = x[i][0] - center[0];
    dx[i][1] = x[i][1] - center[1];
    dx[i][2] = x[i][2] - center[2];
  }

  molradius = 0.0;
  for (int i = 0; i < natoms; i++) {
    double rad = MathExtra::len3(dx[i]);
    if (radiusflag) rad += radius[i];
    molradius = MAX(molradius,rad);
  }
}

/* ----------------------------------------------------------------------
   compute masstotal = total mass of molecule
   could have been set by user, otherwise calculate it
------------------------------------------------------------------------- */

void Molecule::compute_mass()
{
  if (massflag) return;
  massflag = 1;

  atom->check_mass(FLERR);

  masstotal = 0.0;
  for (int i = 0; i < natoms; i++) {
    if (rmassflag) masstotal += rmass[i];
    else masstotal += atom->mass[type[i]];
  }
}

/* ----------------------------------------------------------------------
   compute com = center of mass of molecule
   could have been set by user, otherwise calculate it
   works for finite size particles assuming no overlap
   also compute:
     dxcom = displacement of each atom from COM
     comatom = which atom (1-Natom) is nearest the COM
     maxextent = furthest any atom in molecule is from comatom (not COM)
------------------------------------------------------------------------- */

void Molecule::compute_com()
{
  if (!comflag) {
    comflag = 1;

    atom->check_mass(FLERR);

    double onemass;
    com[0] = com[1] = com[2] = 0.0;
    for (int i = 0; i < natoms; i++) {
      if (rmassflag) onemass = rmass[i];
      else onemass = atom->mass[type[i]];
      com[0] += x[i][0]*onemass;
      com[1] += x[i][1]*onemass;
      com[2] += x[i][2]*onemass;
    }
    if (masstotal > 0.0) {
      com[0] /= masstotal;
      com[1] /= masstotal;
      com[2] /= masstotal;
    }
  }

  memory->destroy(dxcom);
  memory->create(dxcom,natoms,3,"molecule:dxcom");

  for (int i = 0; i < natoms; i++) {
    dxcom[i][0] = x[i][0] - com[0];
    dxcom[i][1] = x[i][1] - com[1];
    dxcom[i][2] = x[i][2] - com[2];
  }

  double rsqmin = BIG;
  for (int i = 0; i < natoms; i++) {
    double rsq = MathExtra::lensq3(dxcom[i]);
    if (rsq < rsqmin) {
      comatom = i;
      rsqmin = rsq;
    }
  }

  double rsqmax = 0.0;
  for (int i = 0; i < natoms; i++) {
    double dx = x[comatom][0] - x[i][0];
    double dy = x[comatom][1] - x[i][1];
    double dz = x[comatom][2] - x[i][2];
    double rsq = dx*dx + dy*dy + dz*dz;
    rsqmax = MAX(rsqmax,rsq);
  }

  comatom++;
  maxextent = sqrt(rsqmax);
}

/* ----------------------------------------------------------------------
   compute itensor = 6 moments of inertia of molecule around xyz axes
   could have been set by user, otherwise calculate it
   accounts for finite size spheres, assuming no overlap
   also compute:
     inertia = 3 principal components of inertia
     ex,ey,ez = principal axes in space coords
     quat = quaternion for orientation of molecule
     dxbody = displacement of each atom from COM in body frame
------------------------------------------------------------------------- */

void Molecule::compute_inertia()
{
  if (!inertiaflag) {
    inertiaflag = 1;

    atom->check_mass(FLERR);

    double onemass,dx,dy,dz;
    for (int i = 0; i < 6; i++) itensor[i] = 0.0;
    for (int i = 0; i < natoms; i++) {
      if (rmassflag) onemass = rmass[i];
      else onemass = atom->mass[type[i]];
      dx = dxcom[i][0];
      dy = dxcom[i][1];
      dz = dxcom[i][2];
      itensor[0] += onemass * (dy*dy + dz*dz);
      itensor[1] += onemass * (dx*dx + dz*dz);
      itensor[2] += onemass * (dx*dx + dy*dy);
      itensor[3] -= onemass * dy*dz;
      itensor[4] -= onemass * dx*dz;
      itensor[5] -= onemass * dx*dy;
    }

    if (radiusflag) {
      for (int i = 0; i < natoms; i++) {
        if (rmassflag) onemass = rmass[i];
        else onemass = atom->mass[type[i]];
        itensor[0] += SINERTIA*onemass * radius[i]*radius[i];
        itensor[1] += SINERTIA*onemass * radius[i]*radius[i];
        itensor[2] += SINERTIA*onemass * radius[i]*radius[i];
      }
    }
  }

  // diagonalize inertia tensor for each body via Jacobi rotations
  // inertia = 3 eigenvalues = principal moments of inertia
  // evectors and exzy = 3 evectors = principal axes of rigid body

  double cross[3];
  double tensor[3][3],evectors[3][3];

  tensor[0][0] = itensor[0];
  tensor[1][1] = itensor[1];
  tensor[2][2] = itensor[2];
  tensor[1][2] = tensor[2][1] = itensor[3];
  tensor[0][2] = tensor[2][0] = itensor[4];
  tensor[0][1] = tensor[1][0] = itensor[5];

  if (MathExtra::jacobi(tensor,inertia,evectors))
    error->all(FLERR,"Insufficient Jacobi rotations for rigid molecule");

  ex[0] = evectors[0][0];
  ex[1] = evectors[1][0];
  ex[2] = evectors[2][0];
  ey[0] = evectors[0][1];
  ey[1] = evectors[1][1];
  ey[2] = evectors[2][1];
  ez[0] = evectors[0][2];
  ez[1] = evectors[1][2];
  ez[2] = evectors[2][2];

  // if any principal moment < scaled EPSILON, set to 0.0

  double max;
  max = MAX(inertia[0],inertia[1]);
  max = MAX(max,inertia[2]);

  if (inertia[0] < EPSILON*max) inertia[0] = 0.0;
  if (inertia[1] < EPSILON*max) inertia[1] = 0.0;
  if (inertia[2] < EPSILON*max) inertia[2] = 0.0;

  // enforce 3 evectors as a right-handed coordinate system
  // flip 3rd vector if needed

  MathExtra::cross3(ex,ey,cross);
  if (MathExtra::dot3(cross,ez) < 0.0) MathExtra::negate3(ez);

  // create quaternion

  MathExtra::exyz_to_q(ex,ey,ez,quat);

  // compute displacements in body frame defined by quat

  memory->destroy(dxbody);
  memory->create(dxbody,natoms,3,"molecule:dxbody");

  for (int i = 0; i < natoms; i++)
    MathExtra::transpose_matvec(ex,ey,ez,dxcom[i],dxbody[i]);
}

/* ----------------------------------------------------------------------
   read molecule info from file
   flag = 0, just scan for sizes of fields
   flag = 1, read and store fields
------------------------------------------------------------------------- */

void Molecule::read(int flag)
{
  char line[MAXLINE],keyword[MAXLINE];
  char *eof,*ptr;

  // skip 1st line of file

  if (me == 0) {
    eof = fgets(line,MAXLINE,fp);
    if (eof == NULL) error->one(FLERR,"Unexpected end of molecule file");
  }

  // read header lines
  // skip blank lines or lines that start with "#"
  // stop when read an unrecognized line

  while (1) {

    readline(line);

    // trim anything from '#' onward
    // if line is blank, continue

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    if (strspn(line," \t\n\r") == strlen(line)) continue;

    // search line for header keywords and set corresponding variable

    int nmatch = 0;
    int nwant = 0;
    if (strstr(line,"atoms")) {
      nmatch = sscanf(line,"%d",&natoms);
      nwant = 1;
    } else if (strstr(line,"bonds")) {
      nmatch = sscanf(line,"%d",&nbonds);
      nwant = 1;
    } else if (strstr(line,"angles")) {
      nmatch = sscanf(line,"%d",&nangles);
      nwant = 1;
    } else if (strstr(line,"dihedrals")) {
      nmatch = sscanf(line,"%d",&ndihedrals);
      nwant = 1;
    } else if (strstr(line,"impropers")) {
      nmatch = sscanf(line,"%d",&nimpropers);
      nwant = 1;
    } else if (strstr(line,"mass")) {
      massflag = 1;
      nmatch = sscanf(line,"%lg",&masstotal);
      nwant = 1;
      masstotal *= sizescale*sizescale*sizescale;
    } else if (strstr(line,"com")) {
      comflag = 1;
      nmatch = sscanf(line,"%lg %lg %lg",&com[0],&com[1],&com[2]);
      nwant = 3;
      com[0] *= sizescale;
      com[1] *= sizescale;
      com[2] *= sizescale;
      if (domain->dimension == 2 && com[2] != 0.0)
        error->all(FLERR,"Molecule file z center-of-mass must be 0.0 for 2d");
    } else if (strstr(line,"inertia")) {
      inertiaflag = 1;
      nmatch = sscanf(line,"%lg %lg %lg %lg %lg %lg",
                      &itensor[0],&itensor[1],&itensor[2],
                      &itensor[3],&itensor[4],&itensor[5]);
      nwant = 6;
      const double scale5 = sizescale*sizescale*sizescale*sizescale*sizescale;
      itensor[0] *= scale5;
      itensor[1] *= scale5;
      itensor[2] *= scale5;
      itensor[3] *= scale5;
      itensor[4] *= scale5;
      itensor[5] *= scale5;
    } else if (strstr(line,"body")) {
      bodyflag = 1;
      avec_body = (AtomVecBody *) atom->style_match("body");
      if (!avec_body)
        error->all(FLERR,"Molecule file requires atom style body");
      nmatch = sscanf(line,"%d %d",&nibody,&ndbody);
      nwant = 2;
    } else break;

    if (nmatch != nwant)
      error->all(FLERR,"Invalid header in molecule file");
  }

  // error checks

  if (natoms < 1)
    error->all(FLERR,"No count or invalid atom count in molecule file");
  if (nbonds < 0) error->all(FLERR,"Invalid bond count in molecule file");
  if (nangles < 0) error->all(FLERR,"Invalid angle count in molecule file");
  if (ndihedrals < 0)
    error->all(FLERR,"Invalid dihedral count in molecule file");
  if (nimpropers < 0)
    error->all(FLERR,"Invalid improper count in molecule file");

  // count = vector for tallying bonds,angles,etc per atom

  if (flag == 0) memory->create(count,natoms,"molecule:count");
  else count = NULL;

  // grab keyword and skip next line

  parse_keyword(0,line,keyword);
  readline(line);

  // loop over sections of molecule file

  while (strlen(keyword) > 0) {
    if (strcmp(keyword,"Coords") == 0) {
      xflag = 1;
      if (flag) coords(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Types") == 0) {
      typeflag = 1;
      if (flag) types(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Charges") == 0) {
      qflag = 1;
      if (flag) charges(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Diameters") == 0) {
      radiusflag = 1;
      if (flag) diameters(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Masses") == 0) {
      rmassflag = 1;
      if (flag) masses(line);
      else skip_lines(natoms,line);

    } else if (strcmp(keyword,"Bonds") == 0) {
      if (nbonds == 0)
        error->all(FLERR,"Molecule file has bonds but no nbonds setting");
      bondflag = tag_require = 1;
      bonds(flag,line);
    } else if (strcmp(keyword,"Angles") == 0) {
      if (nangles == 0)
        error->all(FLERR,"Molecule file has angles but no nangles setting");
      angleflag = tag_require = 1;
      angles(flag,line);
    } else if (strcmp(keyword,"Dihedrals") == 0) {
      if (ndihedrals == 0) error->all(FLERR,"Molecule file has dihedrals "
                                      "but no ndihedrals setting");
      dihedralflag = tag_require = 1;
      dihedrals(flag,line);
    } else if (strcmp(keyword,"Impropers") == 0) {
      if (nimpropers == 0) error->all(FLERR,"Molecule file has impropers "
                                      "but no nimpropers setting");
      improperflag = tag_require = 1;
      impropers(flag,line);

    } else if (strcmp(keyword,"Special Bond Counts") == 0) {
      nspecialflag = 1;
      nspecial_read(flag,line);
    } else if (strcmp(keyword,"Special Bonds") == 0) {
      specialflag = tag_require = 1;
      if (flag) special_read(line);
      else skip_lines(natoms,line);

    } else if (strcmp(keyword,"Shake Flags") == 0) {
      shakeflagflag = 1;
      if (flag) shakeflag_read(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Shake Atoms") == 0) {
      shakeatomflag = tag_require = 1;
      if (shaketypeflag) shakeflag = 1;
      if (!shakeflagflag)
        error->all(FLERR,"Molecule file shake flags not before shake atoms");
      if (flag) shakeatom_read(line);
      else skip_lines(natoms,line);
    } else if (strcmp(keyword,"Shake Bond Types") == 0) {
      shaketypeflag = 1;
      if (shakeatomflag) shakeflag = 1;
      if (!shakeflagflag)
        error->all(FLERR,"Molecule file shake flags not before shake bonds");
      if (flag) shaketype_read(line);
      else skip_lines(natoms,line);

    } else if (strcmp(keyword,"Body Integers") == 0) {
      if (bodyflag == 0 || nibody == 0)
        error->all(FLERR,"Molecule file has body params "
                   "but no setting for them");
      ibodyflag = 1;
      body(flag,0,line);
    } else if (strcmp(keyword,"Body Doubles") == 0) {
      if (bodyflag == 0 || ndbody == 0)
        error->all(FLERR,"Molecule file has body params "
                   "but no setting for them");
      dbodyflag = 1;
      body(flag,1,line);

    } else error->one(FLERR,"Unknown section in molecule file");

    parse_keyword(1,line,keyword);
  }

  // clean up

  memory->destroy(count);

  // error check

  if (flag == 0) {
    if ((nspecialflag && !specialflag) || (!nspecialflag && specialflag))
      error->all(FLERR,"Molecule file needs both Special Bond sections");
    if (specialflag && !bondflag)
      error->all(FLERR,"Molecule file has special flags but no bonds");
    if ((shakeflagflag || shakeatomflag || shaketypeflag) && !shakeflag)
      error->all(FLERR,"Molecule file shake info is incomplete");
    if (bodyflag && nibody && ibodyflag == 0)
      error->all(FLERR,"Molecule file has no Body Integers section");
    if (bodyflag && ndbody && dbodyflag == 0)
      error->all(FLERR,"Molecule file has no Body Doubles section");
  }

  // auto-generate special bonds if needed and not in file

  if (bondflag && specialflag == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR,"Cannot auto-generate special bonds before "
                       "simulation box is defined");

    if (flag) {
      special_generate();
      specialflag = 1;
      nspecialflag = 1;
    }
  }

  // body particle must have natom = 1
  // set radius by having body class compute its own radius

  if (bodyflag) {
    radiusflag = 1;
    if (natoms != 1)
      error->all(FLERR,"Molecule natoms must be 1 for body particle");
    if (sizescale != 1.0)
      error->all(FLERR,"Molecule sizescale must be 1.0 for body particle");
    if (flag) {
      radius[0] = avec_body->radius_body(nibody,ndbody,ibodyparams,dbodyparams);
      maxradius = radius[0];
    }
  }
}

/* ----------------------------------------------------------------------
   read coords from file
------------------------------------------------------------------------- */

void Molecule::coords(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (4 != sscanf(line,"%d %lg %lg %lg",&tmp,&x[i][0],&x[i][1],&x[i][2]))
      error->all(FLERR,"Invalid Coords section in molecule file");

    x[i][0] *= sizescale;
    x[i][1] *= sizescale;
    x[i][2] *= sizescale;
  }

  if (domain->dimension == 2) {
    for (int i = 0; i < natoms; i++)
      if (x[i][2] != 0.0)
        error->all(FLERR,"Molecule file z coord must be 0.0 for 2d");
  }
}

/* ----------------------------------------------------------------------
   read types from file
   set ntypes = max of any atom type
------------------------------------------------------------------------- */

void Molecule::types(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (2 != sscanf(line,"%d %d",&tmp,&type[i]))
      error->all(FLERR,"Invalid Types section in molecule file");
    type[i] += toffset;
  }

  for (int i = 0; i < natoms; i++)
    if (type[i] <= 0 || type[i] > atom->ntypes)
      error->all(FLERR,"Invalid atom type in molecule file");

  for (int i = 0; i < natoms; i++)
    ntypes = MAX(ntypes,type[i]);
}

/* ----------------------------------------------------------------------
   read charges from file
------------------------------------------------------------------------- */

void Molecule::charges(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (2 != sscanf(line,"%d %lg",&tmp,&q[i]))
      error->all(FLERR,"Invalid Charges section in molecule file");
  }
}

/* ----------------------------------------------------------------------
   read diameters from file and set radii
------------------------------------------------------------------------- */

void Molecule::diameters(char *line)
{
  int tmp;
  maxradius = 0.0;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (2 != sscanf(line,"%d %lg",&tmp,&radius[i]))
      error->all(FLERR,"Invalid Diameters section in molecule file");
    radius[i] *= sizescale;
    radius[i] *= 0.5;
    maxradius = MAX(maxradius,radius[i]);
  }

  for (int i = 0; i < natoms; i++)
    if (radius[i] < 0.0)
      error->all(FLERR,"Invalid atom diameter in molecule file");
}

/* ----------------------------------------------------------------------
   read masses from file
------------------------------------------------------------------------- */

void Molecule::masses(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (2 != sscanf(line,"%d %lg",&tmp,&rmass[i]))
      error->all(FLERR,"Invalid Masses section in molecule file");
    rmass[i] *= sizescale*sizescale*sizescale;
  }

  for (int i = 0; i < natoms; i++)
    if (rmass[i] <= 0.0) error->all(FLERR,"Invalid atom mass in molecule file");
}

/* ----------------------------------------------------------------------
   read bonds from file
   set nbondtypes = max type of any bond
   store each with both atoms if newton_bond = 0
   if flag = 0, just count bonds/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::bonds(int flag, char *line)
{
  int tmp,itype;
  tagint m,atom1,atom2;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_bond[i] = 0;

  for (int i = 0; i < nbonds; i++) {
    readline(line);
    if (4 != sscanf(line,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT,
           &tmp,&itype,&atom1,&atom2))
      error->all(FLERR,"Invalid Bonds section in molecule file");
    itype += boffset;

    if ((atom1 <= 0) || (atom1 > natoms) ||
        (atom2 <= 0) || (atom2 > natoms) || (atom1 == atom2))
      error->one(FLERR,"Invalid atom ID in Bonds section of molecule file");
    if (itype <= 0 || itype > atom->nbondtypes)
      error->one(FLERR,"Invalid bond type in Bonds section of molecule file");

    if (flag) {
      m = atom1-1;
      nbondtypes = MAX(nbondtypes,itype);
      bond_type[m][num_bond[m]] = itype;
      bond_atom[m][num_bond[m]] = atom2;
      num_bond[m]++;
      if (newton_bond == 0) {
        m = atom2-1;
        bond_type[m][num_bond[m]] = itype;
        bond_atom[m][num_bond[m]] = atom1;
        num_bond[m]++;
      }
    } else {
      count[atom1-1]++;
      if (newton_bond == 0) count[atom2-1]++;
    }
  }

  // bond_per_atom = max of count vector

  if (flag == 0) {
    bond_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      bond_per_atom = MAX(bond_per_atom,count[i]);
  }
}

/* ----------------------------------------------------------------------
   read angles from file
   store each with all 3 atoms if newton_bond = 0
   if flag = 0, just count angles/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::angles(int flag, char *line)
{
  int tmp,itype;
  tagint m,atom1,atom2,atom3;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_angle[i] = 0;

  for (int i = 0; i < nangles; i++) {
    readline(line);
    if (5 != sscanf(line,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT,
           &tmp,&itype,&atom1,&atom2,&atom3))
      error->all(FLERR,"Invalid Angles section in molecule file");
    itype += aoffset;

    if ((atom1 <= 0) || (atom1 > natoms) ||
        (atom2 <= 0) || (atom2 > natoms) ||
        (atom3 <= 0) || (atom3 > natoms) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom2 == atom3))
      error->one(FLERR,"Invalid atom ID in Angles section of molecule file");
    if (itype <= 0 || itype > atom->nangletypes)
      error->one(FLERR,"Invalid angle type in Angles section of molecule file");

    if (flag) {
      m = atom2-1;
      nangletypes = MAX(nangletypes,itype);
      angle_type[m][num_angle[m]] = itype;
      angle_atom1[m][num_angle[m]] = atom1;
      angle_atom2[m][num_angle[m]] = atom2;
      angle_atom3[m][num_angle[m]] = atom3;
      num_angle[m]++;
      if (newton_bond == 0) {
        m = atom1-1;
        angle_type[m][num_angle[m]] = itype;
        angle_atom1[m][num_angle[m]] = atom1;
        angle_atom2[m][num_angle[m]] = atom2;
        angle_atom3[m][num_angle[m]] = atom3;
        num_angle[m]++;
        m = atom3-1;
        angle_type[m][num_angle[m]] = itype;
        angle_atom1[m][num_angle[m]] = atom1;
        angle_atom2[m][num_angle[m]] = atom2;
        angle_atom3[m][num_angle[m]] = atom3;
        num_angle[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
        count[atom1-1]++;
        count[atom3-1]++;
      }
    }
  }

  // angle_per_atom = max of count vector

  if (flag == 0) {
    angle_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      angle_per_atom = MAX(angle_per_atom,count[i]);
  }
}

/* ----------------------------------------------------------------------
   read dihedrals from file
   store each with all 4 atoms if newton_bond = 0
   if flag = 0, just count dihedrals/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::dihedrals(int flag, char *line)
{
  int tmp,itype;
  tagint m,atom1,atom2,atom3,atom4;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;

  for (int i = 0; i < ndihedrals; i++) {
    readline(line);
    if (6 != sscanf(line,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " "
           TAGINT_FORMAT " " TAGINT_FORMAT " ",
           &tmp,&itype,&atom1,&atom2,&atom3,&atom4))
      error->all(FLERR,"Invalid Dihedrals section in molecule file");
    itype += doffset;

    if ((atom1 <= 0) || (atom1 > natoms) ||
        (atom2 <= 0) || (atom2 > natoms) ||
        (atom3 <= 0) || (atom3 > natoms) ||
        (atom4 <= 0) || (atom4 > natoms) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) ||
        (atom2 == atom3) || (atom2 == atom4) || (atom3 == atom4))
      error->one(FLERR,
                 "Invalid atom ID in dihedrals section of molecule file");
    if (itype <= 0 || itype > atom->ndihedraltypes)
      error->one(FLERR,
                 "Invalid dihedral type in dihedrals section of molecule file");

    if (flag) {
      m = atom2-1;
      ndihedraltypes = MAX(ndihedraltypes,itype);
      dihedral_type[m][num_dihedral[m]] = itype;
      dihedral_atom1[m][num_dihedral[m]] = atom1;
      dihedral_atom2[m][num_dihedral[m]] = atom2;
      dihedral_atom3[m][num_dihedral[m]] = atom3;
      dihedral_atom4[m][num_dihedral[m]] = atom4;
      num_dihedral[m]++;
      if (newton_bond == 0) {
        m = atom1-1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
        m = atom3-1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
        m = atom4-1;
        dihedral_type[m][num_dihedral[m]] = itype;
        dihedral_atom1[m][num_dihedral[m]] = atom1;
        dihedral_atom2[m][num_dihedral[m]] = atom2;
        dihedral_atom3[m][num_dihedral[m]] = atom3;
        dihedral_atom4[m][num_dihedral[m]] = atom4;
        num_dihedral[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
        count[atom1-1]++;
        count[atom3-1]++;
        count[atom4-1]++;
      }
    }
  }

  // dihedral_per_atom = max of count vector

  if (flag == 0) {
    dihedral_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      dihedral_per_atom = MAX(dihedral_per_atom,count[i]);
  }
}

/* ----------------------------------------------------------------------
   read impropers from file
   store each with all 4 atoms if newton_bond = 0
   if flag = 0, just count impropers/atom
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::impropers(int flag, char *line)
{
  int tmp,itype;
  tagint m,atom1,atom2,atom3,atom4;
  int newton_bond = force->newton_bond;

  if (flag == 0)
    for (int i = 0; i < natoms; i++) count[i] = 0;
  else
    for (int i = 0; i < natoms; i++) num_improper[i] = 0;

  for (int i = 0; i < nimpropers; i++) {
    readline(line);
    if (6 != sscanf(line,"%d %d " TAGINT_FORMAT " " TAGINT_FORMAT " "
           TAGINT_FORMAT " " TAGINT_FORMAT " ",
           &tmp,&itype,&atom1,&atom2,&atom3,&atom4))
      error->all(FLERR,"Invalid Impropers section in molecule file");
    itype += ioffset;

    if ((atom1 <= 0) || (atom1 > natoms) ||
        (atom2 <= 0) || (atom2 > natoms) ||
        (atom3 <= 0) || (atom3 > natoms) ||
        (atom4 <= 0) || (atom4 > natoms) ||
        (atom1 == atom2) || (atom1 == atom3) || (atom1 == atom4) ||
        (atom2 == atom3) || (atom2 == atom4) || (atom3 == atom4))
      error->one(FLERR,
                 "Invalid atom ID in impropers section of molecule file");
    if (itype <= 0 || itype > atom->nimpropertypes)
      error->one(FLERR,
                 "Invalid improper type in impropers section of molecule file");

    if (flag) {
      m = atom2-1;
      nimpropertypes = MAX(nimpropertypes,itype);
      improper_type[m][num_improper[m]] = itype;
      improper_atom1[m][num_improper[m]] = atom1;
      improper_atom2[m][num_improper[m]] = atom2;
      improper_atom3[m][num_improper[m]] = atom3;
      improper_atom4[m][num_improper[m]] = atom4;
      num_improper[m]++;
      if (newton_bond == 0) {
        m = atom1-1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
        m = atom3-1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
        m = atom4-1;
        improper_type[m][num_improper[m]] = itype;
        improper_atom1[m][num_improper[m]] = atom1;
        improper_atom2[m][num_improper[m]] = atom2;
        improper_atom3[m][num_improper[m]] = atom3;
        improper_atom4[m][num_improper[m]] = atom4;
        num_improper[m]++;
      }
    } else {
      count[atom2-1]++;
      if (newton_bond == 0) {
        count[atom1-1]++;
        count[atom3-1]++;
        count[atom4-1]++;
      }
    }
  }

  // improper_per_atom = max of count vector

  if (flag == 0) {
    improper_per_atom = 0;
    for (int i = 0; i < natoms; i++)
      improper_per_atom = MAX(improper_per_atom,count[i]);
  }
}

/* ----------------------------------------------------------------------
   read 3 special bonds counts from file
   if flag = 0, just tally maxspecial
   if flag = 1, store them with atoms
------------------------------------------------------------------------- */

void Molecule::nspecial_read(int flag, char *line)
{
  int tmp,c1,c2,c3;

  if (flag == 0) maxspecial = 0;

  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (4 != sscanf(line,"%d %d %d %d",&tmp,&c1,&c2,&c3))
      error->all(FLERR,"Invalid Special Bond Counts section in "
                 "molecule file");

    if (flag) {
      nspecial[i][0] = c1;
      nspecial[i][1] = c1+c2;
      nspecial[i][2] = c1+c2+c3;
    } else maxspecial = MAX(maxspecial,c1+c2+c3);
  }
}

/* ----------------------------------------------------------------------
   read special bond indices from file
------------------------------------------------------------------------- */

void Molecule::special_read(char *line)
{
  int m,nwords;
  char **words = new char*[maxspecial+1];

  for (int i = 0; i < natoms; i++) {
    readline(line);
    nwords = parse(line,words,maxspecial+1);
    if (nwords != nspecial[i][2]+1)
      error->all(FLERR,"Molecule file special list "
                 "does not match special count");

    for (m = 1; m < nwords; m++) {
      special[i][m-1] = ATOTAGINT(words[m]);
      if (special[i][m-1] <= 0 || special[i][m-1] > natoms ||
          special[i][m-1] == i+1)
        error->all(FLERR,"Invalid special atom index in molecule file");
    }
  }

  delete [] words;
}

/* ----------------------------------------------------------------------
   auto generate special bond info
------------------------------------------------------------------------- */

void Molecule::special_generate()
{
  int newton_bond = force->newton_bond;
  tagint atom1,atom2;
  int *count = new int[natoms];

  // temporary array for special atoms

  tagint **tmpspecial;
  memory->create(tmpspecial,natoms,atom->maxspecial,"molecule:tmpspecial");
  memset(&tmpspecial[0][0],0,sizeof(tagint)*natoms*atom->maxspecial);

  for (int i = 0; i < natoms; i++) count[i] = 0;

  // 1-2 neighbors

  if (newton_bond) {
    for (int i = 0; i < natoms; i++) {
      for (int j = 0; j < num_bond[i]; j++) {
        atom1 = i;
        atom2 = bond_atom[i][j]-1;
        nspecial[i][0]++;
        nspecial[atom2][0]++;
        if (count[i] >= atom->maxspecial || count[atom2] >= atom->maxspecial)
          error->one(FLERR,"Molecule auto special bond generation overflow");
        tmpspecial[i][count[i]++] = atom2 + 1;
        tmpspecial[atom2][count[atom2]++] = i + 1;
      }
    }
  } else {
    for (int i = 0; i < natoms; i++) {
      nspecial[i][0] = num_bond[i];
      for (int j = 0; j < num_bond[i]; j++) {
        atom1 = i;
        atom2 = bond_atom[i][j];
        if (count[atom1] >= atom->maxspecial)
          error->one(FLERR,"Molecule auto special bond generation overflow");
        tmpspecial[i][count[atom1]++] = atom2;
      }
    }
  }

  // 1-3 neighbors with no duplicates

  for (int i = 0; i < natoms; i++) nspecial[i][1] = nspecial[i][0];

  int dedup;
  for (int i = 0; i < natoms; i++) {
    for (int m = 0; m < nspecial[i][0]; m++) {
      for (int j = 0; j < nspecial[tmpspecial[i][m]-1][0]; j++) {
        dedup = 0;
        for (int k =0; k < count[i]; k++) {
          if (tmpspecial[tmpspecial[i][m]-1][j] == tmpspecial[i][k] ||
              tmpspecial[tmpspecial[i][m]-1][j] == i+1) {
            dedup = 1;
          }
        }
        if (!dedup) {
          if (count[i] >= atom->maxspecial)
            error->one(FLERR,"Molecule auto special bond generation overflow");
          tmpspecial[i][count[i]++] = tmpspecial[tmpspecial[i][m]-1][j];
          nspecial[i][1]++;
        }
      }
    }
  }

  // 1-4 neighbors with no duplicates

  for (int i = 0; i < natoms; i++) nspecial[i][2] = nspecial[i][1];

  for (int i = 0; i < natoms; i++) {
    for (int m = nspecial[i][0]; m < nspecial[i][1]; m++) {
      for (int j = 0; j < nspecial[tmpspecial[i][m]-1][0]; j++) {
        dedup = 0;
        for (int k =0; k < count[i]; k++) {
          if (tmpspecial[tmpspecial[i][m]-1][j] == tmpspecial[i][k] ||
              tmpspecial[tmpspecial[i][m]-1][j] == i+1) {
            dedup = 1;
          }
        }
        if (!dedup) {
          if (count[i] >= atom->maxspecial)
            error->one(FLERR,"Molecule auto special bond generation overflow");
          tmpspecial[i][count[i]++] = tmpspecial[tmpspecial[i][m]-1][j];
          nspecial[i][2]++;
        }
      }
    }
  }
  delete[] count;

  maxspecial = 0;
  for (int i = 0; i < natoms; i++)
    maxspecial = MAX(maxspecial,nspecial[i][2]);

  memory->create(special,natoms,maxspecial,"molecule:special");
  for (int i = 0; i < natoms; i++)
    for (int j = 0; j < nspecial[i][2]; j++)
      special[i][j] = tmpspecial[i][j];

  memory->destroy(tmpspecial);
}

/* ----------------------------------------------------------------------
   read SHAKE flags from file
------------------------------------------------------------------------- */

void Molecule::shakeflag_read(char *line)
{
  int tmp;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (2 != sscanf(line,"%d %d",&tmp,&shake_flag[i]))
      error->all(FLERR,"Invalid Shake Flags section in molecule file");
  }

  for (int i = 0; i < natoms; i++)
    if (shake_flag[i] < 0 || shake_flag[i] > 4)
      error->all(FLERR,"Invalid shake flag in molecule file");
}

/* ----------------------------------------------------------------------
   read SHAKE atom info from file
------------------------------------------------------------------------- */

void Molecule::shakeatom_read(char *line)
{
  int tmp, nmatch=0, nwant=0;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (shake_flag[i] == 1) {
      nmatch = sscanf(line,"%d " TAGINT_FORMAT " " TAGINT_FORMAT
                      " " TAGINT_FORMAT,&tmp,&shake_atom[i][0],
                      &shake_atom[i][1],&shake_atom[i][2]);
      nwant = 4;
    } else if (shake_flag[i] == 2) {
      nmatch = sscanf(line,"%d " TAGINT_FORMAT " " TAGINT_FORMAT,
                      &tmp,&shake_atom[i][0],&shake_atom[i][1]);
      nwant = 3;
    } else if (shake_flag[i] == 3) {
      nmatch = sscanf(line,"%d " TAGINT_FORMAT " " TAGINT_FORMAT
                      " " TAGINT_FORMAT,&tmp,&shake_atom[i][0],
                      &shake_atom[i][1],&shake_atom[i][2]);
      nwant = 4;
    } else if (shake_flag[i] == 4) {
      nmatch = sscanf(line,"%d " TAGINT_FORMAT " " TAGINT_FORMAT " "
                      TAGINT_FORMAT " " TAGINT_FORMAT,
                      &tmp,&shake_atom[i][0],&shake_atom[i][1],
                      &shake_atom[i][2],&shake_atom[i][3]);
      nwant = 5;
    }
    if (nmatch != nwant)
      error->all(FLERR,"Invalid shake atom in molecule file");
  }

  for (int i = 0; i < natoms; i++) {
    int m = shake_flag[i];
    if (m == 1) m = 3;
    for (int j = 0; j < m; j++)
      if (shake_atom[i][j] <= 0 || shake_atom[i][j] > natoms)
        error->all(FLERR,"Invalid shake atom in molecule file");
  }
}

/* ----------------------------------------------------------------------
   read SHAKE bond type info from file
------------------------------------------------------------------------- */

void Molecule::shaketype_read(char *line)
{
  int tmp, nmatch=0, nwant=0;
  for (int i = 0; i < natoms; i++) {
    readline(line);
    if (shake_flag[i] == 1) {
      nmatch = sscanf(line,"%d %d %d %d",&tmp,&shake_type[i][0],
                      &shake_type[i][1],&shake_type[i][2]);
      nwant = 4;
    } else if (shake_flag[i] == 2) {
      nmatch = sscanf(line,"%d %d",&tmp,&shake_type[i][0]);
      nwant = 2;
    } else if (shake_flag[i] == 3) {
      nmatch = sscanf(line,"%d %d %d",&tmp,&shake_type[i][0],
                      &shake_type[i][1]);
      nwant = 3;
    } else if (shake_flag[i] == 4) {
      nmatch = sscanf(line,"%d %d %d %d",&tmp,&shake_type[i][0],
                      &shake_type[i][1],&shake_type[i][2]);
      nwant = 4;
    }
    if (nmatch != nwant)
      error->all(FLERR,"Invalid shake type data in molecule file");
  }

  for (int i = 0; i < natoms; i++) {
    int m = shake_flag[i];
    if (m == 1) m = 3;
    for (int j = 0; j < m-1; j++)
      if (shake_type[i][j] <= 0)
        error->all(FLERR,"Invalid shake bond type in molecule file");
    if (shake_flag[i] == 1)
      if (shake_type[i][2] <= 0)
        error->all(FLERR,"Invalid shake angle type in molecule file");
  }
}

/* ----------------------------------------------------------------------
   read body params from file
   pflag = 0/1 for integer/double params
------------------------------------------------------------------------- */

void Molecule::body(int flag, int pflag, char *line)
{
  int i,ncount;

  int nparam = nibody;
  if (pflag) nparam = ndbody;

  int nword = 0;
  while (nword < nparam) {
    readline(line);

    ncount = atom->count_words(line);
    if (ncount == 0)
      error->one(FLERR,"Too few values in body section of molecule file");
    if (nword+ncount > nparam)
      error->all(FLERR,"Too many values in body section of molecule file");

    if (flag) {
      if (pflag == 0) {
        ibodyparams[nword++] = force->inumeric(FLERR,strtok(line," \t\n\r\f"));
        for (i = 1; i < ncount; i++)
          ibodyparams[nword++] =
            force->inumeric(FLERR,strtok(NULL," \t\n\r\f"));
      } else {
        dbodyparams[nword++] = force->numeric(FLERR,strtok(line," \t\n\r\f"));
        for (i = 1; i < ncount; i++)
          dbodyparams[nword++] =
            force->numeric(FLERR,strtok(NULL," \t\n\r\f"));
      }
    } else nword += ncount;
  }
}

/* ----------------------------------------------------------------------
   error check molecule attributes and topology against system settings
   flag = 0, just check this molecule
   flag = 1, check all molecules in set, this is 1st molecule in set
------------------------------------------------------------------------- */

void Molecule::check_attributes(int flag)
{
  int n = 1;
  if (flag) n = nset;
  int imol = atom->find_molecule(id);

  for (int i = imol; i < imol+n; i++) {
    Molecule *onemol = atom->molecules[imol];

    // check per-atom attributes of molecule
    // warn if not a match

    int mismatch = 0;
    if (onemol->qflag && !atom->q_flag) mismatch = 1;
    if (onemol->radiusflag && !atom->radius_flag) mismatch = 1;
    if (onemol->rmassflag && !atom->rmass_flag) mismatch = 1;

    if (mismatch && me == 0)
      error->warning(FLERR,
                     "Molecule attributes do not match system attributes");

    // for all atom styles, check nbondtype,etc

    mismatch = 0;
    if (atom->nbondtypes < onemol->nbondtypes) mismatch = 1;
    if (atom->nangletypes < onemol->nangletypes) mismatch = 1;
    if (atom->ndihedraltypes < onemol->ndihedraltypes) mismatch = 1;
    if (atom->nimpropertypes < onemol->nimpropertypes) mismatch = 1;

    if (mismatch)
      error->all(FLERR,"Molecule topology type exceeds system topology type");

    // for molecular atom styles, check bond_per_atom,etc + maxspecial
    // do not check for atom style template, since nothing stored per atom

    if (atom->molecular == 1) {
      if (atom->avec->bonds_allow &&
          atom->bond_per_atom < onemol->bond_per_atom) mismatch = 1;
      if (atom->avec->angles_allow &&
          atom->angle_per_atom < onemol->angle_per_atom) mismatch = 1;
      if (atom->avec->dihedrals_allow &&
          atom->dihedral_per_atom < onemol->dihedral_per_atom) mismatch = 1;
      if (atom->avec->impropers_allow &&
          atom->improper_per_atom < onemol->improper_per_atom) mismatch = 1;
      if (atom->maxspecial < onemol->maxspecial) mismatch = 1;

      if (mismatch)
        error->all(FLERR,"Molecule topology/atom exceeds system topology/atom");
    }

    // warn if molecule topology defined but no special settings

    if (onemol->bondflag && !onemol->specialflag)
      if (me == 0) error->warning(FLERR,"Molecule has bond topology "
                                  "but no special bond settings");
  }
}

/* ----------------------------------------------------------------------
   init all data structures to empty
------------------------------------------------------------------------- */

void Molecule::initialize()
{
  natoms = 0;
  nbonds = nangles = ndihedrals = nimpropers = 0;
  ntypes = 0;
  nbondtypes = nangletypes = ndihedraltypes = nimpropertypes = 0;
  nibody = ndbody = 0;

  bond_per_atom = angle_per_atom = dihedral_per_atom = improper_per_atom = 0;
  maxspecial = 0;

  xflag = typeflag = qflag = radiusflag = rmassflag = 0;
  bondflag = angleflag = dihedralflag = improperflag = 0;
  nspecialflag = specialflag = 0;
  shakeflag = shakeflagflag = shakeatomflag = shaketypeflag = 0;
  bodyflag = ibodyflag = dbodyflag = 0;

  centerflag = massflag = comflag = inertiaflag = 0;
  tag_require = 0;

  x = NULL;
  type = NULL;
  q = NULL;
  radius = NULL;
  rmass = NULL;

  num_bond = NULL;
  bond_type = NULL;
  bond_atom = NULL;

  num_angle = NULL;
  angle_type = NULL;
  angle_atom1 = angle_atom2 = angle_atom3 = NULL;

  num_dihedral = NULL;
  dihedral_type = NULL;
  dihedral_atom1 = dihedral_atom2 = dihedral_atom3 = dihedral_atom4 = NULL;

  num_improper = NULL;
  improper_type = NULL;
  improper_atom1 = improper_atom2 = improper_atom3 = improper_atom4 = NULL;

  nspecial = NULL;
  special = NULL;

  shake_flag = NULL;
  shake_atom = NULL;
  shake_type = NULL;

  ibodyparams = NULL;
  dbodyparams = NULL;

  dx = NULL;
  dxcom = NULL;
  dxbody = NULL;
}

/* ----------------------------------------------------------------------
   allocate all data structures
   also initialize values for data structures that are always allocated
------------------------------------------------------------------------- */

void Molecule::allocate()
{
  if (xflag) memory->create(x,natoms,3,"molecule:x");
  if (typeflag) memory->create(type,natoms,"molecule:type");
  if (qflag) memory->create(q,natoms,"molecule:q");
  if (radiusflag) memory->create(radius,natoms,"molecule:radius");
  if (rmassflag) memory->create(rmass,natoms,"molecule:rmass");

  // always allocate num_bond,num_angle,etc and nspecial
  // even if not in molecule file, initialize to 0
  // this is so methods that use these arrays don't have to check they exist

  memory->create(num_bond,natoms,"molecule:num_bond");
  for (int i = 0; i < natoms; i++) num_bond[i] = 0;
  memory->create(num_angle,natoms,"molecule:num_angle");
  for (int i = 0; i < natoms; i++) num_angle[i] = 0;
  memory->create(num_dihedral,natoms,"molecule:num_dihedral");
  for (int i = 0; i < natoms; i++) num_dihedral[i] = 0;
  memory->create(num_improper,natoms,"molecule:num_improper");
  for (int i = 0; i < natoms; i++) num_improper[i] = 0;
  memory->create(nspecial,natoms,3,"molecule:nspecial");
  for (int i = 0; i < natoms; i++)
    nspecial[i][0] = nspecial[i][1] = nspecial[i][2] = 0;

  if (specialflag)
    memory->create(special,natoms,maxspecial,"molecule:special");

  if (bondflag) {
    memory->create(bond_type,natoms,bond_per_atom,
                   "molecule:bond_type");
    memory->create(bond_atom,natoms,bond_per_atom,
                   "molecule:bond_atom");
  }

  if (angleflag) {
    memory->create(angle_type,natoms,angle_per_atom,
                   "molecule:angle_type");
    memory->create(angle_atom1,natoms,angle_per_atom,
                   "molecule:angle_atom1");
    memory->create(angle_atom2,natoms,angle_per_atom,
                   "molecule:angle_atom2");
    memory->create(angle_atom3,natoms,angle_per_atom,
                   "molecule:angle_atom3");
  }

  if (dihedralflag) {
    memory->create(dihedral_type,natoms,dihedral_per_atom,
                   "molecule:dihedral_type");
    memory->create(dihedral_atom1,natoms,dihedral_per_atom,
                   "molecule:dihedral_atom1");
    memory->create(dihedral_atom2,natoms,dihedral_per_atom,
                   "molecule:dihedral_atom2");
    memory->create(dihedral_atom3,natoms,dihedral_per_atom,
                   "molecule:dihedral_atom3");
    memory->create(dihedral_atom4,natoms,dihedral_per_atom,
                   "molecule:dihedral_atom4");
  }

  if (improperflag) {
    memory->create(improper_type,natoms,improper_per_atom,
                   "molecule:improper_type");
    memory->create(improper_atom1,natoms,improper_per_atom,
                   "molecule:improper_atom1");
    memory->create(improper_atom2,natoms,improper_per_atom,
                   "molecule:improper_atom2");
    memory->create(improper_atom3,natoms,improper_per_atom,
                   "molecule:improper_atom3");
    memory->create(improper_atom4,natoms,improper_per_atom,
                   "molecule:improper_atom4");
  }

  if (shakeflag) {
    memory->create(shake_flag,natoms,"molecule:shake_flag");
    memory->create(shake_atom,natoms,4,"molecule:shake_flag");
    memory->create(shake_type,natoms,3,"molecule:shake_flag");
  }

  if (bodyflag) {
    if (nibody) memory->create(ibodyparams,nibody,"molecule:ibodyparams");
    if (ndbody) memory->create(dbodyparams,ndbody,"molecule:dbodyparams");
  }
}

/* ----------------------------------------------------------------------
   deallocate all data structures
------------------------------------------------------------------------- */

void Molecule::deallocate()
{
  memory->destroy(x);
  memory->destroy(type);
  memory->destroy(q);
  memory->destroy(radius);
  memory->destroy(rmass);

  memory->destroy(num_bond);
  memory->destroy(bond_type);
  memory->destroy(bond_atom);

  memory->destroy(num_angle);
  memory->destroy(angle_type);
  memory->destroy(angle_atom1);
  memory->destroy(angle_atom2);
  memory->destroy(angle_atom3);

  memory->destroy(num_dihedral);
  memory->destroy(dihedral_type);
  memory->destroy(dihedral_atom1);
  memory->destroy(dihedral_atom2);
  memory->destroy(dihedral_atom3);
  memory->destroy(dihedral_atom4);

  memory->destroy(num_improper);
  memory->destroy(improper_type);
  memory->destroy(improper_atom1);
  memory->destroy(improper_atom2);
  memory->destroy(improper_atom3);
  memory->destroy(improper_atom4);

  memory->destroy(nspecial);
  memory->destroy(special);

  memory->destroy(shake_flag);
  memory->destroy(shake_atom);
  memory->destroy(shake_type);

  memory->destroy(dx);
  memory->destroy(dxcom);
  memory->destroy(dxbody);

  memory->destroy(ibodyparams);
  memory->destroy(dbodyparams);
}

/* ----------------------------------------------------------------------
   open molecule file
------------------------------------------------------------------------- */

void Molecule::open(char *file)
{
  fp = fopen(file,"r");
  if (fp == NULL) {
    char str[128];
    snprintf(str,128,"Cannot open molecule file %s",file);
    error->one(FLERR,str);
  }
}

/* ----------------------------------------------------------------------
   read and bcast a line
------------------------------------------------------------------------- */

void Molecule::readline(char *line)
{
  int n;
  if (me == 0) {
    if (fgets(line,MAXLINE,fp) == NULL) n = 0;
    else n = strlen(line) + 1;
  }
  MPI_Bcast(&n,1,MPI_INT,0,world);
  if (n == 0) error->all(FLERR,"Unexpected end of molecule file");
  MPI_Bcast(line,n,MPI_CHAR,0,world);
}

/* ----------------------------------------------------------------------
   extract keyword from line
   flag = 0, read and bcast line
   flag = 1, line has already been read
------------------------------------------------------------------------- */

void Molecule::parse_keyword(int flag, char *line, char *keyword)
{
  if (flag) {

    // read upto non-blank line plus 1 following line
    // eof is set to 1 if any read hits end-of-file

    int eof = 0;
    if (me == 0) {
      if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      while (eof == 0 && strspn(line," \t\n\r") == strlen(line)) {
        if (fgets(line,MAXLINE,fp) == NULL) eof = 1;
      }
      if (fgets(keyword,MAXLINE,fp) == NULL) eof = 1;
    }

    // if eof, set keyword empty and return

    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) {
      keyword[0] = '\0';
      return;
    }

    // bcast keyword line to all procs

    int n;
    if (me == 0) n = strlen(line) + 1;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);
  }

  // copy non-whitespace portion of line into keyword

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t'
         || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';
  strcpy(keyword,&line[start]);
}

/* ----------------------------------------------------------------------
   skip N lines of file
------------------------------------------------------------------------- */

void Molecule::skip_lines(int n, char *line)
{
  for (int i = 0; i < n; i++) readline(line);
}

/* ----------------------------------------------------------------------
   parse line into words separated by whitespace
   return # of words
   max = max pointers storable in words
------------------------------------------------------------------------- */

int Molecule::parse(char *line, char **words, int max)
{
  char *ptr;

  int nwords = 0;
  words[nwords++] = strtok(line," \t\n\r\f");

  while ((ptr = strtok(NULL," \t\n\r\f"))) {
    if (nwords < max) words[nwords] = ptr;
    nwords++;
  }

  return nwords;
}

/* ----------------------------------------------------------------------
   proc 0 prints molecule params
------------------------------------------------------------------------- */

/*

void Molecule::print()
{
  printf("MOLECULE %s\n",id);
  printf("  %d natoms\n",natoms);
  if (nbonds) printf("  %d nbonds\n",nbonds);
  if (nangles) printf("  %d nangles\n",nangles);
  if (ndihedrals) printf("  %d ndihedrals\n",ndihedrals);
  if (nimpropers) printf("  %d nimpropers\n",nimpropers);

  if (xflag) {
    printf(  "Coords:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %g %g %g\n",i+1,x[i][0],x[i][1],x[i][2]);
  }
  if (typeflag) {
    printf(  "Types:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %d\n",i+1,type[i]);
  }
  if (qflag) {
    printf(  "Charges:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %g\n",i+1,q[i]);
  }
  if (radiusflag) {
    printf(  "Radii:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %g\n",i+1,radius[i]);
  }
  if (rmassflag) {
    printf(  "Masses:\n");
    for (int i = 0; i < natoms; i++)
      printf("    %d %g\n",i+1,rmass[i]);
  }

  if (bondflag) {
    printf(  "Bonds:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d\n",i+1,num_bond[i]);
      for (int j = 0; j < num_bond[i]; j++)
        printf("      %d %d %d %d\n",j+1,bond_type[i][j],i+1,bond_atom[i][j]);
    }
  }
  if (angleflag) {
    printf(  "Angles:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d\n",i+1,num_angle[i]);
      for (int j = 0; j < num_angle[i]; j++)
        printf("      %d %d %d %d %d\n",
               j+1,angle_type[i][j],
               angle_atom1[i][j],angle_atom2[i][j],angle_atom3[i][j]);
    }
  }
  if (dihedralflag) {
    printf(  "Dihedrals:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d\n",i+1,num_dihedral[i]);
      for (int j = 0; j < num_dihedral[i]; j++)
        printf("      %d %d %d %d %d %d\n",
               j+1,dihedral_type[i][j],
               dihedral_atom1[i][j],dihedral_atom2[i][j],
               dihedral_atom3[i][j],dihedral_atom4[i][j]);
    }
  }
  if (improperflag) {
    printf(  "Impropers:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d\n",i+1,num_improper[i]);
      for (int j = 0; j < num_improper[i]; j++)
        printf("      %d %d %d %d %d %d\n",
               j+1,improper_type[i][j],
               improper_atom1[i][j],improper_atom2[i][j],
               improper_atom3[i][j],improper_atom4[i][j]);
    }
  }

  if (specialflag) {
    printf(  "Special neighs:\n");
    for (int i = 0; i < natoms; i++) {
      printf("    %d %d %d %d\n",i+1,
             nspecial[i][0],nspecial[i][1]-nspecial[i][0],
             nspecial[i][2]-nspecial[i][1]);
      printf("      ");
      for (int j = 0; j < nspecial[i][2]; j++)
        printf(" %d",special[i][j]);
      printf("\n");
    }
  }
}

*/
