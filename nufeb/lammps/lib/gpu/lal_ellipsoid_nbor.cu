// **************************************************************************
//                              ellipsoid_nbor.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for Ellipsoid neighbor routines
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_preprocessor.h"
#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
#else
texture<int4,1> pos_tex;
#endif
#else
#define pos_tex x_
#endif

// ---------------------------------------------------------------------------
// Unpack neighbors from dev_ij array into dev_nbor matrix for coalesced access
// -- Only unpack neighbors matching the specified inclusive range of forms
// -- Only unpack neighbors within cutoff
// ---------------------------------------------------------------------------
__kernel void kernel_nbor(const __global numtyp4 *restrict x_,
                          const __global numtyp2 *restrict cut_form,
                          const int ntypes,
                          __global int *dev_nbor,
                          const int nbor_pitch, const int start, const int inum,
                          const __global int *dev_ij,
                          const int form_low, const int form_high) {

  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X+start;

  if (ii<inum) {
    int i=dev_ij[ii];
    int nbor=ii+nbor_pitch;
    int numj=dev_ij[nbor];
    nbor+=nbor_pitch;
    int nbor_end=nbor+fast_mul(numj,nbor_pitch);
    int packed=ii+nbor_pitch+nbor_pitch;

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul(iw,ntypes);
    int newj=0;
    for ( ; nbor<nbor_end; nbor+=nbor_pitch) {
      int j=dev_ij[nbor];
      j &= NEIGHMASK;
      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      int mtype=itype+jtype;
      numtyp2 cf=cut_form[mtype];
      if (cf.y>=form_low && cf.y<=form_high) {
        // Compute r12;
        numtyp rsq=jx.x-ix.x;
        rsq*=rsq;
        numtyp t=jx.y-ix.y;
        rsq+=t*t;
        t=jx.z-ix.z;
        rsq+=t*t;

        if (rsq<cf.x) {
          dev_nbor[packed]=j;
          packed+=nbor_pitch;
          newj++;
        }
      }
    }
    dev_nbor[ii+nbor_pitch]=newj;
  }
}

// ---------------------------------------------------------------------------
// Unpack neighbors from dev_ij array into dev_nbor matrix for coalesced access
// -- Only unpack neighbors matching the specified inclusive range of forms
// -- Only unpack neighbors within cutoff
// -- Fast version of routine that uses shared memory for LJ constants
// ---------------------------------------------------------------------------
__kernel void kernel_nbor_fast(const __global numtyp4 *restrict x_,
                               const __global numtyp2 *restrict cut_form,
                               __global int *dev_nbor,
                               const int nbor_pitch, const int start,
                               const int inum,
                               const __global int *dev_ij,
                               const int form_low, const int form_high) {

  int ii=THREAD_ID_X;
  __local int form[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (ii<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    cutsq[ii]=cut_form[ii].x;
    form[ii]=cut_form[ii].y;
  }
  ii+=fast_mul((int)BLOCK_SIZE_X,(int)BLOCK_ID_X)+start;
  __syncthreads();

  if (ii<inum) {
    int i=dev_ij[ii];
    int nbor=ii+nbor_pitch;
    int numj=dev_ij[nbor];
    nbor+=nbor_pitch;
    int nbor_end=nbor+fast_mul(numj,nbor_pitch);
    int packed=ii+nbor_pitch+nbor_pitch;

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    int newj=0;
    for ( ; nbor<nbor_end; nbor+=nbor_pitch) {
      int j=dev_ij[nbor];
      j &= NEIGHMASK;
      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      int mtype=itype+jtype;

      if (form[mtype]>=form_low && form[mtype]<=form_high) {
        // Compute r12;
        numtyp rsq=jx.x-ix.x;
        rsq*=rsq;
        numtyp t=jx.y-ix.y;
        rsq+=t*t;
        t=jx.z-ix.z;
        rsq+=t*t;

        if (rsq<cutsq[mtype]) {
          dev_nbor[packed]=j;
          packed+=nbor_pitch;
          newj++;
        }
      }
    }
    dev_nbor[ii+nbor_pitch]=newj;
  }
}
