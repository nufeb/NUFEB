// **************************************************************************
//                                 tersoff.cu
//                             -------------------
//                              Trung Dac Nguyen
//
//  Device code for acceleration of the tersoff pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//       begin                : Thu April 17, 2014
//       email                : ndactrung@gmail.com
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_tersoff_extra.h"

#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
texture<float4> ts1_tex;
texture<float4> ts2_tex;
texture<float4> ts3_tex;
texture<float4> ts4_tex;
texture<float4> ts5_tex;
#else
texture<int4,1> pos_tex;
texture<int4> ts1_tex;
texture<int4> ts2_tex;
texture<int4> ts3_tex;
texture<int4> ts4_tex;
texture<int4> ts5_tex;
#endif

#else
#define pos_tex x_
#define ts1_tex ts1
#define ts2_tex ts2
#define ts3_tex ts3
#define ts4_tex ts4
#define ts5_tex ts5
#endif

//#define THREE_CONCURRENT

#define TWOTHIRD (numtyp)0.66666666666666666667

#define zeta_idx(nbor_mem, packed_mem, nbor_pitch, n_stride, t_per_atom,    \
                 i, nbor_j, offset_j, idx)                                  \
  if (nbor_mem==packed_mem) {                                               \
    int jj = (nbor_j-offset_j-2*nbor_pitch)/n_stride;                       \
    idx = jj*n_stride + i*t_per_atom + offset_j;                            \
  } else {                                                                  \
    idx = nbor_j;                                                           \
  }

#if (ARCH < 300)

#define store_answers_p(f, energy, virial, ii, inum, tid, t_per_atom,       \
                        offset, eflag, vflag, ans, engv)                    \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[6][BLOCK_PAIR];                                  \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=energy;                                                 \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<4; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    energy=red_acc[3][tid];                                                 \
    if (vflag>0) {                                                          \
      for (int r=0; r<6; r++)                                               \
        red_acc[r][tid]=virial[r];                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
        if (offset < s) {                                                   \
          for (int r=0; r<6; r++)                                           \
            red_acc[r][tid] += red_acc[r][tid+s];                           \
        }                                                                   \
      }                                                                     \
      for (int r=0; r<6; r++)                                               \
        virial[r]=red_acc[r][tid];                                          \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]+=energy*(acctyp)0.5;                                         \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]+=virial[i]*(acctyp)0.5;                                    \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    acctyp4 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }

#define acc_zeta(z, tid, t_per_atom, offset)                                \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[BLOCK_PAIR];                                     \
    red_acc[tid]=z;                                                         \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        red_acc[tid] += red_acc[tid+s];                                     \
      }                                                                     \
    }                                                                       \
    z=red_acc[tid];                                                         \
  }

#else

#define store_answers_p(f, energy, virial, ii, inum, tid, t_per_atom,       \
                        offset, eflag, vflag, ans, engv)                    \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      f.x += shfl_xor(f.x, s, t_per_atom);                                  \
      f.y += shfl_xor(f.y, s, t_per_atom);                                  \
      f.z += shfl_xor(f.z, s, t_per_atom);                                  \
      energy += shfl_xor(energy, s, t_per_atom);                            \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
          for (int r=0; r<6; r++)                                           \
            virial[r] += shfl_xor(virial[r], s, t_per_atom);                \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]+=energy*(acctyp)0.5;                                         \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]+=virial[i]*(acctyp)0.5;                                    \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    acctyp4 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }

#define acc_zeta(z, tid, t_per_atom, offset)                                \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      z += shfl_xor(z, s, t_per_atom);                                      \
    }                                                                       \
  }

#endif

__kernel void k_tersoff_short_nbor(const __global numtyp4 *restrict x_,
                                   const __global numtyp *restrict cutsq,
                                   const __global int *restrict map,
                                   const __global int *restrict elem2param,
                                   const int nelements, const int nparams,
                                   const __global int * dev_nbor,
                                   const __global int * dev_packed,
                                   __global int * dev_short_nbor,
                                   const int inum, const int nbor_pitch,
                                   const int t_per_atom) {
  __local int n_stride;
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    int ncount = 0;
    int m = nbor;
    dev_short_nbor[m] = 0;
    int nbor_short = nbor+n_stride;

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      int nj = j;
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq[ijparam]) {
        dev_short_nbor[nbor_short] = nj;
        nbor_short += n_stride;
        ncount++;
      }
    } // for nbor

    // store the number of neighbors for each thread
    dev_short_nbor[m] = ncount;

  } // if ii
}

// Tersoff is currently used for 3 elements at most: 3*3*3 = 27 entries
// while the block size should never be less than 32.
// SHARED_SIZE = 32 for now to reduce the pressure on the shared memory per block
// must be increased if there will be more than 3 elements in the future.

#define SHARED_SIZE 32

__kernel void k_tersoff_zeta(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict ts1_in,
                             const __global numtyp4 *restrict ts2_in,
                             const __global numtyp4 *restrict ts3_in,
                             const __global numtyp4 *restrict ts4_in,
                             const __global numtyp4 *restrict ts5_in,
                             const __global numtyp *restrict cutsq,
                             const __global int *restrict map,
                             const __global int *restrict elem2param,
                             const int nelements, const int nparams,
                             __global acctyp4 * zetaij,
                             const __global int * dev_nbor,
                             const __global int * dev_packed,
                             const __global int * dev_short_nbor,
                             const int eflag, const int inum,
                             const int nbor_pitch, const int t_per_atom) {
  __local int tpa_sq,n_stride;
  tpa_sq = fast_mul(t_per_atom,t_per_atom);

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  // must be increased if there will be more than 3 elements in the future.
  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  __local numtyp4 ts3[SHARED_SIZE];
  __local numtyp4 ts4[SHARED_SIZE];
  __local numtyp4 ts5[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
    ts3[tid]=ts3_in[tid];
    ts4[tid]=ts4_in[tid];
    ts5[tid]=ts5_in[tid];
  }

  acctyp z = (acctyp)0;

  __syncthreads();

  if (ii<inum) {
    int nbor_j, nbor_end, i, numj;
    const __global int* nbor_mem=dev_packed;
    int offset_j=offset/t_per_atom;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
              n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor_j];
      nbor_j += n_stride;
      nbor_end = nbor_j+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }
    int nborj_start = nbor_j;

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=nbor_mem[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute rij
      numtyp4 delr1, delr2;
      delr1.x = jx.x-ix.x;
      delr1.y = jx.y-ix.y;
      delr1.z = jx.z-ix.z;
      numtyp rsq1 = delr1.x*delr1.x+delr1.y*delr1.y+delr1.z*delr1.z;

      // compute zeta_ij
      z = (acctyp)0;

      int nbor_k = nborj_start-offset_j+offset_k;
      int k_end = nbor_end;
      if (dev_packed==dev_nbor) {
        int numk = dev_short_nbor[nbor_k-n_stride];
        k_end = nbor_k+fast_mul(numk,n_stride);
      }

      for ( ; nbor_k < k_end; nbor_k+=n_stride) {
        int k=nbor_mem[nbor_k];
        k &= NEIGHMASK;

        if (k == j) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex); //x_[k];
        int ktype=kx.w;
        ktype=map[ktype];
        int ijkparam=elem2param[itype*nelements*nelements+jtype*nelements+ktype];

        // Compute rik
        delr2.x = kx.x-ix.x;
        delr2.y = kx.y-ix.y;
        delr2.z = kx.z-ix.z;
        numtyp rsq2 = delr2.x*delr2.x+delr2.y*delr2.y+delr2.z*delr2.z;

        if (rsq2 > cutsq[ijkparam]) continue;

        numtyp4 ts1_ijkparam = ts1[ijkparam]; //fetch4(ts1_ijkparam,ijkparam,ts1_tex);
        numtyp ijkparam_lam3 = ts1_ijkparam.z;
        numtyp ijkparam_powermint = ts1_ijkparam.w;
        numtyp4 ts2_ijkparam = ts2[ijkparam]; //fetch4(ts2_ijkparam,ijkparam,ts2_tex);
        numtyp ijkparam_bigr = ts2_ijkparam.z;
        numtyp ijkparam_bigd = ts2_ijkparam.w;
        numtyp4 ts4_ijkparam = ts4[ijkparam]; //fetch4(ts4_ijkparam,ijkparam,ts4_tex);
        numtyp ijkparam_c = ts4_ijkparam.x;
        numtyp ijkparam_d = ts4_ijkparam.y;
        numtyp ijkparam_h = ts4_ijkparam.z;
        numtyp ijkparam_gamma = ts4_ijkparam.w;
        z += zeta(ijkparam_powermint, ijkparam_lam3, ijkparam_bigr, ijkparam_bigd,
                  ijkparam_c, ijkparam_d, ijkparam_h, ijkparam_gamma,
                  rsq1, rsq2, delr1, delr2);
      }

      // idx to zetaij is shifted by n_stride relative to nbor_j in dev_short_nbor
      int idx = nbor_j;
      if (dev_packed==dev_nbor) idx -= n_stride;
      acc_zeta(z, tid, t_per_atom, offset_k);

      numtyp4 ts1_ijparam = ts1[ijparam]; //fetch4(ts1_ijparam,ijparam,ts1_tex);
      numtyp ijparam_lam2 = ts1_ijparam.y;
      numtyp4 ts2_ijparam = ts2[ijparam]; //fetch4(ts2_ijparam,ijparam,ts2_tex);
      numtyp ijparam_bigb = ts2_ijparam.y;
      numtyp ijparam_bigr = ts2_ijparam.z;
      numtyp ijparam_bigd = ts2_ijparam.w;
      numtyp4 ts3_ijparam = ts3[ijparam]; //fetch4(ts3_ijparam,ijparam,ts3_tex);
      numtyp ijparam_c1 = ts3_ijparam.x;
      numtyp ijparam_c2 = ts3_ijparam.y;
      numtyp ijparam_c3 = ts3_ijparam.z;
      numtyp ijparam_c4 = ts3_ijparam.w;
      numtyp4 ts5_ijparam = ts5[ijparam]; //fetch4(ts5_ijparam,ijparam,ts5_tex);
      numtyp ijparam_beta = ts5_ijparam.x;
      numtyp ijparam_powern = ts5_ijparam.y;

      if (offset_k == 0) {
        numtyp fpfeng[4];
        force_zeta(ijparam_bigb, ijparam_bigr, ijparam_bigd, ijparam_lam2,
                   ijparam_beta, ijparam_powern, ijparam_c1, ijparam_c2, ijparam_c3,
                   ijparam_c4, rsq1, z, eflag, fpfeng);
        acctyp4 zij;
        zij.x = fpfeng[0];
        zij.y = fpfeng[1];
        zij.z = fpfeng[2];
        zij.w = z;
        zetaij[idx] = zij;
      }

    } // for nbor
  } // if ii
}

__kernel void k_tersoff_repulsive(const __global numtyp4 *restrict x_,
                                  const __global numtyp4 *restrict ts1_in,
                                  const __global numtyp4 *restrict ts2_in,
                                  const __global numtyp *restrict cutsq,
                                  const __global int *restrict map,
                                  const __global int *restrict elem2param,
                                  const int nelements, const int nparams,
                                  const __global int * dev_nbor,
                                  const __global int * dev_packed,
                                  const __global int * dev_short_nbor,
                                  __global acctyp4 *restrict ans,
                                  __global acctyp *restrict engv,
                                  const int eflag, const int vflag,
                                  const int inum, const int nbor_pitch,
                                  const int t_per_atom) {
  __local int n_stride;
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
  }

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();

  if (ii<inum) {
    int nbor, nbor_end, i, numj;
    const __global int* nbor_mem=dev_packed;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor];
      nbor += n_stride;
      nbor_end = nbor+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=nbor_mem[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute r12

      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      // rsq<cutsq[ijparam]
      numtyp feng[2];
      numtyp ijparam_lam1 = ts1[ijparam].x;
      numtyp4 ts2_ijparam = ts2[ijparam];
      numtyp ijparam_biga = ts2_ijparam.x;
      numtyp ijparam_bigr = ts2_ijparam.z;
      numtyp ijparam_bigd = ts2_ijparam.w;

      repulsive(ijparam_bigr, ijparam_bigd, ijparam_lam1, ijparam_biga,
                rsq, eflag, feng);

      numtyp force = feng[0];
      f.x+=delx*force;
      f.y+=dely*force;
      f.z+=delz*force;

      if (eflag>0)
        energy+=feng[1];
      if (vflag>0) {
        virial[0] += delx*delx*force;
        virial[1] += dely*dely*force;
        virial[2] += delz*delz*force;
        virial[3] += delx*dely*force;
        virial[4] += delx*delz*force;
        virial[5] += dely*delz*force;
      }
    } // for nbor

    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii

}

__kernel void k_tersoff_three_center(const __global numtyp4 *restrict x_,
                                     const __global numtyp4 *restrict ts1_in,
                                     const __global numtyp4 *restrict ts2_in,
                                     const __global numtyp4 *restrict ts4_in,
                                     const __global numtyp *restrict cutsq,
                                     const __global int *restrict map,
                                     const __global int *restrict elem2param,
                                     const int nelements, const int nparams,
                                     const __global acctyp4 *restrict zetaij,
                                     const __global int * dev_nbor,
                                     const __global int * dev_packed,
                                     const __global int * dev_short_nbor,
                                     __global acctyp4 *restrict ans,
                                     __global acctyp *restrict engv,
                                     const int eflag, const int vflag,
                                     const int inum,  const int nbor_pitch,
                                     const int t_per_atom, const int evatom) {
  __local int tpa_sq, n_stride;
  tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp lam3, powermint, bigr, bigd, c, d, h, gamma;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset); // offset ranges from 0 to tpa_sq-1

  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  __local numtyp4 ts4[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
    ts4[tid]=ts4_in[tid];
  }

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  numtyp tpainv = ucl_recip((numtyp)t_per_atom);

  __syncthreads();

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end;
    const __global int* nbor_mem=dev_packed;
    int offset_j=offset/t_per_atom;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
              n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor_j];
      nbor_j += n_stride;
      nbor_end = nbor_j+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }
    int nborj_start = nbor_j;

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=nbor_mem[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute r12
      numtyp delr1[3];
      delr1[0] = jx.x-ix.x;
      delr1[1] = jx.y-ix.y;
      delr1[2] = jx.z-ix.z;
      numtyp rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      numtyp r1 = ucl_sqrt(rsq1);
      numtyp r1inv = ucl_rsqrt(rsq1);

      // look up for zeta_ij
      // idx to zetaij is shifted by n_stride relative to nbor_j in dev_short_nbor
      int idx = nbor_j;
      if (dev_packed==dev_nbor) idx -= n_stride;
      acctyp4 zeta_ij = zetaij[idx]; // fetch(zeta_ij,idx,zeta_tex);
      numtyp force = zeta_ij.x*tpainv;
      numtyp prefactor = zeta_ij.y;
      f.x += delr1[0]*force;
      f.y += delr1[1]*force;
      f.z += delr1[2]*force;

      if (eflag>0) {
        energy+=zeta_ij.z*tpainv;
      }
      if (vflag>0) {
        numtyp mforce = -force;
        virial[0] += delr1[0]*delr1[0]*mforce;
        virial[1] += delr1[1]*delr1[1]*mforce;
        virial[2] += delr1[2]*delr1[2]*mforce;
        virial[3] += delr1[0]*delr1[1]*mforce;
        virial[4] += delr1[0]*delr1[2]*mforce;
        virial[5] += delr1[1]*delr1[2]*mforce;
      }

      int nbor_k = nborj_start-offset_j+offset_k;
      int k_end = nbor_end;
      if (dev_packed==dev_nbor) {
        int numk = dev_short_nbor[nbor_k-n_stride];
        k_end = nbor_k+fast_mul(numk,n_stride);
      }

      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=nbor_mem[nbor_k];
        k &= NEIGHMASK;

        if (j == k) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int ijkparam=elem2param[itype*nelements*nelements+jtype*nelements+ktype];

        numtyp delr2[3];
        delr2[0] = kx.x-ix.x;
        delr2[1] = kx.y-ix.y;
        delr2[2] = kx.z-ix.z;
        numtyp rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        if (rsq2 > cutsq[ijkparam]) continue;
        numtyp r2 = ucl_sqrt(rsq2);
        numtyp r2inv = ucl_rsqrt(rsq2);

        numtyp fi[3], fj[3], fk[3];
        numtyp4 ts1_ijkparam = ts1[ijkparam]; //fetch4(ts1_ijkparam,ijkparam,ts1_tex);
        lam3 = ts1_ijkparam.z;
        powermint = ts1_ijkparam.w;
        numtyp4 ts2_ijkparam = ts2[ijkparam]; //fetch4(ts2_ijkparam,ijkparam,ts2_tex);
        bigr = ts2_ijkparam.z;
        bigd = ts2_ijkparam.w;
        numtyp4 ts4_ijkparam = ts4[ijkparam]; //fetch4(ts4_ijkparam,ijkparam,ts4_tex);
        c = ts4_ijkparam.x;
        d = ts4_ijkparam.y;
        h = ts4_ijkparam.z;
        gamma = ts4_ijkparam.w;
        if (vflag>0)
          attractive(bigr, bigd, powermint, lam3, c, d, h, gamma,
                     prefactor, r1, r1inv, r2, r2inv, delr1, delr2, fi, fj, fk);
        else
          attractive_fi(bigr, bigd, powermint, lam3, c, d, h, gamma,
                        prefactor, r1, r1inv, r2, r2inv, delr1, delr2, fi);
        f.x += fi[0];
        f.y += fi[1];
        f.z += fi[2];

        if (vflag>0) {
          acctyp v[6];
          numtyp pre = (numtyp)2.0;
          if (evatom==1) pre = TWOTHIRD;
          v[0] = pre*(delr1[0]*fj[0] + delr2[0]*fk[0]);
          v[1] = pre*(delr1[1]*fj[1] + delr2[1]*fk[1]);
          v[2] = pre*(delr1[2]*fj[2] + delr2[2]*fk[2]);
          v[3] = pre*(delr1[0]*fj[1] + delr2[0]*fk[1]);
          v[4] = pre*(delr1[0]*fj[2] + delr2[0]*fk[2]);
          v[5] = pre*(delr1[1]*fj[2] + delr2[1]*fk[2]);

          virial[0] += v[0]; virial[1] += v[1]; virial[2] += v[2];
          virial[3] += v[3]; virial[4] += v[4]; virial[5] += v[5];
        }
      } // nbor_k
    } // for nbor_j

    store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,
                     offset,eflag,vflag,ans,engv);
  } // if ii
}

__kernel void k_tersoff_three_end(const __global numtyp4 *restrict x_,
                                  const __global numtyp4 *restrict ts1_in,
                                  const __global numtyp4 *restrict ts2_in,
                                  const __global numtyp4 *restrict ts4_in,
                                  const __global numtyp *restrict cutsq,
                                  const __global int *restrict map,
                                  const __global int *restrict elem2param,
                                  const int nelements, const int nparams,
                                  const __global acctyp4 *restrict zetaij,
                                  const __global int * dev_nbor,
                                  const __global int * dev_packed,
                                  const __global int * dev_ilist,
                                  const __global int * dev_short_nbor,
                                  __global acctyp4 *restrict ans,
                                  __global acctyp *restrict engv,
                                  const int eflag, const int vflag,
                                  const int inum,  const int nbor_pitch,
                                  const int t_per_atom, const int gpu_nbor) {
  __local int tpa_sq, n_stride;
  tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp lam3, powermint, bigr, bigd, c, d, h, gamma;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  __local numtyp4 ts4[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
    ts4[tid]=ts4_in[tid];
  }

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __local int red_acc[2*BLOCK_PAIR];

  __syncthreads();

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end, k_end;
    const __global int* nbor_mem=dev_packed;
    int offset_j=offset/t_per_atom;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
              n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    numtyp tpainv = ucl_recip((numtyp)t_per_atom);

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor_j];
      nbor_j += n_stride;
      nbor_end = nbor_j+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=nbor_mem[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute r12
      numtyp delr1[3];
      delr1[0] = jx.x-ix.x;
      delr1[1] = jx.y-ix.y;
      delr1[2] = jx.z-ix.z;
      numtyp rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      numtyp mdelr1[3];
      mdelr1[0] = -delr1[0];
      mdelr1[1] = -delr1[1];
      mdelr1[2] = -delr1[2];

      int nbor_k,numk;
      if (dev_nbor==dev_packed) {
        if (gpu_nbor) nbor_k=j+nbor_pitch;
        else nbor_k=dev_ilist[j]+nbor_pitch;
        numk=dev_nbor[nbor_k];
        nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
        k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk & (t_per_atom-1));
        nbor_k+=offset_k;
      } else {
        nbor_k=dev_ilist[j]+nbor_pitch;
        numk=dev_nbor[nbor_k];
        nbor_k+=nbor_pitch;
        nbor_k=dev_nbor[nbor_k];
        k_end=nbor_k+numk;
        nbor_k+=offset_k;
      }

      // recalculate numk and k_end for the use of short neighbor list
      if (dev_packed==dev_nbor) {
        numk = dev_short_nbor[nbor_k];
        nbor_k += n_stride;
        k_end = nbor_k+fast_mul(numk,n_stride);
      }
      int nbork_start = nbor_k;

      // look up for zeta_ji: find i in the j's neighbor list
      int m = tid / t_per_atom;
      int ijnum = -1;
      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=nbor_mem[nbor_k];
        k &= NEIGHMASK;
        if (k == i) {
          ijnum = nbor_k;
          red_acc[2*m+0] = ijnum;
          red_acc[2*m+1] = offset_k;
          break;
        }
      }

      numtyp r1 = ucl_sqrt(rsq1);
      numtyp r1inv = ucl_rsqrt(rsq1);
      int offset_kf;
      if (ijnum >= 0) {
        offset_kf = offset_k;
      } else {
        ijnum = red_acc[2*m+0];
        offset_kf = red_acc[2*m+1];
      }

      // idx to zetaij is shifted by n_stride relative to ijnum in dev_short_nbor
      int idx = ijnum;
      if (dev_packed==dev_nbor) idx -= n_stride;
      acctyp4 zeta_ji = zetaij[idx]; // fetch(zeta_ji,idx,zeta_tex);
      numtyp force = zeta_ji.x*tpainv;
      numtyp prefactor_ji = zeta_ji.y;
      f.x += delr1[0]*force;
      f.y += delr1[1]*force;
      f.z += delr1[2]*force;

      if (eflag>0) {
        energy+=zeta_ji.z*tpainv;
      }
      if (vflag>0) {
        numtyp mforce = -force;
        virial[0] += mdelr1[0]*mdelr1[0]*mforce;
        virial[1] += mdelr1[1]*mdelr1[1]*mforce;
        virial[2] += mdelr1[2]*mdelr1[2]*mforce;
        virial[3] += mdelr1[0]*mdelr1[1]*mforce;
        virial[4] += mdelr1[0]*mdelr1[2]*mforce;
        virial[5] += mdelr1[1]*mdelr1[2]*mforce;
      }

      // attractive forces
      for (nbor_k = nbork_start ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=nbor_mem[nbor_k];
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int jikparam=elem2param[jtype*nelements*nelements+itype*nelements+ktype];

        numtyp delr2[3];
        delr2[0] = kx.x-jx.x;
        delr2[1] = kx.y-jx.y;
        delr2[2] = kx.z-jx.z;
        numtyp rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        if (rsq2 > cutsq[jikparam]) continue;
        numtyp r2 = ucl_sqrt(rsq2);
        numtyp r2inv = ucl_rsqrt(rsq2);
        numtyp4 ts1_param, ts2_param, ts4_param;
        numtyp fi[3];

        ts1_param = ts1[jikparam]; //fetch4(ts1_jikparam,jikparam,ts1_tex);
        lam3 = ts1_param.z;
        powermint = ts1_param.w;
        ts2_param = ts2[jikparam]; //fetch4(ts2_jikparam,jikparam,ts2_tex);
        bigr = ts2_param.z;
        bigd = ts2_param.w;
        ts4_param = ts4[jikparam]; //fetch4(ts4_jikparam,jikparam,ts4_tex);
        c = ts4_param.x;
        d = ts4_param.y;
        h = ts4_param.z;
        gamma = ts4_param.w;
        attractive_fj(bigr, bigd, powermint, lam3, c, d, h, gamma,
                      prefactor_ji, r1, r1inv, r2, r2inv, mdelr1, delr2, fi);
        f.x += fi[0];
        f.y += fi[1];
        f.z += fi[2];

        // idx to zetaij is shifted by n_stride relative to nbor_k in dev_short_nbor
        int idx = nbor_k;
        if (dev_packed==dev_nbor) idx -= n_stride;

        acctyp4 zeta_jk = zetaij[idx]; // fetch(zeta_jk,idx,zeta_tex);
        numtyp prefactor_jk = zeta_jk.y;
        int jkiparam=elem2param[jtype*nelements*nelements+ktype*nelements+itype];
        ts1_param = ts1[jkiparam]; //fetch4(ts1_jkiparam,jkiparam,ts1_tex);
        lam3 = ts1_param.z;
        powermint = ts1_param.w;
        ts2_param = ts2[jkiparam]; //fetch4(ts2_jkiparam,jkiparam,ts2_tex);
        bigr = ts2_param.z;
        bigd = ts2_param.w;
        ts4_param = ts4[jkiparam]; //fetch4(ts4_jkiparam,jkiparam,ts4_tex);
        c = ts4_param.x;
        d = ts4_param.y;
        h = ts4_param.z;
        gamma = ts4_param.w;
        attractive_fk(bigr, bigd, powermint, lam3, c, d, h, gamma,
                      prefactor_jk, r2, r2inv, r1, r1inv, delr2, mdelr1, fi);
        f.x += fi[0];
        f.y += fi[1];
        f.z += fi[2];
      } // for nbor_k
    } // for nbor_j

    #ifdef THREE_CONCURRENT
    store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                  eflag,vflag,ans,engv);
    #else
    store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                    eflag,vflag,ans,engv);
    #endif
  } // if ii
}

__kernel void k_tersoff_three_end_vatom(const __global numtyp4 *restrict x_,
                                        const __global numtyp4 *restrict ts1_in,
                                        const __global numtyp4 *restrict ts2_in,
                                        const __global numtyp4 *restrict ts4_in,
                                        const __global numtyp *restrict cutsq,
                                        const __global int *restrict map,
                                        const __global int *restrict elem2param,
                                        const int nelements, const int nparams,
                                        const __global acctyp4 *restrict zetaij,
                                        const __global int * dev_nbor,
                                        const __global int * dev_packed,
                                        const __global int * dev_ilist,
                                        const __global int * dev_short_nbor,
                                        __global acctyp4 *restrict ans,
                                        __global acctyp *restrict engv,
                                        const int eflag, const int vflag,
                                        const int inum,  const int nbor_pitch,
                                        const int t_per_atom, const int gpu_nbor) {
  __local int tpa_sq, n_stride;
  tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp lam3, powermint, bigr, bigd, c, d, h, gamma;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  __local numtyp4 ts4[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
    ts4[tid]=ts4_in[tid];
  }

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __local int red_acc[2*BLOCK_PAIR];

  __syncthreads();

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end, k_end;
    const __global int* nbor_mem = dev_packed;
    int offset_j=offset/t_per_atom;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
              n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    numtyp tpainv = ucl_recip((numtyp)t_per_atom);

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor_j];
      nbor_j += n_stride;
      nbor_end = nbor_j+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=nbor_mem[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute r12
      numtyp delr1[3];
      delr1[0] = jx.x-ix.x;
      delr1[1] = jx.y-ix.y;
      delr1[2] = jx.z-ix.z;
      numtyp rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      numtyp mdelr1[3];
      mdelr1[0] = -delr1[0];
      mdelr1[1] = -delr1[1];
      mdelr1[2] = -delr1[2];

      int nbor_k,numk;
      if (dev_nbor==dev_packed) {
        if (gpu_nbor) nbor_k=j+nbor_pitch;
        else nbor_k=dev_ilist[j]+nbor_pitch;
        numk=dev_nbor[nbor_k];
        nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
        k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk & (t_per_atom-1));
        nbor_k+=offset_k;
      } else {
        nbor_k=dev_ilist[j]+nbor_pitch;
        numk=dev_nbor[nbor_k];
        nbor_k+=nbor_pitch;
        nbor_k=dev_nbor[nbor_k];
        k_end=nbor_k+numk;
        nbor_k+=offset_k;
      }

      // recalculate numk and k_end for the use of short neighbor list
      if (dev_packed==dev_nbor) {
        numk = dev_short_nbor[nbor_k];
        nbor_k += n_stride;
        k_end = nbor_k+fast_mul(numk,n_stride);
      }
      int nbork_start = nbor_k;

      // look up for zeta_ji
      int m = tid / t_per_atom;
      int ijnum = -1;
      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=nbor_mem[nbor_k];
        k &= NEIGHMASK;
        if (k == i) {
          ijnum = nbor_k;
          red_acc[2*m+0] = ijnum;
          red_acc[2*m+1] = offset_k;
          break;
        }
      }

      numtyp r1 = ucl_sqrt(rsq1);
      numtyp r1inv = ucl_rsqrt(rsq1);
      int offset_kf;
      if (ijnum >= 0) {
        offset_kf = offset_k;
      } else {
        ijnum = red_acc[2*m+0];
        offset_kf = red_acc[2*m+1];
      }

      // idx to zetaij is shifted by n_stride relative to ijnum in dev_short_nbor
      int idx = ijnum;
      if (dev_packed==dev_nbor) idx -= n_stride;
      acctyp4 zeta_ji = zetaij[idx]; //  fetch(zeta_ji,idx,zeta_tex);
      numtyp force = zeta_ji.x*tpainv;
      numtyp prefactor_ji = zeta_ji.y;
      f.x += delr1[0]*force;
      f.y += delr1[1]*force;
      f.z += delr1[2]*force;

      if (eflag>0) {
        energy+=zeta_ji.z*tpainv;
      }
      if (vflag>0) {
        numtyp mforce = -force;
        virial[0] += mdelr1[0]*mdelr1[0]*mforce;
        virial[1] += mdelr1[1]*mdelr1[1]*mforce;
        virial[2] += mdelr1[2]*mdelr1[2]*mforce;
        virial[3] += mdelr1[0]*mdelr1[1]*mforce;
        virial[4] += mdelr1[0]*mdelr1[2]*mforce;
        virial[5] += mdelr1[1]*mdelr1[2]*mforce;
      }

      // attractive forces
      for (nbor_k = nbork_start; nbor_k<k_end; nbor_k+=n_stride) {
        int k=nbor_mem[nbor_k];
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int jikparam=elem2param[jtype*nelements*nelements+itype*nelements+ktype];

        numtyp delr2[3];
        delr2[0] = kx.x-jx.x;
        delr2[1] = kx.y-jx.y;
        delr2[2] = kx.z-jx.z;
        numtyp rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        if (rsq2 > cutsq[jikparam]) continue;
        numtyp r2 = ucl_sqrt(rsq2);
        numtyp r2inv = ucl_rsqrt(rsq2);

        numtyp fi[3], fj[3], fk[3];
        numtyp4 ts1_param, ts2_param, ts4_param;
        ts1_param = ts1[jikparam]; //fetch4(ts1_jikparam,jikparam,ts1_tex);
        lam3 = ts1_param.z;
        powermint = ts1_param.w;
        ts2_param = ts2[jikparam]; //fetch4(ts2_jikparam,jikparam,ts2_tex);
        bigr = ts2_param.z;
        bigd = ts2_param.w;
        ts4_param = ts4[jikparam]; //fetch4(ts4_jikparam,jikparam,ts4_tex);
        c = ts4_param.x;
        d = ts4_param.y;
        h = ts4_param.z;
        gamma = ts4_param.w;
        attractive(bigr, bigd, powermint, lam3, c, d, h, gamma,
                   prefactor_ji, r1, r1inv, r2, r2inv, mdelr1, delr2, fi, fj, fk);
        f.x += fj[0];
        f.y += fj[1];
        f.z += fj[2];

        virial[0] += TWOTHIRD*(mdelr1[0]*fj[0] + delr2[0]*fk[0]);
        virial[1] += TWOTHIRD*(mdelr1[1]*fj[1] + delr2[1]*fk[1]);
        virial[2] += TWOTHIRD*(mdelr1[2]*fj[2] + delr2[2]*fk[2]);
        virial[3] += TWOTHIRD*(mdelr1[0]*fj[1] + delr2[0]*fk[1]);
        virial[4] += TWOTHIRD*(mdelr1[0]*fj[2] + delr2[0]*fk[2]);
        virial[5] += TWOTHIRD*(mdelr1[1]*fj[2] + delr2[1]*fk[2]);

        // idx to zetaij is shifted by n_stride relative to nbor_k in dev_short_nbor
        int idx = nbor_k;
        if (dev_packed==dev_nbor) idx -= n_stride;
        acctyp4 zeta_jk = zetaij[idx]; // fetch(zeta_jk,idx,zeta_tex);
        numtyp prefactor_jk = zeta_jk.y;

        int jkiparam=elem2param[jtype*nelements*nelements+ktype*nelements+itype];
        ts1_param = ts1[jkiparam]; //fetch4(ts1_jkiparam,jkiparam,ts1_tex);
        lam3 = ts1_param.z;
        powermint = ts1_param.w;
        ts2_param = ts2[jkiparam]; //fetch4(ts2_jkiparam,jkiparam,ts2_tex);
        bigr = ts2_param.z;
        bigd = ts2_param.w;
        ts4_param = ts4[jkiparam]; //fetch4(ts4_jkiparam,jkiparam,ts4_tex);
        c = ts4_param.x;
        d = ts4_param.y;
        h = ts4_param.z;
        gamma = ts4_param.w;
        attractive(bigr, bigd, powermint, lam3, c, d, h, gamma,
                   prefactor_jk, r2, r2inv, r1, r1inv, delr2, mdelr1, fi, fj, fk);
        f.x += fk[0];
        f.y += fk[1];
        f.z += fk[2];

        virial[0] += TWOTHIRD*(delr2[0]*fj[0] + mdelr1[0]*fk[0]);
        virial[1] += TWOTHIRD*(delr2[1]*fj[1] + mdelr1[1]*fk[1]);
        virial[2] += TWOTHIRD*(delr2[2]*fj[2] + mdelr1[2]*fk[2]);
        virial[3] += TWOTHIRD*(delr2[0]*fj[1] + mdelr1[0]*fk[1]);
        virial[4] += TWOTHIRD*(delr2[0]*fj[2] + mdelr1[0]*fk[2]);
        virial[5] += TWOTHIRD*(delr2[1]*fj[2] + mdelr1[1]*fk[2]);
      }
    } // for nbor

    #ifdef THREE_CONCURRENT
    store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                  eflag,vflag,ans,engv);
    #else
    store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                    eflag,vflag,ans,engv);
    #endif
  } // if ii
}

