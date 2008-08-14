      module mod_xc
c
#ifdef CCSMCOUPLED
      use data_mct, only: mpicom_mct
#endif
c
      implicit none
c
c --- HYCOM communication interface.
c --- see README.src.mod_xc for more details.
c
      include 'dimensions.h'
c
c --- tile dimensions and tile numbers (counting from 1), see xcspmd
      integer, public, save      :: ipr,  jpr,  ijpr,
     &                              mproc,nproc,mnproc
c
c --- region type (-1==unknown,
c ---               0==  closed/closed,
c ---               1==periodic/closed,
c ---               2==periodic/arctic,
c ---               3==periodic/periodic
c ---               4==  closed/periodic)
      integer, public, save      :: nreg
c
c --- timers on, usually and default .true.
      logical, public, save      :: timer_on=.true.
c
c --- fill value for land, usually 0.0
      real,    public, save      :: vland
c
c --- xctilr halo options
      integer, public, parameter :: halo_ps=1, halo_pv=11,
     &                              halo_qs=2, halo_qv=12,
     &                              halo_us=3, halo_uv=13,
     &                              halo_vs=4, halo_vv=14
c
c --- xcsync stdout flushing options
      logical, public, parameter :: flush_lp=.true.,
     &                              no_flush=.false.
c
c --- generic subroutine names
      interface xcmaxr
         module procedure xcmaxr_0  ! rank 0 array (i.e. scalar)
         module procedure xcmaxr_1  ! rank 1 array
      end interface

      interface xcminr
         module procedure xcminr_0  ! rank 0 array (i.e. scalar)
         module procedure xcminr_1  ! rank 1 array
      end interface
c
c --- private timer variables, see xctmri
      character*6, private, dimension(97), save :: cc
      integer,     private,                save :: nxc
      integer,     private, dimension(97), save :: nc
      real*8,      private, dimension(97), save :: tc,t0
      real*8,      private, dimension(2),  save :: tcxc,tcxl
#if defined(MPI) || defined(SHMEM)
c
c --- private message passing data structures, see xcspmd
      integer, private, save :: idproc( 0: iqr+1,0:jqr+1),
     &                          idproc1(0:ijqr+1),idhalo(2),
     &                          i0_pe(iqr,jqr),ii_pe(iqr,jqr),
     &                          j0_pe(iqr,jqr),jj_pe(iqr,jqr),
     &                          mpe_1(     jqr),
     &                          mpe_e(     jqr),
     &                          mpe_i(itdm,jqr),npe_j(jtdm)
      integer, private, save :: i1sum(iqr,jqr),iisum(iqr,jqr)
      integer, private, save :: m0_top,i0_st(iqr+1),ii_st(iqr+1),
     &                          mm_top,i0_gt(iqr+1),ii_gt(iqr+1),
     &                          m0_bot,i0_sb(iqr+1),ii_sb(iqr+1),
     &                          mm_bot,i0_gb(iqr+1),ii_gb(iqr+1)
      integer, private, save :: null_tile
#endif
c
c --- actual module subroutines
      contains
#if defined(MPI) || defined(SHMEM)
# include "mod_xc_mp.inc"
#else
# include "mod_xc_sm.inc"
#endif
      end module mod_xc
