      module mod_xc
c
      use dimensions
#ifdef CCSMCOUPLED
      use data_mct, only: mpicom_mct
#endif
c
      implicit none
c
c --- HYCOM communication interface.
c --- see README.src.mod_xc for more details.
c
c --- mxthrd= maximum number of OpenMP threads
      integer    mxthrd
      parameter (mxthrd=8)  ! NOMP=0,1
c
c --- halo size
      integer    nbdy
      parameter (nbdy=3)
c
c --- OpenMP will allocate jblk rows to each thread in turn
      integer    jblk
      parameter (jblk=(jdm+2*nbdy+mxthrd-1)/mxthrd)
c
c --- how far out the halo is valid (margin<=nbdy)
      integer, save :: margin
c
c --- actual extent of this tile is (i0+1:i0+ii,j0+1:j0+jj,1:kk)
      integer, save :: i0,j0,ii,jj
      integer      kk
      parameter   (kk=kdm)
c
c --- ms-1  = max. number of interruptions of any tile row or column by land
      integer    ms
      parameter (ms=100)  ! should be enough for any region
c
c --- information in /gindex/ keeps do loops from running into land
      integer, save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     &         ip,iu,iv,iq,ips,ipwocn,cplmsk
      integer, save, dimension (1-nbdy:jdm+nbdy,ms) :: 
     &         ifp,ilp,ifq,ilq,ifu,ilu,ifv,ilv
      integer, save, dimension (1-nbdy:idm+nbdy,ms) :: 
     &         jfp,jlp,jfq,jlq,jfu,jlu,jfv,jlv
      integer, save, dimension (1-nbdy:jdm+nbdy) :: 
     &         isp,isq,isu,isv
      integer, save, dimension (1-nbdy:idm+nbdy) :: 
     &         jsp,jsq,jsu,jsv
c
c --- line printer unit (stdout) and file unit
      integer, save :: lp,nfu
c-----------------------------------------------------------------------------
c
c --- tile dimensions and tile numbers (counting from 1), see xcspmd
      integer, public, save      :: ipr,  jpr,  ijpr,
     &                              mproc,nproc,mnproc
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
