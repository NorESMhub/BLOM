! ------------------------------------------------------------------------------
! Copyright (C) 2005 HYCOM Consortium and contributors
! Copyright (C) 2006-2024 Lars Inge Enstad, Mats Bentsen, Alok Kumar Gupta,
!                         Mariana Vertenstein
!
! This file is part of BLOM.
!
! BLOM is free software: you can redistribute it and/or modify it under the
! terms of the GNU Lesser General Public License as published by the Free
! Software Foundation, either version 3 of the License, or (at your option)
! any later version.
!
! BLOM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with BLOM. If not, see <https://www.gnu.org/licenses/>.
! ------------------------------------------------------------------------------

module mod_xc

  ! Note that dimensions.F is auto-generated
  use dimensions, only: idm,jdm,kdm,itdm,jtdm,iqr,jqr,ijqr,&
                        ii_pe,jj_pe,i0_pe,j0_pe,nreg
  use mod_wtime,  only: wtime
  use mod_ifdefs, only: use_TIMER, use_DEBUG_TIMER, use_DEBUG_TIMER_ALL, &
                        use_DEBUG_ALL, use_ARCTIC

  implicit none
  public

  ! HYCOM communication interface.
  ! see README.src.mod_xc for more details.

  ! Externally set mpi communicator
  integer :: mpicom_external = -1

  ! mxthrd= maximum number of OpenMP threads
  integer, parameter :: mxthrd=8  ! NOMP = 0,1

  ! halo size
  integer, parameter :: nbdy = 3

  ! OpenMP will allocate jblk rows to each thread in turn
  integer, parameter :: jblk = (jdm+2*nbdy+mxthrd-1)/mxthrd

  ! how far out the halo is valid (margin<=nbdy)
  integer :: margin

  ! actual extent of this tile is (i0+1:i0+ii,j0+1:j0+jj,1:kk)
  integer :: i0,j0,ii,jj
  integer, parameter :: kk = kdm

  ! ms-1  = max. number of interruptions of any tile row or column by land
  integer, parameter :: ms = 100  ! should be enough for any region

  ! information in /gindex/ keeps do loops from running into land
  ! Set in mod_bigrid.F90
  integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ip
  integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iu
  integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iv
  integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: iq
  integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ips
  integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ipwocn
  integer, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: cplmsk
  integer, dimension (1-nbdy:jdm+nbdy,ms) :: ifp
  integer, dimension (1-nbdy:jdm+nbdy,ms) :: ilp
  integer, dimension (1-nbdy:jdm+nbdy,ms) :: ifq
  integer, dimension (1-nbdy:jdm+nbdy,ms) :: ilq
  integer, dimension (1-nbdy:jdm+nbdy,ms) :: ifu
  integer, dimension (1-nbdy:jdm+nbdy,ms) :: ilu
  integer, dimension (1-nbdy:jdm+nbdy,ms) :: ifv
  integer, dimension (1-nbdy:jdm+nbdy,ms) :: ilv
  integer, dimension (1-nbdy:idm+nbdy,ms) :: jfp
  integer, dimension (1-nbdy:idm+nbdy,ms) :: jlp
  integer, dimension (1-nbdy:idm+nbdy,ms) :: jfq
  integer, dimension (1-nbdy:idm+nbdy,ms) :: jlq
  integer, dimension (1-nbdy:idm+nbdy,ms) :: jfu
  integer, dimension (1-nbdy:idm+nbdy,ms) :: jlu
  integer, dimension (1-nbdy:idm+nbdy,ms) :: jfv
  integer, dimension (1-nbdy:idm+nbdy,ms) :: jlv
  integer, dimension (1-nbdy:jdm+nbdy) :: isp
  integer, dimension (1-nbdy:jdm+nbdy) :: isq
  integer, dimension (1-nbdy:jdm+nbdy) :: isu
  integer, dimension (1-nbdy:jdm+nbdy) :: isv
  integer, dimension (1-nbdy:idm+nbdy) :: jsp
  integer, dimension (1-nbdy:idm+nbdy) :: jsq
  integer, dimension (1-nbdy:idm+nbdy) :: jsu
  integer, dimension (1-nbdy:idm+nbdy) :: jsv

  ! line printer unit (stdout) and file unit with default values 6 and
  ! 12, respectively
  integer :: lp=6

  ! tile dimensions and tile numbers (counting from 1), see xcspmd
  integer, public :: ipr, jpr, ijpr, mproc, nproc, mnproc

  ! timers on, usually and default .true.
  logical, public :: timer_on = .true.

  ! fill value for land, usually 0.0
  real, public :: vland

  ! xctilr halo options
  integer, public, parameter :: halo_ps=1, halo_pv = 11
  integer, public, parameter :: halo_qs=2, halo_qv = 12
  integer, public, parameter :: halo_us=3, halo_uv = 13
  integer, public, parameter :: halo_vs=4, halo_vv = 14

  ! xcsync stdout flushing options
  logical, public, parameter :: flush_lp = .true.
  logical, public, parameter :: no_flush = .false.

  ! generic subroutine names
  interface xcmax
    module procedure xcmax_i0  ! rank 0 integer array (i.e. scalar)
    module procedure xcmax_i1  ! rank 1 integer array
    module procedure xcmax_r0  ! rank 0 real array (i.e. scalar)
    module procedure xcmax_r1  ! rank 1 real array
  end interface xcmax

  interface xcmin
    module procedure xcmin_i0  ! rank 0 integer array (i.e. scalar)
    module procedure xcmin_i1  ! rank 1 integer array
    module procedure xcmin_r0  ! rank 0 real array (i.e. scalar)
    module procedure xcmin_r1  ! rank 1 real array
  end interface xcmin

  interface xcbcst
    module procedure xcbcst_i0 ! rank 0 integer array (i.e. scalar)
    module procedure xcbcst_i1 ! rank 1 integer array
    module procedure xcbcst_r0 ! rank 0 real array (i.e. scalar)
    module procedure xcbcst_r1 ! rank 1 real array
    module procedure xcbcst_l0 ! rank 0 logical array
    module procedure xcbcst_l1 ! rank 1 logical array
    module procedure xcbcst_c  ! rank 1 character array
  end interface xcbcst

  ! private timer variables, see xctmri
  character*6, private, dimension(97) :: cc
  integer,     private                :: nxc
  integer,     private, dimension(97) :: nc
  real(8),     private, dimension(97) :: tc,t0
  real(8),     private, dimension(2)  :: tcxc,tcxl

  ! The following is only used for MPI
  integer, public   :: idproc( 0: iqr+1,0:jqr+1)
  integer, public   :: idproc1(0:ijqr+1),idhalo(2)
  integer, public   :: mpe_1(jqr)
  integer, public   :: mpe_e(jqr)
  integer, public   :: mpe_i(itdm,jqr),npe_j(jtdm)

#if defined(MPI)
  include 'mpif.h'

  ! set in xcspmd
  integer, protected :: mpicomm
  integer, protected :: mpierr
  integer, protected :: mpireq(4)
  integer, protected :: mpistat(mpi_status_size,4*max(iqr,jqr))

  ! private message passing data structures, see xcspmd
  integer, private  :: i1sum(iqr,jqr),iisum(iqr,jqr)
  integer, private  :: m0_top,i0_st(iqr+1),ii_st(iqr+1)
  integer, private  :: mm_top,i0_gt(iqr+1),ii_gt(iqr+1)
  integer, private  :: m0_bot,i0_sb(iqr+1),ii_sb(iqr+1)
  integer, private  :: mm_bot,i0_gb(iqr+1),ii_gb(iqr+1)
  integer, private  :: null_tile

  real, private :: al(itdm,jdm),ala(itdm,jdm,jqr),at(idm*jdm),ata(idm*jdm,iqr)
#endif

contains

!**************************************************************************************
#if defined(MPI)
!**************************************************************************************

  ! MTYPER        mpi type for real
  ! MTYPED        mpi type for real*8
  ! MPI_SEND      either mpi_send  or mpi_ssend
  ! MPI_ISEND     either mpi_isend or mpi_issend

#if defined(NOMPIR8)
#if defined(REAL4)
# define MTYPER mpi_real
# define MTYPED mpi_double_precision
#else ! REAL8
# define MTYPER mpi_double_precision
# define MTYPED mpi_double_precision
#endif

#else ! most MPI's allow mpi_real[48]

#if defined(REAL4)
# define MTYPER mpi_real4
# define MTYPED mpi_real8
#else ! REAL8
# define MTYPER mpi_real8
# define MTYPED mpi_real8
#endif
#endif

#if defined(SSEND)
# define MPI_SEND mpi_ssend
# define MPI_ISEND mpi_issend
#else
# define MPI_SEND mpi_send
# define MPI_ISEND mpi_isend
#endif ! SSEND:else

  !-----------------------------------------------------------------------
  ! auxillary routines that involve off-processor communication.
  ! message passing version, contained in module mod_xc.
  ! author:  Alan J. Wallcraft,  NRL.
  !-----------------------------------------------------------------------

  subroutine xcaget(aa, a, mnflg)
    real,    intent(out)   :: aa(itdm,jtdm) ! non-tiled target array
    real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! tiled source array
    integer, intent(in)    :: mnflg ! node return flag

    !-----------
    ! convert an entire 2-D array from tiled to non-tiled layout.
    ! mnflg selects which nodes must return the array
    ! = 0; all nodes
    ! = n; node number n (mnproc=n)
    !-----------

    integer :: mpireqa(jpr),mpireqb(ipr)

    integer :: i,j,l,mp,np,mnp
    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 1)
        nxc = 1
      end if
    end if

    ! gather each row of tiles onto the first tile in the row.

    if (mproc == mpe_1(nproc)) then
      do j= 1,jj
        do i= 1,i0
          al(i,j) = vland
        end do
        do i= 1,ii
          al(i0+i,j) = a(i,j)
        end do
        do i= i0+ii+1,itdm
          al(i,j) = vland
        end do
      end do
      l = 0
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = l + 1
        call MPI_IRECV(ata(1,l),ii_pe(mp,nproc)*jj,MTYPER, &
             idproc(mp,nproc), 9941, &
             mpicomm, mpireqb(l), mpierr)
      end do
      call mpi_barrier(mpicomm, mpierr)
      l = 0
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = l + 1
        call MPI_WAIT(mpireqb(l), mpistat, mpierr)
        do j= 1,jj
          do i= 1,ii_pe(mp,nproc)
            al(i0_pe(mp,nproc)+i,j) = ata(i+(j-1)*ii_pe(mp,nproc),l)
          end do
        end do
      end do
    else  !mproc>1
      do j= 1,jj
        do i= 1,ii
          at(i+(j-1)*ii) = a(i,j)
        end do
      end do
      call mpi_barrier(mpicomm, mpierr)
      call MPI_SEND(at,ii*jj,MTYPER, &
           idproc(mpe_1(nproc),nproc), 9941, &
           mpicomm, mpierr)
    end if

    ! gather each row of tiles onto the output array.

    mnp = max(mnflg,1)

    if (mnproc == mnp) then
      if (mproc == mpe_1(nproc)) then
        do j= 1,jj
          do i= 1,itdm
            aa(i,j+j0) = al(i,j)
          end do
        end do
      end if
      l = 0
      do np= 1,jpr
        mp = mpe_1(np)
        if (idproc(mp,np) /= idproc(mproc,nproc)) then
          l = l + 1
          call MPI_IRECV(ala(1,1,l),itdm*jj_pe(mp,np),MTYPER, &
               idproc(mp,np), 9942,mpicomm, mpireqa(l), mpierr)
        end if
      end do
      call mpi_barrier(mpicomm, mpierr)
      l = 0
      do np= 1,jpr
        mp = mpe_1(np)
        if (idproc(mp,np) /= idproc(mproc,nproc)) then
          l = l + 1
          call MPI_WAIT(mpireqa(l), mpistat, mpierr)
          do j= 1,jj_pe(mp,np)
            do i= 1,itdm
              aa(i,j+j0_pe(mp,np)) = ala(i,j,l)
            end do
          end do
        end if
      end do
    else if (mproc == mpe_1(nproc)) then
      call mpi_barrier(mpicomm, mpierr)
      call MPI_SEND(al,itdm*jj,MTYPER,idproc1(mnp), 9942,mpicomm, mpierr)
    else
      call mpi_barrier(mpicomm, mpierr)
    end if

    if (mnflg == 0) then
      call mpi_bcast(aa,itdm*jtdm,MTYPER,idproc1(1),mpicomm,mpierr)
    end if

    if (use_TIMER) then
      if (nxc ==  1) then
        call xctmr1( 1)
        nxc = 0
      end if
    end if
  end subroutine xcaget

  !-----------------------------------------------------------------------
  subroutine xcaput(aa, a, mnflg)
    real,    intent(inout) :: aa(itdm,jtdm)
    real,    intent(out)   :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
    integer, intent(in)    :: mnflg

    !-----------
    !  1) convert an entire 2-D array from non-tiled to tiled layout.

    !  3) mnflg selects which nodes must contain the non-tiled array
    !        = 0; all nodes
    !        = n; node number n (mnproc=n)
    !     if mnflg.ne.0 the array aa may be broadcast to all nodes,
    !      so aa must exist and be overwritable on all nodes.

    !  4) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aa              real           input     non-tiled source array
    !    a               real           output    tiled target array
    !    mnflg           integer        input     node source flag
    !-----------

    integer :: mpireqa(jpr),mpireqb(ipr)
    integer :: i,j,np,mp,mnp

    if (use_TIMER) then
      if (nxc == 0) then
        call mpi_barrier(mpicomm, mpierr)
        call xctmr0( 4)
        nxc = 4
      end if
    end if

    ! use xclput for now,
    ! this is slow for mnflg.ne.0, but easy to implement.

    if (mnflg /= 0) then
      ! "broadcast" row sections of aa to all processors in the row.
      if (mnproc /= mnflg) then
        aa(:,:) = vland
      end if

      if (mnproc == mnflg) then
        j = 0
        do np= 1,jpr
          mp = mpe_1(np)
          if (np == nproc .and. mp == mproc) then
            al(:,1:jj) = aa(:,j0+1:j0+jj)
          else
            j = j + 1
            call MPI_ISEND(aa(1,j0_pe(mp,np)+1),itdm*jj_pe(mp,np),MTYPER,idproc(mp,np), 9951, &
                 mpicomm, mpireqa(j), mpierr)
          end if
        end do
        call mpi_waitall( j, mpireqa, mpistat, mpierr)
      else if (mproc == mpe_1(nproc)) then
        call MPI_RECV(al,itdm*jj,MTYPER,idproc1(mnflg), 9951,mpicomm, mpistat, mpierr)
      end if

      if (mproc == mpe_1(nproc)) then
        i = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          i = i + 1
          call MPI_ISEND(al,itdm*jj,MTYPER,idproc(mp,nproc), 9952,mpicomm, mpireqb(i), mpierr)
        end do
        call mpi_waitall( i, mpireqb, mpistat, mpierr)
      else
        call MPI_RECV(al,itdm*jj,MTYPER,idproc(mpe_1(nproc),nproc), 9952,mpicomm, mpistat, mpierr)
      end if

      aa(:,j0+1:j0+jj) = al(:,1:jj)
    end if

    do j= 1,jtdm
      call xclput(aa(1,j),itdm, a, 1,j,1,0)
    end do

    if (use_TIMER) then
      if (nxc ==  4) then
        call xctmr1( 4)
        nxc = 0
      end if
    end if
  end subroutine xcaput

  !-----------------------------------------------------------------------
  subroutine xceget(aelem, a, ia,ja)
    real,    intent(out)   :: aelem
    real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
    integer, intent(in)    :: ia,ja

    !-----------
    !  1) find the value of a(ia,ja) on the non-tiled 2-D grid.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aelem           real           output    required element
    !    a               real           input     source array
    !    ia              integer        input     1st index into a
    !    ja              integer        input     2nd index into a
    !  3) the global variable vland is returned when outside active tiles.
    !-----------

    ! double buffer to reduce the number of required barriers.
    real :: elem(0:1)
    integer, save :: kdb = 0
    integer :: i,j,mp,np

    kdb = mod(kdb+1,2)  ! the least recently used of the two buffers

    ! find the host tile.
    np = npe_j(ja)
    mp = mpe_i(ia,np)

    if (mp <= 0) then
      ! no tile.
      elem(kdb) = vland
    else if (mp == mproc .and. np == nproc) then
      ! this tile.
      i = ia - i0
      j = ja - j0
      elem(kdb) = a(i,j)
      call mpi_bcast(elem(kdb),1,MTYPER,idproc(mp,np),mpicomm,mpierr)
    else
      ! another tile.
      call mpi_bcast(elem(kdb),1,MTYPER, idproc(mp,np),mpicomm,mpierr)
    end if
    aelem = elem(kdb)
  end subroutine xceget

  !-----------------------------------------------------------------------
  subroutine xceput(aelem, a, ia,ja)
    integer, intent(in)    :: ia,ja
    real,    intent(in)    :: aelem
    real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

    !-----------
    !  1) fill a single element in the non-tiled 2-D grid.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aelem           real           input     element value
    !    a               real           in/out    target array
    !    ia              integer        input     1st index into a
    !    ja              integer        input     2nd index into a
    !-----------

    integer :: mp,np

    if (i0 < ia .and. ia <= i0+ii .and. j0 < ja .and. ja <= j0+jj) then
      ! this tile.
      a(ia-i0,ja-j0) = aelem
    end if
  end subroutine xceput

  !-----------------------------------------------------------------------
  subroutine xcgetc(iline)
    integer, intent(inout) :: iline(81)

    !-----------
    !  1) machine specific routine for broadcasting iline.
    !  2) only use in zagetc (hence the name).
    !-----------

    ! broadcast to all processors
    call mpi_bcast(iline,81,mpi_integer,idproc1(1),mpicomm,mpierr)
  end subroutine xcgetc

  !-----------------------------------------------------------------------
  subroutine xchalt(cerror)
    character*(*), intent(in) :: cerror ! error message
    !-----------
    ! stop all processes.
    ! only one processes need call this routine, i.e. it is for
    !  emergency stops.  use 'xcstop' for ordinary stops called
    !  by all processes.
    !-----------

    ! message passing version.
    if (cerror /= ' ') then
      write(lp,*) '**************************************************'
      write(lp,*) cerror
    end if
    write(lp,*) '**************************************************'
    write(lp,*) 'XCHALT CALLED ON PROC = ',mnproc,mproc,nproc
    write(lp,*) '**************************************************'

    if (mpicomm == mpicom_external) then
      call external_abort(cerror)
    else
      call mpi_abort(mpicomm,9,mpierr)
    end if
    error stop '(xchalt)'

  end subroutine xchalt

  !-----------------------------------------------------------------------
  subroutine xclget(aline,nl, a, i1,j1,iinc,jinc, mnflg)
    integer, intent(in)    ::  nl,i1,j1,iinc,jinc,mnflg
    real,    intent(out)   ::  aline(nl)
    real,    intent(in)    ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

    !-----------
    !  1) extract a line of elements from the non-tiled 2-D grid.
    !  2) aline(i) == aa(i1+iinc*(i-1),j1+jinc*(i-1)), for i=1...nl.
    !     where aa is the non-tiled representation of a, and
    !     iinc and jinc can each be -1, 0, or +1.

    !  3) mnflg selects which nodes must return the line
    !        = 0; all nodes
    !        = n; node number n (mnproc=n)

    !  4) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aline           real           output    required line of elements
    !    nl              integer        input     dimension of aline
    !    a               real           input     source array
    !    i1              integer        input     1st index into a
    !    j1              integer        input     2nd index into a
    !    iinc            integer        input     1st index increment
    !    jinc            integer        input     2nd index increment
    !    mnflg           integer        input     node return flag

    !  5) the global variable vland is returned when outside active tiles.
    !-----------

    ! pad al to guard against false sharing.
    ! double buffer to reduce the number of required barriers.
    real :: al(-47:itdm+jtdm+48,0:1)
    integer, save :: kdb = 0
    real :: dummy
    integer :: i1n,iif,iil,j1n,jjf,jjl,l,lx0,nxl,mp,np

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 3)
        nxc = 3
      end if
    end if

    kdb = mod(kdb+1,2)  ! the least recently used of the two buffers

    if ((jinc == 0 .and. iinc == 1) .or. nl == 1) then
      !   horizontal forward line.
      al(1:nl,kdb) = vland
      np  = npe_j(j1)
      do mp= mpe_1(np),mpe_e(np)
        iif = max(i1,     i0_pe(mp,np)+1)
        iil = min(i1+nl-1,i0_pe(mp,np)+ii_pe(mp,np))
        lx0 = iif - i1
        nxl = iil - iif + 1
        if (nxl <= 0) then
          cycle  ! no elements from tile (mp,np)
        end if

        if (mp == mproc .and. np == nproc) then
          ! this tile.
          do l= lx0+1,lx0+nxl
            al(l,kdb) = a(i1+l-1-i0,j1-j0)
          end do
          if (mnflg == 0) then
            call mpi_bcast(al(lx0+1,kdb),nxl,MTYPER,idproc(mp,np),mpicomm,mpierr)
          else if (mnflg /= mnproc) then
            call MPI_SEND(al(lx0+1,kdb),nxl,MTYPER,idproc1(mnflg), 9931,mpicomm, mpierr)
          end if
        else
          ! another tile.
          if (mnflg == 0) then
            call mpi_bcast(al(lx0+1,kdb),nxl,MTYPER,idproc(mp,np),mpicomm,mpierr)
          else if (mnflg == mnproc) then
            call MPI_RECV(al(lx0+1,kdb),nxl,MTYPER,idproc(mp,np), 9931,mpicomm, mpistat, mpierr)
          end if
        end if

        if (lx0+nxl == nl) then
          exit
        end if
      end do  ! np = 1,jpr

      if (mnflg == 0 .or. mnflg == mnproc) then
        aline(1:nl) = al(1:nl,kdb)
      end if

    else if (iinc == 0 .and. jinc == 1) then

      !   vertical forward line.
      al(1:nl,kdb) = vland

      do np= 1,jpr
        mp = mpe_i(i1,np)
        if (mp <= 0) then
          cycle  ! an all-land column
        end if
        jjf = max(j1,     j0_pe(mp,np)+1)
        jjl = min(j1+nl-1,j0_pe(mp,np)+jj_pe(mp,np))
        lx0 = jjf - j1
        nxl = jjl - jjf + 1
        if (nxl <= 0) then
          cycle  ! no elements from tile (mp,np)
        end if

        if (mp == mproc .and. np == nproc) then
          ! this tile.
          do l= lx0+1,lx0+nxl
            al(l,kdb) = a(i1-i0,j1+l-1-j0)
          end do
          if (mnflg == 0) then
            call mpi_bcast(al(lx0+1,kdb),nxl,MTYPER,idproc(mp,np),mpicomm,mpierr)
          else if (mnflg /= mnproc) then
            call MPI_SEND(al(lx0+1,kdb),nxl,MTYPER,idproc1(mnflg), 9931,mpicomm, mpierr)
          end if
        else
          ! another tile.
          if (mnflg == 0) then
            call mpi_bcast(al(lx0+1,kdb),nxl,MTYPER,idproc(mp,np),mpicomm,mpierr)
          else if (mnflg == mnproc) then
            call MPI_RECV(al(lx0+1,kdb),nxl,MTYPER,idproc(mp,np), 9931,mpicomm, mpistat, mpierr)
          end if
        end if
        if (lx0+nxl == nl) then
          exit
        end if
      end do  ! np = 1,jpr
      if (mnflg == 0 .or. mnflg == mnproc) then
        aline(1:nl) = al(1:nl,kdb)
      end if

    else

      ! diagonal and reversing lines - repeatedly call xceget.
      ! this always works, but is very slow.
      do l= 1,nl
        if (mnflg == 0 .or. mnflg == mnproc) then
          call xceget(aline(l), a, i1+iinc*(l-1),j1+jinc*(l-1))
        else
          call xceget(dummy,    a, i1+iinc*(l-1),j1+jinc*(l-1))
        end if
      end do
    end if

    if (use_TIMER) then
      if (nxc ==  3) then
        call xctmr1( 3)
        nxc = 0
      end if
    end if
  end subroutine xclget

  !-----------------------------------------------------------------------
  subroutine xclput(aline,nl, a, i1,j1,iinc,jinc)
    integer, intent(in)    ::  nl,i1,j1,iinc,jinc
    real,    intent(in)    ::  aline(nl)
    real,    intent(inout) ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

    !-----------
    !  1) fill a line of elements in the non-tiled 2-D grid.
    !  2) aline(i) == aa(i1+i1*(i-1),j1+j1*(i-1)), for i=1...nl.
    !     where aa is the non-tiled representation of a, and
    !     one of iinc and jinc must be 0, and the other must be 1.
    !     also updates the halo.

    !  3) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aline           real           input     line of element values
    !    nl              integer        input     dimension of aline
    !    a               real           in/out    target array
    !    i1              integer        input     1st index into a
    !    j1              integer        input     2nd index into a
    !    iinc            integer        input     1st index increment
    !    jinc            integer        input     2nd index increment
    !-----------

    integer :: i,j
    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 5)
        nxc = 5
      end if
    end if

    if (jinc == 0) then
      if (j1-j0 >= 1-nbdy .and. j1-j0 <= jj+nbdy) then
        do i= max(1-nbdy,i1-i0),min(i1-i0+nl-1,ii+nbdy)
          a(i,j1-j0) = aline(i+i0-i1+1)
        end do
        if (nreg /= 0 .and. &
             i0 == 0 .and. i1+nl-1 >= itdm-nbdy+1) then  ! periodic
          do i= max(itdm-nbdy+1,i1),i1+nl-1
            a(i-itdm,j1-j0) = aline(i)
          end do
        end if
        if (nreg /= 0 .and. &
             i0+ii == itdm .and. i1 <= nbdy) then        ! periodic
          do i= i1,min(nbdy,i1+nl-1)
            a(ii+i,j1-j0) = aline(i)
          end do
        end if
      end if
    else if (iinc == 0) then
      if (i1-i0 >= 1-nbdy .and. i1-i0 <= ii+nbdy) then
        do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
          a(i1-i0,j) = aline(j+j0-j1+1)
        end do
      end if
      if (nreg /= 0 .and. &
           i0 == 0 .and. i1 >= itdm-nbdy+1) then       ! periodic
        do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
          a(i1-itdm,j) = aline(j+j0-j1+1)
        end do
      end if
      if (nreg /= 0 .and. &
           i0+ii == itdm .and. i1 <= nbdy) then        ! periodic
        do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
          a(ii+i1,j) = aline(j+j0-j1+1)
        end do
      end if
    end if

    if (use_TIMER) then
      if (nxc ==  5) then
        call xctmr1( 5)
        nxc = 0
      end if
    end if
  end subroutine xclput

  !-----------------------------------------------------------------------
  subroutine xclput4(aline,nl, a, i1,j1,iinc,jinc)

    integer, intent(in)    ::  nl,i1,j1,iinc,jinc
    real*4,  intent(in)    ::  aline(nl)
    real*4,  intent(inout) ::  a(ii,jj)

    !-----------
    !  1) fill a line of elements in the non-tiled 2-D grid.
    !     Special version for xcaput4 only.
    !  2) aline(i) == aa(i1+i1*(i-1),j1+j1*(i-1)), for i=1...nl.
    !     where aa is the non-tiled representation of a, and
    !     one of iinc and jinc must be 0, and the other must be 1.
    !     also updates the halo.

    !  3) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aline           real           input     line of element values
    !    nl              integer        input     dimension of aline
    !    a               real           in/out    target array
    !    i1              integer        input     1st index into a
    !    j1              integer        input     2nd index into a
    !    iinc            integer        input     1st index increment
    !    jinc            integer        input     2nd index increment
    !-----------

    integer :: i,j

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 5)
        nxc = 5
      end if
    end if

    if (jinc == 0) then
      if (j1-j0 >= 1 .and. j1-j0 <= jj) then
        do i= max(1,i1-i0),min(i1-i0+nl-1,ii)
          a(i,j1-j0) = aline(i+i0-i1+1)
        end do
        if (nreg /= 0 .and. &
             i0 == 0 .and. i1+nl-1 >= itdm+1) then  ! periodic
          do i= max(itdm+1,i1),i1+nl-1
            a(i-itdm,j1-j0) = aline(i)
          end do
        end if
      end if
    else if (iinc == 0) then
      if (i1-i0 >= 1 .and. i1-i0 <= ii) then
        do j= max(1,j1-j0),min(j1-j0+nl-1,jj)
          a(i1-i0,j) = aline(j+j0-j1+1)
        end do
      end if
      if (nreg /= 0 .and. &
           i0 == 0 .and. i1 >= itdm+1) then       ! periodic
        do j= max(1,j1-j0),min(j1-j0+nl-1,jj)
          a(i1-itdm,j) = aline(j+j0-j1+1)
        end do
      end if
    end if

    if (use_TIMER) then
      if (nxc ==  5) then
        call xctmr1( 5)
        nxc = 0
      end if
    end if
  end subroutine xclput4

  !-----------------------------------------------------------------------
  subroutine xcbcst_i0(ia)
    integer, intent(inout) :: ia ! target variable
    !-----------
    !  1) broadcast integer scalar from tile mnproc = 1 to all tiles
    !  2) parameters:
    !-----------
    integer :: ia1(1)

    ia1(1) = ia
    call xcbcst_i1(ia1)
    ia = ia1(1)
  end subroutine xcbcst_i0

  !-----------------------------------------------------------------------
  subroutine xcbcst_i1(ia)
    integer, intent(inout) :: ia(:) ! target array

    !-----------
    !  1) broadcast integer array from tile mnproc = 1 to all tiles
    !-----------

    integer, parameter :: nmax=1024
    integer :: ib(nmax)
    integer :: i,is0,isl,n,nn

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 8)
        nxc = 8
      end if
    end if

    ! stripmine ia.
    n = size(ia)
    do is0= 0,n-1,nmax
      isl = min(is0+nmax,n)
      nn = isl - is0
      do i= 1,nn
        ib(i) = ia(is0+i)
      end do
      call mpi_bcast(ib,nn,mpi_integer,idproc1(1),mpicomm,mpierr)
      do i= 1,nn
        ia(is0+i) = ib(i)
      end do
    end do  ! stripmine loop

    if (use_TIMER) then
      if (nxc == 8) then
        call xctmr1( 8)
        nxc = 0
      end if
    end if
  end subroutine xcbcst_i1

  !-----------------------------------------------------------------------
  subroutine xcbcst_r0(ra)

    real, intent(inout) :: ra

    !-----------
    !  1) broadcast real scalar from tile mnproc = 1 to all tiles
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ra              real           in/out    target variable
    !-----------

    real :: ra1(1)

    ra1(1) = ra
    call xcbcst_r1(ra1)
    ra = ra1(1)
  end subroutine xcbcst_r0

  !-----------------------------------------------------------------------
  subroutine xcbcst_r1(ra)

    real, intent(inout) :: ra(:)

    !-----------
    !  1) broadcast real array from tile mnproc = 1 to all tiles
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ra              real           in/out    target array
    !-----------

    integer, parameter :: nmax=1024
    real :: rb(nmax)
    integer :: i,is0,isl,n,nn

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 8)
        nxc = 8
      end if
    end if

    ! stripmine ra.
    n = size(ra)
    do is0= 0,n-1,nmax
      isl = min(is0+nmax,n)
      nn = isl - is0
      do i= 1,nn
        rb(i) = ra(is0+i)
      end do
      call mpi_bcast(rb,nn,MTYPED,idproc1(1),mpicomm,mpierr)
      do i= 1,nn
        ra(is0+i) = rb(i)
      end do
    end do  ! stripmine loop

    if (use_TIMER) then
      if (nxc == 8) then
        call xctmr1( 8)
        nxc = 0
      end if
    end if
  end subroutine xcbcst_r1

  !-----------------------------------------------------------------------
  subroutine xcbcst_l0(la)
    logical, intent(inout) :: la
    !-----------
    !  1) broadcast logical from tile mnproc = 1 to all tiles
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    la              logical        in/out    target variable
    !-----------

    logical :: la1(1)

    la1(1) = la
    call xcbcst_l1(la1)
    la = la1(1)
  end subroutine xcbcst_l0

  !-----------------------------------------------------------------------
  subroutine xcbcst_l1(la)

    logical, intent(inout) :: la(:)

    !-----------
    !  1) broadcast logical array from tile mnproc = 1 to all tiles
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    la              logical        in/out    target array
    !-----------

    integer, parameter :: nmax=1024
    logical :: lb(nmax)
    integer :: i,is0,isl,n,nn

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 8)
        nxc = 8
      end if
    end if

    ! stripmine la.
    n = size(la)
    do is0= 0,n-1,nmax
      isl = min(is0+nmax,n)
      nn = isl - is0
      do i= 1,nn
        lb(i) = la(is0+i)
      end do
      call mpi_bcast(lb,nn,mpi_logical,idproc1(1),mpicomm,mpierr)
      do i= 1,nn
        la(is0+i) = lb(i)
      end do
    end do  ! stripmine loop

    if (use_TIMER) then
      if (nxc == 8) then
        call xctmr1( 8)
        nxc = 0
      end if
    end if
  end subroutine xcbcst_l1

  !-----------------------------------------------------------------------
  subroutine xcbcst_c(ca)
    character*(*), intent(inout) :: ca
    !-----------
    !  1) broadcast character string from tile mnproc = 1 to all tiles
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ca              character      in/out    target array
    !-----------

    integer, parameter :: nmax=1024
    character :: cb*(nmax)
    integer :: i,is0,isl,n,nn

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 8)
        nxc = 8
      end if
    end if

    ! stripmine ca.
    n = len(ca)
    do is0= 0,n-1,nmax
      isl = min(is0+nmax,n)
      nn = isl - is0
      cb(1:nn) = ca(is0+1:is0+nn)
      call mpi_bcast(cb,nn,mpi_character,idproc1(1),mpicomm,mpierr)
      ca(is0+1:is0+nn) = cb(1:nn)
    end do  ! stripmine loop

    if (use_TIMER) then
      if (nxc == 8) then
        call xctmr1( 8)
        nxc = 0
      end if
    end if
  end subroutine xcbcst_c

  !-----------------------------------------------------------------------
  subroutine xcmax_i0(ia)

    integer, intent(inout) :: ia

    !-----------
    !  1) replace integer scalar a with its element-wise maximum over all tiles.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ia              integer        in/out    target variable
    !-----------

    integer :: ia1(1)

    ia1(1) = ia
    call xcmax_i1(ia1)
    ia = ia1(1)
  end subroutine xcmax_i0

  !-----------------------------------------------------------------------
  subroutine xcmax_i1(ia)

    integer, intent(inout) :: ia(:)

    !-----------
    !  1) replace integer array a with its element-wise maximum over all tiles.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ia              integer        in/out    target array
    !-----------

    integer, parameter :: nmax=1024
    integer :: ib(nmax),ic(nmax)
    integer :: i,is0,isl,n,nn

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 9)
        nxc = 9
      end if
    end if

    !     stripmine ia.
    n = size(ia)
    do is0= 0,n-1,nmax
      isl = min(is0+nmax,n)
      nn = isl - is0
      do i= 1,nn
        ib(i) = ia(is0+i)
      end do
      call mpi_allreduce(ib,ic,nn,mpi_integer,mpi_max,mpicomm,mpierr)
      do i= 1,nn
        ia(is0+i) = ic(i)
      end do
    end do  ! stripmine loop

    if (use_TIMER) then
      if (nxc == 9) then
        call xctmr1( 9)
        nxc = 0
      end if
    end if
  end subroutine xcmax_i1

  !-----------------------------------------------------------------------
  subroutine xcmax_r0(ra)

    real, intent(inout) :: ra

    !-----------
    !  1) replace scalar a with its element-wise maximum over all tiles.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ra              real           in/out    target variable
    !-----------

    real :: ra1(1)

    ra1(1) = ra
    call xcmax_r1(ra1)
    ra = ra1(1)
  end subroutine xcmax_r0

  !-----------------------------------------------------------------------
  subroutine xcmax_r1(ra)
    real, intent(inout) :: ra(:)

    !-----------
    !  1) replace array a with its element-wise maximum over all tiles.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ra              real           in/out    target array
    !-----------

    integer, parameter :: nmax=1024
    real :: rb(nmax),rc(nmax)
    integer :: i,is0,isl,n,nn

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 9)
        nxc = 9
      end if
    end if

    !     stripmine ra.
    n = size(ra)

    do is0= 0,n-1,nmax
      isl = min(is0+nmax,n)
      nn = isl - is0
      do i= 1,nn
        rb(i) = ra(is0+i)
      end do
      call mpi_allreduce(rb,rc,nn,MTYPER,mpi_max,mpicomm,mpierr)
      do i= 1,nn
        ra(is0+i) = rc(i)
      end do
    end do  ! stripmine loop

    if (use_TIMER) then
      if (nxc == 9) then
        call xctmr1( 9)
        nxc = 0
      end if
    end if
  end subroutine xcmax_r1

  !-----------------------------------------------------------------------
  subroutine xcmin_i0(ia)
    integer, intent(inout) :: ia

    !-----------
    !  1) replace integer scalar a with its element-wise minimum over all tiles.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ia              integer        in/out    target variable
    !-----------

    integer :: ia1(1)

    ia1(1) = ia
    call xcmin_i1(ia1)
    ia = ia1(1)
  end subroutine xcmin_i0

  !-----------------------------------------------------------------------
  subroutine xcmin_i1(ia)
    integer, intent(inout) :: ia(:)

    !-----------
    !  1) replace integer array a with its element-wise minimum over all tiles.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ia              integer        in/out    target array
    !-----------

    integer, parameter :: nmax=1024
    integer :: ib(nmax),ic(nmax)
    integer :: i,is0,isl,n,nn

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0(10)
        nxc = 10
      end if
    end if

    !     stripmine ia.
    n = size(ia)
    do is0= 0,n-1,nmax
      isl = min(is0+nmax,n)
      nn = isl - is0
      do i= 1,nn
        ib(i) = ia(is0+i)
      end do
      call mpi_allreduce(ib,ic,nn,mpi_integer,mpi_min,mpicomm,mpierr)
      do i= 1,nn
        ia(is0+i) = ic(i)
      end do
    end do

    if (use_TIMER) then
      if (nxc == 10) then
        call xctmr1(10)
        nxc = 0
      end if
    end if
  end subroutine xcmin_i1

  !-----------------------------------------------------------------------
  subroutine xcmin_r0(ra)
    real, intent(inout) :: ra
    !-----------
    !  1) replace real scalar a with its element-wise minimum over all tiles.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ra              real           in/out    target variable
    !-----------

    real :: ra1(1)

    ra1(1) = ra
    call xcmin_r1(ra1)
    ra = ra1(1)
  end subroutine xcmin_r0

  !-----------------------------------------------------------------------
  subroutine xcmin_r1(ra)
    real, intent(inout) :: ra(:)


    !-----------
    !  1) replace real array a with its element-wise minimum over all tiles.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    ra              real           in/out    target array
    !-----------

    integer, parameter :: nmax=1024
    real    :: rb(nmax),rc(nmax)
    integer :: i,is0,isl,mn,n,nn

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0(10)
        nxc = 10
      end if
    end if

    !     stripmine ra.
    n = size(ra)

    do is0= 0,n-1,nmax
      isl = min(is0+nmax,n)
      nn = isl - is0
      do i= 1,nn
        rb(i) = ra(is0+i)
      end do
      call mpi_allreduce(rb,rc,nn,MTYPER,mpi_min,mpicomm,mpierr)
      do i= 1,nn
        ra(is0+i) = rc(i)
      end do
    end do

    if (use_TIMER) then
      if (nxc == 10) then
        call xctmr1(10)
        nxc = 0
      end if
    end if
  end subroutine xcmin_r1

  !-----------------------------------------------------------------------
  subroutine xcspmd

    !-----------
    ! 1) initialize data structures that identify the tiles.
    ! 2) data structures:
    ! ipr     - 1st 2-D node dimension (<=iqr)
    ! jpr     - 2nd 2-D node dimension (<=jqr)
    ! ijpr    -     1-D node dimension (ipr*jpr)
    ! mproc   - 1st 2-D node index
    ! nproc   - 2nd 2-D node index
    ! mnproc  -     1-D node index
    ! i0      -     i0_pe(mproc,nproc)
    ! ii      -     ii_pe(mproc,nproc)
    ! j0      -     j0_pe(mproc,nproc)
    ! jj      -     jj_pe(mproc,nproc)
    ! margin  -     how much of the halo is currently valid
    ! nreg    -     region type
    ! vland   -     fill value for land (standard value 0.0)
    ! idproc  -     2-D node addresses, with periodic wrap
    ! idproc1 -     1-D node addresses, with periodic wrap
    ! idhalo  -     left and right halo target nodes
    ! i0_pe   - 1st dimension tile offsets
    ! ii_pe   - 1st dimension tile extents (<=idm)
    ! j0_pe   - 2nd dimension tile offsets
    ! jj_pe   - 2nd dimension tile extents (<=jdm)
    ! mpe_1   - 1st node in each row of 2-D nodes
    ! mpe_e   - end node in each row of 2-D nodes
    ! mpe_i   - mapping from 1st global dimension to 2-D nodes
    ! npe_j   - mapping from 2nd global dimension to 2-D nodes
    ! i1sum   - local index of 1st partial sum on each tile
    ! iisum   - number of partial sums on each tile
    ! m0_top  - tile offset:       top neighbors (0:jpr-1)
    ! mm_top  - tile extent:       top neighbors (<=jpr)
    ! i0_st   - halo offsets: send top neighbors
    ! ii_st   - halo lengths: send top neighbors
    ! i0_gt   - halo offsets:  get top neighbors
    ! ii_gt   - halo lengths:  get top neighbors
    ! m0_bot  - tile offset:       bot neighbors (0:jpr-1)
    ! mm_bot  - tile extent:       bot neighbors (<=jpr)
    ! i0_sb   - halo offsets: send bot neighbors
    ! ii_sb   - halo lengths: send bot neighbors
    ! i0_gb   - halo offsets:  get bot neighbors
    ! ii_gb   - halo lengths:  get bot neighbors

    ! 3) all data structures are based on the processor number and
    ! the patch distribution file, 'patch.input'.
    !-----------

    integer :: i,j,l,m,mm,mn,mypei,n,npesi
    logical :: linit

    ! standard mpi (message passing) version.

    call mpi_initialized(linit,mpierr)
    if (linit) then
      mpicomm = mpicom_external
    else
      call mpi_init(mpierr)
      mpicomm = mpi_comm_world
    end if

    call mpi_comm_rank(mpicomm, mypei, mpierr)
    call mpi_comm_size(mpicomm, npesi, mpierr)

    mnproc = mypei + 1  ! mnproc counts from 1

    call setlogunit

    if (use_DEBUG_ALL) then
       write(lp,'(a,i5)') 'mnproc =',mnproc
       call xcsync(flush_lp)
    end if

    vland = 0.0

    ! read in the tile locations and sizes.
    ! patch distibution file on unit 21 (fort.21).

    ! here is an example of a patch file, with a leading "!" added:

    ! npes   npe   mpe   idm   jdm  ibig  jbig  nreg
    !    16     4     4    57    52    19    15     0

    !ispt(  1) =     22    35    42    49
    !iipe(  1) =     13     7     7     7
    !ispt(  2) =      1    15    25    34
    !iipe(  2) =     14    10     9     9
    !ispt(  3) =      2    21    29    38
    !iipe(  3) =     19     8     9    10
    !ispt(  4) =     18    28    35    42
    !iipe(  4) =     10     7     7     9

    !jspt(  1) =      1    15    26    38
    !jjpe(  1) =     14    11    12    15

    ! ispt(j) is the 1st i-point on each tile in the j-th row
    ! iipe(j) is the i-extent    of each tile in the j-th row
    ! jspt(1) is the 1st j-point on each tile in all columns
    ! jjpe(1) is the j-extent    of each tile in all columns

    ! note that each tile row can have a different tile layout,
    ! but each tile column is identical.  So a tile can have at
    ! most one nieghbor in the E and W directions, but several
    ! nieghbors in the N and S directions.

    ! iipe can be zero, indicating an empty tile (not included in the
    ! active tile count, npes) and therefore at least a halo widths
    ! land separation between active tiles to its east and west.  In
    ! periodic cases, the last non-empty tile can periodic wrap to the
    ! first tile in the row (i.e. trailing "empty" tiles can be null
    ! tiles, rather than all-land tiles).

    ijpr = ijqr
    ipr = iqr
    jpr = jqr

    if (use_ARCTIC) then
       if (nreg /= 2) then  ! not arctic
          if (mnproc == 1) then
             write(lp,'(a,i5)') 'input: nreg =',nreg
          end if
          call xcstop('xcspmd: patch.input must be for arctic')
          stop '(xcspmd)'
       end if
    else
       if (nreg == 2) then  ! arctic
          if (mnproc == 1) then
             write(lp,'(a,i5)') 'input: nreg =',nreg
          end if
          call xcstop('xcspmd: patch.input for arctic but ARCTIC undef.')
          stop '(xcspmd)'
       else if (nreg < 0) then  ! nreg = -1 reserved for type = one/omp
          if (mnproc == 1) then
             write(lp,'(a,i5)') 'input: nreg =',nreg
          end if
          call xcstop('xcspmd: patch.input for wrong nreg')
          stop '(xcspmd)'
       end if
    end if

    ! individual tile rows.
    do n= 1,jpr
      if (maxval(ii_pe(1:ipr,n)) <= 0) then
        call xcstop('xcspmd: decomposition has an empty row')
        stop '(xcspmd)'
      end if
    end do

    if (use_ARCTIC) then
       ! all arctic patch tiles must be the same size or empty,
       ! and empty tiles must be "twinned" across the top boundary.

       if (ipr > 1) then
          do m= 1,ipr
             if (ii_pe(m,jpr) == 0) then
                if (ii_pe(ipr+1-m,jpr) /= 0) then
                   if (mnproc == 1) then
                      write(lp,'(a,i3,a,i3,a)') &
                           'error - tile',m,',jpr is empty but tile', &
                           ipr+1-m,',jpr is not'
                   end if
                   call xcstop('xcspmd: arctic empty tiles are not twins')
                   stop '(xcspmd)'
                end if
             else if (ii_pe(m,jpr) /= itdm/ipr) then
                if (mnproc == 1) then
                   write(lp,'(a,i5)') &
                        'error - arctic patch tiles should have ii =',itdm/ipr
                end if
                call xcstop('xcspmd: arctic tiles are not the right size')
                stop '(xcspmd)'
             end if
          end do !m
       end if
    end if

    ! the generic tile column (must cover entire column).

    if (j0_pe(1,1) /= 0) then
      call xcstop('xcspmd: decomposition has wrong j0_pe')
      stop '(xcspmd)'
    end if
    do n= 2,jpr
      if (j0_pe(1,n) /= j0_pe(1,n-1)+jj_pe(1,n-1)) then
        call xcstop('xcspmd: decomposition is non-contiguous')
        stop '(xcspmd)'
      end if
    end do
    if (j0_pe(1,jpr)+jj_pe(1,jpr) /= jtdm) then
      call xcstop('xcspmd: decomposition has wrong jj_pe')
      stop '(xcspmd)'
    end if

    ! do we have the right number of pes?

    if (npesi /= ijpr) then
      if (mnproc == 1) then
        write(lp,*)
        write(lp,*) '***** ERROR - WRONG MPI SIZE *****'
        write(lp,*)
        write(lp,*) 'NPES    = ',npesi
        write(lp,*) '   IJPR = ',    ijpr
        write(lp,*) 'IPR,JPR = ',ipr,jpr
        write(lp,*)
      end if
      call xcstop('Error in xcspmd')
      stop
    end if

    ! mpi messages are sent and received by pe number (0:ijpr-1).

    null_tile = mpi_proc_null

    mn = 0
    do n= 1,jpr
      mpe_1(n) = 0
      do m= 1,ipr
        if (ii_pe(m,n) == 0) then
          idproc(m,n)   = null_tile
        else
          idproc1(mn+1) = mn
          idproc(m,n)   = mn
          mn = mn + 1
          if (mnproc == mn) then
            mproc = m
            nproc = n
          end if
          mpe_e(n) = m
          if (mpe_1(n) == 0) then
            mpe_1(n) = m
          end if
        end if
      end do
    end do
    if (mn /= ijpr) then
      if (mnproc == 1) then
        write(lp,'(a,i5)') 'input: ijpr =',ijpr
        write(lp,'(a,i5)') 'calc:  ijpr =',mn
      end if
      call xcstop('xcspmd: wrong number of sea tiles')
      stop '(xcspmd)'
    end if

    if (nreg == 0.or.nreg == 4) then

      ! longitudinal tile dimension is closed (not periodic)

      do n= 1,jpr
        idproc(    0,n) = null_tile
        idproc(ipr+1,n) = null_tile
      end do
    else

      ! longitudinal tile dimension is potentially periodic.

      do n= 1,jpr
        idproc(    0,n) = null_tile
        idproc(ipr+1,n) = null_tile

        i = maxval((i0_pe(1:ipr,n)+ii_pe(1:ipr,n)))
        if (i0_pe(1,n) == 0 .and. i == itdm) then
          idproc(         0,n) = idproc(mpe_e(n),n)
          idproc(mpe_e(n)+1,n) = idproc(       1,n)
        end if
      end do
    end if

    if (use_ARCTIC) then
       ! must have ipr even or 1 for arctic boundary case.
       if (ipr > 1 .and. mod(ipr,2) /= 0) then
          call xcstop('Error in xcspmd (arctic) - ipr must be even')
          stop '(xcspmd)'
       end if

       ! latitudinal tile dimension is closed/arctic.

       do m= 1,ipr
          idproc(m,    0) = null_tile
          idproc(m,jpr+1) = idproc(ipr+1-m,jpr)  !arctic tile mapping
       end do
       idproc(    0,    0) = null_tile
       idproc(ipr+1,    0) = null_tile
       idproc(    0,jpr+1) = idproc(ipr,jpr+1)
       idproc(ipr+1,jpr+1) = idproc(1,  jpr+1)
    else
       if (nreg >= 3) then
          ! latitudinal tile dimension is periodic
          do m= 0,ipr+1
             idproc(m,    0) = idproc(m,jpr)
             idproc(m,jpr+1) = idproc(m,  1)
          end do
       else
          ! latitudinal tile dimension is closed
          do m= 0,ipr+1
             idproc(m,    0) = null_tile
             idproc(m,jpr+1) = null_tile
          end do
       end if
    end if  ! end if use_ARCTIC

    ! 1-d tiling logic is easier if assumed periodic.

    idproc1(     0) = idproc1(ijpr)
    idproc1(ijpr+1) = idproc1(   1)

    ! mapping from global i,j to mp,np.
    ! ia,ja is on tile mpe_i(ia,npe_j(ja)),npe_j(ja),
    ! or on no tile if mpe_i(ia,npe_j(ja)) is 0 or -1.

    do n= 1,jpr
      mpe_i(1:itdm,n) = 0  ! default is an empty tile
      do m= 1,ipr  ! i-map potentially varies with n
        if (ii_pe(m,n) > 0) then
          do i= i0_pe(m,n)+1,i0_pe(m,n)+ii_pe(m,n)
            mpe_i(i,n) = m
          end do
          if (m /= ipr) then
            if (ii_pe(m+1,n) > 0) then
              do i= i0_pe(m,n)+ii_pe(m,n)+1,i0_pe(m+1,n)
                mpe_i(i,n) = -1  ! gap between tiles
              end do
            end if
          end if
        end if
      end do
      m = 1  ! only one j-map
      do j= j0_pe(m,n)+1,j0_pe(m,n)+jj_pe(m,n)
        npe_j(j)   = n
      end do
    end do

    ! do each partial sum on the tile that owns its center point.
    ! i1sum - local index of 1st partial sum on each tile
    ! iisum - number of partial sums on each tile
    ! see xcsum for how i1sum and iisum are used.

    do n= 1,jpr
      do m= 1,ipr
        if (ii_pe(m,n) <= 0) then
          i1sum(m,n) =  0
          iisum(m,n) =  0
        else
          idhalo(1) = idproc(m-1,n)
          idhalo(2) = idproc(m+1,n)
          if (idhalo(1) /= null_tile .and. m /= 1) then
            if (i0_pe(m,n) /= i0_pe(m-1,n)+ii_pe(m-1,n)) then
              idhalo(1) = null_tile
            end if
          end if
          if (idhalo(2) /= null_tile .and. m /= mpe_e(n)) then
            if (i0_pe(m,n)+ii_pe(m,n) /= i0_pe(m+1,n)) then
              idhalo(2) = null_tile
            end if
          end if
          i1sum(m,n) = -99
          iisum(m,n) =  0
          do i= 1+nbdy,itdm+nbdy,2*nbdy+1
            if (i0_pe(m,n) < i .and. &
                 i <= i0_pe(m,n)+ii_pe(m,n)) then
              iisum(m,n) = iisum(m,n) + 1
              if (iisum(m,n) == 1) then
                i1sum(m,n) = i - nbdy - i0_pe(m,n)
              end if
            else if (idhalo(1) == null_tile .and. &
                 i > i0_pe(m,n)-nbdy   .and. &
                 i <= i0_pe(m,n)             ) then
              iisum(m,n) = iisum(m,n) + 1
              if (iisum(m,n) == 1) then
                i1sum(m,n) = i - nbdy - i0_pe(m,n)
              end if
            else if ((idhalo(2) == null_tile.or.m == mpe_e(n)) .and. &
                 i > i0_pe(m,n)+ii_pe(m,n)                .and. &
                 i <= i0_pe(m,n)+ii_pe(m,n)+nbdy      ) then
              iisum(m,n) = iisum(m,n) + 1
              if (iisum(m,n) == 1) then
                i1sum(m,n) = i - nbdy - i0_pe(m,n)
              end if
            end if
          end do
        end if
      end do
    end do

    ! local tile extents.

    i0 = i0_pe(mproc,nproc)
    ii = ii_pe(mproc,nproc)
    j0 = j0_pe(mproc,nproc)
    jj = jj_pe(mproc,nproc)

    margin = 0

    ! left and right halo targets

    idhalo(1) = idproc(mproc-1,nproc)
    idhalo(2) = idproc(mproc+1,nproc)

    if (idhalo(1) /= null_tile .and. mproc /= 1) then

      ! is the left tile touching this one?

      if (i0 /= i0_pe(mproc-1,nproc)+ii_pe(mproc-1,nproc)) then
        idhalo(1) = null_tile
      end if
    end if

    if (idhalo(2) /= null_tile .and. mproc /= mpe_e(nproc)) then

      ! is the right tile touching this one?

      if (i0+ii /= i0_pe(mproc+1,nproc)) then
        idhalo(2) = null_tile
      end if
    end if

    ! local halo exchange data structures

    ! m0_top - tile offset:       top neighbors
    ! mm_top - tile extent:       top neighbors (<=jpr)
    ! i0_st  - halo offsets: send top neighbors
    ! ii_st  - halo lengths: send top neighbors
    ! i0_gt  - halo offsets:  get top neighbors
    ! ii_gt  - halo lengths:  get top neighbors
    ! m0_bot - tile offset:       bot neighbors
    ! mm_bot - tile extent:       bot neighbors (<=jpr)
    ! i0_sb  - halo offsets: send bot neighbors
    ! ii_sb  - halo lengths: send bot neighbors
    ! i0_gb  - halo offsets:  get bot neighbors
    ! ii_gb  - halo lengths:  get bot neighbors

    ! note that send is also receive, and is w.r.t. the local  tile.
    ! similarly get  is also put,     and is w.r.t. the remote tile.

    if (nproc == jpr) then
       if (use_ARCTIC) then
          ! single, same size, top arctic nieghbor
          m0_top   = mproc - 1
          mm_top   =  1
          i0_st(1) =  0
          i0_gt(1) =  0
          ii_st(1) = ii
          ii_gt(1) = ii
       else
          m0_top = 0
          mm_top = 0
          if (nreg >= 3) then
             ! latitudinal tile dimension is open
             n = 1
             m = 0
             do i= 1,ii
                if (mpe_i(i0+i,n) /= m) then
                   if (mm_top == 0) then
                      m0_top = mpe_i(i0+i,n) - 1
                   else if (m /= -1) then
                      ii_st(mm_top) = i-1 - i0_st(mm_top)
                      ii_gt(mm_top) = ii_st(mm_top)
                   end if
                   m = mpe_i(i0+i,n)
                   if (m > 0) then
                      mm_top        = mm_top + 1
                      i0_st(mm_top) = i-1
                      i0_gt(mm_top) = i-1 + i0-i0_pe(m,n)
                   else if (m == 0) then
                      mm_top        = mm_top + 1
                      i0_st(mm_top) = i-1
                      i0_gt(mm_top) = i0_gt(mm_top-1) + ii_gt(mm_top-1)
                      ! elseif (m.eq.-1) then  !do nothing
                   end if
                end if
             end do
             if (mm_top > 0) then
                if (m > 0) then
                   ii_st(mm_top) = ii - i0_st(mm_top)
                   ii_gt(mm_top) = ii_st(mm_top)
                else if (m == 0) then
                   mm_top = mm_top-1
                   ! elseif (m.eq.-1) then  !do nothing
                end if
             end if
          end if
       end if ! end if use_ARCTIC
   else
      n = nproc + 1
      m0_top = 0
      mm_top = 0
      m      = 0
      do i= 1,ii
        if (mpe_i(i0+i,n) /= m) then
          if (mm_top == 0) then
            m0_top = mpe_i(i0+i,n) - 1
          else if (m /= -1) then
            ii_st(mm_top) = i-1 - i0_st(mm_top)
            ii_gt(mm_top) = ii_st(mm_top)
          end if
          m = mpe_i(i0+i,n)
          if (m > 0) then
            mm_top        = mm_top + 1
            i0_st(mm_top) = i-1
            i0_gt(mm_top) = i-1 + i0-i0_pe(m,n)
          else if (m == 0) then
            mm_top        = mm_top + 1
            i0_st(mm_top) = i-1
            i0_gt(mm_top) = i0_gt(mm_top-1) + ii_gt(mm_top-1)
            ! elseif (m.eq.-1) then  !do nothing
          end if
        end if
      end do
      if (mm_top > 0) then
        if (m > 0) then
          ii_st(mm_top) = ii - i0_st(mm_top)
          ii_gt(mm_top) = ii_st(mm_top)
        else if (m == 0) then
          mm_top = mm_top-1
          ! elseif (m.eq.-1) then  !do nothing
        end if
      end if
    end if !nproc == 1:else

    if (nproc == 1) then
      ! no bottom nieghbor (closed boundary)
      m0_bot = 0
      mm_bot = 0
      if (nreg >= 3) then
        ! latitudinal tile dimension is open
        n = jpr
        m = 0
        do i= 1,ii
          if (mpe_i(i0+i,n) /= m) then
            if (mm_bot == 0) then
              m0_bot = mpe_i(i0+i,n) - 1
            else if (m /= -1) then
              ii_sb(mm_bot) = i-1 - i0_sb(mm_bot)
              ii_gb(mm_bot) = ii_sb(mm_bot)
            end if
            m = mpe_i(i0+i,n)
            if (m > 0) then
              mm_bot        = mm_bot + 1
              i0_sb(mm_bot) = i-1
              i0_gb(mm_bot) = i-1 + i0-i0_pe(m,n)
            else if (m == 0) then
              mm_bot        = mm_bot + 1
              i0_sb(mm_bot) = i-1
              i0_gb(mm_bot) = i0_gb(mm_bot-1) + ii_gb(mm_bot-1)
              ! elseif (m.eq.-1) then  !do nothing
            end if
          end if
        end do
        if (mm_bot > 0) then
          if (m > 0) then
            ii_sb(mm_bot) = ii - i0_sb(mm_bot)
            ii_gb(mm_bot) = ii_sb(mm_bot)
          else if (m == 0) then
            mm_bot = mm_bot-1
            ! elseif (m.eq.-1) then  !do nothing
          end if
        end if
      end if
    else
      n = nproc - 1
      m0_bot = 0
      mm_bot = 0
      m      = 0
      do i= 1,ii
        if (mpe_i(i0+i,n) /= m) then
          if (mm_bot == 0) then
            m0_bot = mpe_i(i0+i,n) - 1
          else if (m /= -1) then
            ii_sb(mm_bot) = i-1 - i0_sb(mm_bot)
            ii_gb(mm_bot) = ii_sb(mm_bot)
          end if
          m = mpe_i(i0+i,n)
          if (m > 0) then
            mm_bot        = mm_bot + 1
            i0_sb(mm_bot) = i-1
            i0_gb(mm_bot) = i-1 + i0-i0_pe(m,n)
          else if (m == 0) then
            mm_bot        = mm_bot + 1
            i0_sb(mm_bot) = i-1
            i0_gb(mm_bot) = i0_gb(mm_bot-1) + ii_gb(mm_bot-1)
            ! elseif (m.eq.-1) then  !do nothing
          end if
        end if
      end do
      if (mm_bot > 0) then
        if (m > 0) then
          ii_sb(mm_bot) = ii - i0_sb(mm_bot)
          ii_gb(mm_bot) = ii_sb(mm_bot)
        else if (m == 0) then
          mm_bot = mm_bot-1
          ! elseif (m.eq.-1) then  !do nothing
        end if
      end if
    end if !nproc == 1:else

    ! printout the tile data structures.

    if (mnproc == 1) then
       write(lp,'(/a)')'mnproc mproc nproc     i0   ii     j0   jj  i1sum iisum'
       mn = 0
       do n= 1,jpr
          do m= 1,ipr
             if (ii_pe(m,n) /= 0) then
                mn= mn + 1
                write(lp,'(i6,2i6,i7,i5,i7,i5,i7,i6)') &
                     mn,m,n, &
                     i0_pe(m,n),ii_pe(m,n), &
                     j0_pe(m,n),jj_pe(m,n), &
                     i1sum(m,n),iisum(m,n)
             end if
          end do !m
       end do !n
       write(lp,*)
       if (use_ARCTIC) then
          write(lp,'(a)')'mnproc mproc nproc mnarct'
          mn = 0
          do n= 1,jpr
             do m= 1,ipr
                if (ii_pe(m,n) /= 0) then
                   mn= mn + 1
                   if (n == jpr) then
                      write(lp,'(i6,2i6,i7)') mn,m,n,idproc(m,n+1)
                   end if
                end if
             end do !m
          end do !n
          write(lp,*)
       end if
       if (use_DEBUG_ALL) then
          do n= 1,jpr
             write(lp,*) 'mpe_1,mpe_e = ',mpe_1(n),mpe_e(n)
          end do
          do n= 1,jpr
             write(lp,*) 'mpe_i = ',mpe_i(:,n)
          end do
          write(lp,*)
          write(lp,*) 'npe_j = ',npe_j(:)
          write(lp,*)
       end if
    end if ! if (mnproc == 1)
    call xcsync(flush_lp)

    if (use_DEBUG_ALL) then
       do n= 1,jpr
          do m= 1,ipr
             if (mproc == m .and. nproc == n) then
                write(lp,'(a,2i3,i4,i3,16(i5,i4))') &
                     'm,n,_top,_st = ', &
                     m,n,m0_top,mm_top,(i0_st(l),ii_st(l), l= 1,mm_top)
                if (m == ipr) then
                   write(lp,*) !blank line
                end if
             end if
             call xcsync(flush_lp)
          end do !m
          do m= 1,ipr
             if (mproc == m .and. nproc == n) then
                write(lp,'(a,2i3,i4,i3,16(i5,i4))') &
                     'm,n,_bot,_sb = ', &
                     m,n,m0_bot,mm_bot,(i0_sb(l),ii_sb(l), l= 1,mm_bot)
                if (m == ipr) then
                   write(lp,*) !blank line
                end if
             end if
             call xcsync(flush_lp)
          end do !m
          do m= 1,ipr
             if (mproc == m .and. nproc == n) then
                write(lp,'(a,2i3,3i5)') &
                     'm,n,id,idhalo   = ', &
                     m,n,idproc(m,n),idhalo
                if (m == ipr) then
                   write(lp,*) !blank line
                end if
             end if
             call xcsync(flush_lp)
          end do !m
       end do !n
    end if

    ! initialize timers.
    call xctmri
    if (use_TIMER) then
      call xctmrn( 1,'xcaget')
      call xctmrn( 2,'xceget')
      call xctmrn( 3,'xclget')
      call xctmrn( 4,'xcaput')
      call xctmrn( 5,'xcXput')
      call xctmrn( 6,'xcsum ')
      call xctmrn( 8,'xcbcst')
      call xctmrn( 9,'xcmax ')
      call xctmrn(10,'xcmin ')
      call xctmrn(12,'xctilr')
    end if
  end subroutine xcspmd

  !-----------------------------------------------------------------------
  subroutine xcstop(cerror)
    character*(*), intent(in) :: cerror

    !-----------
    !  1) stop all processes.
    !  2) all processes must call this routine.
    !     use 'xchalt' for emergency stops.

    !  3) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    cerror          char*(*)       input     error message
    !-----------

    !     print active timers.
    call xctmrp

    !     message passing version, set barrier and then stop everything.
    !     use xcsync as the barrier so that stdout is flushed.
    !     note that the system will hang unless all processes call xcstop.

    call xcsync(flush_lp)
    if (mnproc == 1 .and. cerror /= ' ') then
      write(lp,*) '**************************************************'
      write(lp,*) cerror
      write(lp,*) '**************************************************'
      write(lp,*)
    end if
    call xcsync(flush_lp)

    if (mpicomm == mpicom_external) then
      call external_abort(cerror)
    else
      write(lp,*) 'mpi_finalize called on processor ',mnproc
      call xcsync(flush_lp)
      call mpi_finalize(mpierr)
    end if

    stop '(xcstop)'
  end subroutine xcstop

  !-----------------------------------------------------------------------
  subroutine xcsum(sum, a,mask)
    real(8),  intent(out)   :: sum
    real,    intent(inout) :: a(   1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
    integer, intent(in)    :: mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

    !-----------
    !  1) sum a 2-d array, where mask==1
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    sum             real*8         output    sum of a
    !    a               real           input     source array
    !    mask            integer        input     mask array

    !  3) sum is bit for bit reproducable for the same halo size, nbdy.
    !-----------

    real(8), parameter :: zero8 = 0.0
    integer, parameter :: mxsum = (idm+3*nbdy)/(2*nbdy+1)
    real(8) :: sum8t(mxsum*jdm),sum8j(jdm),sum8s
    real(8) :: sum8
    real    :: vsave
    integer :: i,i1,j,l,mp,np

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0( 6)
        nxc = 6
      end if
    end if

    ! halo update so that 2*nbdy+1 wide strips are on chip.

    vsave = vland
    vland = 0.0
    call xctilr(a,1,1, nbdy,0, halo_ps)
    vland = vsave

    ! row sums in 2*nbdy+1 wide strips.

    !$omp parallel do private(j,i1,i,l,sum8) &
    !$omp schedule(static,jblk)
    do j = 1,jj
      do l= 1,iisum(mproc,nproc)
        i1   = i1sum(mproc,nproc) + (l-1)*(2*nbdy+1)
        sum8 = zero8
        do i= max(i1,1-nbdy),min(i1+2*nbdy,ii+nbdy,itdm-i0)
          if (mask(i,j) == 1) then
            sum8 = sum8 + a(i,j)
          end if
        end do
        sum8t(l + (j-1)*iisum(mproc,nproc)) = sum8
      end do
    end do
    !$omp end parallel do

    ! complete row sums on first processor in each row.
    if (mproc == mpe_1(nproc)) then
      do j = 1,jj
        sum8j(j) = zero8
        do l= 1,iisum(mproc,nproc)
          sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mproc,nproc))
        end do
      end do

      ! remote sums.
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = iisum(mp,nproc)*jj
        if (l > 0) then
          call MPI_RECV(sum8t,l,MTYPED,idproc(mp,nproc), 9900,mpicomm, mpistat, mpierr)
          do j = 1,jj
            do l= 1,iisum(mp,nproc)
              sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mp,nproc))
            end do
          end do
        end if
      end do
    else
      l = iisum(mproc,nproc)*jj
      if (l > 0) then
        call MPI_SEND(sum8t,l,MTYPED,idproc(mpe_1(nproc),nproc), 9900,mpicomm, mpierr)
      end if
    end if

    ! sum of row sums, on first processor.
    if (mnproc == 1) then
      sum8 = zero8
      do j= 1,jj
        sum8 = sum8 + sum8j(j)
      end do
      do np= 2,jpr
        mp = mpe_1(np)
        call MPI_RECV(sum8j,jj_pe(mp,np),MTYPED,idproc(mp,np), 9901,mpicomm, mpistat, mpierr)
        do j= 1,jj_pe(mp,np)
          sum8 = sum8 + sum8j(j)
        end do
      end do
      sum8s = sum8
    else if (mproc == mpe_1(nproc)) then
      call MPI_SEND(sum8j,jj,MTYPED,idproc1(1), 9901,mpicomm, mpierr)
    end if

    ! broadcast result to all processors.
    call mpi_bcast(sum8s,1,MTYPED,idproc1(1),mpicomm,mpierr)
    sum = sum8s

    if (use_TIMER) then
      if (nxc ==  6) then
        call xctmr1( 6)
        nxc = 0
      end if
    end if
  end subroutine xcsum

  !-----------------------------------------------------------------------
  subroutine xcsync(lflush)
    logical, intent(in) :: lflush

    !-----------
    !  1) barrier, no processor exits until all arrive (and flush stdout).
    !  2) some MPI implementations only flush stdout as a collective
    !     operation, and hence the lflush=.true. option to flush stdout.
    !-----------

    if (lflush) then
      call mpi_barrier(mpicomm, mpierr)
    else
      call mpi_barrier(mpicomm, mpierr)
    end if
  end subroutine xcsync

  !-----------------------------------------------------------------------
  subroutine xctilr(a,l1,ld,mh,nh,itype)
     integer, intent(in)    :: l1,ld,mh,nh,itype
     real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)

     if (use_ARCTIC) then
        call xctilr_arctic(a,l1,ld,mh,nh,itype)
     else
        call xctilr_nonarctic(a,l1,ld,mh,nh,itype)
     end if
  end subroutine xctilr

  !-----------------------------------------------------------------------
  recursive subroutine xctilr_arctic(a,l1,ld,mh,nh,itype)
    integer, intent(in)    :: l1,ld,mh,nh,itype
    real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)

    !-----------
    !  1) update the tile overlap halo of a real array.

    !     this version of arctic bi-polar patch only
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    a               real           in/out    target array
    !    l1              integer        input     3rd dim. start index
    !    ld              integer        input     3rd dimension of a
    !    mh              integer        input     1st (EW) update halo size
    !    nh              integer        input     2nd (NS) update halo size
    !    itype           integer        input     grid and field type

    !  3) itype selects both the grid and field type
    !        itype= 1; p-grid, scalar field
    !        itype= 2; q-grid, scalar field
    !        itype= 3; u-grid, scalar field
    !        itype= 4; v-grid, scalar field
    !        itype=11; p-grid, vector field
    !        itype=12; q-grid, vector field
    !        itype=13; u-grid, vector field
    !        itype=14; v-grid, vector field

    !  4) the global variable vland is returned by halos over land.
    !-----------

    integer, parameter :: lsize = kdm
    integer, parameter :: ilen= idm*lsize*(nbdy+1)+64
    integer, parameter :: jlen=(jdm+2*nbdy)*lsize*nbdy+64

    ! halo buffer
    real, save :: ai(ilen,4),aj(jlen,4)
    real, save :: aia(lsize*(nbdy+1)+64,2)
    integer    :: i,io,j,k,l,lts,lbs,lbr,ltr,las,lar,lg0,ls0,lm,m,mhl,nhl
    real       :: sarc
    integer, save :: mpireqa(4*iqr),mpireqb(4) ! communication handles.
    integer, save :: ilold,ltsold,lbsold,lbrold,ltrold,lasold,larold,nreqa ! communication handles.
    data ilold,ltsold,lbsold,lbrold,ltrold,lasold,larold &
         / 0,0,0,0,0,0,0 /

    ! split large requests into smaller pieces

    if (ld-l1+1 > lsize) then
      do k= l1,ld,lsize
        l = min(k+lsize-1,ld)
        call xctilr(a,k,l,mh,nh,itype)
      end do
      return
    end if

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0(12)
        nxc = 12
      end if
    end if

    mhl = max(0,min(mh,nbdy))
    nhl = max(0,min(nh,nbdy))

    if (itype < 10) then
      sarc =  1.0
    else
      sarc = -1.0
    end if

    if (ipr == 1 .and. jpr == 1) then
      do k= l1,ld
        do j= 1,nhl
          do i= 1,ii
            a(i,1-j,k) = vland
          end do
        end do
        if (itype == 1 .or. itype == 11) then
          do j= 0,nhl
            do i= 1,ii
              io = ii-mod(i-1,ii)
              if (a(io,jj-1-j,k) /= vland) then
                a(i,jj+j,k) = sarc*a(io,jj-1-j,k)
              else
                a(i,jj+j,k) = vland
              end if
            end do !i
          end do !j
        else if (itype == 2 .or. itype == 12) then
          do i= ii/2+1,ii
            io = mod(ii-(i-1),ii)+1
            if (a(io,jj,k) /= vland) then
              a(i,jj,k) = sarc*a(io,jj,k)
            else
              a(i,jj,k) = vland
            end if
          end do !i
          do j= 1,nhl
            do i= 1,ii
              io = mod(ii-(i-1),ii)+1
              if (a(io,jj-j,k) /= vland) then
                a(i,jj+j,k) = sarc*a(io,jj-j,k)
              else
                a(i,jj+j,k) = vland
              end if
            end do !i
          end do !j
        else if (itype == 3 .or. itype == 13) then
          do j= 0,nhl
            do i= 1,ii
              io = mod(ii-(i-1),ii)+1
              if (a(io,jj-1-j,k) /= vland) then
                a(i,jj+j,k) = sarc*a(io,jj-1-j,k)
              else
                a(i,jj+j,k) = vland
              end if
            end do !i
          end do !j
        else if (itype == 4 .or. itype == 14) then
          do i= ii/2+1,ii
            io = ii-mod(i-1,ii)
            if (a(io,jj,k) /= vland) then
              a(i,jj,k) = sarc*a(io,jj,k)
            else
              a(i,jj,k) = vland
            end if
          end do !i
          do j= 1,nhl
            do i= 1,ii
              io = ii-mod(i-1,ii)
              if (a(io,jj-j,k) /= vland) then
                a(i,jj+j,k) = sarc*a(io,jj-j,k)
              else
                a(i,jj+j,k) = vland
              end if
            end do !i
          end do !j
        end if !itype
      end do !k
    else
      lts = 0
      lbs = 0
      lbr = 0
      ltr = 0
      las = 0
      lar = 0
      if (nproc /= jpr) then
        if (nhl > 0) then
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            do k= l1,ld
              do j= 1,nhl
                lts = lts + 1
                ai(lts,1) = a(i,jj+1-j,k)
                ai(lts,2) = a(i,     j,k)
                ai(lts,3) = vland
                ai(lts,4) = vland
              end do !j
            end do !k
          end do !i
          lbs = lts
          lbr = lts
          ltr = lts
        end if
      else if (itype == 1 .or. itype == 11) then !p-grid
        do i= 1,ii  ! outer loop to simplify multiple neighbor case
          io = ii+1-i !ii:1:-1
          do k= l1,ld
            do j= 0,nhl
              lts = lts + 1
              if (a(io,jj-1-j,k) /= vland) then
                ai(lts,1) = sarc*a(io,jj-1-j,k)
              else
                ai(lts,1) = vland
              end if
              ai(lts,4) = vland
            end do !j
            do j= 1,nhl
              lbs = lbs + 1
              ai(lbs,2) = a(i,j,k)
              ai(lbs,3) = vland
            end do !j
          end do !k
        end do !i
        lbr = lbs
        ltr = lts
      else if (itype == 3 .or. itype == 13) then  !u-grid
        do i= 1,ii  ! outer loop to simplify multiple neighbor case
          io = ii+2-i !ii+1:2:-1
          do k= l1,ld
            do j= 0,nhl
              lts = lts + 1
              if (a(io,jj-1-j,k) /= vland) then
                ai(lts,1) = sarc*a(io,jj-1-j,k)
              else
                ai(lts,1) = vland
              end if
              ai(lts,4) = vland
            end do !j
            do j= 1,nhl
              lbs = lbs + 1
              ai(lbs,2) = a(i,j,k)
              ai(lbs,3) = vland
            end do !j
          end do !k
        end do !i
        lbr = lbs
        ltr = lts
        do k= l1,ld
          do j= 0,nhl
            las = las + 1
            if (a(1,jj-1-j,k) /= vland) then
              aia(las,1) = sarc*a(1,jj-1-j,k)
            else
              aia(las,1) = vland
            end if
            aia(las,2) = vland
          end do !j
        end do !k
        lar = las
      else if (itype == 2 .or. itype == 12) then  !q-grid
        if (mproc > ipr/2) then
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            io = ii+2-i !ii+1:2:-1
            do k= l1,ld
              do j= 1,nhl
                lts = lts + 1
                if (a(io,jj-j,k) /= vland) then
                  ai(lts,1) = sarc*a(io,jj-j,k)
                else
                  ai(lts,1) = vland
                end if
                ai(lts,2) = a(i,j,k)
                ai(lts,3) = vland
              end do !j
              do j= 0,nhl
                ltr = ltr + 1
                ai(ltr,4) = vland
              end do !j
            end do !k
          end do !i
          lbs = lts
          lbr = lts
          do k= l1,ld
            do j= 0,nhl
              lar = lar + 1
              aia(lar,2) = vland
            end do !j
          end do !k
        else
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            io = ii+2-i !ii+1:2:-1
            do k= l1,ld
              do j= 0,nhl
                lts = lts + 1
                if (a(io,jj-j,k) /= vland) then
                  ai(lts,1) = sarc*a(io,jj-j,k)
                else
                  ai(lts,1) = vland
                end if
              end do !j
              do j= 1,nhl
                lbs = lbs + 1
                ai(lbs,2) = a(i,j,k)
                ai(lbs,3) = vland
                ai(lbs,4) = vland
              end do !j
            end do !k
          end do !i
          lbr = lbs
          ltr = lbs
          do k= l1,ld
            do j= 1,nhl
              lar = lar + 1
              aia(lar,2) = vland
            end do !j
          end do !k
        end if
        if (mod(ipr+1-mproc,ipr)+1 <= ipr/2) then
          do k= l1,ld
            do j= 1,nhl
              las = las + 1
              if (a(1,jj-j,k) /= vland) then
                aia(las,1) = sarc*a(1,jj-j,k)
              else
                aia(las,1) = vland
              end if
            end do !j
          end do !k
        else
          do k= l1,ld
            do j= 0,nhl
              las = las + 1
              if (a(1,jj-j,k) /= vland) then
                aia(las,1) = sarc*a(1,jj-j,k)
              else
                aia(las,1) = vland
              end if
            end do !j
          end do !k
        end if
      else  !v-grid
        if (mproc > ipr/2) then
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            io = ii+1-i  !ii:1:-1
            do k= l1,ld
              do j= 1,nhl
                lts = lts + 1
                if (a(io,jj-j,k) /= vland) then
                  ai(lts,1) = sarc*a(io,jj-j,k)
                else
                  ai(lts,1) = vland
                end if
                ai(lts,2) = a(i,j,k)
                ai(lts,3) = vland
              end do !j
              do j= 0,nhl
                ltr = ltr + 1
                ai(ltr,4) = vland
              end do !j
            end do !k
          end do !i
          lbs = lts
          lbr = lts
        else
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            io = ii+1-i  !ii:1:-1
            do k= l1,ld
              do j= 0,nhl
                lts = lts + 1
                if (a(io,jj-j,k) /= vland) then
                  ai(lts,1) = sarc*a(io,jj-j,k)
                else
                  ai(lts,1) = vland
                end if
              end do !j
              do j= 1,nhl
                lbs = lbs + 1
                ai(lbs,2) = a(i,j,k)
                ai(lbs,3) = vland
                ai(lbs,4) = vland
              end do !j
            end do !k
          end do !i
          lbr = lbs
          ltr = lbs
        end if
      end if !itype

      if (lts /= ltsold .or. lbs /= lbsold .or. lbr /= lbrold .or. &
           ltr /= ltrold .or. las /= lasold .or. lar /= larold) then
        if (ltsold + lbsold + lbrold + &
             ltrold + lasold + larold /= 0) then
          do i= 1,nreqa
            call mpi_request_free(mpireqa(i), mpierr)
          end do
        end if
        ltsold = lts
        lbsold = lbs
        lbrold = lbr
        ltrold = ltr
        lasold = las
        larold = lar

        !         loop through all neigboring tiles.

        l = 0
        if (lts > 0) then
          do m= 1,mm_top
            l = l + 1
            if (nproc /= jpr) then
              ls0 = i0_st(m)*nhl*(ld-l1+1)
              lm  = ii_st(m)*nhl*(ld-l1+1)
              call mpi_send_init( &
                   ai(ls0+1,1),lm,MTYPER, &
                   idproc(m0_top+m,nproc+1), 9905, &
                   mpicomm, mpireqa(l), mpierr)
            else !arctic
              call mpi_send_init( &
                   ai(1,1),lts,MTYPER, &
                   idproc(m0_top+m,nproc+1), 99051, &
                   mpicomm, mpireqa(l), mpierr)
            end if
          end do
        end if ! lts > 0
        if (las > 0) then
          l = l + 1
          call mpi_send_init( &
               aia(1,1),las,MTYPER, &
               idproc(mod(ipr+1-mproc,ipr)+1,nproc), 99052, &
               mpicomm, mpireqa(l), mpierr)
        end if ! las > 0
        if (lbs > 0) then
          do m= 1,mm_bot
            l   = l + 1
            ls0 = i0_sb(m)*nhl*(ld-l1+1)
            lm  = ii_sb(m)*nhl*(ld-l1+1)
            call mpi_send_init( &
                 ai(ls0+1,2),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9906, &
                 mpicomm, mpireqa(l), mpierr)
          end do
        end if ! lbs > 0
        if (ltr > 0) then
          do m= 1,mm_top
            l   = l + 1
            if (nproc /= jpr) then
              ls0 = i0_st(m)*nhl*(ld-l1+1)
              lm  = ii_st(m)*nhl*(ld-l1+1)
              call mpi_recv_init( &
                   ai(ls0+1,4),lm,MTYPER, &
                   idproc(m0_top+m,nproc+1), 9906, &
                   mpicomm, mpireqa(l), mpierr)
            else !arctic
              call mpi_recv_init( &
                   ai(1,4),ltr,MTYPER, &
                   idproc(m0_top+m,nproc+1), 99051, &
                   mpicomm, mpireqa(l), mpierr)
            end if
          end do
        end if ! ltr > 0
        if (lar > 0) then
          l = l + 1
          call mpi_recv_init( &
               aia(1,2),lar,MTYPER, &
               idproc(mod(ipr+1-mproc,ipr)+1,nproc), 99052, &
               mpicomm, mpireqa(l), mpierr)
        end if ! lar > 0
        if (lbr > 0) then
          do m= 1,mm_bot
            l   = l + 1
            ls0 = i0_sb(m)*nhl*(ld-l1+1)
            lm  = ii_sb(m)*nhl*(ld-l1+1)
            call mpi_recv_init( &
                 ai(ls0+1,3),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9905, &
                 mpicomm, mpireqa(l), mpierr)
          end do
        end if ! lbr > 0
        nreqa = l
      end if
      if (nreqa > 0) then
        call mpi_startall(nreqa, mpireqa,          mpierr)
        call mpi_waitall( nreqa, mpireqa, mpistat, mpierr)
      end if

      if (nproc /= jpr) then  !arctic
        if (lbr > 0) then
          lbr = 0
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            do k= l1,ld
              do j= 1,nhl
                lbr = lbr + 1
                a(i, 1-j,k) = ai(lbr,3)
                a(i,jj+j,k) = ai(lbr,4)
              end do
            end do
          end do
        end if ! lbr > 0
      else !arctic
        if (lar > 0) then
          lar  = 0
          if (itype == 3 .or. itype == 13 .or. mproc > ipr/2) then
            do k= l1,ld
              do j= 0,nhl
                lar  = lar + 1
                ai(lar,4) = aia(lar,2)
              end do !j
            end do !k
          else
            do k= l1,ld
              do j= 1,nhl
                lar  = lar + 1
                ai(lar,4) = aia(lar,2)
              end do !j
            end do !k
          end if
        end if ! lar > 0
        if (ltr > 0) then
          if (lbr == 0) then
            ltr = 0
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              do k= l1,ld
                do j= 0,nhl
                  ltr = ltr + 1
                  a(i,jj+j,k) = ai(ltr,4)
                end do
              end do
            end do
          else if (ltr > lbr) then
            lbr = 0
            ltr = 0
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              do k= l1,ld
                do j= 0,nhl
                  ltr = ltr + 1
                  a(i,jj+j,k) = ai(ltr,4)
                end do
                do j= 1,nhl
                  lbr = lbr + 1
                  a(i, 1-j,k) = ai(lbr,3)
                end do
              end do
            end do
          else
            lbr = 0
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              do k= l1,ld
                do j= 1,nhl
                  lbr = lbr + 1
                  a(i, 1-j,k) = ai(lbr,3)
                  a(i,jj+j,k) = ai(lbr,4)
                end do
              end do
            end do
          end if
        end if ! ltr > 0
      end if !arctic
    end if ! ipr == 1 .and. jpr == 1

    if (mhl > 0) then
      if (ipr == 1) then
        if (nreg == 0.or.nreg == 4) then
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                a( 1-i,j,k) = vland
                a(ii+i,j,k) = vland
              end do
            end do
          end do
        else
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                a( 1-i,j,k) = a(ii+1-i,j,k)
                a(ii+i,j,k) = a(     i,j,k)
              end do
            end do
          end do
        end if
      else
        l = 0
        do k= l1,ld
          do j= 1-nhl,jj+nhl
            do i= 1,mhl
              l = l + 1
              aj(l,1) = a(ii+1-i,j,k)
              aj(l,2) = a(     i,j,k)
              aj(l,3) = vland
              aj(l,4) = vland
            end do
          end do
        end do

        if (ilold /= l) then
          if (ilold /= 0) then
            call mpi_request_free(mpireqb(1), mpierr)
            call mpi_request_free(mpireqb(2), mpierr)
            call mpi_request_free(mpireqb(3), mpierr)
            call mpi_request_free(mpireqb(4), mpierr)
          end if
          ilold = l
          call mpi_send_init( &
               aj(1,1),l,MTYPER,idhalo(2), 9907, &
               mpicomm, mpireqb(1), mpierr)
          call mpi_send_init( &
               aj(1,2),l,MTYPER,idhalo(1), 9908, &
               mpicomm, mpireqb(2), mpierr)
          call mpi_recv_init( &
               aj(1,3),l,MTYPER,idhalo(1), 9907, &
               mpicomm, mpireqb(3), mpierr)
          call mpi_recv_init( &
               aj(1,4),l,MTYPER,idhalo(2), 9908, &
               mpicomm, mpireqb(4), mpierr)
        end if
        call mpi_startall(4, mpireqb,          mpierr)
        call mpi_waitall( 4, mpireqb, mpistat, mpierr)

        l = 0
        do k= l1,ld
          do j= 1-nhl,jj+nhl
            do i= 1,mhl
              l = l + 1
              a( 1-i,j,k) = aj(l,3)
              a(ii+i,j,k) = aj(l,4)
            end do
          end do
        end do
      end if ! ipr == 1:else
    end if ! mhl > 0

    if (use_TIMER) then
      if (nxc == 12) then
        call xctmr1(12)
        nxc = 0
      end if
    end if

  end subroutine xctilr_arctic

  !-----------------------------------------------------------------------
  recursive subroutine xctilr_nonarctic(a,l1,ld,mh,nh,itype)
    integer, intent(in)    :: l1,ld,mh,nh,itype
    real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)

    !-----------
    !  1) update the tile overlap halo of a real array.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    a               real           in/out    target array
    !    l1              integer        input     3rd dim. start index
    !    ld              integer        input     3rd dimension of a
    !    mh              integer        input     1st (EW) update halo size
    !    nh              integer        input     2nd (NS) update halo size
    !    itype           integer        input     grid and field type

    !  3) itype selects both the grid and field type
    !        itype= 1; p-grid, scalar field
    !        itype= 2; q-grid, scalar field
    !        itype= 3; u-grid, scalar field
    !        itype= 4; v-grid, scalar field
    !        itype=11; p-grid, vector field
    !        itype=12; q-grid, vector field
    !        itype=13; u-grid, vector field
    !        itype=14; v-grid, vector field
    !     it is ignored here because all types are the same unless
    !      the grid includes the arctic ocean

    !  4) the global variable vland is returned by halos over land.
    !-----------

    integer, parameter :: lsize = kdm
    integer, parameter :: ilen  = idm*lsize*nbdy+64
    integer, parameter :: jlen  = (jdm+2*nbdy)*lsize*nbdy+64

    !  halo buffer
    real   , save :: ai(ilen,4),aj(jlen,4)
    integer       :: i,j,k,l,ls0,lm,m,mhl,nhl
    integer, save :: mpireqa(4*iqr),mpireqb(4),ilold,jlold,nreqa ! persistent communication handles.
    data ilold,jlold / 0,0 /

    ! split large requests into smaller pieces
    if (ld-l1+1 > lsize) then
      do k= l1,ld,lsize
        l = min(k+lsize-1,ld)
        call xctilr(a,k,l,mh,nh,itype)
      end do
      return
    end if

    if (use_TIMER) then
      if (nxc == 0) then
        call xctmr0(12)
        nxc = 12
      end if
    end if

    mhl = max(0,min(mh,nbdy))
    nhl = max(0,min(nh,nbdy))

    if (nhl > 0) then
      if (jpr == 1) then
        if (nreg == 0.or.nreg == 1) then
          do k= l1,ld
            do j= 1,nhl
              do i= 1,ii
                a(i, 1-j,k) = vland
                a(i,jj+j,k) = vland
              end do
            end do
          end do
        else
          do k= l1,ld
            do j= 1,nhl
              do i= 1,ii
                a(i, 1-j,k) = a(i,jj+1-j,k)
                a(i,jj+j,k) = a(i,     j,k)
              end do
            end do
          end do
        end if
      else
        l = 0
        do i= 1,ii  ! outer loop to simplify multiple neighbor case
          do k= l1,ld
            do j= 1,nhl
              l = l + 1
              ai(l,1) = a(i,jj+1-j,k)
              ai(l,2) = a(i,     j,k)
              ai(l,3) = vland
              ai(l,4) = vland
            end do
          end do
        end do
        if (jlold /= l) then
          if (jlold /= 0) then
            do i= 1,nreqa
              call mpi_request_free(mpireqa(i), mpierr)
            end do
          end if
          jlold = l
          !           loop through all neigboring tiles.
          l = 0
          do m= 1,mm_top
            l   = l + 1
            ls0 = i0_st(m)*nhl*(ld-l1+1)
            lm  = ii_st(m)*nhl*(ld-l1+1)
            call mpi_send_init( &
                 ai(ls0+1,1),lm,MTYPER,idproc(m0_top+m,nproc+1), 9905, &
                 mpicomm, mpireqa(l), mpierr)
          end do
          do m= 1,mm_bot
            l   = l + 1
            ls0 = i0_sb(m)*nhl*(ld-l1+1)
            lm  = ii_sb(m)*nhl*(ld-l1+1)
            call mpi_send_init( &
                 ai(ls0+1,2),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9906, &
                 mpicomm, mpireqa(l), mpierr)
          end do
          do m= 1,mm_top
            l   = l + 1
            ls0 = i0_st(m)*nhl*(ld-l1+1)
            lm  = ii_st(m)*nhl*(ld-l1+1)
            call mpi_recv_init( &
                 ai(ls0+1,4),lm,MTYPER,idproc(m0_top+m,nproc+1), 9906, &
                 mpicomm, mpireqa(l), mpierr)
          end do
          do m= 1,mm_bot
            l   = l + 1
            ls0 = i0_sb(m)*nhl*(ld-l1+1)
            lm  = ii_sb(m)*nhl*(ld-l1+1)
            call mpi_recv_init( &
                 ai(ls0+1,3),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9905, &
                 mpicomm, mpireqa(l), mpierr)
          end do
          nreqa = l
        end if
        if (nreqa > 0) then
          call mpi_startall(nreqa, mpireqa,          mpierr)
          call mpi_waitall( nreqa, mpireqa, mpistat, mpierr)
        end if

        l = 0
        do i= 1,ii  ! outer loop to simplify multiple neighbor case
          do k= l1,ld
            do j= 1,nhl
              l = l + 1
              a(i, 1-j,k) = ai(l,3)
              a(i,jj+j,k) = ai(l,4)
            end do
          end do
        end do
      end if ! jpr == 1:else
    end if ! nhl > 0

    if (mhl > 0) then
      if (ipr == 1) then
        if (nreg == 0.or.nreg == 4) then
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                a( 1-i,j,k) = vland
                a(ii+i,j,k) = vland
              end do
            end do
          end do
        else
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                a( 1-i,j,k) = a(ii+1-i,j,k)
                a(ii+i,j,k) = a(     i,j,k)
              end do
            end do
          end do
        end if
      else
        l = 0
        do k= l1,ld
          do j= 1-nhl,jj+nhl
            do i= 1,mhl
              l = l + 1
              aj(l,1) = a(ii+1-i,j,k)
              aj(l,2) = a(     i,j,k)
              aj(l,3) = vland
              aj(l,4) = vland
            end do
          end do
        end do

        if (ilold /= l) then
          if (ilold /= 0) then
            call mpi_request_free(mpireqb(1), mpierr)
            call mpi_request_free(mpireqb(2), mpierr)
            call mpi_request_free(mpireqb(3), mpierr)
            call mpi_request_free(mpireqb(4), mpierr)
          end if
          ilold = l
          call mpi_send_init( &
               aj(1,1),l,MTYPER,idhalo(2), 9907, &
               mpicomm, mpireqb(1), mpierr)
          call mpi_send_init( &
               aj(1,2),l,MTYPER,idhalo(1), 9908, &
               mpicomm, mpireqb(2), mpierr)
          call mpi_recv_init( &
               aj(1,3),l,MTYPER,idhalo(1), 9907, &
               mpicomm, mpireqb(3), mpierr)
          call mpi_recv_init( &
               aj(1,4),l,MTYPER,idhalo(2), 9908, &
               mpicomm, mpireqb(4), mpierr)
        end if
        call mpi_startall(4, mpireqb,          mpierr)
        call mpi_waitall( 4, mpireqb, mpistat, mpierr)

        l = 0
        do k= l1,ld
          do j= 1-nhl,jj+nhl
            do i= 1,mhl
              l = l + 1
              a( 1-i,j,k) = aj(l,3)
              a(ii+i,j,k) = aj(l,4)
            end do
          end do
        end do
      end if ! ipr == 1:else
    end if ! mhl > 0

    if (use_TIMER) then
      if (nxc == 12) then
        call xctmr1(12)
        nxc = 0
      end if
    end if

 end subroutine xctilr_nonarctic

  !-----------------------------------------------------------------------
  subroutine xctmri()
    !-----------
    !  1) initialize timers.
    !  2) timers  1:32 are for message passing routines,
    !     timers 33:80 are for general hycom routines,
    !     timers 81:96 are for user selected routines.
    !     timer     97 is the total time.
    !  3) call xctmri    to initialize timers (called in xcspmd),
    !     call xctmr0(n) to start timer n,
    !     call xctmr1(n) to stop  timer n and add event to timer sum,
    !     call xctnrn(n,cname) to register a name for timer n,
    !     call xctmrp to printout timer statistics (called by xcstop).
    !  4) time every 50-th event above 1,000.
    !-----------
    integer :: i
    real(8), parameter :: zero8 = 0.0

    nxc = 0
    do i= 1,97
      cc(i) = '      '
      nc(i) = 0
      tc(i) = zero8
    end do

    call xctmrn(97,'total ')
    call xctmr0(97)
  end subroutine xctmri

  !-----------------------------------------------------------------------
  subroutine xctmr0(n)
    integer, intent(in) :: n ! timer number
    !-----------
    !  start timer n.
    !  time every 50-th event above 1,000.
    !-----------
    if (use_DEBUG_TIMER_ALL) then
      if (cc(n) /= '      ') then
        write(lp,'(i5,2x,a,a)') mnproc,'call ',cc(n)
      end if
    end if
    if (use_DEBUG_TIMER) then
      if (n > 32 .and. cc(n) /= '      ') then
        if (mnproc == 1) then
          write(lp,*) 'call ',cc(n)
        end if
      end if
    end if
    if (timer_on) then
      if (mod(nc(n),50) == 0 .or. nc(n) <= 1000) then
        t0(n) = wtime()
      end if
    end if !timer_on
  end subroutine xctmr0

  !-----------------------------------------------------------------------
  subroutine xctmr1(n)
    integer, intent(in) :: n ! timer number
    !-----------
    !  add time since call to xctim0 to timer n.
    !  time every 50-th event above 1,000.
    !-----------
    if (timer_on) then
      if (nc(n) > 1000) then
        if (mod(nc(n),50) == 0) then
          tc(n) = tc(n) + 50.0*(wtime() - t0(n))
        end if
      else
        tc(n) = tc(n) + (wtime() - t0(n))
      end if
      nc(n) = nc(n) + 1
    end if !timer_on
    if (use_DEBUG_TIMER_ALL) then
      if (cc(n) /= '      ') then
        write(lp,'(i5,2x,a,a)') mnproc,'exit ',cc(n)
      end if
    end if
    if (use_DEBUG_TIMER) then
      if (n > 32 .and. cc(n) /= '      ') then
        if (mnproc == 1) then
          write(lp,*) 'exit ',cc(n)
        end if
      end if
    end if
  end subroutine xctmr1

  !-----------------------------------------------------------------------
  subroutine xctmrn(n,cname)
    character(len=6), intent(in) :: cname ! timer name
    integer,          intent(in) :: n     !  timer number
    !-----------
    ! register name of timer n.
    !-----------
    cc(n) = cname
  end subroutine xctmrn

  !-----------------------------------------------------------------------
  subroutine xctmrp
    !-----------
    !  print all active timers.
    !  on exit all timers are reset to zero.
    !-----------
    integer :: i,mnloc
    real(8), parameter :: zero8 = 0.0

    ! get total time.
    call xctmr1(97)

    ! report time on the processor with the least communication overhead
    tcxc(2) = mnproc
    tcxc(1) = zero8
    do i= 1,32
      if (nc(i) /= 0 .and. cc(i)(1:2) == 'xc') then
        tcxc(1) = tcxc(1) + tc(i)  !communication overhead
      end if
    end do !i

    if (ijpr /= 1) then
      tcxl(1) = tcxc(1)
      tcxl(2) = tcxc(2)
      call mpi_allreduce(tcxl,tcxc,1, &
           mpi_2double_precision,mpi_minloc, &
           mpicomm,mpierr)
      mnloc = tcxc(2)  !processor with the least comm. overhead
      if (mnproc == 1) then
        if (mnloc /= 1) then
          call MPI_RECV(tc,97,MTYPED, &
               idproc1(mnloc), 9949, &
               mpicomm, mpistat, mpierr)
        end if
      else if (mnproc == mnloc) then
        call MPI_SEND(tc,97,MTYPED, &
             idproc1(1), 9949, &
             mpicomm, mpierr)
      end if
    end if

    call xcsync(flush_lp)
    if (mnproc == 1) then
      write(lp,6000) mnloc,ijpr
      do i= 1,32
        if (nc(i) /= 0) then
          if (cc(i) /= '      ') then
            write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
          else
            write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
          end if
        end if
      end do !i
      write(lp,6100) 'xc****',1,tcxc(1),tcxc(1)
      do i= 33,97
        if (nc(i) /= 0) then
          if (cc(i) /= '      ') then
            write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
          else
            write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
          end if
        end if
      end do !i
      write(lp,6200)
    end if !mnproc == 1
    call xcsync(flush_lp)

    ! reset timers to zero.

    do i= 1,97
      nc(i) = 0
      tc(i) = zero8
    end do
    tcxc(1) = zero8

    ! start a new total time measurement.

    call xctmr0(97)

6000 format(/ / &
         3x,' timer statistics, processor',i5,' out of',i5 / &
         3x,'-----------------------------------------------' /)
6100 format(3x,a6, &
         '   calls =',i9, &
         '   time =',f11.5, &
         '   time/call =',f14.8)
6150 format(3x,'   #',i2, &
         '   calls =',i9, &
         '   time =',f11.5, &
         '   time/call =',f14.8)
6200 format(/ /)
  end subroutine xctmrp

  !-----------------------------------------------------------------------
  subroutine xcgetrow(outm, inm, kt)
    integer, intent(in)    :: kt
    real,    intent(out)   :: outm(itdm,jj,kt)
    real,    intent(in)    :: inm(ii,jj,kt)

    !-----------
    !  convert an entire 2-D array from tiled to non-tiled layout.
    !-----------

    integer :: mpireqb(ipr)
    real    :: at(idm*jdm*kt),ata(idm*jdm*kt,iqr)
    integer :: i,j,k,l,mp

    ! gather each row of tiles onto the first tile in the row.
    if (mproc == mpe_1(nproc)) then
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            outm(i0+i,j,k) = inm(i,j,k)
          end do
        end do
      end do
      l = 0
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = l + 1
        call MPI_IRECV(ata(1,l),ii_pe(mp,nproc)*jj*kt,MTYPER, &
             idproc(mp,nproc), 9941, &
             mpicomm, mpireqb(l), mpierr)
      end do
      call mpi_barrier(mpicomm, mpierr)
      l = 0
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = l + 1
        call MPI_WAIT(mpireqb(l), mpistat, mpierr)
        do k= 1,kt
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
              outm(i0_pe(mp,nproc)+i,j,k) = ata(i+(j-1)*ii_pe(mp,nproc) &
                   +(k-1)*jj*ii_pe(mp,nproc),l)
            end do
          end do
        end do
      end do
    else  !mproc>1
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            at(i+(j-1)*ii+(k-1)*ii*jj) = inm(i,j,k)
          end do
        end do
      end do
      call mpi_barrier(mpicomm, mpierr)
      call MPI_SEND(at,ii*jj*kt,MTYPER, &
           idproc(mpe_1(nproc),nproc), 9941, &
           mpicomm, mpierr)
    end if

  end subroutine xcgetrow

  !-----------------------------------------------------------------------
  subroutine xcgetrow4(outm,inm, kt)
    integer, intent(in)      :: kt
    real*4,    intent(out)   :: outm(itdm,jj,kt)
    real*4,    intent(in)    :: inm(ii,jj,kt)

    !-----------
    !  convert an entire 2-D array from tiled to non-tiled layout.
    !-----------

    integer :: mpireqb(ipr)
    real*4  :: at(idm*jdm*2*kk),ata(idm*jdm*2*kk,iqr)
    integer :: i,j,k,l,mp

    !  gather each row of tiles onto the first tile in the row.
    if (mproc == mpe_1(nproc)) then
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            outm(i0+i,j,k) = inm(i,j,k)
          end do
        end do
      end do
      l = 0
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = l + 1
        call MPI_IRECV(ata(1,l),ii_pe(mp,nproc)*jj*kt,mpi_real, &
             idproc(mp,nproc), 9941, &
             mpicomm, mpireqb(l), mpierr)
      end do
      call mpi_barrier(mpicomm, mpierr)
      l = 0
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = l + 1
        call MPI_WAIT(mpireqb(l), mpistat, mpierr)
        do k= 1,kt
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
              outm(i0_pe(mp,nproc)+i,j,k) = ata(i+(j-1)*ii_pe(mp,nproc) &
                   +(k-1)*jj*ii_pe(mp,nproc),l)
            end do
          end do
        end do
      end do
    else  !mproc>1
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            at(i+(j-1)*ii+(k-1)*ii*jj) = inm(i,j,k)
          end do
        end do
      end do
      call mpi_barrier(mpicomm, mpierr)
      call MPI_SEND(at,ii*jj*kt,mpi_real, &
           idproc(mpe_1(nproc),nproc), 9941, &
           mpicomm, mpierr)
    end if

  end subroutine xcgetrow4

  !-----------------------------------------------------------------------
  subroutine xcgetrowint2(outm,inm, kt)
    integer,   intent(in)  :: kt
    integer*2, intent(out) :: outm(itdm,jj,kt)
    integer*2, intent(in)  :: inm(ii,jj,kt)

    !-----------
    ! convert an entire 2-D array from tiled to non-tiled layout.
    !-----------

    integer   :: mpireqb(ipr)
    integer*2 :: at(idm*jdm*kt),ata(idm*jdm*kt,iqr)
    integer   :: i,j,k,l,mp,np

    !     gather each row of tiles onto the first tile in the row.

    if (mproc == mpe_1(nproc)) then
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            outm(i0+i,j,k) = inm(i,j,k)
         end do
        end do
      end do
      l = 0
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = l + 1
        call MPI_IRECV(ata(1,l),ii_pe(mp,nproc)*jj*kt,mpi_integer2, &
             idproc(mp,nproc), 9941, &
             mpicomm, mpireqb(l), mpierr)
      end do
      call mpi_barrier(mpicomm, mpierr)
      l = 0
      do mp= mpe_1(nproc)+1,mpe_e(nproc)
        l = l + 1
        call MPI_WAIT(mpireqb(l), mpistat, mpierr)
        do k= 1,kt
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
              outm(i0_pe(mp,nproc)+i,j,k) = ata(i+(j-1)*ii_pe(mp,nproc) &
                   +(k-1)*jj*ii_pe(mp,nproc),l)
            end do
          end do
        end do
      end do
    else  !mproc>1
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            at(i+(j-1)*ii+(k-1)*ii*jj) = inm(i,j,k)
          end do
        end do
      end do
      call mpi_barrier(mpicomm, mpierr)
      call MPI_SEND(at,ii*jj*kt,mpi_integer2, &
           idproc(mpe_1(nproc),nproc), 9941, &
           mpicomm, mpierr)
    end if

  end subroutine xcgetrowint2

!**************************************************************************************
#else
!**************************************************************************************

  !-----------------------------------------------------------------------
  !     auxillary routines that involve off-processor communication.
  !     shared memory version, contained in module mod_xc.
  !     author:  Alan J. Wallcraft,  NRL.
  !     The following was previously contained in  mod_xc_sm.inc
  !-----------------------------------------------------------------------

  subroutine xcaget(aa, a, mnflg)

    real,    intent(out)   :: aa(itdm,jtdm)                      ! non-tiled target array
    real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! tiled source array
    integer, intent(in)    :: mnflg ! node return flag
                                    ! = 0; all nodes
                                    ! = n; node number n (mnproc=n)

    !-----------
    !  onvert an entire 2-D array from tiled to non-tiled layout.
    !-----------

    integer :: j

    if (use_TIMER) then
      ! call xctmr0( 1)
    end if
    do j= 1,jtdm
      call xclget(aa(1,j),itdm, a, 1,j,1,0, mnflg)
    end do
    if (use_TIMER) then
      !     call xctmr1( 1)
    end if
  end subroutine xcaget

  !-----------------------------------------------------------------------
  subroutine xcaput(aa, a, mnflg)
    real,    intent(in)  :: aa(itdm,jtdm)                      ! non-tiled target array
    real,    intent(out) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! tiled source array
    integer, intent(in)  :: mnflg ! node source flag
                                    ! = 0; all nodes
                                    ! = n; node number n (mnproc=n)
    !-----------
    !  convert an entire 2-D array from non-tiled to tiled layout.
    !-----------

    integer :: j

    if (use_TIMER) then
      !  call xctmr0( 4)
    end if
    do j= 1,jtdm
      call xclput(aa(1,j),itdm, a, 1,j,1,0)
    end do
    if (use_TIMER) then
      !     call xctmr1( 4)
    end if
  end subroutine xcaput

  !-----------------------------------------------------------------------
  subroutine xceget(aelem, a, ia,ja)
    real,    intent(out) :: aelem
    real,    intent(in)  :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
    integer, intent(in)  :: ia,ja

    !-----------
    !  1) find the value of a(ia,ja) on the non-tiled 2-D grid.
    !  2) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aelem           real           output    required element
    !    a               real           input     source array
    !    ia              integer        input     1st index into a
    !    ja              integer        input     2nd index into a
    !-----------
    if (use_TIMER) then
      call xctmr0( 2)
    end if
    ! single node version - trivial indexing.
    aelem = a(ia,ja)
    if (use_TIMER) then
      call xctmr1( 2)
    end if
  end subroutine xceget

  !-----------------------------------------------------------------------
  subroutine xceput(aelem, a, ia,ja)
    integer, intent(in)    :: ia    ! 1st index into a
    integer, intent(in)    :: ja    ! 2nd index into a
    real,    intent(in)    :: aelem ! element value
    real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! target array
    !-----------
    ! fill a single element in the non-tiled 2-D grid.
    !-----------
    if (use_TIMER) then
      call xctmr0( 4)
    end if
    ! single node version - trivial indexing.
    a(ia,ja) = aelem
    if (use_TIMER) then
      call xctmr1( 4)
    end if
  end subroutine xceput

  !-----------------------------------------------------------------------
  subroutine xchalt(cerror)
    character*(*), intent(in) :: cerror

    !-----------
    !  1) stop all processes.
    !  2) only one processes need call this routine, i.e. it is for
    !      emergency stops.  use 'xcstop' for ordinary stops called
    !      by all processes.
    !  3) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    cerror          char*(*)       input     error message
    !-----------

    !     shared memory version, just stop.

    if (cerror /= ' ') then
      write(lp,*) '**************************************************'
      write(lp,*) cerror
      write(lp,*) '**************************************************'
    end if
    error stop '(xchalt)'
  end subroutine xchalt

  !-----------------------------------------------------------------------
  subroutine xclget(aline,nl, a, i1,j1,iinc,jinc, mnflg)
    integer, intent(in)    ::  nl,i1,j1,iinc,jinc,mnflg
    real,    intent(out)   ::  aline(nl)
    real,    intent(in)    ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)

    !-----------
    !  1) extract a line of elements from the non-tiled 2-D grid.
    !  2) aline(i) = a(i1+iinc*(i-1),j1+jinc*(i-1)), for i=1...nl.
    !     iinc and jinc can each be -1, 0, or +1.

    !     if jinc=0, j1 can be between jtdm+1 and jtdm+nbdy to return
    !     values from the top halo.  This is for debugging the arctic
    !     patch halo exchange only.

    !  3) mnflg selects which nodes must return the line
    !        =-n; node number n (mnproc=n), nl,i1,j1 only on node n
    !        = 0; all nodes
    !        = n; node number n (mnproc=n)
    !     normally all integer arguments must be identical on all nodes,
    !     but a negative mnflg indicates that only the target node is
    !     providing the nl,i1,j1 values.  These are broadcast to all other
    !     nodes, and returned in nl,i1,j1 by all of them.

    !     mnflg is ignored here (only have a single node).

    !  4) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aline           real           output    required line of elements
    !    nl              integer        in/out    dimension of aline
    !    a               real           input     source array
    !    i1              integer        in/out    1st index into a
    !    j1              integer        in/out    2nd index into a
    !    iinc            integer        input     1st index increment
    !    jinc            integer        input     2nd index increment
    !    mnflg           integer        input     node return flag

    !    nl,i1,j1 are input only unless mnflg is negative.
    !-----------

    integer :: i
    if (use_TIMER) then
      call xctmr0( 3)
    end if

    !     single node version - trivial indexing and no error checking.

    if (jinc == 0) then
      do i= 1,nl
        aline(i) = a(i1+iinc*(i-1),j1)
      end do
    else if (iinc == 0) then
      do i= 1,nl
        aline(i) = a(i1,j1+jinc*(i-1))
      end do
    else
      do i= 1,nl
        aline(i) = a(i1+iinc*(i-1),j1+jinc*(i-1))
      end do
    end if
    if (use_TIMER) then
      call xctmr1( 3)
    end if
  end subroutine xclget

  !-----------------------------------------------------------------------
  subroutine xclput(aline,nl, a, i1,j1,iinc,jinc)
    real,    intent(in)  ::  aline(nl)
    integer, intent(in)  ::  nl
    real,    intent(out) ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
    integer, intent(in)  ::  i1,j1,iinc,jinc

    !-----------
    !  1) fill a line of elements in the non-tiled 2-D grid.
    !  2) aline(i) = a(i1+i1*(i-1),j1+j1*(i-1)), for i=1...nl.
    !     one of iinc and jinc must be 0, and the other must be 1.
    !  3) parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    aline           real           input     line of element values
    !    nl              integer        input     dimension of aline
    !    a               real           in/out    target array
    !    i1              integer        input     1st index into a
    !    j1              integer        input     2nd index into a
    !    iinc            integer        input     1st index increment
    !    jinc            integer        input     2nd index increment
    !-----------

    integer :: i
    if (use_TIMER) then
      call xctmr0( 4)
    end if

    ! single node version - trivial indexing.

    if (jinc == 0) then
      do i= 1,nl
        a(i1+i-1,j1) = aline(i)
      end do
    else if (iinc == 0) then
      do i= 1,nl
        a(i1,j1+i-1) = aline(i)
      end do
    end if
    if (use_TIMER) then
      call xctmr1( 4)
    end if
  end subroutine xclput

  !-----------------------------------------------------------------------
  subroutine xcbcst_i0(ia)
    integer, intent(inout) :: ia  ! target variable
    !-----------
    ! broadcast integer scalar from tile mnproc = 1 to all tiles
    !-----------
    if (use_TIMER) then
      call xctmr0( 8)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 8)
    end if
  end subroutine xcbcst_i0

  !-----------------------------------------------------------------------
  subroutine xcbcst_i1(ia)
    integer, intent(inout) :: ia(:) ! target array
    !-----------
    !  broadcast integer array from tile mnproc = 1 to all tiles
    !-----------
    if (use_TIMER) then
      call xctmr0( 8)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 8)
    end if
  end subroutine xcbcst_i1

  !-----------------------------------------------------------------------
  subroutine xcbcst_r0(ra)
    real, intent(inout) :: ra ! target variable
    !-----------
    !  broadcast real scalar from tile mnproc = 1 to all tiles
    !-----------
    if (use_TIMER) then
      call xctmr0( 8)
    end if
    !     single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 8)
    end if
  end subroutine xcbcst_r0

  !-----------------------------------------------------------------------
  subroutine xcbcst_r1(ra)
    real, intent(inout) :: ra(:) ! target array
    !-----------
    !  broadcast real array from tile mnproc = 1 to all tiles
    !-----------
    if (use_TIMER) then
      call xctmr0( 8)
    end if
    !     single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 8)
    end if
  end subroutine xcbcst_r1

  !-----------------------------------------------------------------------
  subroutine xcbcst_l0(la)
    logical, intent(inout) :: la ! target variable
    !-----------
    !  broadcast logical from tile mnproc = 1 to all tiles
    !-----------

    if (use_TIMER) then
      call xctmr0( 8)
    end if
    !     single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 8)
    end if
  end subroutine xcbcst_l0

  !-----------------------------------------------------------------------
  subroutine xcbcst_l1(la)
    logical, intent(inout) :: la(:) ! target array
    !-----------
    !  broadcast logical array from tile mnproc = 1 to all tiles
    !-----------
    if (use_TIMER) then
      call xctmr0( 8)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 8)
    end if
  end subroutine xcbcst_l1

  !-----------------------------------------------------------------------
  subroutine xcbcst_c(ca)
    character*(*), intent(inout) :: ca ! target array
    !-----------
    !  broadcast character string from tile mnproc = 1 to all tiles
    !-----------

    if (use_TIMER) then
      call xctmr0( 8)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 8)
    end if
  end subroutine xcbcst_c

  !-----------------------------------------------------------------------
  subroutine xcmax_i0(ia)
    integer, intent(inout) :: ia ! target variable
    !-----------
    !  replace scalar a with its element-wise maximum over all tiles.
    !-----------
    if (use_TIMER) then
      call xctmr0( 9)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 9)
    end if
  end subroutine xcmax_i0

  !-----------------------------------------------------------------------
  subroutine xcmax_i1(ia)
    integer, intent(inout) :: ia(:) ! target array
    !-----------
    !  replace array a with its element-wise maximum over all tiles.
    !-----------
    if (use_TIMER) then
      call xctmr0( 9)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 9)
    end if
  end subroutine xcmax_i1

  !-----------------------------------------------------------------------
  subroutine xcmax_r0(ra)
    real, intent(inout) :: ra ! target variable
    !-----------
    !  replace scalar a with its element-wise maximum over all tiles.
    !-----------
    if (use_TIMER) then
      call xctmr0( 9)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 9)
    end if
  end subroutine xcmax_r0

  !-----------------------------------------------------------------------
  subroutine xcmax_r1(ra)
    real, intent(inout) :: ra(:) ! target array
    !-----------
    !  replace array a with its element-wise maximum over all tiles.
    !-----------
    if (use_TIMER) then
      call xctmr0( 9)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1( 9)
    end if
  end subroutine xcmax_r1

  !-----------------------------------------------------------------------
  subroutine xcmin_i0(ia)
    integer, intent(inout) :: ia ! target variable
    !-----------
    !  replace scalar a with its element-wise minimum over all tiles.
    !-----------
    if (use_TIMER) then
      call xctmr0(10)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1(10)
    end if
  end subroutine xcmin_i0

  !-----------------------------------------------------------------------
  subroutine xcmin_i1(ia)
    integer, intent(inout) :: ia(:) ! target array
    !-----------
    !  replace array a with its element-wise minimum over all tiles.
    !-----------
    if (use_TIMER) then
      call xctmr0(10)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1(10)
    end if
  end subroutine xcmin_i1

  !-----------------------------------------------------------------------
  subroutine xcmin_r0(ra)
    real, intent(inout) :: ra ! target variable
    !-----------
    !  replace scalar a with its element-wise minimum over all tiles.
    !-----------
    if (use_TIMER) then
      call xctmr0(10)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1(10)
    end if
  end subroutine xcmin_r0

  !-----------------------------------------------------------------------
  subroutine xcmin_r1(ra)
    real, intent(inout) :: ra(:) ! target array
    !-----------
    !  replace array a with its element-wise minimum over all tiles.
    !-----------
    if (use_TIMER) then
      call xctmr0(10)
    end if
    ! single node version - do nothing.
    if (use_TIMER) then
      call xctmr1(10)
    end if
  end subroutine xcmin_r1

  !-----------------------------------------------------------------------
  subroutine xcspmd()
    !-----------
    !  1) initialize data structures that identify the tiles.
    !  2) data structures:
    !      ipr     - 1st 2-D node dimension
    !      jpr     - 2nd 2-D node dimension
    !      ijpr    -     1-D node dimension (ipr*jpr)
    !      mproc   - 1st 2-D node index
    !      nproc   - 2nd 2-D node index
    !      mnproc  -     1-D node index
    !      i0      -     1st dimension tile offset
    !      ii      -     1st dimension tile extent
    !      j0      -     2nd dimension tile offset
    !      jj      -     2nd dimension tile extent
    !      margin  -     how much of the halo is currently valid
    !      nreg    -     region type
    !      vland   -     fill value for land (standard value 0.0)
    !  3) ipr,jpr,ijpr are global (tile independent) values.
    !     all other values depend on the processor number,
    !     but in this case there is only one processor.
    !-----------

    ! shared memory version, mproc=nproc=1.

    if (iqr /= 1 .or. jqr /= 1 .or. ijqr /= 1) then
      call xcstop('Error in xcspmd: must have iqr=jqr=ijqr = 1')
      stop '(xcspmd)'
    end if

    ipr    = 1
    jpr    = 1
    ijpr   = 1
    mnproc = 1
    mproc  = 1
    nproc  = 1

    i0  = 0
    ii  = itdm
    j0  = 0
    jj  = jtdm

    margin = 0

    if (use_ARCTIC) then
      nreg   =  2  ! arctic patch region type
    else
      nreg   = -1  ! unknown region type
    end if

    vland  = 0.0

    ! initialize timers.
    call xctmri
    if (use_TIMER) then
      call xctmrn( 1,'xcaget')
      call xctmrn( 2,'xceget')
      call xctmrn( 3,'xclget')
      call xctmrn( 4,'xcXput')
      call xctmrn( 5,'xcsum ')
      call xctmrn( 8,'xcbcst')
      call xctmrn( 9,'xcmax ')
      call xctmrn(10,'xcmin ')
      call xctmrn(12,'xctilr')
    end if
  end subroutine xcspmd

  !-----------------------------------------------------------------------
  subroutine xcstop(cerror)
    character*(*), intent(in) :: cerror ! error message
    !-----------
    !  1) stop all processes.
    !  2) all processes must call this routine.
    !     use 'xchalt' for emergency stops.
    !-----------
    ! print active timers.
    call xctmrp

    ! shared memory version, just stop.
    if (cerror /= ' ') then
      write(lp,*) '**************************************************'
      write(lp,*) cerror
      write(lp,*) '**************************************************'
    end if
    stop '(xcstop)'
  end subroutine xcstop

  !-----------------------------------------------------------------------
  subroutine xcsum(sum, a,mask)
    real(8),  intent(out)   :: sum ! sum of a
    real,     intent(inout) :: a(   1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! source array
    integer,  intent(in)    :: mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ! mask array

    !-----------
    ! sum a 2-d array, where mask==1
    ! sum is bit for bit reproducable for the same halo size, nbdy.
    !-----------
    ! Local variables
    real(8) :: zero8
    parameter (zero8 = 0.0)
    real(8) :: sum8,sum8p,sum8j(jdm)
    integer :: i,i1,j
    if (use_TIMER) then
      call xctmr0( 5)
    end if

    ! row sums in 2*nbdy+1 wide strips.
    !$OMP PARALLEL DO PRIVATE(j,i1,i,sum8,sum8p) &
    !$OMP SCHEDULE(STATIC,jblk)
    do j = 1,jdm
      sum8 = zero8
      do i1 = 1,idm,2*nbdy+1
        sum8p = zero8
        do i= i1,min(i1+2*nbdy,idm)
          if (mask(i,j) == 1) then
            sum8p = sum8p + a(i,j)
          end if
        end do
        sum8 = sum8 + sum8p
      end do
      sum8j(j) = sum8  ! use of sum8 minimizes false sharing of sum8j
    end do
    !$OMP END PARALLEL DO

    ! serial sum of rwo-sum loop.
    sum8 = sum8j(1)
    do j = 2,jdm
      sum8 = sum8 + sum8j(j)
    end do
    sum = sum8
    if (use_TIMER) then
      call xctmr1( 5)
    end if
  end subroutine xcsum

  !-----------------------------------------------------------------------
  subroutine xcsync(lflush)
    logical, intent(in) :: lflush
    !-----------
    !  1) barrier, no processor exits until all arrive (and flush stdout).
    !  2) some MPI implementations only flush stdout as a collective
    !     operation, and hence the lflush=.true. option to flush stdout.
    !  3) Only one processor, so the barrier is a no-op in this case.
    !-----------
    if (lflush) then
      call flush(lp)
    end if
  end subroutine xcsync

  !-----------------------------------------------------------------------
  subroutine xctilr(a,l1,ld,mh,nh,itype)
    integer, intent(in)    :: l1,ld,mh,nh,itype
    real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)

    !-----------
    !  update the tile overlap halo of a real array.
    !  parameters:
    !       name            type         usage            description
    !    ----------      ----------     -------  ----------------------------
    !    a               real           in/out    target array
    !    l1              integer        input     3rd dim. start index
    !    ld              integer        input     3rd dimension of a
    !    mh              integer        input     1st (EW) update halo size
    !    nh              integer        input     2nd (NS) update halo size
    !    itype           integer        input     grid and field type

    !   itype selects both the grid and field type
    !        itype= 1; p-grid, scalar field
    !        itype= 2; q-grid, scalar field
    !        itype= 3; u-grid, scalar field
    !        itype= 4; v-grid, scalar field
    !        itype=11; p-grid, vector field
    !        itype=12; q-grid, vector field
    !        itype=13; u-grid, vector field
    !        itype=14; v-grid, vector field
    !
    !  itype is ignored if region does not include artic ocean since
    !  it is ignored here because all types are the same unless
    !    the grid includes the arctic ocean
    !-----------

    integer :: i,io,j,k,mhl,nhl
    if (use_TIMER) then
      call xctmr0(12)
    end if

    mhl = max(0,min(mh,nbdy))
    nhl = max(0,min(nh,nbdy))

    ! ---------------------------
    if (use_ARCTIC) then
    ! ---------------------------

      do k= l1,ld
        ! southern boundary is closed.
        do j= 1,nhl
          do i= 1,ii
            a(i,1-j,k) = vland
          end do
        end do

        if (itype < 10) then
          ! scalar field
          if (itype ==  1) then
            ! p-grid
            do j= 0,nhl
              do i= 1,ii
                io = ii-mod(i-1,ii)
                a(i,jj+j,k) = a(io,jj-1-j,k)
              end do
            end do
          else if (itype ==  2) then
            ! q-grid
            do i= ii/2+1,ii
              io = mod(ii-(i-1),ii)+1
              a(i,jj,k) = a(io,jj,k)
            end do
            do j= 1,nhl
              do i= 1,ii
                io = mod(ii-(i-1),ii)+1
                a(i,jj+j,k) = a(io,jj-j,k)
              end do
            end do
          else if (itype ==  3) then
            ! u-grid
            do j= 0,nhl
              do i= 1,ii
                io = mod(ii-(i-1),ii)+1
                a(i,jj+j,k) = a(io,jj-1-j,k)
              end do
            end do
          else
            ! v-grid
            do i= ii/2+1,ii
              io = ii-mod(i-1,ii)
              a(i,jj,k) = a(io,jj,k)
            end do
            do j= 1,nhl
              do i= 1,ii
                io = ii-mod(i-1,ii)
                a(i,jj+j,k) = a(io,jj-j,k)
              end do
            end do
          end if
        else
          ! vector field, swap sign
          if (itype == 11) then
            ! p-grid
            do j= 0,nhl
              do i= 1,ii
                io = ii-mod(i-1,ii)
                a(i,jj+j,k) = -a(io,jj-1-j,k)
              end do
            end do
          else if (itype == 12) then
            ! q-grid
            do i= ii/2+1,ii
              io = mod(ii-(i-1),ii)+1
              a(i,jj,k) = -a(io,jj,k)
            end do
            do j= 1,nhl
              do i= 1,ii
                io = mod(ii-(i-1),ii)+1
                a(i,jj+j,k) = -a(io,jj-j,k)
              end do
            end do
          else if (itype == 13) then
            ! u-grid
            do j= 0,nhl
              do i= 1,ii
                io = mod(ii-(i-1),ii)+1
                a(i,jj+j,k) = -a(io,jj-1-j,k)
              end do
            end do
          else
            ! v-grid
            do i= ii/2+1,ii
              io = ii-mod(i-1,ii)
              a(i,jj,k) = -a(io,jj,k)
            end do
            do j= 1,nhl
              do i= 1,ii
                io = ii-mod(i-1,ii)
                a(i,jj+j,k) = -a(io,jj-j,k)
              end do
            end do
          end if
        end if
      end do

      if (mhl > 0) then
        do k= 1,ld
          do j= 1-nhl,jj+nhl
            do i= 1,mhl
              a( 1-i,j,k) = a(ii+1-i,j,k)
              a(ii+i,j,k) = a(     i,j,k)
            end do
          end do
        end do
      end if

    ! ---------------------------
    else  ! NOT use_ARCTIC
    ! ---------------------------

      if (nhl > 0) then
        if (nreg <= 2) then  ! closed in latitude
          do k= l1,ld
            do j= 1,nhl
              do i= 1,ii
                a(i, 1-j,k) = vland
                a(i,jj+j,k) = vland
              end do
            end do
          end do
        else  ! periodic (f-plane) in latitude
          do k= l1,ld
            do j= 1,nhl
              do i= 1,ii
                a(i, 1-j,k) = a(i,jj+1-j,k)
                a(i,jj+j,k) = a(i,     j,k)
              end do
            end do
          end do
        end if
      end if

      if (mhl > 0) then
        if (nreg == 0 .or. nreg == 4) then  ! closed in longitude
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                a( 1-i,j,k) = vland
                a(ii+i,j,k) = vland
              end do
            end do
          end do
        else  ! periodic in longitude
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                a( 1-i,j,k) = a(ii+1-i,j,k)
                a(ii+i,j,k) = a(     i,j,k)
              end do
            end do
          end do
        end if
      end if

    ! ---------------------------
    end if ! if use_ARCTIC
    ! ---------------------------

    if (use_TIMER) then
      call xctmr1(12)
    end if
  end subroutine xctilr

  subroutine xctmri
    !-----------
    !  1) initialize timers.
    !  2) timers  1:32 are for message passing routines,
    !     timers 33:80 are for general hycom routines,
    !     timers 81:96 are for user selected routines.
    !     timer     97 is the total time.
    !  3) call xctmri    to initialize timers (called in xcspmd),
    !     call xctmr0(n) to start timer n,
    !     call xctmr1(n) to stop  timer n and add event to timer sum,
    !     call xctnrn(n,cname) to register a name for timer n,
    !     call xctmrp to printout timer statistics (called by xcstop).
    !-----------

    integer :: i
    real(8) :: zero8
    parameter (zero8 = 0.0)

    do i= 1,97
      cc(i) = '      '
      nc(i) = 0
      tc(i) = zero8
    end do

    call xctmrn(97,'total ')
    call xctmr0(97)
  end subroutine xctmri

  subroutine xctmr0(n)
    integer, intent(in) :: n !  timer number
    !-----------
    ! start timer n.
    ! time every 50-th event above 1,000.
    !-----------
    if (use_DEBUG_TIMER) then
      if (n > 24 .and. cc(n) /= '      ') then
        write(lp,*) 'call ',cc(n)
      end if
    end if
    if (timer_on) then
      if (mod(nc(n),50) == 0 .or. nc(n) <= 1000) then
        t0(n) = wtime()
      end if
    end if !timer_on
  end subroutine xctmr0

  subroutine xctmr1(n)
    integer, intent(in) :: n !  timer number
    !-----------
    ! add time since call to xctim0 to timer n.
    ! time every 50-th event above 1,000.
    !-----------
    if (timer_on) then
      if (nc(n) > 1000) then
        if (mod(nc(n),50) == 0) then
          tc(n) = tc(n) + 50.0*(wtime() - t0(n))
        end if
      else
        tc(n) = tc(n) + (wtime() - t0(n))
      end if
      nc(n) = nc(n) + 1
    end if !timer_on
    if (use_DEBUG_TIMER) then
      if (n > 24 .and. cc(n) /= '      ') then
        write(lp,*) 'exit ',cc(n)
      end if
    end if
  end subroutine xctmr1

  subroutine xctmrn(n,cname)
    character*6, intent(in) :: cname  ! timer name
    integer,     intent(in) :: n      ! timer number
    !-----------
    !  register name of timer n.
    !-----------
    cc(n) = cname
  end subroutine xctmrn

  subroutine xctmrp
    !-----------
    ! print all active timers.
    ! on exit all timers are reset to zero.
    !-----------
    integer :: i
    real(8) :: zero8
    parameter (zero8 = 0.0)

    ! get total time.
    call xctmr1(97)

    ! print timers.
    write(lp,6000)
    do i= 1,97
      if (nc(i) /= 0) then
        if (cc(i) /= '      ') then
          write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
        else
          write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
        end if
      end if
    end do
    write(lp,6200)

    ! reset timers to zero.
    do i= 1,97
      nc(i) = 0
      tc(i) = zero8
    end do

    ! start a new total time measurement.
    call xctmr0(97)

6000 format(/ / &
         4x,' timer statistics ' / &
         4x,'------------------' /)
6100 format(5x,a6, &
         '   calls =',i9, &
         '   time =',f11.5, &
         '   time/call =',f14.8)
6150 format(5x,'   #',i2, &
         '   calls =',i9, &
         '   time =',f11.5, &
         '   time/call =',f14.8)
6200 format(/ /)
  end subroutine xctmrp

!**************************************************************************************
#endif
!**************************************************************************************

end module mod_xc
