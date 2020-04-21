! ------------------------------------------------------------------------------
! Copyright (C) 2005 HYCOM Consortium and contributors
! Copyright (C) 2006-2015 Lars Inge Enstad, Mats Bentsen, Alok Kumar Gupta
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

c/* BARRIER       set a barrier; for SPMD versions      */
#if defined(MPI)
# define BARRIER call mpi_barrier(mpicomm,mpierr)
#elif defined(SHMEM)
# define BARRIER call shmem_barrier_all()
#endif

c/* BARRIER_MP    halo synchronization; SHMEM only      */
c/* BARRIER_NP    halo synchronization; SHMEM only      */
#if defined(RINGB)
#define BARRIER_MP call xctbar(idproc(mproc-1,nproc),idproc(mproc+1,nproc))
#define BARRIER_NP call xctbar(idproc(mproc,nproc-1),idproc(mproc,nproc+1))
#else
#define BARRIER_MP BARRIER
#define BARRIER_NP BARRIER
#endif

#if defined(MPI)
c/* MTYPE4        mpi type for real*4                   */
c/* MTYPER        mpi type for real                     */
c/* MTYPED        mpi type for real*8                   */
c/* MTYPEI        mpi type for integer                  */
c/* MTYPEL        mpi type for logical                  */
c/* MTYPEC        mpi type for character                */
c/* MPI_SEND      either mpi_send  or mpi_ssend         */
c/* MPI_ISEND     either mpi_isend or mpi_issend        */
#if defined(NOMPIR8) /* LAM does not support mpi_real[48] */
#if defined(REAL4)
# define MTYPE4 mpi_real
# define MTYPER mpi_real
# define MTYPED mpi_double_precision
# define MTYPEI mpi_integer
# define MTYPEL mpi_logical
# define MTYPEC mpi_character
# define MTYPEI2 mpi_integer2
#else /* REAL8 */
# define MTYPE4 mpi_real
# define MTYPER mpi_double_precision
# define MTYPED mpi_double_precision
# define MTYPEI mpi_integer
# define MTYPEL mpi_logical
# define MTYPEC mpi_character
# define MTYPEI2 mpi_integer2
#endif
#else /* most MPI's allow mpi_real[48] */
#if defined(REAL4)
# define MTYPE4 mpi_real4
# define MTYPER mpi_real4
# define MTYPED mpi_real8
# define MTYPEI mpi_integer
# define MTYPEL mpi_logical
# define MTYPEC mpi_character
# define MTYPEI2 mpi_integer2
#else /* REAL8 */
# define MTYPE4 mpi_real4
# define MTYPER mpi_real8
# define MTYPED mpi_real8
# define MTYPEI mpi_integer
# define MTYPEL mpi_logical
# define MTYPEC mpi_character
# define MTYPEI2 mpi_integer2
#endif
#endif
#if defined(SSEND)
# define MPI_SEND mpi_ssend
# define MPI_ISEND mpi_issend
#else
# define MPI_SEND mpi_send
# define MPI_ISEND mpi_isend
#endif /* SSEND:else */
#endif /* MPI */

#if defined(SHMEM)
c/* SHMEM_GET4    get real*4  variables                 */
c/* SHMEM_GETR    get real    variables                 */
c/* SHMEM_GETD    get real*8  variables                 */
c/* SHMEM_GETI    get integer variables                 */
c/* SHMEM_GETL    get logical variables                 */
c/* SHMEM_GETC    get character variables               */
c/* SHMEM_MYPE    return number of this PE (0...npes-1) */
c/* SHMEM_NPES    return number of PEs                  */
#if defined(REAL4)
# define SHMEM_GET4 shmem_get32
# define SHMEM_GETR shmem_get32
# define SHMEM_GETD shmem_get64
# define SHMEM_GETI shmem_integer_get
# define SHMEM_GETL shmem_logical_get
# define SHMEM_GETC shmem_character_get
#else /* REAL8 */
# define SHMEM_GET4 shmem_get32
# define SHMEM_GETR shmem_get64
# define SHMEM_GETD shmem_get64
# define SHMEM_GETI shmem_integer_get
# define SHMEM_GETL shmem_logical_get
# define SHMEM_GETC shmem_character_get
#endif
# define SHMEM_MYPE shmem_my_pe
# define SHMEM_NPES shmem_n_pes
#endif /* SHMEM */

c
c-----------------------------------------------------------------------
c
c     auxillary routines that involve off-processor communication.
c     message passing version, contained in module mod_xc.
c
c     author:  Alan J. Wallcraft,  NRL.
c
c-----------------------------------------------------------------------
c
      subroutine xcaget(aa, a, mnflg)
      implicit none
c
      real,    intent(out)   :: aa(itdm,jtdm)
      real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) convert an entire 2-D array from tiled to non-tiled layout.
c
c  3) mnflg selects which nodes must return the array
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  4) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aa              real           output    non-tiled target array
c    a               real           input     tiled source array
c    mnflg           integer        input     node return flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
c
      integer        mpireqa(jpr),mpireqb(ipr)
#endif
      real            al,ala,at,ata
      common/xcagetr/ al(itdm,jdm),ala(itdm,jdm,jqr),
     &                at(idm*jdm),ata(idm*jdm,iqr)
      save  /xcagetr/
c
      integer i,j,l,mp,np,mnp
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 1)
        nxc = 1
      endif
#endif
c
c     gather each row of tiles onto the first tile in the row.
c
#if defined(MPI)
      if     (mproc.eq.mpe_1(nproc)) then
        do j= 1,jj
          do i= 1,i0
            al(i,j) = vland
          enddo
          do i= 1,ii
            al(i0+i,j) = a(i,j)
          enddo
          do i= i0+ii+1,itdm
            al(i,j) = vland
          enddo
        enddo
        l = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = l + 1
          call MPI_IRECV(ata(1,l),ii_pe(mp,nproc)*jj,MTYPER,
     &                   idproc(mp,nproc), 9941,
     &                   mpicomm, mpireqb(l), mpierr)
        enddo
        BARRIER
        l = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = l + 1
          call MPI_WAIT(mpireqb(l), mpistat, mpierr)
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
              al(i0_pe(mp,nproc)+i,j) = ata(i+(j-1)*ii_pe(mp,nproc),l)
            enddo
          enddo
        enddo
      else  !mproc>1
        do j= 1,jj
          do i= 1,ii
            at(i+(j-1)*ii) = a(i,j)
          enddo
        enddo
        BARRIER
        call MPI_SEND(at,ii*jj,MTYPER,
     &                idproc(mpe_1(nproc),nproc), 9941,
     &                mpicomm, mpierr)
      endif
#elif defined(SHMEM)
      do j= 1,jj
        do i= 1,ii
          at(i+(j-1)*ii) = a(i,j)
        enddo
      enddo
      BARRIER
      if     (mproc.eq.mpe_1(nproc)) then
        do j= 1,jj
          do i= 1,i0
            al(i,j) = vland
          enddo
          do i= 1,ii
            al(i0+i,j) = a(i,j)
          enddo
          do i= i0+ii+1,itdm
            al(i,j) = vland
          enddo
        enddo
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          call SHMEM_GETR(at,
     &                    at,ii_pe(mp,nproc)*jj, idproc(mp,nproc))
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
              al(i0_pe(mp,nproc)+i,j) = at(i+(j-1)*ii_pe(mp,nproc))
            enddo
          enddo
        enddo
      endif
      BARRIER
#endif /* MPI:SHMEM */
c
c     gather each row of tiles onto the output array.
c
#if defined(MPI)
      mnp = max(mnflg,1)
c
      if     (mnproc.eq.mnp) then
        if (mproc.eq.mpe_1(nproc)) then
          do j= 1,jj
            do i= 1,itdm
              aa(i,j+j0) = al(i,j)
            enddo
          enddo
        endif
        l = 0
        do np= 1,jpr
          mp = mpe_1(np)
          if     (idproc(mp,np).ne.idproc(mproc,nproc)) then
            l = l + 1
            call MPI_IRECV(ala(1,1,l),itdm*jj_pe(mp,np),MTYPER,
     &                     idproc(mp,np), 9942,
     &                     mpicomm, mpireqa(l), mpierr)
          endif
        enddo
        BARRIER
        l = 0
        do np= 1,jpr
          mp = mpe_1(np)
          if     (idproc(mp,np).ne.idproc(mproc,nproc)) then
            l = l + 1
            call MPI_WAIT(mpireqa(l), mpistat, mpierr)
            do j= 1,jj_pe(mp,np)
              do i= 1,itdm
                aa(i,j+j0_pe(mp,np)) = ala(i,j,l)
              enddo
            enddo
          endif
        enddo
      elseif (mproc.eq.mpe_1(nproc)) then
        BARRIER
        call MPI_SEND(al,itdm*jj,MTYPER,
     &                idproc1(mnp), 9942,
     &                mpicomm, mpierr)
      else
        BARRIER
      endif
c
      if     (mnflg.eq.0) then
        call mpi_bcast(aa,itdm*jtdm,MTYPER,
     &                 idproc1(1),mpicomm,mpierr)
      endif
#elif defined(SHMEM)
      if     (mnflg.eq.0 .or. mnproc.eq.mnflg) then
        if (mproc.eq.mpe_1(nproc)) then
          do j= 1,jj
            do i= 1,itdm
              aa(i,j+j0) = al(i,j)
            enddo
          enddo
        endif
        do np= 1,jpr
          mp = mpe_1(np)
          if     (idproc(mp,np).ne.idproc(mproc,nproc)) then
            call SHMEM_GETR(al,
     &                      al,itdm*jj_pe(mp,np), idproc(mp,np))
            do j= 1,jj_pe(mp,np)
              do i= 1,itdm
                aa(i,j+j0_pe(mp,np)) = al(i,j)
              enddo
            enddo
          endif
        enddo
      endif
      ! no barrier needed here because of double buffering (at and then al)
#endif /* MPI:SHMEM */
#if defined(TIMER)
c
      if     (nxc.eq. 1) then
        call xctmr1( 1)
        nxc = 0
      endif
#endif
      return
      end subroutine xcaget

      subroutine xcaget4(aa, a, mnflg)
      implicit none
c
      real*4,  intent(out)   :: aa(itdm,jtdm)
      real*4,  intent(in)    :: a(ii,jj)
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) convert an entire 2-D array from tiled to non-tiled layout.
c
c  2) Special version for zaiord and zaiowr only.
c     arrays are real*4 and tiled array has no halo.
c
c  3) mnflg selects which nodes must return the array
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  4) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aa              real           output    non-tiled target array
c    a               real           input     tiled source array
c    mnflg           integer        input     node return flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
      real            al8,ala8,at8,ata8
      common/xcagetr/ al8(itdm,jdm),ala8(itdm,jdm,jqr),
     &                at8(idm*jdm),ata8(idm*jdm,iqr)
      save  /xcagetr/
c
      real*4          al4(itdm,jdm),       at4(idm*jdm)
      equivalence    (al4(1,1),al8(1,1)), (at4(1),at8(1))
c
      integer i,j,mp,np,mnp
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 1)
        nxc = 1
      endif
#endif
c
c     gather each row of tiles onto the first tile in the row.
c
#if defined(MPI)
      if     (mproc.eq.mpe_1(nproc)) then
        do j= 1,jj
          do i= 1,i0
            al4(i,j) = vland
          enddo
          do i= 1,ii
            al4(i0+i,j) = a(i,j)
          enddo
          do i= i0+ii+1,itdm
            al4(i,j) = vland
          enddo
        enddo
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          call MPI_RECV(at4,ii_pe(mp,nproc)*jj,MTYPE4,
     &                  idproc(mp,nproc), 9941,
     &                  mpicomm, mpistat, mpierr)
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
              al4(i0_pe(mp,nproc)+i,j) = at4(i+(j-1)*ii_pe(mp,nproc))
            enddo
          enddo
        enddo
      else  !mproc>1
        call MPI_SEND(a,ii*jj,MTYPE4,
     &                idproc(mpe_1(nproc),nproc), 9941,
     &                mpicomm, mpierr)
      endif
#elif defined(SHMEM)
      do j= 1,jj
        do i= 1,ii
          at4(i+(j-1)*ii) = a(i,j)
        enddo
      enddo
      BARRIER
      if     (mproc.eq.mpe_1(nproc)) then
        do j= 1,jj
          do i= 1,i0
            al4(i,j) = vland
          enddo
          do i= 1,ii
            al4(i0+i,j) = a(i,j)
          enddo
          do i= i0+ii+1,itdm
            al4(i,j) = vland
          enddo
        enddo
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          call SHMEM_GET4(at4,
     &                    at4,ii_pe(mp,nproc)*jj, idproc(mp,nproc))
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
              al4(i0_pe(mp,nproc)+i,j) = at4(i+(j-1)*ii_pe(mp,nproc))
            enddo
          enddo
        enddo
      endif
      BARRIER
#endif /* MPI:SHMEM */
c
c     gather each row of tiles onto the output array.
c
#if defined(MPI)
      mnp = max(mnflg,1)
c
      if     (mnproc.eq.mnp) then
        if (mproc.eq.mpe_1(nproc)) then
          do j= 1,jj
            do i= 1,itdm
              aa(i,j+j0) = al4(i,j)
            enddo
          enddo
        endif
        do np= 1,jpr
          mp = mpe_1(np)
          if     (idproc(mp,np).ne.idproc(mproc,nproc)) then
            call MPI_RECV(al4,itdm*jj_pe(mp,np),MTYPE4,
     &                    idproc(mp,np), 9942,
     &                    mpicomm, mpistat, mpierr)
            do j= 1,jj_pe(mp,np)
              do i= 1,itdm
                aa(i,j+j0_pe(mp,np)) = al4(i,j)
              enddo
            enddo
          endif
        enddo
      elseif (mproc.eq.mpe_1(nproc)) then
        call MPI_SEND(al4,itdm*jj,MTYPE4,
     &                idproc1(mnp), 9942,
     &                mpicomm, mpierr)
      endif
c
      if     (mnflg.eq.0) then
        call mpi_bcast(aa,itdm*jtdm,MTYPE4,
     &                 idproc1(1),mpicomm,mpierr)
      endif
#elif defined(SHMEM)
      if     (mnflg.eq.0 .or. mnproc.eq.mnflg) then
        if (mproc.eq.mpe_1(nproc)) then
          do j= 1,jj
            do i= 1,itdm
              aa(i,j+j0) = al4(i,j)
            enddo
          enddo
        endif
        do np= 1,jpr
          mp = mpe_1(np)
          if     (idproc(mp,np).ne.idproc(mproc,nproc)) then
            call SHMEM_GET4(al4,
     &                      al4,itdm*jj_pe(mp,np), idproc(mp,np))
            do j= 1,jj_pe(mp,np)
              do i= 1,itdm
                aa(i,j+j0_pe(mp,np)) = al4(i,j)
              enddo
            enddo
          endif
        enddo
      endif
      ! no barrier needed here because of double buffering (at and then al)
#endif /* MPI:SHMEM */
#if defined(TIMER)
c
      if     (nxc.eq. 1) then
        call xctmr1( 1)
        nxc = 0
      endif
#endif
      return
      end subroutine xcaget4

      subroutine xcaput(aa, a, mnflg)
      implicit none
c
      real,    intent(inout) :: aa(itdm,jtdm)
      real,    intent(out)   :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) convert an entire 2-D array from non-tiled to tiled layout.
c
c  3) mnflg selects which nodes must contain the non-tiled array
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c     if mnflg.ne.0 the array aa may be broadcast to all nodes,
c      so aa must exist and be overwritable on all nodes.
c
c  4) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aa              real           input     non-tiled source array
c    a               real           output    tiled target array
c    mnflg           integer        input     node source flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
c
      integer        mpireqa(jpr),mpireqb(ipr)
#endif
      real            al,ala,at,ata
      common/xcagetr/ al(itdm,jdm),ala(itdm,jdm,jqr),
     &                at(idm*jdm),ata(idm*jdm,iqr)
      save  /xcagetr/
c
      integer i,j,mp,np,mnp
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        BARRIER
        call xctmr0( 4)
        nxc = 4
      endif
#endif
*
*     if     (mnproc.eq.1) then
*       write(lp,'(a,g20.10)') 'xcaput - vland =',vland
*     endif
c
c     use xclput for now,
c     this is slow for mnflg.ne.0, but easy to implement.
c
      if     (mnflg.ne.0) then
c       "broadcast" row sections of aa to all processors in the row.
        if     (mnproc.ne.mnflg) then
          aa(:,:) = vland
        endif
#if defined(MPI)
        if     (mnproc.eq.mnflg) then
          j = 0
          do np= 1,jpr
            mp = mpe_1(np)
            if     (np.eq.nproc .and. mp.eq.mproc) then
              al(:,1:jj) = aa(:,j0+1:j0+jj)
            else
              j = j + 1
              call MPI_ISEND(aa(1,j0_pe(mp,np)+1),
     &                      itdm*jj_pe(mp,np),MTYPER,
     &                      idproc(mp,np), 9951,
     &                      mpicomm, mpireqa(j), mpierr)
            endif
          enddo
          call mpi_waitall( j, mpireqa, mpistat, mpierr)
        elseif (mproc.eq.mpe_1(nproc)) then
          call MPI_RECV(al,itdm*jj,MTYPER,
     &                  idproc1(mnflg), 9951,
     &                  mpicomm, mpistat, mpierr)
        endif
c
        if     (mproc.eq.mpe_1(nproc)) then
          i = 0
          do mp= mpe_1(nproc)+1,mpe_e(nproc)
            i = i + 1
            call MPI_ISEND(al,itdm*jj,MTYPER,
     &                    idproc(mp,nproc), 9952,
     &                    mpicomm, mpireqb(i), mpierr)
          enddo
          call mpi_waitall( i, mpireqb, mpistat, mpierr)
        else
          call MPI_RECV(al,itdm*jj,MTYPER,
     &                  idproc(mpe_1(nproc),nproc), 9952,
     &                  mpicomm, mpistat, mpierr)
        endif
c
        aa(:,j0+1:j0+jj) = al(:,1:jj)
#elif defined(SHMEM)
c       assume aa is in common.
        BARRIER
        if     (mnproc.ne.mnflg) then
          do j= 1,jj
            call SHMEM_GETR(aa(1,j0+j),
     &                      aa(1,j0+j),itdm,idproc1(mnflg))
          enddo
        endif
        BARRIER
#endif
      endif
      do j= 1,jtdm
        call xclput(aa(1,j),itdm, a, 1,j,1,0)
      enddo
#if defined(TIMER)
c
      if     (nxc.eq. 4) then
        call xctmr1( 4)
        nxc = 0
      endif
#endif
      return
      end subroutine xcaput

      subroutine xcaput4(aa, a, mnflg)
      implicit none
c
      real*4,  intent(inout) :: aa(itdm,jtdm)
      real*4,  intent(out)   :: a(ii,jj)
      integer, intent(in)    :: mnflg
c
c**********
c*
c  1) convert an entire 2-D array from non-tiled to tiled layout.
c
c  2) Special version for zaiord and zaiowr only.
c     arrays are real*4 and tiled array has no halo.
c
c  3) mnflg selects which nodes must contain the non-tiled array
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c     if mnflg.ne.0 the array aa may be broadcast to all nodes,
c      so aa must exist and be overwritable on all nodes.
c
c  4) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aa              real           input     non-tiled source array
c    a               real           output    tiled target array
c    mnflg           integer        input     node source flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
c
      integer        mpireqa(jpr),mpireqb(ipr)
#endif
      real            al8,ala8,at8,ata8
      common/xcagetr/ al8(itdm,jdm),ala8(itdm,jdm,jqr),
     &                at8(idm*jdm),ata8(idm*jdm,iqr)
      save  /xcagetr/
c                                               
      real*4          al4(itdm,jdm),       at4(idm*jdm)
      equivalence    (al4(1,1),al8(1,1)), (at4(1),at8(1))
c
      integer i,j,mp,np,mnp
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        BARRIER
        call xctmr0( 4)
        nxc = 4
      endif
#endif
*
*     if     (mnproc.eq.1) then
*       write(lp,'(a,g20.10)') 'xcaput - vland =',vland
*     endif
c
c     use xclput for now,
c     this is slow for mnflg.ne.0, but easy to implement.
c
      if     (mnflg.ne.0) then
c       "broadcast" row sections of aa to all processors in the row.
        if     (mnproc.ne.mnflg) then
          aa(:,:) = vland
        endif
#if defined(MPI)
        if     (mnproc.eq.mnflg) then
          j = 0
          do np= 1,jpr
            mp = mpe_1(np)
            if     (np.eq.nproc .and. mp.eq.mproc) then
              al4(:,1:jj) = aa(:,j0+1:j0+jj)
            else
              j = j + 1
              call MPI_ISEND(aa(1,j0_pe(mp,np)+1),
     &                      itdm*jj_pe(mp,np),MTYPE4,
     &                      idproc(mp,np), 9951,
     &                      mpicomm, mpireqa(j), mpierr)
            endif
          enddo
          call mpi_waitall( j, mpireqa, mpistat, mpierr)
        elseif (mproc.eq.mpe_1(nproc)) then
          call MPI_RECV(al4,itdm*jj,MTYPE4,
     &                  idproc1(mnflg), 9951,
     &                  mpicomm, mpistat, mpierr)
        endif
c
        if     (mproc.eq.mpe_1(nproc)) then
          i = 0
          do mp= mpe_1(nproc)+1,mpe_e(nproc)
            i = i + 1
            call MPI_ISEND(al4,itdm*jj,MTYPE4,
     &                    idproc(mp,nproc), 9952,
     &                    mpicomm, mpireqb(i), mpierr)
          enddo
          call mpi_waitall( i, mpireqb, mpistat, mpierr)
        else
          call MPI_RECV(al4,itdm*jj,MTYPE4,
     &                  idproc(mpe_1(nproc),nproc), 9952,
     &                  mpicomm, mpistat, mpierr)
        endif
c
        aa(:,j0+1:j0+jj) = al4(:,1:jj)
#elif defined(SHMEM)
c       assume aa is in common.
        BARRIER
        if     (mnproc.ne.mnflg) then
          do j= 1,jj
            call SHMEM_GET4(aa(1,j0+j),
     &                      aa(1,j0+j),itdm,idproc1(mnflg))
          enddo
        endif
        BARRIER
#endif
      endif
      do j= 1,jtdm
        call xclput4(aa(1,j),itdm, a, 1,j,1,0)
      enddo
#if defined(TIMER)
c
      if     (nxc.eq. 4) then
        call xctmr1( 4)
        nxc = 0
      endif
#endif
      return
      end subroutine xcaput4

      subroutine xceget(aelem, a, ia,ja)
      implicit none
c
      real,    intent(out)   :: aelem
      real,    intent(in)    :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: ia,ja
c
c**********
c*
c  1) find the value of a(ia,ja) on the non-tiled 2-D grid.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aelem           real           output    required element
c    a               real           input     source array
c    ia              integer        input     1st index into a
c    ja              integer        input     2nd index into a
c
c  3) the global variable vland is returned when outside active tiles.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
c     double buffer to reduce the number of required barriers.
c
      real            elem(0:1)
      common/xcegetr/ elem
      save  /xcegetr/
c
      integer, save :: kdb = 0
c
      integer i,j,mp,np
#if defined(TIMER)
c
*     if     (nxc.eq.0) then
*       call xctmr0( 2)
*       nxc = 2
*     endif
#endif
c
      kdb = mod(kdb+1,2)  ! the least recently used of the two buffers
c
c     find the host tile.
c
      np = npe_j(ja)
      mp = mpe_i(ia,np)
c
      if      (mp.le.0) then
c
c       no tile.
c
        elem(kdb) = vland
      elseif  (mp.eq.mproc .and. np.eq.nproc) then
c
c       this tile.
c
        i = ia - i0
        j = ja - j0
c
        elem(kdb) = a(i,j)
#if defined(MPI)
        call mpi_bcast(elem(kdb),1,MTYPER,
     &                 idproc(mp,np),mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
#endif
      else
c
c       another tile.
c
#if defined(MPI)
        call mpi_bcast(elem(kdb),1,MTYPER,
     &                 idproc(mp,np),mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        call SHMEM_GETR(elem(kdb),
     &                  elem(kdb),1,idproc(mp,np))
        ! no barrier needed here because of double buffering
#endif
      endif
      aelem = elem(kdb)
#if defined(TIMER)
c
*     if     (nxc.eq. 2) then
*       call xctmr1( 2)
*       nxc = 0
*     endif
#endif
      return
      end subroutine xceget

      subroutine xceput(aelem, a, ia,ja)
      implicit none
c
      integer, intent(in)    :: ia,ja
      real,    intent(in)    :: aelem
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) fill a single element in the non-tiled 2-D grid.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aelem           real           input     element value
c    a               real           in/out    target array
c    ia              integer        input     1st index into a
c    ja              integer        input     2nd index into a
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer mp,np
#if defined(TIMER)
c
*     if     (nxc.eq.0) then
*       call xctmr0( 5)
*       nxc = 5
*     endif
#endif
      if     (i0.lt.ia .and. ia.le.i0+ii .and.
     &        j0.lt.ja .and. ja.le.j0+jj      ) then
c
c       this tile.
c
        a(ia-i0,ja-j0) = aelem
      endif
#if defined(TIMER)
c
*     if     (nxc.eq. 5) then
*       call xctmr1( 5)
*       nxc = 0
*     endif
#endif
      return
      end subroutine xceput

      subroutine xcgetc(iline)
      implicit none
c
      integer, intent(inout) :: iline(81)
c
c**********
c*
c  1) machine specific routine for broadcasting iline.
c
c  2) only use in zagetc (hence the name).
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
c     broadcast to all processors
c
#if defined(MPI)
      call mpi_bcast(iline,81,MTYPEI,
     &               idproc1(1),mpicomm,mpierr)
#elif defined(SHMEM)
      BARRIER
      if     (mnproc.ne.1) then
        call SHMEM_GETI(iline,
     &                  iline,81, idproc1(1))
      endif
      ! no barrier needed here because zagetc is using two buffers
#endif
      return
      end subroutine xcgetc

      subroutine xchalt(cerror)
      implicit none
c
      character*(*), intent(in) :: cerror
c
c**********
c*
c  1) stop all processes.
c
c  2) only one processes need call this routine, i.e. it is for
c      emergency stops.  use 'xcstop' for ordinary stops called
c      by all processes.
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    cerror          char*(*)       input     error message
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
c     message passing version.
c
      if     (cerror.ne.' ') then
        write(lp,*) '**************************************************'
        write(lp,*) cerror
      endif
      write(lp,*) '**************************************************'
      write(lp,*) 'XCHALT CALLED ON PROC = ',mnproc,mproc,nproc
      write(lp,*) '**************************************************'
      call flush(lp)
c
#if defined(MPI)
      if (mpicomm.eq.mpicom_external) then
        call external_abort(cerror)
      else
        call mpi_abort(mpicomm,9,mpierr)
      endif
#else
      call abort()
#endif
      stop '(xchalt)'
      end subroutine xchalt

      subroutine xclget(aline,nl, a, i1,j1,iinc,jinc, mnflg)
      implicit none
c
      integer, intent(in)    ::  nl,i1,j1,iinc,jinc,mnflg
      real,    intent(out)   ::  aline(nl)
      real,    intent(in)    ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) extract a line of elements from the non-tiled 2-D grid.
c
c  2) aline(i) == aa(i1+iinc*(i-1),j1+jinc*(i-1)), for i=1...nl.
c     where aa is the non-tiled representation of a, and
c     iinc and jinc can each be -1, 0, or +1.
c
c  3) mnflg selects which nodes must return the line
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  4) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aline           real           output    required line of elements
c    nl              integer        input     dimension of aline
c    a               real           input     source array
c    i1              integer        input     1st index into a
c    j1              integer        input     2nd index into a
c    iinc            integer        input     1st index increment
c    jinc            integer        input     2nd index increment
c    mnflg           integer        input     node return flag
c
c  5) the global variable vland is returned when outside active tiles.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
c     pad al to guard against false sharing.
c     double buffer to reduce the number of required barriers.
c
      real            al
      common/xclgetr/ al(-47:itdm+jtdm+48,0:1)
      save  /xclgetr/
c
      integer, save :: kdb = 0
c
      real    dummy
      integer i1n,iif,iil,j1n,jjf,jjl,l,lx0,nxl,mp,np
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 3)
        nxc = 3
      endif
#endif
c
      kdb = mod(kdb+1,2)  ! the least recently used of the two buffers
c
      if     ((jinc.eq.0 .and. iinc.eq.1) .or. nl.eq.1) then
c
c ---   horizontal forward line.
c
        al(1:nl,kdb) = vland
c
        np  = npe_j(j1)
        do mp= mpe_1(np),mpe_e(np)
          iif = max(i1,     i0_pe(mp,np)+1)
          iil = min(i1+nl-1,i0_pe(mp,np)+ii_pe(mp,np))
          lx0 = iif - i1
          nxl = iil - iif + 1
          if     (nxl.le.0) then
            cycle  ! no elements from tile (mp,np)
          endif
c
          if      (mp.eq.mproc .and. np.eq.nproc) then
c
c           this tile.
c
            do l= lx0+1,lx0+nxl
              al(l,kdb) = a(i1+l-1-i0,j1-j0)
            enddo
#if defined(MPI)
            if     (mnflg.eq.0) then
              call mpi_bcast(al(lx0+1,kdb),nxl,MTYPER,
     &                       idproc(mp,np),mpicomm,mpierr)
            elseif (mnflg.ne.mnproc) then
              call MPI_SEND(al(lx0+1,kdb),nxl,MTYPER,
     &                      idproc1(mnflg), 9931,
     &                      mpicomm, mpierr)
            endif
          else
c
c           another tile.
c
            if     (mnflg.eq.0) then
              call mpi_bcast(al(lx0+1,kdb),nxl,MTYPER,
     &                       idproc(mp,np),mpicomm,mpierr)
            elseif (mnflg.eq.mnproc) then
              call MPI_RECV(al(lx0+1,kdb),nxl,MTYPER,
     &                      idproc(mp,np), 9931,
     &                      mpicomm, mpistat, mpierr)
            endif
#endif /* MPI */
          endif
c
          if     (lx0+nxl.eq.nl) then
            exit
          endif
        enddo  ! np=1,jpr
#if defined(SHMEM)
c
c       spliting process into two phases saves on barriers.
c
        BARRIER
        do mp= mpe_1(np),mpe_e(np)
          iif = max(i1,     i0_pe(mp,np)+1)
          iil = min(i1+nl-1,i0_pe(mp,np)+ii_pe(mp,np))
          lx0 = iif - i1
          nxl = iil - iif + 1
          if     (nxl.le.0) then
            cycle  ! no elements from tile (mp,np)
          endif
c
          if      (mp.eq.mproc .and. np.eq.nproc) then
c
c           nothing to do here (see 1st phase, above).
c
          else
c
c           another tile.
c
            if     (mnflg.eq.0 .or. mnflg.eq.mnproc) then
              call SHMEM_GETR(al(lx0+1,kdb),
     &                        al(lx0+1,kdb),nxl,idproc(mp,np))
            endif
          endif
c
          if     (lx0+nxl.eq.nl) then
            exit
          endif
        enddo  ! np=1,jpr
        ! no barrier needed here because of double buffering
#endif /* SHMEM */
c
        if     (mnflg.eq.0 .or. mnflg.eq.mnproc) then
          aline(1:nl) = al(1:nl,kdb)
        endif
      elseif (iinc.eq.0 .and. jinc.eq.1) then
c
c ---   vertical forward line.
c
        al(1:nl,kdb) = vland
c
        do np= 1,jpr
          mp = mpe_i(i1,np)
          if     (mp.le.0) then
            cycle  ! an all-land column
          endif
          jjf = max(j1,     j0_pe(mp,np)+1)
          jjl = min(j1+nl-1,j0_pe(mp,np)+jj_pe(mp,np))
          lx0 = jjf - j1
          nxl = jjl - jjf + 1
          if     (nxl.le.0) then
            cycle  ! no elements from tile (mp,np)
          endif
c
          if      (mp.eq.mproc .and. np.eq.nproc) then
c
c           this tile.
c
            do l= lx0+1,lx0+nxl
              al(l,kdb) = a(i1-i0,j1+l-1-j0)
            enddo
#if defined(MPI)
            if     (mnflg.eq.0) then
              call mpi_bcast(al(lx0+1,kdb),nxl,MTYPER,
     &                       idproc(mp,np),mpicomm,mpierr)
            elseif (mnflg.ne.mnproc) then
              call MPI_SEND(al(lx0+1,kdb),nxl,MTYPER,
     &                      idproc1(mnflg), 9931,
     &                      mpicomm, mpierr)
            endif
          else
c
c           another tile.
c
            if     (mnflg.eq.0) then
              call mpi_bcast(al(lx0+1,kdb),nxl,MTYPER,
     &                       idproc(mp,np),mpicomm,mpierr)
            elseif (mnflg.eq.mnproc) then
              call MPI_RECV(al(lx0+1,kdb),nxl,MTYPER,
     &                      idproc(mp,np), 9931,
     &                      mpicomm, mpistat, mpierr)
            endif
#endif /* MPI */
          endif
c
          if     (lx0+nxl.eq.nl) then
            exit
          endif
        enddo  ! np=1,jpr
#if defined(SHMEM)
c
c       spliting process into two phases saves on barriers.
c
        BARRIER
        do np= 1,jpr
          mp = mpe_i(i1,np)
          if     (mp.le.0) then
            cycle  ! an all-land column
          endif
          jjf = max(j1,     j0_pe(mp,np)+1)
          jjl = min(j1+nl-1,j0_pe(mp,np)+jj_pe(mp,np))
          lx0 = jjf - j1
          nxl = jjl - jjf + 1
          if     (nxl.le.0) then
            cycle  ! no elements from tile (mp,np)
          endif
c
          if      (mp.eq.mproc .and. np.eq.nproc) then
c
c           nothing to do here (see 1st phase, above).
c
          else
c
c           another tile.
c
            if     (mnflg.eq.0 .or. mnflg.eq.mnproc) then
              call SHMEM_GETR(al(lx0+1,kdb),
     &                        al(lx0+1,kdb),nxl,idproc(mp,np))
            endif
          endif
c
          if     (lx0+nxl.eq.nl) then
            exit
          endif
        enddo  ! np=1,jpr
        ! no barrier needed here because of double buffering
#endif /* SHMEM */
c
        if     (mnflg.eq.0 .or. mnflg.eq.mnproc) then
          aline(1:nl) = al(1:nl,kdb)
        endif
      else
c
c       diagonal and reversing lines - repeatedly call xceget.
c       this always works, but is very slow.
c
        do l= 1,nl
          if     (mnflg.eq.0 .or. mnflg.eq.mnproc) then
            call xceget(aline(l), a, i1+iinc*(l-1),j1+jinc*(l-1))
          else
            call xceget(dummy,    a, i1+iinc*(l-1),j1+jinc*(l-1))
          endif
        enddo
      endif
#if defined(TIMER)
c
      if     (nxc.eq. 3) then
        call xctmr1( 3)
        nxc = 0
      endif
#endif
      return
      end subroutine xclget

      subroutine xclput(aline,nl, a, i1,j1,iinc,jinc)
      implicit none
c
      integer, intent(in)    ::  nl,i1,j1,iinc,jinc
      real,    intent(in)    ::  aline(nl)
      real,    intent(inout) ::  a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) fill a line of elements in the non-tiled 2-D grid.
c
c  2) aline(i) == aa(i1+i1*(i-1),j1+j1*(i-1)), for i=1...nl.
c     where aa is the non-tiled representation of a, and
c     one of iinc and jinc must be 0, and the other must be 1.
c     also updates the halo.
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aline           real           input     line of element values
c    nl              integer        input     dimension of aline
c    a               real           in/out    target array
c    i1              integer        input     1st index into a
c    j1              integer        input     2nd index into a
c    iinc            integer        input     1st index increment
c    jinc            integer        input     2nd index increment
c*
c**********
c
      integer i,j
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 5)
        nxc = 5
      endif
#endif
c
      if     (jinc.eq.0) then
        if     (j1-j0.ge.1-nbdy .and. j1-j0.le.jj+nbdy) then
          do i= max(1-nbdy,i1-i0),min(i1-i0+nl-1,ii+nbdy)
            a(i,j1-j0) = aline(i+i0-i1+1)
          enddo
          if     (nreg.ne.0 .and.
     &            i0.eq.0 .and. i1+nl-1.ge.itdm-nbdy+1) then  ! periodic
            do i= max(itdm-nbdy+1,i1),i1+nl-1
              a(i-itdm,j1-j0) = aline(i)
            enddo
          endif
          if     (nreg.ne.0 .and.
     &            i0+ii.eq.itdm .and. i1.le.nbdy) then        ! periodic
            do i= i1,min(nbdy,i1+nl-1)
              a(ii+i,j1-j0) = aline(i)
            enddo
          endif
        endif
      elseif (iinc.eq.0) then
        if     (i1-i0.ge.1-nbdy .and. i1-i0.le.ii+nbdy) then
          do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
            a(i1-i0,j) = aline(j+j0-j1+1)
          enddo
        endif
        if     (nreg.ne.0 .and.
     &          i0.eq.0 .and. i1.ge.itdm-nbdy+1) then       ! periodic
          do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
            a(i1-itdm,j) = aline(j+j0-j1+1)
          enddo
        endif
        if     (nreg.ne.0 .and.
     &          i0+ii.eq.itdm .and. i1.le.nbdy) then        ! periodic
          do j= max(1-nbdy,j1-j0),min(j1-j0+nl-1,jj+nbdy)
            a(ii+i1,j) = aline(j+j0-j1+1)
          enddo
        endif
      endif
#if defined(TIMER)
c
      if     (nxc.eq. 5) then
        call xctmr1( 5)
        nxc = 0
      endif
#endif
      return
      end subroutine xclput

      subroutine xclput4(aline,nl, a, i1,j1,iinc,jinc)
      implicit none
c
      integer, intent(in)    ::  nl,i1,j1,iinc,jinc
      real*4,  intent(in)    ::  aline(nl)
      real*4,  intent(inout) ::  a(ii,jj)
c
c**********
c*
c  1) fill a line of elements in the non-tiled 2-D grid.
c     Special version for xcaput4 only.
c
c  2) aline(i) == aa(i1+i1*(i-1),j1+j1*(i-1)), for i=1...nl.
c     where aa is the non-tiled representation of a, and
c     one of iinc and jinc must be 0, and the other must be 1.
c     also updates the halo.
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aline           real           input     line of element values
c    nl              integer        input     dimension of aline
c    a               real           in/out    target array
c    i1              integer        input     1st index into a
c    j1              integer        input     2nd index into a
c    iinc            integer        input     1st index increment
c    jinc            integer        input     2nd index increment
c*
c**********
c
      integer i,j
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 5)
        nxc = 5
      endif
#endif
c
      if     (jinc.eq.0) then
        if     (j1-j0.ge.1 .and. j1-j0.le.jj) then
          do i= max(1,i1-i0),min(i1-i0+nl-1,ii)
            a(i,j1-j0) = aline(i+i0-i1+1)
          enddo
          if     (nreg.ne.0 .and.
     &            i0.eq.0 .and. i1+nl-1.ge.itdm+1) then  ! periodic
            do i= max(itdm+1,i1),i1+nl-1
              a(i-itdm,j1-j0) = aline(i)
            enddo
          endif
        endif
      elseif (iinc.eq.0) then
        if     (i1-i0.ge.1 .and. i1-i0.le.ii) then
          do j= max(1,j1-j0),min(j1-j0+nl-1,jj)
            a(i1-i0,j) = aline(j+j0-j1+1)
          enddo
        endif
        if     (nreg.ne.0 .and.
     &          i0.eq.0 .and. i1.ge.itdm+1) then       ! periodic
          do j= max(1,j1-j0),min(j1-j0+nl-1,jj)
            a(i1-itdm,j) = aline(j+j0-j1+1)
          enddo
        endif
      endif
#if defined(TIMER)
c
      if     (nxc.eq. 5) then
        call xctmr1( 5)
        nxc = 0
      endif
#endif
      return
      end subroutine xclput4

      subroutine xcbcst_i0(ia)
      implicit none
c
      integer, intent(inout) :: ia
c
c**********
c*
c  1) broadcast integer scalar from tile mnproc = 1 to all tiles
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ia              integer        in/out    target variable
c*
c**********
c
      integer ia1(1)
c
      ia1(1) = ia
      call xcbcst_i1(ia1)
      ia = ia1(1)
      return
      end subroutine xcbcst_i0

      subroutine xcbcst_i1(ia)
      implicit none
c
      integer, intent(inout) :: ia(:)
c
c**********
c*
c  1) broadcast integer array from tile mnproc = 1 to all tiles
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ia              integer        in/out    target array
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      integer        ib
      common/xcbcsi/ ib(nmax)
      save  /xcbcsi/
c
      integer i,is0,isl,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 8)
        nxc = 8
      endif
#endif
c
c     stripmine ia.
c
      n = size(ia)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          ib(i) = ia(is0+i)
        enddo
c
#if defined(MPI)
        call mpi_bcast(ib,nn,MTYPEI,
     &                 idproc1(1),mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.ne.1) then
          call SHMEM_GETI(ib,ib,nn,idproc1(1))
        endif
        BARRIER
#endif
        do i= 1,nn
          ia(is0+i) = ib(i)
        enddo
      enddo  ! stripmine loop
#if defined(TIMER)
c
      if     (nxc.eq.8) then
        call xctmr1( 8)
        nxc = 0
      endif
#endif
      return
      end subroutine xcbcst_i1

      subroutine xcbcst_r0(ra)
      implicit none
c
      real, intent(inout) :: ra
c
c**********
c*
c  1) broadcast real scalar from tile mnproc = 1 to all tiles
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ra              real           in/out    target variable
c*
c**********
c
      real ra1(1)
c
      ra1(1) = ra
      call xcbcst_r1(ra1)
      ra = ra1(1)
      return
      end subroutine xcbcst_r0

      subroutine xcbcst_r1(ra)
      implicit none
c
      real, intent(inout) :: ra(:)
c
c**********
c*
c  1) broadcast real array from tile mnproc = 1 to all tiles
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ra              real           in/out    target array
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      real           rb
      common/xcbcsr/ rb(nmax)
      save  /xcbcsr/
c
      integer i,is0,isl,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 8)
        nxc = 8
      endif
#endif
c
c     stripmine ra.
c
      n = size(ra)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          rb(i) = ra(is0+i)
        enddo
c
#if defined(MPI)
        call mpi_bcast(rb,nn,MTYPED,
     &                 idproc1(1),mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.ne.1) then
          call SHMEM_GETD(rb,rb,nn,idproc1(1))
        endif
        BARRIER
#endif
        do i= 1,nn
          ra(is0+i) = rb(i)
        enddo
      enddo  ! stripmine loop
#if defined(TIMER)
c
      if     (nxc.eq.8) then
        call xctmr1( 8)
        nxc = 0
      endif
#endif
      return
      end subroutine xcbcst_r1

      subroutine xcbcst_l0(la)
      implicit none
c
      logical, intent(inout) :: la
c
c**********
c*
c  1) broadcast logical from tile mnproc = 1 to all tiles
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    la              logical        in/out    target variable
c*
c**********
c
      logical la1(1)
c
      la1(1) = la
      call xcbcst_l1(la1)
      la = la1(1)
      return
      end subroutine xcbcst_l0

      subroutine xcbcst_l1(la)
      implicit none
c
      logical, intent(inout) :: la(:)
c
c**********
c*
c  1) broadcast logical array from tile mnproc = 1 to all tiles
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    la              logical        in/out    target array
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      logical        lb
      common/xcbcsl/ lb(nmax)
      save  /xcbcsl/
c
      integer i,is0,isl,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 8)
        nxc = 8
      endif
#endif
c
c     stripmine la.
c
      n = size(la)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          lb(i) = la(is0+i)
        enddo
c
#if defined(MPI)
        call mpi_bcast(lb,nn,MTYPEL,
     &                 idproc1(1),mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.ne.1) then
          call SHMEM_GETL(lb,lb,nn,idproc1(1))
        endif
        BARRIER
#endif
        do i= 1,nn
          la(is0+i) = lb(i)
        enddo
      enddo  ! stripmine loop
#if defined(TIMER)
c
      if     (nxc.eq.8) then
        call xctmr1( 8)
        nxc = 0
      endif
#endif
      return
      end subroutine xcbcst_l1

      subroutine xcbcst_c(ca)
      implicit none
c
      character*(*), intent(inout) :: ca
c
c**********
c*
c  1) broadcast character string from tile mnproc = 1 to all tiles
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ca              character      in/out    target array
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      character cb*(nmax)
      common/xcbcsc/ cb
      save  /xcbcsc/
c
      integer i,is0,isl,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 8)
        nxc = 8
      endif
#endif
c
c     stripmine ca.
c
      n = len(ca)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        cb(1:nn) = ca(is0+1:is0+nn)
c
#if defined(MPI)
        call mpi_bcast(cb,nn,MTYPEC,
     &                 idproc1(1),mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.ne.1) then
          call SHMEM_GETC(cb,cb,nn,idproc1(1))
        endif
        BARRIER
#endif
        ca(is0+1:is0+nn) = cb(1:nn)
      enddo  ! stripmine loop
#if defined(TIMER)
c
      if     (nxc.eq.8) then
        call xctmr1( 8)
        nxc = 0
      endif
#endif
      return
      end subroutine xcbcst_c

      subroutine xcmax_i0(ia)
      implicit none
c
      integer, intent(inout) :: ia
c
c**********
c*
c  1) replace integer scalar a with its element-wise maximum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ia              integer        in/out    target variable
c*
c**********
c
      integer ia1(1)
c
      ia1(1) = ia
      call xcmax_i1(ia1)
      ia = ia1(1)
      return
      end subroutine xcmax_i0

      subroutine xcmax_i1(ia)
      implicit none
c
      integer, intent(inout) :: ia(:)
c
c**********
c*
c  1) replace integer array a with its element-wise maximum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ia              integer        in/out    target array
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      integer        ib,ic
      common/xcmaxi/ ib(nmax),ic(nmax)
      save  /xcmaxi/
c
      integer i,is0,isl,mn,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 9)
        nxc = 9
      endif
#endif
c
c     stripmine ia.
c
      n = size(ia)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          ib(i) = ia(is0+i)
        enddo
c
#if defined(MPI)
        call mpi_allreduce(ib,ic,nn,MTYPEI,mpi_max,
     &                     mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.eq.1) then
c         form global maximum on 1st processor
          do i= 1,nn
            ic(i) = ib(i)
          enddo
          do mn= 2,ijpr
            call SHMEM_GETI(ib,ib,nn, idproc1(mn))
            do i= 1,nn
              ic(i) = max(ic(i),ib(i))
            enddo
          enddo
          BARRIER
        else
c         get global maximum from 1st processor
          BARRIER
          call SHMEM_GETI(ic,ic,nn, idproc1(1))
        endif
        ! no barrier needed here because using two buffers (ib and ic)
#endif
        do i= 1,nn
          ia(is0+i) = ic(i)
        enddo
      enddo  ! stripmine loop
#if defined(TIMER)
c
      if     (nxc.eq.9) then
        call xctmr1( 9)
        nxc = 0
      endif
#endif
      return
      end subroutine xcmax_i1

      subroutine xcmax_r0(ra)
      implicit none
c
      real, intent(inout) :: ra
c
c**********
c*
c  1) replace scalar a with its element-wise maximum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ra              real           in/out    target variable
c*
c**********
c
      real ra1(1)
c
      ra1(1) = ra
      call xcmax_r1(ra1)
      ra = ra1(1)
      return
      end subroutine xcmax_r0

      subroutine xcmax_r1(ra)
      implicit none
c
      real, intent(inout) :: ra(:)
c
c**********
c*
c  1) replace array a with its element-wise maximum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ra              real           in/out    target array
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      real           rb,rc
      common/xcmaxr/ rb(nmax),rc(nmax)
      save  /xcmaxr/
c
      integer i,is0,isl,mn,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 9)
        nxc = 9
      endif
#endif
c
c     stripmine ra.
c
      n = size(ra)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          rb(i) = ra(is0+i)
        enddo
c
#if defined(MPI)
        call mpi_allreduce(rb,rc,nn,MTYPER,mpi_max,
     &                     mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.eq.1) then
c         form global maximum on 1st processor
          do i= 1,nn
            rc(i) = rb(i)
          enddo
          do mn= 2,ijpr
            call SHMEM_GETR(rb,rb,nn, idproc1(mn))
            do i= 1,nn
              rc(i) = max(rc(i),rb(i))
            enddo
          enddo
          BARRIER
        else
c         get global maximum from 1st processor
          BARRIER
          call SHMEM_GETR(rc,rc,nn, idproc1(1))
        endif
        ! no barrier needed here because using two buffers (b and c)
#endif
        do i= 1,nn
          ra(is0+i) = rc(i)
        enddo
      enddo  ! stripmine loop
#if defined(TIMER)
c
      if     (nxc.eq.9) then
        call xctmr1( 9)
        nxc = 0
      endif
#endif
      return
      end subroutine xcmax_r1

      subroutine xcmin_i0(ia)
      implicit none
c
      integer, intent(inout) :: ia
c
c**********
c*
c  1) replace integer scalar a with its element-wise minimum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ia              integer        in/out    target variable
c*
c**********
c
      integer ia1(1)
c
      ia1(1) = ia
      call xcmin_i1(ia1)
      ia = ia1(1)
      return
      end subroutine xcmin_i0

      subroutine xcmin_i1(ia)
      implicit none
c
      integer, intent(inout) :: ia(:)
c
c**********
c*
c  1) replace integer array a with its element-wise minimum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ia              integer        in/out    target array
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      integer        ib,ic
      common/xcmaxi/ ib(nmax),ic(nmax)
      save  /xcmaxi/
c
      integer i,is0,isl,mn,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0(10)
        nxc = 10
      endif
#endif
c
c     stripmine ia.
c
      n = size(ia)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          ib(i) = ia(is0+i)
        enddo
c
#if defined(MPI)
        call mpi_allreduce(ib,ic,nn,MTYPEI,mpi_min,
     &                     mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.eq.1) then
c         form global minimum on 1st processor
          do i= 1,nn
            ic(i) = ib(i)
          enddo
          do mn= 2,ijpr
            call SHMEM_GETI(ib,ib,nn, idproc1(mn))
            do i= 1,nn
              ic(i) = min(ic(i),ib(i))
            enddo
          enddo
          BARRIER
        else
c         get global minimum from 1st processor
          BARRIER
          call SHMEM_GETI(ic,ic,nn, idproc1(1))
        endif
        ! no barrier needed here because using two buffers (b and c)
#endif
        do i= 1,nn
          ia(is0+i) = ic(i)
        enddo
      enddo
#if defined(TIMER)
c
      if     (nxc.eq.10) then
        call xctmr1(10)
        nxc = 0
      endif
#endif
      return
      end subroutine xcmin_i1

      subroutine xcmin_r0(ra)
      implicit none
c
      real, intent(inout) :: ra
c
c**********
c*
c  1) replace real scalar a with its element-wise minimum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ra              real           in/out    target variable
c*
c**********
c
      real ra1(1)
c
      ra1(1) = ra
      call xcmin_r1(ra1)
      ra = ra1(1)
      return
      end subroutine xcmin_r0

      subroutine xcmin_r1(ra)
      implicit none
c
      real, intent(inout) :: ra(:)
c
c**********
c*
c  1) replace real array a with its element-wise minimum over all tiles.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    ra              real           in/out    target array
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    nmax
      parameter (nmax=1024)
c
      real           rb,rc
      common/xcmaxr/ rb(nmax),rc(nmax)
      save  /xcmaxr/
c
      integer i,is0,isl,mn,n,nn
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0(10)
        nxc = 10
      endif
#endif
c
c     stripmine ra.
c
      n = size(ra)
c
      do is0= 0,n-1,nmax
        isl = min(is0+nmax,n)
        nn = isl - is0
        do i= 1,nn
          rb(i) = ra(is0+i)
        enddo
c
#if defined(MPI)
        call mpi_allreduce(rb,rc,nn,MTYPER,mpi_min,
     &                     mpicomm,mpierr)
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.eq.1) then
c         form global minimum on 1st processor
          do i= 1,nn
            rc(i) = rb(i)
          enddo
          do mn= 2,ijpr
            call SHMEM_GETR(rb,rb,nn, idproc1(mn))
            do i= 1,nn
              rc(i) = min(rc(i),rb(i))
            enddo
          enddo
          BARRIER
        else
c         get global minimum from 1st processor
          BARRIER
          call SHMEM_GETR(rc,rc,nn, idproc1(1))
        endif
        ! no barrier needed here because using two buffers (b and c)
#endif
        do i= 1,nn
          ra(is0+i) = rc(i)
        enddo
      enddo
#if defined(TIMER)
c
      if     (nxc.eq.10) then
        call xctmr1(10)
        nxc = 0
      endif
#endif
      return
      end subroutine xcmin_r1

      subroutine xcspmd
      implicit none
c
c**********
c*
c  1) initialize data structures that identify the tiles.
c
c  2) data structures:
c      ipr     - 1st 2-D node dimension (<=iqr)
c      jpr     - 2nd 2-D node dimension (<=jqr)
c      ijpr    -     1-D node dimension (ipr*jpr)
c      mproc   - 1st 2-D node index
c      nproc   - 2nd 2-D node index
c      mnproc  -     1-D node index
c      i0      -     i0_pe(mproc,nproc)
c      ii      -     ii_pe(mproc,nproc)
c      j0      -     j0_pe(mproc,nproc)
c      jj      -     jj_pe(mproc,nproc)
c      margin  -     how much of the halo is currently valid
c      nreg    -     region type
c      vland   -     fill value for land (standard value 0.0)
c      idproc  -     2-D node addresses, with periodic wrap
c      idproc1 -     1-D node addresses, with periodic wrap
c      idhalo  -     left and right halo target nodes
c      i0_pe   - 1st dimension tile offsets
c      ii_pe   - 1st dimension tile extents (<=idm)
c      j0_pe   - 2nd dimension tile offsets
c      jj_pe   - 2nd dimension tile extents (<=jdm)
c      mpe_1   - 1st node in each row of 2-D nodes
c      mpe_e   - end node in each row of 2-D nodes
c      mpe_i   - mapping from 1st global dimension to 2-D nodes
c      npe_j   - mapping from 2nd global dimension to 2-D nodes
c      i1sum   - local index of 1st partial sum on each tile
c      iisum   - number of partial sums on each tile
c      m0_top  - tile offset:       top neighbors (0:jpr-1)
c      mm_top  - tile extent:       top neighbors (<=jpr)
c      i0_st   - halo offsets: send top neighbors
c      ii_st   - halo lengths: send top neighbors
c      i0_gt   - halo offsets:  get top neighbors
c      ii_gt   - halo lengths:  get top neighbors
c      m0_bot  - tile offset:       bot neighbors (0:jpr-1)
c      mm_bot  - tile extent:       bot neighbors (<=jpr)
c      i0_sb   - halo offsets: send bot neighbors
c      ii_sb   - halo lengths: send bot neighbors
c      i0_gb   - halo offsets:  get bot neighbors
c      ii_gb   - halo lengths:  get bot neighbors
c
c  3) all data structures are based on the processor number and
c     the patch distribution file, 'patch.input'.
c*
c**********
c
      integer i,j,l,m,mm,mn,mypei,n,npesi
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
      logical linit
c
c     standard mpi (message passing) version.
c
      call mpi_initialized(linit,mpierr)
      if (linit) then
        mpicomm = mpicom_external
      else
        call mpi_init(mpierr)
        mpicomm = mpi_comm_world
      endif
c
      call mpi_comm_rank(mpicomm, mypei, mpierr)
      call mpi_comm_size(mpicomm, npesi, mpierr)
c
      mnproc = mypei + 1  ! mnproc counts from 1
c
      call setlogunit
c
#  if defined(DEBUG_ALL)
      write(lp,'(a,i5)') 'mnproc =',mnproc
      call xcsync(flush_lp)
#  endif
#elif defined(SHMEM)
      integer SHMEM_MYPE,SHMEM_NPES
c
c     shmem version.
c
      call start_pes(0)
c
      mypei = SHMEM_MYPE()
      npesi = SHMEM_NPES()
c
      mnproc = mypei + 1  ! mnproc counts from 1
#else
c
      write(lp,*)
      write(lp,*) '***** ERROR - UNDEFINED SPMD MACHINE TYPE *****'
      write(lp,*)
      call flush(lp)
      stop '(xcspmd)'
#endif
c
      vland = 0.0
c
#if defined(T3E)
      open(unit=lp,delim='none')  ! forces open on stdout
#endif
c
c     read in the tile locations and sizes.
c     patch distibution file on unit 21 (fort.21).
c
c     here is an example of a patch file, with a leading "!" added:
c
!  npes   npe   mpe   idm   jdm  ibig  jbig  nreg
!    16     4     4    57    52    19    15     0
!
!ispt(  1) =     22    35    42    49
!iipe(  1) =     13     7     7     7
!ispt(  2) =      1    15    25    34
!iipe(  2) =     14    10     9     9
!ispt(  3) =      2    21    29    38
!iipe(  3) =     19     8     9    10
!ispt(  4) =     18    28    35    42
!iipe(  4) =     10     7     7     9
! 
!jspt(  1) =      1    15    26    38
!jjpe(  1) =     14    11    12    15
c
c     ispt(j) is the 1st i-point on each tile in the j-th row
c     iipe(j) is the i-extent    of each tile in the j-th row
c     jspt(1) is the 1st j-point on each tile in all columns
c     jjpe(1) is the j-extent    of each tile in all columns
c
c     note that each tile row can have a different tile layout, 
c     but each tile column is identical.  So a tile can have at
c     most one nieghbor in the E and W directions, but several
c     nieghbors in the N and S directions.
c
c     iipe can be zero, indicating an empty tile (not included in the
c     active tile count, npes) and therefore at least a halo widths
c     land separation between active tiles to its east and west.  In
c     periodic cases, the last non-empty tile can periodic wrap to the
c     first tile in the row (i.e. trailing "empty" tiles can be null
c     tiles, rather than all-land tiles).
c
      ijpr=ijqr
      ipr=iqr
      jpr=jqr
#if defined(ARCTIC)
      if     (nreg.ne.2) then  ! not arctic
        if     (mnproc.eq.1) then
          write(lp,'(a,i5)') 'input: nreg =',nreg
          call flush(lp)
        endif
        call xcstop('xcspmd: patch.input must be for arctic')
        stop '(xcspmd)'
      endif
#else
      if     (nreg.eq.2) then  ! arctic
        if     (mnproc.eq.1) then
          write(lp,'(a,i5)') 'input: nreg =',nreg
          call flush(lp)
        endif
        call xcstop('xcspmd: patch.input for arctic but ARCTIC undef.')
        stop '(xcspmd)'
      elseif (nreg.lt.0) then  ! nreg = -1 reserved for TYPE=one/omp 
        if     (mnproc.eq.1) then
          write(lp,'(a,i5)') 'input: nreg =',nreg
          call flush(lp)
        endif
        call xcstop('xcspmd: patch.input for wrong nreg')
        stop '(xcspmd)'
      endif
#endif /* ARCTIC:else */
c
c     individual tile rows.
c
      do n= 1,jpr
        if (maxval(ii_pe(1:ipr,n)).le.0) then
          call xcstop('xcspmd: decomposition has an empty row')
          stop '(xcspmd)'
        endif
      enddo
#if defined(ARCTIC)
c
c --- all arctic patch tiles must be the same size or empty,
c --- and empty tiles must be "twinned" across the top boundary.
c
      if     (ipr.gt.1) then
        do m= 1,ipr
          if     (ii_pe(m,jpr).eq.0) then
            if     (ii_pe(ipr+1-m,jpr).ne.0) then
              if     (mnproc.eq.1) then
                write(lp,'(a,i3,a,i3,a)') 
     &           'error - tile',m,',jpr is empty but tile',
     &                    ipr+1-m,',jpr is not'
                call flush(lp)
              endif
              call xcstop('xcspmd: arctic empty tiles are not twins')
              stop '(xcspmd)'
            endif
          elseif (ii_pe(m,jpr).ne.itdm/ipr) then
            if     (mnproc.eq.1) then
              write(lp,'(a,i5)') 
     &         'error - arctic patch tiles should have ii =',itdm/ipr
              call flush(lp)
            endif
            call xcstop('xcspmd: arctic tiles are not the right size')
            stop '(xcspmd)'
          endif
        enddo !m
      endif
#endif /* ARCTIC */
c
c     the generic tile column (must cover entire column).
c
      if     (j0_pe(1,1).ne.0) then
        call xcstop('xcspmd: decomposition has wrong j0_pe')
        stop '(xcspmd)'
      endif 
      do n= 2,jpr
        if (j0_pe(1,n).ne.j0_pe(1,n-1)+jj_pe(1,n-1)) then
          call xcstop('xcspmd: decomposition is non-contiguous')
          stop '(xcspmd)'
        endif
      enddo
      if     (j0_pe(1,jpr)+jj_pe(1,jpr).ne.jtdm) then
        call xcstop('xcspmd: decomposition has wrong jj_pe')
        stop '(xcspmd)'
      endif 
c
#if defined(MPI)
c
c     do we have the right number of pes?
c
      if     (npesi.ne.ijpr) then
        if     (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) '***** ERROR - WRONG MPI SIZE *****'
          write(lp,*)
          write(lp,*) 'NPES    = ',npesi
          write(lp,*) '   IJPR = ',    ijpr
          write(lp,*) 'IPR,JPR = ',ipr,jpr
          write(lp,*)
          call flush(lp)
        endif
        call xcstop('Error in xcspmd')
        stop
      endif
c
c     mpi messages are sent and received by pe number (0:ijpr-1).
c
      null_tile = mpi_proc_null
c
      mn = 0
      do n= 1,jpr
        mpe_1(n) = 0
        do m= 1,ipr
          if     (ii_pe(m,n).eq.0) then
            idproc(m,n)   = null_tile
          else
            idproc1(mn+1) = mn
            idproc(m,n)   = mn
            mn = mn + 1
            if     (mnproc.eq.mn) then
              mproc = m
              nproc = n
            endif
            mpe_e(n) = m
            if     (mpe_1(n).eq.0) then
              mpe_1(n) = m
            endif
          endif
        enddo
      enddo
#elif defined(SHMEM)
c
c     do we have the right number of pes?
c
      if     (npesi.ne.ijpr) then
        if     (mnproc.eq.1) then
          write(lp,*)
          write(lp,*) '***** ERROR - WRONG SHMEM SIZE *****'
          write(lp,*)
          write(lp,*) 'NPES    = ',npesi
          write(lp,*) '   IJPR = ',    ijpr
          write(lp,*) 'IPR,JPR = ',ipr,jpr
          write(lp,*)
          call flush(lp)
        endif
        call xcstop('Error in xcspmd')
        stop
      endif
c
c     shmem messages are sent and received by pe number (0:ijpr-1).
c
      null_tile = -1
c
      mn = 0
      do n= 1,jpr
        mpe_1(n) = 0
        do m= 1,ipr
          if     (ii_pe(m,n).eq.0) then
            idproc(m,n)   = null_tile
          else
            idproc1(mn+1) = mn
            idproc(m,n)   = mn
            mn = mn + 1
            if     (mnproc.eq.mn) then
              mproc = m
              nproc = n
            endif
            mpe_e(n) = m
            if     (mpe_1(n).eq.0) then
              mpe_1(n) = m
            endif
          endif
        enddo
      enddo
#endif
c
      if     (mn.ne.ijpr) then
        if     (mnproc.eq.1) then
          write(lp,'(a,i5)') 'input: ijpr =',ijpr
          write(lp,'(a,i5)') 'calc:  ijpr =',mn
          call flush(lp)
        endif
        call xcstop('xcspmd: wrong number of sea tiles')
        stop '(xcspmd)'
      endif
c
      if     (nreg.eq.0.or.nreg.eq.4) then
c
c       longitudinal tile dimension is closed (not periodic)
c
        do n= 1,jpr
          idproc(    0,n) = null_tile
          idproc(ipr+1,n) = null_tile
        enddo
      else
c
c       longitudinal tile dimension is potentially periodic.
c
        do n= 1,jpr
          idproc(    0,n) = null_tile
          idproc(ipr+1,n) = null_tile
c
          i = maxval((i0_pe(1:ipr,n)+ii_pe(1:ipr,n)))
          if     (i0_pe(1,n).eq.0 .and. i.eq.itdm) then
            idproc(         0,n) = idproc(mpe_e(n),n)
            idproc(mpe_e(n)+1,n) = idproc(       1,n)
          endif
        enddo
      endif
#if defined(ARCTIC)
c
c     must have ipr even or 1 for arctic boundary case.
c
      if     (ipr.gt.1 .and. mod(ipr,2).ne.0) then
        call xcstop('Error in xcspmd (arctic) - ipr must be even')
        stop '(xcspmd)'
      endif
c
c     latitudinal tile dimension is closed/arctic.
c
      do m= 1,ipr
        idproc(m,    0) = null_tile
        idproc(m,jpr+1) = idproc(ipr+1-m,jpr)  !arctic tile mapping
      enddo
      idproc(    0,    0) = null_tile
      idproc(ipr+1,    0) = null_tile
      idproc(    0,jpr+1) = idproc(ipr,jpr+1)
      idproc(ipr+1,jpr+1) = idproc(1,  jpr+1)
#else
c
c
      if (nreg.ge.3) then
c       latitudinal tile dimension is periodic
        do m= 0,ipr+1
          idproc(m,    0) = idproc(m,jpr)
          idproc(m,jpr+1) = idproc(m,  1)
        enddo
      else
c       latitudinal tile dimension is closed
        do m= 0,ipr+1
          idproc(m,    0) = null_tile
          idproc(m,jpr+1) = null_tile
        enddo
      endif
#endif /* ARCTIC:else */
c
c     1-d tiling logic is easier if assumed periodic.
c
      idproc1(     0) = idproc1(ijpr)
      idproc1(ijpr+1) = idproc1(   1)
c
c     mapping from global i,j to mp,np.
c     ia,ja is on tile mpe_i(ia,npe_j(ja)),npe_j(ja), 
c     or on no tile if mpe_i(ia,npe_j(ja)) is 0 or -1.
c
      do n= 1,jpr
        mpe_i(1:itdm,n) = 0  ! default is an empty tile
        do m= 1,ipr  ! i-map potentially varies with n
          if     (ii_pe(m,n).gt.0) then
            do i= i0_pe(m,n)+1,i0_pe(m,n)+ii_pe(m,n)
              mpe_i(i,n) = m
            enddo
            if     (m.ne.ipr) then
              if     (ii_pe(m+1,n).gt.0) then
                do i= i0_pe(m,n)+ii_pe(m,n)+1,i0_pe(m+1,n)
                  mpe_i(i,n) = -1  ! gap between tiles
                enddo
              endif
            endif
          endif
        enddo
        m = 1  ! only one j-map
          do j= j0_pe(m,n)+1,j0_pe(m,n)+jj_pe(m,n)
            npe_j(j)   = n
          enddo
      enddo
c
c     do each partial sum on the tile that owns its center point.
c       i1sum - local index of 1st partial sum on each tile
c       iisum - number of partial sums on each tile
c     see xcsum for how i1sum and iisum are used.
c
      do n= 1,jpr
        do m= 1,ipr
          if     (ii_pe(m,n).le.0) then
            i1sum(m,n) =  0
            iisum(m,n) =  0
          else
            idhalo(1) = idproc(m-1,n)
            idhalo(2) = idproc(m+1,n)
            if     (idhalo(1).ne.null_tile .and. m.ne.1) then
              if (i0_pe(m,n).ne.i0_pe(m-1,n)+ii_pe(m-1,n)) then
                idhalo(1) = null_tile
              endif
            endif
            if     (idhalo(2).ne.null_tile .and. m.ne.mpe_e(n)) then
              if     (i0_pe(m,n)+ii_pe(m,n).ne.i0_pe(m+1,n)) then
                idhalo(2) = null_tile
              endif
            endif
            i1sum(m,n) = -99
            iisum(m,n) =  0
            do i= 1+nbdy,itdm+nbdy,2*nbdy+1
              if     (i0_pe(m,n).lt.i .and.
     &                              i.le.i0_pe(m,n)+ii_pe(m,n)) then
                iisum(m,n) = iisum(m,n) + 1
                if     (iisum(m,n).eq.1) then
                  i1sum(m,n) = i - nbdy - i0_pe(m,n)
                endif
              elseif (idhalo(1).eq.null_tile .and.
     &                i.gt.i0_pe(m,n)-nbdy   .and.
     &                i.le.i0_pe(m,n)             ) then
                iisum(m,n) = iisum(m,n) + 1
                if     (iisum(m,n).eq.1) then
                  i1sum(m,n) = i - nbdy - i0_pe(m,n)
                endif
              elseif ((idhalo(2).eq.null_tile.or.m.eq.mpe_e(n)) .and.
     &                i.gt.i0_pe(m,n)+ii_pe(m,n)                .and.
     &                i.le.i0_pe(m,n)+ii_pe(m,n)+nbdy      ) then
                iisum(m,n) = iisum(m,n) + 1
                if     (iisum(m,n).eq.1) then
                  i1sum(m,n) = i - nbdy - i0_pe(m,n)
                endif
              endif
            enddo
          endif
        enddo
      enddo
c
c     local tile extents.
c
      i0 = i0_pe(mproc,nproc)
      ii = ii_pe(mproc,nproc)
      j0 = j0_pe(mproc,nproc)
      jj = jj_pe(mproc,nproc)
c
      margin = 0
c
c     left and right halo targets
c
      idhalo(1) = idproc(mproc-1,nproc)
      idhalo(2) = idproc(mproc+1,nproc)
c
      if     (idhalo(1).ne.null_tile .and. mproc.ne.1) then
c
c       is the left tile touching this one?
c
        if (i0.ne.i0_pe(mproc-1,nproc)+ii_pe(mproc-1,nproc)) then
          idhalo(1) = null_tile
        endif
      endif
c
      if     (idhalo(2).ne.null_tile .and. mproc.ne.mpe_e(nproc)) then
c
c       is the right tile touching this one?
c
        if     (i0+ii.ne.i0_pe(mproc+1,nproc)) then
          idhalo(2) = null_tile
        endif
      endif
c
c     local halo exchange data structures
c
c     m0_top - tile offset:       top neighbors
c     mm_top - tile extent:       top neighbors (<=jpr)
c     i0_st  - halo offsets: send top neighbors
c     ii_st  - halo lengths: send top neighbors
c     i0_gt  - halo offsets:  get top neighbors
c     ii_gt  - halo lengths:  get top neighbors
c     m0_bot - tile offset:       bot neighbors
c     mm_bot - tile extent:       bot neighbors (<=jpr)
c     i0_sb  - halo offsets: send bot neighbors
c     ii_sb  - halo lengths: send bot neighbors
c     i0_gb  - halo offsets:  get bot neighbors
c     ii_gb  - halo lengths:  get bot neighbors
c
c     note that send is also receive, and is w.r.t. the local  tile.
c     similarly get  is also put,     and is w.r.t. the remote tile.
c
      if     (nproc.eq.jpr) then
#if defined(ARCTIC)
c       single, same size, top arctic nieghbor
        m0_top   = mproc - 1
        mm_top   =  1
        i0_st(1) =  0
        i0_gt(1) =  0
        ii_st(1) = ii
        ii_gt(1) = ii
#else
        m0_top = 0
        mm_top = 0
        if (nreg.ge.3) then
c         latitudinal tile dimension is open
          n = 1
          m = 0
          do i= 1,ii
            if     (mpe_i(i0+i,n).ne.m) then
              if     (mm_top.eq.0) then
                m0_top = mpe_i(i0+i,n) - 1
              elseif (m.ne.-1) then
                ii_st(mm_top) = i-1 - i0_st(mm_top) 
                ii_gt(mm_top) = ii_st(mm_top)
              endif
              m = mpe_i(i0+i,n)
              if     (m.gt.0) then
                mm_top        = mm_top + 1
                i0_st(mm_top) = i-1
                i0_gt(mm_top) = i-1 + i0-i0_pe(m,n)
              elseif (m.eq.0) then
                mm_top        = mm_top + 1
                i0_st(mm_top) = i-1
                i0_gt(mm_top) = i0_gt(mm_top-1) + ii_gt(mm_top-1)
*             elseif (m.eq.-1) then  !do nothing
              endif
            endif
          enddo
          if     (mm_top.gt.0) then
            if     (m.gt.0) then
              ii_st(mm_top) = ii - i0_st(mm_top) 
              ii_gt(mm_top) = ii_st(mm_top)
            elseif (m.eq.0) then
              mm_top = mm_top-1
*           elseif (m.eq.-1) then  !do nothing
            endif
          endif
        endif
#endif /* ARCTIC:else */
      else
        n = nproc + 1
        m0_top = 0
        mm_top = 0
        m      = 0
        do i= 1,ii
          if     (mpe_i(i0+i,n).ne.m) then
            if     (mm_top.eq.0) then
              m0_top = mpe_i(i0+i,n) - 1
            elseif (m.ne.-1) then
              ii_st(mm_top) = i-1 - i0_st(mm_top) 
              ii_gt(mm_top) = ii_st(mm_top)
            endif
            m = mpe_i(i0+i,n)
            if     (m.gt.0) then
              mm_top        = mm_top + 1
              i0_st(mm_top) = i-1
              i0_gt(mm_top) = i-1 + i0-i0_pe(m,n)
            elseif (m.eq.0) then
              mm_top        = mm_top + 1
              i0_st(mm_top) = i-1
              i0_gt(mm_top) = i0_gt(mm_top-1) + ii_gt(mm_top-1)
*           elseif (m.eq.-1) then  !do nothing
            endif
          endif
        enddo
        if     (mm_top.gt.0) then
          if     (m.gt.0) then
            ii_st(mm_top) = ii - i0_st(mm_top) 
            ii_gt(mm_top) = ii_st(mm_top)
          elseif (m.eq.0) then
            mm_top = mm_top-1
*         elseif (m.eq.-1) then  !do nothing
          endif
        endif
      endif  !nproc.eq.1:else
c
      if     (nproc.eq.1) then
c       no bottom nieghbor (closed boundary)
        m0_bot = 0
        mm_bot = 0
        if (nreg.ge.3) then
c         latitudinal tile dimension is open
          n = jpr
          m = 0
          do i= 1,ii
            if     (mpe_i(i0+i,n).ne.m) then
              if     (mm_bot.eq.0) then
                m0_bot = mpe_i(i0+i,n) - 1
              elseif (m.ne.-1) then
                ii_sb(mm_bot) = i-1 - i0_sb(mm_bot) 
                ii_gb(mm_bot) = ii_sb(mm_bot)
              endif
              m = mpe_i(i0+i,n)
              if     (m.gt.0) then
                mm_bot        = mm_bot + 1
                i0_sb(mm_bot) = i-1
                i0_gb(mm_bot) = i-1 + i0-i0_pe(m,n)
              elseif (m.eq.0) then
                mm_bot        = mm_bot + 1
                i0_sb(mm_bot) = i-1
                i0_gb(mm_bot) = i0_gb(mm_bot-1) + ii_gb(mm_bot-1)
*             elseif (m.eq.-1) then  !do nothing
              endif
            endif
          enddo
          if     (mm_bot.gt.0) then
            if     (m.gt.0) then
              ii_sb(mm_bot) = ii - i0_sb(mm_bot) 
              ii_gb(mm_bot) = ii_sb(mm_bot)
            elseif (m.eq.0) then
              mm_bot = mm_bot-1
*           elseif (m.eq.-1) then  !do nothing
            endif
          endif
        endif
      else
        n = nproc - 1
        m0_bot = 0
        mm_bot = 0
        m      = 0
        do i= 1,ii
          if     (mpe_i(i0+i,n).ne.m) then
            if     (mm_bot.eq.0) then
              m0_bot = mpe_i(i0+i,n) - 1
            elseif (m.ne.-1) then
              ii_sb(mm_bot) = i-1 - i0_sb(mm_bot) 
              ii_gb(mm_bot) = ii_sb(mm_bot)
            endif
            m = mpe_i(i0+i,n)
            if     (m.gt.0) then
              mm_bot        = mm_bot + 1
              i0_sb(mm_bot) = i-1
              i0_gb(mm_bot) = i-1 + i0-i0_pe(m,n)
            elseif (m.eq.0) then
              mm_bot        = mm_bot + 1
              i0_sb(mm_bot) = i-1
              i0_gb(mm_bot) = i0_gb(mm_bot-1) + ii_gb(mm_bot-1)
*           elseif (m.eq.-1) then  !do nothing
            endif
          endif
        enddo
        if     (mm_bot.gt.0) then
          if     (m.gt.0) then
            ii_sb(mm_bot) = ii - i0_sb(mm_bot) 
            ii_gb(mm_bot) = ii_sb(mm_bot)
          elseif (m.eq.0) then
            mm_bot = mm_bot-1
*         elseif (m.eq.-1) then  !do nothing
          endif
        endif
      endif  !nproc.eq.1:else
c
c     printout the tile data structures.
c
      if     (mnproc.eq.1) then
        write(lp,'(/a)')
     &   'mnproc mproc nproc     i0   ii     j0   jj  i1sum iisum'
        mn = 0
        do n= 1,jpr
          do m= 1,ipr
            if     (ii_pe(m,n).ne.0) then
              mn= mn + 1
              write(lp,'(i6,2i6,i7,i5,i7,i5,i7,i6)')
     &           mn,m,n,
     &           i0_pe(m,n),ii_pe(m,n),
     &           j0_pe(m,n),jj_pe(m,n),
     &           i1sum(m,n),iisum(m,n)
            endif
          enddo !m
        enddo !n
        write(lp,*)
#if defined(ARCTIC)
        write(lp,'(a)')
     &   'mnproc mproc nproc mnarct'
        mn = 0
        do n= 1,jpr
          do m= 1,ipr
            if     (ii_pe(m,n).ne.0) then
              mn= mn + 1
              if     (n.eq.jpr) then
                write(lp,'(i6,2i6,i7)')
     &           mn,m,n,idproc(m,n+1)
              endif
            endif
          enddo !m
        enddo !n
        write(lp,*)
#endif /* ARCTIC */
#if defined(DEBUG_ALL)
        do n= 1,jpr
          write(lp,*) 'mpe_1,mpe_e = ',mpe_1(n),mpe_e(n)
        enddo
        do n= 1,jpr
          write(lp,*) 'mpe_i = ',mpe_i(:,n)
        enddo
        write(lp,*)
        write(lp,*) 'npe_j = ',npe_j(:)
        write(lp,*)
#endif
      endif
      call xcsync(flush_lp)
c
#if defined(DEBUG_ALL)
      do n= 1,jpr
        do m= 1,ipr
          if     (mproc.eq.m .and. nproc.eq.n) then
            write(lp,'(a,2i3,i4,i3,16(i5,i4))')
     &         'm,n,_top,_st = ',
     &          m,n,m0_top,mm_top,(i0_st(l),ii_st(l), l= 1,mm_top)
#if defined(SHMEM)
            write(lp,'(a,2i3,i4,i3,16(i5,i4))')
     &         'm,n,_top,_gt = ',
     &          m,n,m0_top,mm_top,(i0_gt(l),ii_gt(l), l= 1,mm_top)
#endif
            if     (m.eq.ipr) then
              write(lp,*) !blank line
            endif
          endif
          call xcsync(flush_lp)
        enddo !m
        do m= 1,ipr
          if     (mproc.eq.m .and. nproc.eq.n) then
            write(lp,'(a,2i3,i4,i3,16(i5,i4))')
     &         'm,n,_bot,_sb = ',
     &          m,n,m0_bot,mm_bot,(i0_sb(l),ii_sb(l), l= 1,mm_bot)
#if defined(SHMEM)
            write(lp,'(a,2i3,i4,i3,16(i5,i4))')
     &         'm,n,_bot,_gb = ',
     &          m,n,m0_bot,mm_bot,(i0_gb(l),ii_gb(l), l= 1,mm_bot)
#endif
            if     (m.eq.ipr) then
              write(lp,*) !blank line
            endif
          endif
          call xcsync(flush_lp)
        enddo !m
        do m= 1,ipr
          if     (mproc.eq.m .and. nproc.eq.n) then
            write(lp,'(a,2i3,3i5)')
     &         'm,n,id,idhalo   = ',
     &          m,n,idproc(m,n),idhalo
            if     (m.eq.ipr) then
              write(lp,*) !blank line
            endif
          endif
          call xcsync(flush_lp)
        enddo !m
      enddo !n
#endif /* DEBUG_ALL */
c
#if defined(MPI) && defined(MPISRJ)
      call xcmax(mm_top)
      call xcmax(mm_bot)
      if     (max(mm_top,mm_bot).gt.1) then
        call xcstop(
     &    'xcspmd: if MPISRJ is defined, tiles must be uniform')
        stop '(xcspmd)'
      endif
#endif
c
c     initialize timers.
c
      call xctmri
#if defined(TIMER)
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
#endif
      return
      end subroutine xcspmd

      subroutine xcstop(cerror)
      implicit none
c
      character*(*), intent(in) :: cerror
c
c**********
c*
c  1) stop all processes.
c
c  2) all processes must call this routine.
c     use 'xchalt' for emergency stops.
c
c  3) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    cerror          char*(*)       input     error message
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
c     print active timers.
c
      call xctmrp
c
c     message passing version, set barrier and then stop everything.
c     use xcsync as the barrier so that stdout is flushed.
c     note that the system will hang unless all processes call xcstop.
c
      call xcsync(flush_lp)
      if     (mnproc.eq.1 .and. cerror.ne.' ') then
        write(lp,*) '**************************************************'
        write(lp,*) cerror
        write(lp,*) '**************************************************'
        write(lp,*)
      endif
      call xcsync(flush_lp)
c
#if defined(MPI)
      if (mpicomm.eq.mpicom_external) then
        call external_abort(cerror)
      else
        write(lp,*) 'mpi_finalize called on processor ',mnproc
        call xcsync(flush_lp)
        call mpi_finalize(mpierr)
      endif
#endif
c
      stop '(xcstop)'
      end subroutine xcstop

      subroutine xcsum(sum, a,mask)
      implicit none
c
      real*8,  intent(out)   :: sum
      real,    intent(inout) :: a(   1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) sum a 2-d array, where mask==1
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    sum             real*8         output    sum of a
c    a               real           input     source array
c    mask            integer        input     mask array
c
c  3) sum is bit for bit reproducable for the same halo size, nbdy.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    mxsum
      parameter (mxsum=(idm+3*nbdy)/(2*nbdy+1))
c
      real*8         sum8t,sum8j,sum8s
      common/xcsum8/ sum8t(mxsum*jdm),sum8j(jdm),sum8s
      save  /xcsum8/
c
      real*8     zero8
      parameter (zero8=0.0)
c
      real*8  sum8
      real    vsave
      integer i,i1,j,l,mp,np
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 6)
        nxc = 6
      endif
#endif
c
c     halo update so that 2*nbdy+1 wide strips are on chip.
c
      vsave = vland
      vland = 0.0
      call xctilr(a,1,1, nbdy,0, halo_ps)
      vland = vsave
c
c     row sums in 2*nbdy+1 wide strips.
c
!$OMP PARALLEL DO PRIVATE(j,i1,i,l,sum8)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        do l= 1,iisum(mproc,nproc)
          i1   = i1sum(mproc,nproc) + (l-1)*(2*nbdy+1)
          sum8 = zero8
          do i= max(i1,1-nbdy),min(i1+2*nbdy,ii+nbdy,itdm-i0)
            if     (mask(i,j).eq.1) then
              sum8 = sum8 + a(i,j)
            endif
          enddo
          sum8t(l + (j-1)*iisum(mproc,nproc)) = sum8
        enddo
      enddo
!$OMP END PARALLEL DO
c
c     complete row sums on first processor in each row.
c
#if defined(SHMEM)
      BARRIER
#endif
      if     (mproc.eq.mpe_1(nproc)) then
        do j=1,jj
          sum8j(j) = zero8
          do l= 1,iisum(mproc,nproc)
            sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mproc,nproc))
          enddo
*         write(lp,'(a,i3,i5,f12.2)') 'xcsum: np,j,sum = ',
*    &                                1,j,sum8j(j)
        enddo
c
c       remote sums.
c
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = iisum(mp,nproc)*jj
          if     (l.gt.0) then
#if defined(MPI)
            call MPI_RECV(sum8t,l,MTYPED,
     &                    idproc(mp,nproc), 9900,
     &                    mpicomm, mpistat, mpierr)
#elif defined(SHMEM)
            call SHMEM_GETD(sum8t,
     &                      sum8t,l,idproc(mp,nproc))
#endif
            do j=1,jj
              do l= 1,iisum(mp,nproc)
                sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mp,nproc))
              enddo
*             write(lp,'(a,i3,i5,f12.2)') 'xcsum: np,j,sum = ',
*    &                                    mp,j,sum8j(j)
            enddo
          endif
        enddo
#if defined(MPI)
      else
        l = iisum(mproc,nproc)*jj
        if     (l.gt.0) then
          call MPI_SEND(sum8t,l,MTYPED,
     &                  idproc(mpe_1(nproc),nproc), 9900,
     &                  mpicomm, mpierr)
        endif
#endif
      endif
c
c     sum of row sums, on first processor.
c
#if defined(SHMEM)
      BARRIER
#endif
      if     (mnproc.eq.1) then
        sum8 = zero8
        do j= 1,jj
          sum8 = sum8 + sum8j(j)
        enddo
*       write(lp,'(a,i5,f12.2)') 'xcsum: jj,sum = ',jj,sum8
c
        do np= 2,jpr
          mp = mpe_1(np)
#if defined(MPI)
          call MPI_RECV(sum8j,jj_pe(mp,np),MTYPED,
     &                  idproc(mp,np), 9901,
     &                  mpicomm, mpistat, mpierr)
#elif defined(SHMEM)
          call SHMEM_GETD(sum8j,
     &                    sum8j,jj_pe(mp,np),idproc(mp,np))
#endif
          do j= 1,jj_pe(mp,np)
            sum8 = sum8 + sum8j(j)
          enddo
*         write(lp,'(a,i5,f12.2)') 'xcsum: jj,sum = ',
*    &                             jj_pe(mp,np),sum8
        enddo
        sum8s = sum8
#if defined(MPI)
      elseif (mproc.eq.mpe_1(nproc)) then
        call MPI_SEND(sum8j,jj,MTYPED,
     &                idproc1(1), 9901,
     &                mpicomm, mpierr)
#endif
      endif
c
c     broadcast result to all processors.
c
#if defined(MPI)
      call mpi_bcast(sum8s,1,MTYPED,
     &               idproc1(1),mpicomm,mpierr)
#elif defined(SHMEM)
      BARRIER
      if     (mnproc.ne.1) then
        call SHMEM_GETD(sum8s,
     &                  sum8s,1,idproc1(1))
      endif
#endif
c
      sum = sum8s
#if defined(TIMER)
c
      if     (nxc.eq. 6) then
        call xctmr1( 6)
        nxc = 0
      endif
#endif
      return
      end subroutine xcsum

      subroutine xcsumj(sumj, a,mask)
      implicit none
c
      real*8,  intent(out)   :: sumj(jtdm)
      real,    intent(inout) :: a(   1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
      integer, intent(in)    :: mask(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)
c
c**********
c*
c  1) rwo-sum of a 2-d array, where mask==1, on first processor only.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    sumj            real*8         output    row-sum of a
c    a               real           input     source array
c    mask            integer        input     mask array
c
c  3) sum is bit for bit reproducable for the same halo size, nbdy.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer    mxsum
      parameter (mxsum=(idm+3*nbdy)/(2*nbdy+1))
c
      real*8         sum8t,sum8j,sum8s
      common/xcsum8/ sum8t(mxsum*jdm),sum8j(jdm),sum8s
      save  /xcsum8/
c
      real*8     zero8
      parameter (zero8=0.0)
c
      real*8  sum8
      real    vsave
      integer i,i1,j,l,mp,np
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0( 6)
        nxc = 6
      endif
#endif
c
c     halo update so that 2*nbdy+1 wide strips are on chip.
c
      vsave = vland
      vland = 0.0
      call xctilr(a,1,1, nbdy,0, halo_ps)
      vland = vsave
c
c     row sums in 2*nbdy+1 wide strips.
c
!$OMP PARALLEL DO PRIVATE(j,i1,i,l,sum8)
!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1,jj
        do l= 1,iisum(mproc,nproc)
          i1   = i1sum(mproc,nproc) + (l-1)*(2*nbdy+1)
          sum8 = zero8
          do i= i1,min(i1+2*nbdy,ii+nbdy,itdm-i0)
            if     (mask(i,j).eq.1) then
              sum8 = sum8 + a(i,j)
            endif
          enddo
          sum8t(l + (j-1)*iisum(mproc,nproc)) = sum8
        enddo
      enddo
!$OMP END PARALLEL DO
c
c     complete row sums on first processor in each row.
c
#if defined(SHMEM)
      BARRIER
#endif
      if     (mproc.eq.mpe_1(nproc)) then
        do j=1,jj
          sum8j(j) = zero8
          do l= 1,iisum(mproc,nproc)
            sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mproc,nproc))
          enddo
*         write(lp,'(a,i3,i5,f12.2)') 'xcsum: np,j,sum = ',
*    &                                1,j,sum8j(j)
        enddo
c
c       remote sums.
c
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = iisum(mp,nproc)*jj
          if     (l.gt.0) then
#if defined(MPI)
            call MPI_RECV(sum8t,l,MTYPED,
     &                    idproc(mp,nproc), 9900,
     &                    mpicomm, mpistat, mpierr)
#elif defined(SHMEM)
            call SHMEM_GETD(sum8t,
     &                      sum8t,l,idproc(mp,nproc))
#endif
            do j=1,jj
              do l= 1,iisum(mp,nproc)
                sum8j(j) = sum8j(j) + sum8t(l + (j-1)*iisum(mp,nproc))
              enddo
*             write(lp,'(a,i3,i5,f12.2)') 'xcsum: np,j,sum = ',
*    &                                    mp,j,sum8j(j)
            enddo
          endif
        enddo
#if defined(MPI)
      else
        l = iisum(mproc,nproc)*jj
        if     (l.gt.0) then
          call MPI_SEND(sum8t,l,MTYPED,
     &                  idproc(mpe_1(nproc),nproc), 9900,
     &                  mpicomm, mpierr)
        endif
#endif
      endif
c
c     send row sums to first processor.
c
#if defined(SHMEM)
      BARRIER
#endif
      if     (mnproc.eq.1) then
        do j= 1,jj
          sumj(j) = sum8j(j)
        enddo
c
        do np= 2,jpr
          mp = mpe_1(np)
#if defined(MPI)
          call MPI_RECV(sum8j,jj_pe(mp,np),MTYPED,
     &                  idproc(mp,np), 9901,
     &                  mpicomm, mpistat, mpierr)
#elif defined(SHMEM)
          call SHMEM_GETD(sum8j,
     &                    sum8j,jj_pe(mp,np),idproc(mp,np))
#endif
          do j= 1,jj_pe(1,np)
            sumj(j+j0_pe(1,np)) = sum8j(j)
          enddo
        enddo
#if defined(MPI)
      elseif (mproc.eq.mpe_1(nproc)) then
        call MPI_SEND(sum8j,jj,MTYPED,
     &                idproc1(1), 9901,
     &                mpicomm, mpierr)
#endif
      endif
#if defined(TIMER)
c
      if     (nxc.eq. 6) then
        call xctmr1( 6)
        nxc = 0
      endif
#endif
      return
      end subroutine xcsumj

      subroutine xcsync(lflush)
      implicit none
c
      logical, intent(in) :: lflush
c
c**********
c*
c  1) barrier, no processor exits until all arrive (and flush stdout).
c
c  2) some MPI implementations only flush stdout as a collective
c     operation, and hence the lflush=.true. option to flush stdout.
c
c  3) typically this is just a wrapper to the "BARRIER" macro.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      if     (lflush) then
        call flush(lp)
        BARRIER
      else
        BARRIER
      endif
      return
      end subroutine xcsync

#if defined(SHMEM) && defined(RINGB)
      subroutine xctbar(ipe1,ipe2)
      implicit none
c
      integer, intent(in) :: ipe1,ipe2
c
c**********
c*
c  1) sync with processors ipe1 and ipe2.
c
c  2) this is a global collective operation, and the calls on ipe1
c     and ipe2 must list this processor as one of the two targets.
c
c  3) this is used in place of a global barrier in halo operations,
c     but it only provides syncronization of one or two processors 
c     with the local processor.
c
c  4) ipe1 and/or ipe2 can be null_tile, to indicate no processor.
c*
c**********
c
      integer    cache_line,ilarge
      parameter (cache_line=32, ilarge=2**30)
c
      integer        ibp
      common/halobp/ ibp(cache_line,-1:ijpr-1)
      save  /halobp/
c
      integer i
c
      integer icount
      save    icount
      data    icount / -1 /
c
      icount = mod(icount+1,ilarge)
      if     (icount.eq.0) then
        call shmem_barrier_all()
        do i= -1,ijpr-1
          ibp(1,i) = -1
        enddo
        call shmem_barrier_all()
      endif
c
      ibp(1,-1) = icount
      if     (ipe1.ne.null_tile) then
        call shmem_integer_put(ibp(1,mnproc-1),icount,1,ipe1)
      endif
      if     (ipe2.ne.null_tile) then
        call shmem_integer_put(ibp(1,mnproc-1),icount,1,ipe2)
      endif
      call shmem_fence
c
      i = -1
      do while (i.lt.icount)
c       this assignment statement must not be optimized away.
cdir$   suppress
        i = min(ibp(1,ipe1),ibp(1,ipe2))
      enddo
      return
      end subroutine xctbar
#endif /* SHMEM && RINGB */

#if defined(ARCTIC)
      recursive subroutine xctilr(a,l1,ld,mh,nh,itype)
      implicit none
c
      integer, intent(in)    :: l1,ld,mh,nh,itype
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)
c
c**********
c*
c  1) update the tile overlap halo of a real array.
c
c     this version of arctic bi-polar patch only
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target array
c    l1              integer        input     3rd dim. start index
c    ld              integer        input     3rd dimension of a
c    mh              integer        input     1st (EW) update halo size
c    nh              integer        input     2nd (NS) update halo size
c    itype           integer        input     grid and field type
c
c  3) itype selects both the grid and field type
c        itype= 1; p-grid, scalar field
c        itype= 2; q-grid, scalar field
c        itype= 3; u-grid, scalar field
c        itype= 4; v-grid, scalar field
c        itype=11; p-grid, vector field
c        itype=12; q-grid, vector field
c        itype=13; u-grid, vector field
c        itype=14; v-grid, vector field
c
c  4) the global variable vland is returned by halos over land.
c*
c**********
c
      integer    lsize,ilen,jlen
      parameter (lsize=kdm,
     &           ilen= idm        *lsize*(nbdy+1)+64,
     &           jlen=(jdm+2*nbdy)*lsize* nbdy   +64)
c
c     halo buffer (in common for enhanced MPI safety).
c
      real            ai,aj
      common/xctilr4/ ai(ilen,4),aj(jlen,4)
      save  /xctilr4/
c
      real            aia
      common/xctilra/ aia(lsize*(nbdy+1)+64,2)
      save  /xctilra/
c
      integer i,io,j,k,l,lts,lbs,lbr,ltr,las,lar,lg0,ls0,lm,m,mhl,nhl
      real    sarc
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
c
c     communication handles.
c
      integer mpireqa(4*iqr),mpireqb(4),
     &        ilold,ltsold,lbsold,lbrold,ltrold,lasold,larold,
     &        nreqa
      save    mpireqa,mpireqb,
     &        ilold,ltsold,lbsold,lbrold,ltrold,lasold,larold,
     &        nreqa
c
      data ilold,ltsold,lbsold,lbrold,ltrold,lasold,larold
     &     / 0,0,0,0,0,0,0 /
#endif /* MPI */
c
c --- split large requests into smaller pieces
c
      if     (ld-l1+1.gt.lsize) then
        do k= l1,ld,lsize
          l = min(k+lsize-1,ld)
          call xctilr(a,k,l,mh,nh,itype)
        enddo
        return
      endif
c
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0(12)
        nxc = 12
      endif
#endif
c
      mhl = max(0,min(mh,nbdy))
      nhl = max(0,min(nh,nbdy))
c
      if     (itype.lt.10) then
        sarc =  1.0
      else
        sarc = -1.0
      endif
c
      if     (ipr.eq.1 .and. jpr.eq.1) then
        do k= l1,ld
          do j= 1,nhl
            do i= 1,ii
              a(i,1-j,k) = vland
            enddo
          enddo
          if     (itype.eq.1 .or. itype.eq.11) then
            do j= 0,nhl
              do i= 1,ii
                io = ii-mod(i-1,ii)
                if     (a(io,jj-1-j,k).ne.vland) then
                  a(i,jj+j,k) = sarc*a(io,jj-1-j,k)
                else
                  a(i,jj+j,k) = vland
                endif
              enddo !i
            enddo !j
          elseif (itype.eq.2 .or. itype.eq.12) then
            do i= ii/2+1,ii
              io = mod(ii-(i-1),ii)+1
              if     (a(io,jj,k).ne.vland) then
                a(i,jj,k) = sarc*a(io,jj,k)
              else
                a(i,jj,k) = vland
              endif
            enddo !i
            do j= 1,nhl
              do i= 1,ii
                io = mod(ii-(i-1),ii)+1
                if     (a(io,jj-j,k).ne.vland) then
                  a(i,jj+j,k) = sarc*a(io,jj-j,k)
                else
                  a(i,jj+j,k) = vland
                endif
              enddo !i
            enddo !j
          elseif (itype.eq.3 .or. itype.eq.13) then
            do j= 0,nhl
              do i= 1,ii
                io = mod(ii-(i-1),ii)+1
                if     (a(io,jj-1-j,k).ne.vland) then
                  a(i,jj+j,k) = sarc*a(io,jj-1-j,k)
                else
                  a(i,jj+j,k) = vland
                endif
              enddo !i
            enddo !j
          elseif (itype.eq.4 .or. itype.eq.14) then
            do i= ii/2+1,ii
              io = ii-mod(i-1,ii)
              if     (a(io,jj,k).ne.vland) then
                a(i,jj,k) = sarc*a(io,jj,k)
              else
                a(i,jj,k) = vland
              endif
            enddo !i
            do j= 1,nhl
              do i= 1,ii
                io = ii-mod(i-1,ii)
                if     (a(io,jj-j,k).ne.vland) then
                  a(i,jj+j,k) = sarc*a(io,jj-j,k)
                else
                  a(i,jj+j,k) = vland
                endif
              enddo !i
            enddo !j
          endif !itype
        enddo !k
      else
        lts = 0
        lbs = 0
        lbr = 0
        ltr = 0
        las = 0
        lar = 0
        if     (nproc.ne.jpr) then
          if     (nhl.gt.0) then
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              do k= l1,ld
                do j= 1,nhl
                  lts = lts + 1
                  ai(lts,1) = a(i,jj+1-j,k)
                  ai(lts,2) = a(i,     j,k)
                  ai(lts,3) = vland
                  ai(lts,4) = vland
                enddo !j
              enddo !k
            enddo !i
            lbs = lts
            lbr = lts
            ltr = lts
          endif
        elseif (itype.eq.1 .or. itype.eq.11) then !p-grid
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            io = ii+1-i !ii:1:-1
            do k= l1,ld
              do j= 0,nhl
                lts = lts + 1
                if     (a(io,jj-1-j,k).ne.vland) then
                  ai(lts,1) = sarc*a(io,jj-1-j,k)
                else
                  ai(lts,1) = vland
                endif
                ai(lts,4) = vland
              enddo !j
              do j= 1,nhl
                lbs = lbs + 1
                ai(lbs,2) = a(i,j,k)
                ai(lbs,3) = vland
              enddo !j
            enddo !k
          enddo !i
          lbr = lbs
          ltr = lts
        elseif (itype.eq.3 .or. itype.eq.13) then  !u-grid
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            io = ii+2-i !ii+1:2:-1
            do k= l1,ld
              do j= 0,nhl
                lts = lts + 1
                if     (a(io,jj-1-j,k).ne.vland) then
                  ai(lts,1) = sarc*a(io,jj-1-j,k)
                else
                  ai(lts,1) = vland
                endif
                ai(lts,4) = vland
              enddo !j
              do j= 1,nhl
                lbs = lbs + 1
                ai(lbs,2) = a(i,j,k)
                ai(lbs,3) = vland
              enddo !j
            enddo !k
          enddo !i
          lbr = lbs
          ltr = lts
          do k= l1,ld
            do j= 0,nhl
              las = las + 1
              if     (a(1,jj-1-j,k).ne.vland) then
                aia(las,1) = sarc*a(1,jj-1-j,k)
              else
                aia(las,1) = vland
              endif
              aia(las,2) = vland
            enddo !j
          enddo !k
          lar = las
        elseif (itype.eq.2 .or. itype.eq.12) then  !q-grid
          if (mproc.gt.ipr/2) then
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              io = ii+2-i !ii+1:2:-1
              do k= l1,ld
                do j= 1,nhl
                  lts = lts + 1
                  if     (a(io,jj-j,k).ne.vland) then
                    ai(lts,1) = sarc*a(io,jj-j,k)
                  else
                    ai(lts,1) = vland
                  endif
                  ai(lts,2) = a(i,j,k)
                  ai(lts,3) = vland
                enddo !j
                do j= 0,nhl
                  ltr = ltr + 1
                  ai(ltr,4) = vland
                enddo !j
              enddo !k
            enddo !i
            lbs = lts
            lbr = lts
            do k= l1,ld
              do j= 0,nhl
                lar = lar + 1
                aia(lar,2) = vland
              enddo !j
            enddo !k
          else
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              io = ii+2-i !ii+1:2:-1
              do k= l1,ld
                do j= 0,nhl
                  lts = lts + 1
                  if     (a(io,jj-j,k).ne.vland) then
                    ai(lts,1) = sarc*a(io,jj-j,k)
                  else
                    ai(lts,1) = vland
                  endif
                enddo !j
                do j= 1,nhl
                  lbs = lbs + 1
                  ai(lbs,2) = a(i,j,k)
                  ai(lbs,3) = vland
                  ai(lbs,4) = vland
                enddo !j
              enddo !k
            enddo !i
            lbr = lbs
            ltr = lbs
            do k= l1,ld
              do j= 1,nhl
                lar = lar + 1
                aia(lar,2) = vland
              enddo !j
            enddo !k
          endif
          if (mod(ipr+1-mproc,ipr)+1.le.ipr/2) then
            do k= l1,ld
              do j= 1,nhl
                las = las + 1
                if     (a(1,jj-j,k).ne.vland) then
                  aia(las,1) = sarc*a(1,jj-j,k)
                else
                  aia(las,1) = vland
                endif
              enddo !j
            enddo !k
          else
            do k= l1,ld
              do j= 0,nhl
                las = las + 1
                if     (a(1,jj-j,k).ne.vland) then
                  aia(las,1) = sarc*a(1,jj-j,k)
                else
                  aia(las,1) = vland
                endif
              enddo !j
            enddo !k
          endif
        else  !v-grid
          if (mproc.gt.ipr/2) then
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              io = ii+1-i  !ii:1:-1
              do k= l1,ld
                do j= 1,nhl
                  lts = lts + 1
                  if     (a(io,jj-j,k).ne.vland) then
                    ai(lts,1) = sarc*a(io,jj-j,k)
                  else
                    ai(lts,1) = vland
                  endif
                  ai(lts,2) = a(i,j,k)
                  ai(lts,3) = vland
                enddo !j
                do j= 0,nhl
                  ltr = ltr + 1
                  ai(ltr,4) = vland
                enddo !j
              enddo !k
            enddo !i
            lbs = lts
            lbr = lts
          else
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              io = ii+1-i  !ii:1:-1
              do k= l1,ld
                do j= 0,nhl
                  lts = lts + 1
                  if     (a(io,jj-j,k).ne.vland) then
                    ai(lts,1) = sarc*a(io,jj-j,k)
                  else
                    ai(lts,1) = vland
                  endif
                enddo !j
                do j= 1,nhl
                  lbs = lbs + 1
                  ai(lbs,2) = a(i,j,k)
                  ai(lbs,3) = vland
                  ai(lbs,4) = vland
                enddo !j
              enddo !k
            enddo !i
            lbr = lbs
            ltr = lbs
          endif
        endif  !itype
c
#if defined(MPI)
        if (lts.ne.ltsold .or. lbs.ne.lbsold .or. lbr.ne.lbrold .or.
     &      ltr.ne.ltrold .or. las.ne.lasold .or. lar.ne.larold) then
          if (ltsold + lbsold + lbrold + 
     &        ltrold + lasold + larold.ne.0) then
            do i= 1,nreqa
                call mpi_request_free(mpireqa(i), mpierr)
            enddo
          endif
          ltsold = lts
          lbsold = lbs
          lbrold = lbr
          ltrold = ltr
          lasold = las
          larold = lar
c
c         loop through all neigboring tiles.
c
          l = 0
          if (lts.gt.0) then
            do m= 1,mm_top
              l = l + 1
              if     (nproc.ne.jpr) then
                ls0 = i0_st(m)*nhl*(ld-l1+1)
                lm  = ii_st(m)*nhl*(ld-l1+1)
                call mpi_send_init(
     &            ai(ls0+1,1),lm,MTYPER,
     &            idproc(m0_top+m,nproc+1), 9905,
     &            mpicomm, mpireqa(l), mpierr)
              else !arctic
                call mpi_send_init(
     &            ai(1,1),lts,MTYPER,
     &            idproc(m0_top+m,nproc+1), 99051,
     &            mpicomm, mpireqa(l), mpierr)
              endif
            enddo
          endif ! lts.gt.0
          if (las.gt.0) then
            l = l + 1
            call mpi_send_init(
     &        aia(1,1),las,MTYPER,
     &        idproc(mod(ipr+1-mproc,ipr)+1,nproc), 99052,
     &        mpicomm, mpireqa(l), mpierr)
          endif ! las.gt.0
          if (lbs.gt.0) then
            do m= 1,mm_bot
              l   = l + 1
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call mpi_send_init(
     &          ai(ls0+1,2),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9906,
     &          mpicomm, mpireqa(l), mpierr)
            enddo
          endif ! lbs.gt.0
          if (ltr.gt.0) then
            do m= 1,mm_top
              l   = l + 1
              if     (nproc.ne.jpr) then
                ls0 = i0_st(m)*nhl*(ld-l1+1)
                lm  = ii_st(m)*nhl*(ld-l1+1)
                call mpi_recv_init(
     &            ai(ls0+1,4),lm,MTYPER,
     &            idproc(m0_top+m,nproc+1), 9906,
     &            mpicomm, mpireqa(l), mpierr)
              else !arctic
                call mpi_recv_init(
     &            ai(1,4),ltr,MTYPER,
     &            idproc(m0_top+m,nproc+1), 99051,
     &            mpicomm, mpireqa(l), mpierr)
              endif
            enddo
          endif ! ltr.gt.0
          if (lar.gt.0) then
            l = l + 1
            call mpi_recv_init(
     &        aia(1,2),lar,MTYPER,
     &        idproc(mod(ipr+1-mproc,ipr)+1,nproc), 99052,
     &        mpicomm, mpireqa(l), mpierr)
          endif ! lar.gt.0
          if (lbr.gt.0) then
            do m= 1,mm_bot
              l   = l + 1
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call mpi_recv_init(
     &          ai(ls0+1,3),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9905,
     &          mpicomm, mpireqa(l), mpierr)
            enddo
          endif ! lbr.gt.0
          nreqa = l
        endif
        if     (nreqa.gt.0) then
          call mpi_startall(nreqa, mpireqa,          mpierr)
          call mpi_waitall( nreqa, mpireqa, mpistat, mpierr)
        endif
#elif defined(SHMEM)
        BARRIER
c
c       loop through all neigboring tiles.
c
        if (ltr.gt.0) then
          do m= 1,mm_top
            if     (idproc(m0_top+m,nproc+1).ne.null_tile) then
              if     (nproc.ne.jpr) then
                lg0 = i0_gt(m)*nhl*(ld-l1+1)
                ls0 = i0_st(m)*nhl*(ld-l1+1)
                lm  = ii_st(m)*nhl*(ld-l1+1)
                call SHMEM_GETR(ai(ls0+1,4),
     &                          ai(lg0+1,2),lm,
     &                          idproc(m0_top+m,nproc+1))
              else !arctic
                call SHMEM_GETR(ai(1,4),
     &                          ai(1,1),ltr,  !buffer 1
     &                          idproc(m0_top+m,nproc+1))
              endif
            endif
          enddo
        endif ! ltr.gt.0
        if (lar.gt.0) then
          call SHMEM_GETR(aia(1,2),
     &                    aia(1,1),lar,
     &                    idproc(mod(ipr+1-mproc,ipr)+1,nproc))
        endif ! lar.gt.0
c
        if (lbr.gt.0) then
          do m= 1,mm_bot
            if     (idproc(m0_bot+m,nproc-1).ne.null_tile) then
              lg0 = i0_gb(m)*nhl*(ld-l1+1)
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call SHMEM_GETR(ai(ls0+1,3),
     &                        ai(lg0+1,1),lm, idproc(m0_bot+m,nproc-1))
            endif
          enddo
        endif ! lbr.gt.0
#endif  /* MPI:SHMEM */
c
        if     (nproc.ne.jpr) then  !arctic
          if (lbr.gt.0) then
            lbr = 0
            do i= 1,ii  ! outer loop to simplify multiple neighbor case
              do k= l1,ld
                do j= 1,nhl
                  lbr = lbr + 1
                  a(i, 1-j,k) = ai(lbr,3)
                  a(i,jj+j,k) = ai(lbr,4)
                enddo
              enddo
            enddo
          endif ! lbr.gt.0
        else !arctic
          if (lar.gt.0) then
            lar  = 0
            if (itype.eq.3 .or. itype.eq.13 .or. mproc.gt.ipr/2) then
              do k= l1,ld
                do j= 0,nhl
                  lar  = lar + 1
                  ai(lar,4) = aia(lar,2)
                enddo !j
              enddo !k
            else
              do k= l1,ld
                do j= 1,nhl
                  lar  = lar + 1
                  ai(lar,4) = aia(lar,2)
                enddo !j
              enddo !k
            endif
          endif ! lar.gt.0
          if (ltr.gt.0) then
            if     (lbr.eq.0) then
              ltr = 0
              do i= 1,ii  ! outer loop to simplify multiple neighbor case
                do k= l1,ld
                  do j= 0,nhl
                    ltr = ltr + 1
                    a(i,jj+j,k) = ai(ltr,4)
                  enddo
                enddo
              enddo
            elseif (ltr.gt.lbr) then
              lbr = 0
              ltr = 0
              do i= 1,ii  ! outer loop to simplify multiple neighbor case
                do k= l1,ld
                  do j= 0,nhl
                    ltr = ltr + 1
                    a(i,jj+j,k) = ai(ltr,4)
                  enddo
                  do j= 1,nhl
                    lbr = lbr + 1
                    a(i, 1-j,k) = ai(lbr,3)
                  enddo
                enddo
              enddo
            else
              lbr = 0
              do i= 1,ii  ! outer loop to simplify multiple neighbor case
                do k= l1,ld
                  do j= 1,nhl
                    lbr = lbr + 1
                    a(i, 1-j,k) = ai(lbr,3)
                    a(i,jj+j,k) = ai(lbr,4)
                  enddo
                enddo
              enddo
            endif
          endif ! ltr.gt.0
        endif !arctic
      endif  ! ipr.eq.1 .and. jpr.eq.1
c
      if     (mhl.gt.0) then
        if     (ipr.eq.1) then
          if     (nreg.eq.0.or.nreg.eq.4) then
            do k= l1,ld
              do j= 1-nhl,jj+nhl
                do i= 1,mhl
                  a( 1-i,j,k) = vland
                  a(ii+i,j,k) = vland
                enddo
              enddo
            enddo
          else
            do k= l1,ld
              do j= 1-nhl,jj+nhl
                do i= 1,mhl
                  a( 1-i,j,k) = a(ii+1-i,j,k)
                  a(ii+i,j,k) = a(     i,j,k)
                enddo
              enddo
            enddo
          endif
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
              enddo
            enddo
          enddo
*         write(lp,'(a,6i6)') 'xctilr - mhl,l1,ld,jj,l,mnproc = ',
*    &                                  mhl,l1,ld,jj,l,mnproc
*         call xcsync(flush_lp)
#if defined(MPI)
#if defined(MPISRI)
          call mpi_sendrecv(
     &          aj(1,1),l,MTYPER,idhalo(2), 9907,
     &          aj(1,3),l,MTYPER,idhalo(1), 9907,
     &          mpicomm, mpistat, mpierr)
          call mpi_sendrecv(
     &          aj(1,2),l,MTYPER,idhalo(1), 9908,
     &          aj(1,4),l,MTYPER,idhalo(2), 9908,
     &          mpicomm, mpistat, mpierr)
#else /* MPISRI */
          if     (ilold.ne.l) then
            if     (ilold.ne.0) then
              call mpi_request_free(mpireqb(1), mpierr)
              call mpi_request_free(mpireqb(2), mpierr)
              call mpi_request_free(mpireqb(3), mpierr)
              call mpi_request_free(mpireqb(4), mpierr)
            endif
            ilold = l
            call mpi_send_init(
     &            aj(1,1),l,MTYPER,idhalo(2), 9907,
     &            mpicomm, mpireqb(1), mpierr)
            call mpi_send_init(
     &            aj(1,2),l,MTYPER,idhalo(1), 9908,
     &            mpicomm, mpireqb(2), mpierr)
            call mpi_recv_init(
     &            aj(1,3),l,MTYPER,idhalo(1), 9907,
     &            mpicomm, mpireqb(3), mpierr)
            call mpi_recv_init(
     &            aj(1,4),l,MTYPER,idhalo(2), 9908,
     &            mpicomm, mpireqb(4), mpierr)
          endif
          call mpi_startall(4, mpireqb,          mpierr)
          call mpi_waitall( 4, mpireqb, mpistat, mpierr)
#endif /* MPISRI */
#elif defined(SHMEM)
          BARRIER_MP
          if     (idhalo(1).ne.null_tile) then
            call SHMEM_GETR(aj(1,3),
     &                      aj(1,1),l,idhalo(1))
          endif
          if     (idhalo(2).ne.null_tile) then
            call SHMEM_GETR(aj(1,4),
     &                      aj(1,2),l,idhalo(2))
          endif
          BARRIER_MP
#endif
          l = 0
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                l = l + 1
                a( 1-i,j,k) = aj(l,3)
                a(ii+i,j,k) = aj(l,4)
              enddo
            enddo
          enddo
        endif  ! ipr.eq.1:else
      endif  ! mhl.gt.0
#if defined(TIMER)
c
      if     (nxc.eq.12) then
        call xctmr1(12)
        nxc = 0
      endif
#endif
      return
      end subroutine xctilr
#else /* !ARCTIC */
      recursive subroutine xctilr(a,l1,ld,mh,nh,itype)
      implicit none
c
      integer, intent(in)    :: l1,ld,mh,nh,itype
      real,    intent(inout) :: a(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ld)
c
c**********
c*
c  1) update the tile overlap halo of a real array.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    a               real           in/out    target array
c    l1              integer        input     3rd dim. start index
c    ld              integer        input     3rd dimension of a
c    mh              integer        input     1st (EW) update halo size
c    nh              integer        input     2nd (NS) update halo size
c    itype           integer        input     grid and field type
c
c  3) itype selects both the grid and field type
c        itype= 1; p-grid, scalar field
c        itype= 2; q-grid, scalar field
c        itype= 3; u-grid, scalar field
c        itype= 4; v-grid, scalar field
c        itype=11; p-grid, vector field
c        itype=12; q-grid, vector field
c        itype=13; u-grid, vector field
c        itype=14; v-grid, vector field
c     it is ignored here because all types are the same unless
c      the grid includes the arctic ocean
c
c  4) the global variable vland is returned by halos over land.
c*
c**********
c
      integer    lsize,ilen,jlen
      parameter (lsize=kdm,
     &           ilen= idm        *lsize*nbdy+64,
     &           jlen=(jdm+2*nbdy)*lsize*nbdy+64)
c
c     halo buffer (in common for enhanced MPI safety).
c
      real            ai,aj
      common/xctilr4/ ai(ilen,4),aj(jlen,4)
      save  /xctilr4/
c
      integer i,j,k,l,lg0,ls0,lst0,lsb0,lm,ltm,lbm,m,mhl,nhl
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
c
c     persistent communication handles.
c
      integer mpireqa(4*iqr),mpireqb(4),ilold,jlold,nreqa
      save    mpireqa,mpireqb,ilold,jlold,nreqa
c
      data ilold,jlold / 0,0 /
#endif /* MPI */
c
c --- split large requests into smaller pieces
c
      if     (ld-l1+1.gt.lsize) then
        do k= l1,ld,lsize
          l = min(k+lsize-1,ld)
          call xctilr(a,k,l,mh,nh,itype)
        enddo
        return
      endif
c
#if defined(TIMER)
c
      if     (nxc.eq.0) then
        call xctmr0(12)
        nxc = 12
      endif
#endif
c
      mhl = max(0,min(mh,nbdy))
      nhl = max(0,min(nh,nbdy))
c
      if     (nhl.gt.0) then
        if     (jpr.eq.1) then
          if     (nreg.eq.0.or.nreg.eq.1) then
            do k= l1,ld
              do j= 1,nhl
                do i= 1,ii
                  a(i, 1-j,k) = vland
                  a(i,jj+j,k) = vland
                enddo
              enddo
            enddo
          else
            do k= l1,ld
              do j= 1,nhl
                do i= 1,ii
                  a(i, 1-j,k) = a(i,jj+1-j,k)
                  a(i,jj+j,k) = a(i,     j,k)
                enddo
              enddo
            enddo
          endif
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
              enddo
            enddo
          enddo
*         write(lp,'(a,6i6)') 'xctilr - nhl,l1,ld,ii,l,mnproc = ',
*    &                                  nhl,l1,ld,ii,l,mnproc
*         call xcsync(flush_lp)
#if defined(MPI)
#if defined(MPISRJ)
c
          lst0 = i0_st(1)*nhl*(ld-l1+1)
          ltm  = ii_st(1)*nhl*(ld-l1+1)
          lsb0 = i0_sb(1)*nhl*(ld-l1+1)
          lbm  = ii_sb(1)*nhl*(ld-l1+1)
          call mpi_sendrecv(
     &      ai(lst0+1,1),ltm,MTYPER,idproc(m0_top+1,nproc+1), 9905,
     &      ai(lsb0+1,3),lbm,MTYPER,idproc(m0_bot+1,nproc-1), 9905,
     &      mpicomm, mpistat, mpierr)
          call mpi_sendrecv(
     &      ai(lsb0+1,2),lbm,MTYPER,idproc(m0_bot+1,nproc-1), 9906,
     &      ai(lst0+1,4),ltm,MTYPER,idproc(m0_top+1,nproc+1), 9906,
     &      mpicomm, mpistat, mpierr)
#else /* MPISRJ */
c
          if     (jlold.ne.l) then
            if     (jlold.ne.0) then
              do i= 1,nreqa
                call mpi_request_free(mpireqa(i), mpierr)
              enddo
            endif
            jlold = l
c
c           loop through all neigboring tiles.
c
            l = 0
            do m= 1,mm_top
              l   = l + 1
              ls0 = i0_st(m)*nhl*(ld-l1+1)
              lm  = ii_st(m)*nhl*(ld-l1+1)
              call mpi_send_init(
     &          ai(ls0+1,1),lm,MTYPER,idproc(m0_top+m,nproc+1), 9905,
     &          mpicomm, mpireqa(l), mpierr)
            enddo
            do m= 1,mm_bot
              l   = l + 1
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call mpi_send_init(
     &          ai(ls0+1,2),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9906,
     &          mpicomm, mpireqa(l), mpierr)
            enddo
            do m= 1,mm_top
              l   = l + 1
              ls0 = i0_st(m)*nhl*(ld-l1+1)
              lm  = ii_st(m)*nhl*(ld-l1+1)
              call mpi_recv_init(
     &          ai(ls0+1,4),lm,MTYPER,idproc(m0_top+m,nproc+1), 9906,
     &          mpicomm, mpireqa(l), mpierr)
            enddo
            do m= 1,mm_bot
              l   = l + 1
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call mpi_recv_init(
     &          ai(ls0+1,3),lm,MTYPER,idproc(m0_bot+m,nproc-1), 9905,
     &          mpicomm, mpireqa(l), mpierr)
            enddo
            nreqa = l
          endif
          if     (nreqa.gt.0) then
            call mpi_startall(nreqa, mpireqa,          mpierr)
            call mpi_waitall( nreqa, mpireqa, mpistat, mpierr)
          endif
#endif /* MPISRJ */
#elif defined(SHMEM)
          BARRIER
c
c         loop through all neigboring tiles.
c
          do m= 1,mm_top
            if     (idproc(m0_top+m,nproc+1).ne.null_tile) then
              lg0 = i0_gt(m)*nhl*(ld-l1+1)
              ls0 = i0_st(m)*nhl*(ld-l1+1)
              lm  = ii_st(m)*nhl*(ld-l1+1)
              call SHMEM_GETR(ai(ls0+1,4),
     &                        ai(lg0+1,2),lm, idproc(m0_top+m,nproc+1))
            endif
          enddo
c
          do m= 1,mm_bot
            if     (idproc(m0_bot+m,nproc-1).ne.null_tile) then
              lg0 = i0_gb(m)*nhl*(ld-l1+1)
              ls0 = i0_sb(m)*nhl*(ld-l1+1)
              lm  = ii_sb(m)*nhl*(ld-l1+1)
              call SHMEM_GETR(ai(ls0+1,3),
     &                        ai(lg0+1,1),lm, idproc(m0_bot+m,nproc-1))
            endif
          enddo
#endif  /* MPI:SHMEM */
c
          l = 0
          do i= 1,ii  ! outer loop to simplify multiple neighbor case
            do k= l1,ld
              do j= 1,nhl
                l = l + 1
                a(i, 1-j,k) = ai(l,3)
                a(i,jj+j,k) = ai(l,4)
              enddo
            enddo
          enddo
        endif  ! jpr.eq.1:else
      endif  ! nhl.gt.0
c
      if     (mhl.gt.0) then
        if     (ipr.eq.1) then
          if     (nreg.eq.0.or.nreg.eq.4) then
            do k= l1,ld
              do j= 1-nhl,jj+nhl
                do i= 1,mhl
                  a( 1-i,j,k) = vland
                  a(ii+i,j,k) = vland
                enddo
              enddo
            enddo
          else
            do k= l1,ld
              do j= 1-nhl,jj+nhl
                do i= 1,mhl
                  a( 1-i,j,k) = a(ii+1-i,j,k)
                  a(ii+i,j,k) = a(     i,j,k)
                enddo
              enddo
            enddo
          endif
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
              enddo
            enddo
          enddo
*         write(lp,'(a,6i6)') 'xctilr - mhl,l1,ld,jj,l,mnproc = ',
*    &                                  mhl,l1,ld,jj,l,mnproc
*         call xcsync(flush_lp)
#if defined(MPI)
#if defined(MPISRI)
          call mpi_sendrecv(
     &          aj(1,1),l,MTYPER,idhalo(2), 9907,
     &          aj(1,3),l,MTYPER,idhalo(1), 9907,
     &          mpicomm, mpistat, mpierr)
          call mpi_sendrecv(
     &          aj(1,2),l,MTYPER,idhalo(1), 9908,
     &          aj(1,4),l,MTYPER,idhalo(2), 9908,
     &          mpicomm, mpistat, mpierr)
#else /* MPISRI */
          if     (ilold.ne.l) then
            if     (ilold.ne.0) then
              call mpi_request_free(mpireqb(1), mpierr)
              call mpi_request_free(mpireqb(2), mpierr)
              call mpi_request_free(mpireqb(3), mpierr)
              call mpi_request_free(mpireqb(4), mpierr)
            endif
            ilold = l
            call mpi_send_init(
     &            aj(1,1),l,MTYPER,idhalo(2), 9907,
     &            mpicomm, mpireqb(1), mpierr)
            call mpi_send_init(
     &            aj(1,2),l,MTYPER,idhalo(1), 9908,
     &            mpicomm, mpireqb(2), mpierr)
            call mpi_recv_init(
     &            aj(1,3),l,MTYPER,idhalo(1), 9907,
     &            mpicomm, mpireqb(3), mpierr)
            call mpi_recv_init(
     &            aj(1,4),l,MTYPER,idhalo(2), 9908,
     &            mpicomm, mpireqb(4), mpierr)
          endif
          call mpi_startall(4, mpireqb,          mpierr)
          call mpi_waitall( 4, mpireqb, mpistat, mpierr)
#endif /* MPISRI */
#elif defined(SHMEM)
          BARRIER_MP
          if     (idhalo(1).ne.null_tile) then
            call SHMEM_GETR(aj(1,3),
     &                      aj(1,1),l,idhalo(1))
          endif
          if     (idhalo(2).ne.null_tile) then
            call SHMEM_GETR(aj(1,4),
     &                      aj(1,2),l,idhalo(2))
          endif
          BARRIER_MP
#endif
          l = 0
          do k= l1,ld
            do j= 1-nhl,jj+nhl
              do i= 1,mhl
                l = l + 1
                a( 1-i,j,k) = aj(l,3)
                a(ii+i,j,k) = aj(l,4)
              enddo
            enddo
          enddo
        endif  ! ipr.eq.1:else
      endif  ! mhl.gt.0
#if defined(TIMER)
c
      if     (nxc.eq.12) then
        call xctmr1(12)
        nxc = 0
      endif
#endif
      return
      end subroutine xctilr
#endif /* ARCTIC:else */

      subroutine xctmri
      implicit none
c
c**********
c*
c  1) initialize timers.
c
c  2) timers  1:32 are for message passing routines,
c     timers 33:80 are for general hycom routines,
c     timers 81:96 are for user selected routines.
c     timer     97 is the total time.
c
c  3) call xctmri    to initialize timers (called in xcspmd),
c     call xctmr0(n) to start timer n,
c     call xctmr1(n) to stop  timer n and add event to timer sum,
c     call xctnrn(n,cname) to register a name for timer n,
c     call xctmrp to printout timer statistics (called by xcstop).
c
c  4) time every 50-th event above 1,000.
c*
c**********
c
      integer i
c
      real*8     zero8
      parameter (zero8=0.0)
c
      nxc = 0
      do i= 1,97
        cc(i) = '      '
        nc(i) = 0
        tc(i) = zero8
      enddo
c
      call xctmrn(97,'total ')
      call xctmr0(97)
      return
      end subroutine xctmri

      subroutine xctmr0(n)
      implicit none
c
      integer, intent(in) :: n
c
c**********
c*
c  1) start timer n.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    n               integer        input     timer number
c
c  3) time every 50-th event above 1,000.
c*
c**********
c
      real*8 wtime
c
#if defined(DEBUG_TIMER_ALL)
      if     (              cc(n).ne.'      ') then
        write(lp,'(i5,2x,a,a)') mnproc,'call ',cc(n)
        call flush(lp)
      endif
#endif
#if defined(DEBUG_TIMER)
      if     (n.gt.32 .and. cc(n).ne.'      ') then
        if     (mnproc.eq.1) then
        write(lp,*) 'call ',cc(n)
        call flush(lp)
        endif
      endif
#endif
      if     (timer_on) then
        if     (mod(nc(n),50).eq.0 .or. nc(n).le.1000) then
          t0(n) = wtime()
        endif
      endif !timer_on
      return
      end subroutine xctmr0

      subroutine xctmr1(n)
      implicit none
c
      integer, intent(in) :: n
c
c**********
c*
c  1) add time since call to xctim0 to timer n.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    n               integer        input     timer number
c
c  3) time every 50-th event above 1,000.
c*
c**********
c
      real*8  wtime
c
      if     (timer_on) then
        if     (nc(n).gt.1000) then
          if     (mod(nc(n),50).eq.0) then
            tc(n) = tc(n) + 50.0*(wtime() - t0(n))
          endif
        else
          tc(n) = tc(n) + (wtime() - t0(n))
        endif
        nc(n) = nc(n) + 1
      endif !timer_on
#if defined(DEBUG_TIMER_ALL)
      if     (              cc(n).ne.'      ') then
        write(lp,'(i5,2x,a,a)') mnproc,'exit ',cc(n)
        call flush(lp)
      endif
#endif
#if defined(DEBUG_TIMER)
      if     (n.gt.32 .and. cc(n).ne.'      ') then
        if     (mnproc.eq.1) then
        write(lp,*) 'exit ',cc(n)
        call flush(lp)
        endif
      endif
#endif
      return
      end subroutine xctmr1

      subroutine xctmrn(n,cname)
      implicit none
c
      character*6, intent(in) :: cname
      integer,     intent(in) :: n
c
c**********
c*
c  1) register name of timer n.
c
c  2) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    n               integer        input     timer number
c    cname           char*(6)       input     timer name
c*
c**********
c
      cc(n) = cname
      return
      end subroutine xctmrn

      subroutine xctmrp
      implicit none
c
c**********
c*
c  1) print all active timers.
c
c  2) on exit all timers are reset to zero.
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
#endif
c
      integer i,mn,mnloc
c
      real*8     zero8
      parameter (zero8=0.0)
c
c     get total time.
c
      call xctmr1(97)
c
c     report time on the processor with the least communication overhead
c
      tcxc(2) = mnproc
      tcxc(1) = zero8
      do i= 1,32
        if     (nc(i).ne.0 .and. cc(i)(1:2).eq.'xc') then
          tcxc(1) = tcxc(1) + tc(i)  !communication overhead
        endif
      enddo !i
c
      if     (ijpr.ne.1) then
#if defined(MPI)
        tcxl(1) = tcxc(1)
        tcxl(2) = tcxc(2)
        call mpi_allreduce(tcxl,tcxc,1,
     &                     mpi_2double_precision,mpi_minloc,
     &                     mpicomm,mpierr)
        mnloc = tcxc(2)  !processor with the least comm. overhead
        if     (mnproc.eq.1) then
          if     (mnloc.ne.1) then
            call MPI_RECV(tc,97,MTYPED,
     &                    idproc1(mnloc), 9949,
     &                    mpicomm, mpistat, mpierr)
          endif
        elseif (mnproc.eq.mnloc) then
          call MPI_SEND(tc,97,MTYPED,
     &                  idproc1(1), 9949,
     &                  mpicomm, mpierr)
        endif
#elif defined(SHMEM)
        BARRIER
        if     (mnproc.eq.1) then
          mnloc = 1
          do mn= 2,ijpr
            call SHMEM_GETD(tcxc(2),tcxc(1),1,idproc1(mn))
            if     (tcxc(2).gt.tcxc(1)) then
              tcxc(1) = tcxc(2)
              mnloc   = mn
            endif
          enddo
          tcxc(2) = mnloc  !processor with the least comm. overhead
        endif
        if     (mnloc.ne.1) then
          call SHMEM_GETD(tc,tc,97,idproc1(mnloc))
        endif
        BARRIER
#endif
      endif
c
      call xcsync(flush_lp)
      if     (mnproc.eq.1) then
        write(lp,6000) mnloc,ijpr
        do i= 1,32
          if     (nc(i).ne.0) then
            if     (cc(i).ne.'      ') then
              write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
            else
              write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
            endif
          endif
        enddo !i
        write(lp,6100) 'xc****',1,tcxc(1),tcxc(1)
        do i= 33,97
          if     (nc(i).ne.0) then
            if     (cc(i).ne.'      ') then
              write(lp,6100) cc(i),nc(i),tc(i),tc(i)/nc(i)
            else
              write(lp,6150)    i, nc(i),tc(i),tc(i)/nc(i)
            endif
          endif
        enddo !i
        write(lp,6200)
      endif !mnproc.eq.1
      call xcsync(flush_lp)
c
c     reset timers to zero.
c
      do i= 1,97
        nc(i) = 0
        tc(i) = zero8
      enddo
      tcxc(1) = zero8
c
c     start a new total time measurement.
c
      call xctmr0(97)
      return
c
 6000 format(/ /
     &    3x,' timer statistics, processor',i5,' out of',i5 /
     &    3x,'-----------------------------------------------' /)
 6100 format(3x,a6,
     &   '   calls =',i9,
     &   '   time =',f11.5,
     &   '   time/call =',f14.8)
 6150 format(3x,'   #',i2,
     &   '   calls =',i9,
     &   '   time =',f11.5,
     &   '   time/call =',f14.8)
 6200 format(/ /)
      end subroutine xctmrp
c-----------------------------------------------------------------------
c
      subroutine xcgetrow(outm, inm, kt)
      implicit none
 
c
      integer, intent(in)    :: kt
      real,    intent(out)   :: outm(itdm,jj,kt)
      real,    intent(in)    :: inm(ii,jj,kt)
c
c**********
c*
c  1) convert an entire 2-D array from tiled to non-tiled layout.
c
c  3) mnflg selects which nodes must return the array
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  4) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aa              real           output    non-tiled target array
c    a               real           input     tiled source array
c    mnflg           integer        input     node return flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
c
      integer        mpireqa(jpr),mpireqb(ipr)
#endif
      real  at(idm*jdm*kt),ata(idm*jdm*kt,iqr)
      integer i,j,k,l,mp,np,mnp
c
c     gather each row of tiles onto the first tile in the row.
c
#if defined(MPI)
      if     (mproc.eq.mpe_1(nproc)) then
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            outm(i0+i,j,k) = inm(i,j,k)
          enddo
        enddo
      enddo
        l = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = l + 1
          call MPI_IRECV(ata(1,l),ii_pe(mp,nproc)*jj*kt,MTYPER,
     &                   idproc(mp,nproc), 9941,
     &                   mpicomm, mpireqb(l), mpierr)
        enddo
        BARRIER
        l = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = l + 1
          call MPI_WAIT(mpireqb(l), mpistat, mpierr)
        do k= 1,kt
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
            outm(i0_pe(mp,nproc)+i,j,k) = ata(i+(j-1)*ii_pe(mp,nproc)
     &                                   +(k-1)*jj*ii_pe(mp,nproc),l)
            enddo
          enddo
       enddo
        enddo
      else  !mproc>1
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            at(i+(j-1)*ii+(k-1)*ii*jj) = inm(i,j,k)
          enddo
        enddo
      enddo
        BARRIER
        call MPI_SEND(at,ii*jj*kt,MTYPER,
     &                idproc(mpe_1(nproc),nproc), 9941,
     &                mpicomm, mpierr)
      endif
#endif
      return
      end subroutine xcgetrow

c-----------------------------------------------------------------------
c
      subroutine xcgetrow4(outm,inm, kt)
      implicit none
c
      integer, intent(in)      :: kt
      real*4,    intent(out)   :: outm(itdm,jj,kt)
      real*4,    intent(in)    :: inm(ii,jj,kt)
c
c**********
c*
c  1) convert an entire 2-D array from tiled to non-tiled layout.
c
c  3) mnflg selects which nodes must return the array
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  4) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aa              real           output    non-tiled target array
c    a               real           input     tiled source array
c    mnflg           integer        input     node return flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
c
      integer        mpireqa(jpr),mpireqb(ipr)
#endif
      real*4  at(idm*jdm*2*kk),ata(idm*jdm*2*kk,iqr)
      integer i,j,k,l,mp,np,mnp
c
c     gather each row of tiles onto the first tile in the row.
c
#if defined(MPI)
      if     (mproc.eq.mpe_1(nproc)) then
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            outm(i0+i,j,k) = inm(i,j,k)
          enddo
        enddo
      enddo
        l = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = l + 1
          call MPI_IRECV(ata(1,l),ii_pe(mp,nproc)*jj*kt,MTYPE4,
     &                   idproc(mp,nproc), 9941,
     &                   mpicomm, mpireqb(l), mpierr)
        enddo
        BARRIER
        l = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = l + 1
          call MPI_WAIT(mpireqb(l), mpistat, mpierr)
        do k= 1,kt
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
            outm(i0_pe(mp,nproc)+i,j,k) = ata(i+(j-1)*ii_pe(mp,nproc)
     &                                   +(k-1)*jj*ii_pe(mp,nproc),l)
            enddo
          enddo
       enddo
        enddo
      else  !mproc>1
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            at(i+(j-1)*ii+(k-1)*ii*jj) = inm(i,j,k)
          enddo
        enddo
      enddo
        BARRIER
        call MPI_SEND(at,ii*jj*kt,MTYPE4,
     &                idproc(mpe_1(nproc),nproc), 9941,
     &                mpicomm, mpierr)
      endif
#endif
      return
      end subroutine xcgetrow4

c-----------------------------------------------------------------------
      subroutine xcgetrowint2(outm,inm, kt)
      implicit none
c
      integer, intent(in)      :: kt
      integer*2,    intent(out)   :: outm(itdm,jj,kt)
      integer*2,    intent(in)    :: inm(ii,jj,kt)
c
c**********
c*
c  1) convert an entire 2-D array from tiled to non-tiled layout.
c
c  3) mnflg selects which nodes must return the array
c        = 0; all nodes
c        = n; node number n (mnproc=n)
c
c  4) parameters:
c       name            type         usage            description
c    ----------      ----------     -------  ----------------------------
c    aa              real           output    non-tiled target array
c    a               real           input     tiled source array
c    mnflg           integer        input     node return flag
c*
c**********
c
#if defined(MPI)
      include 'mpif.h'
      integer        mpicomm,mpierr,mpireq,mpistat
      common/xcmpii/ mpicomm,mpierr,mpireq(4),
     &               mpistat(mpi_status_size,4*max(iqr,jqr))
      save  /xcmpii/
c
      integer        mpireqa(jpr),mpireqb(ipr)
#endif
      integer*2  at(idm*jdm*kt),ata(idm*jdm*kt,iqr)
      integer i,j,k,l,mp,np,mnp
c
c     gather each row of tiles onto the first tile in the row.
c
#if defined(MPI)
      if     (mproc.eq.mpe_1(nproc)) then
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            outm(i0+i,j,k) = inm(i,j,k)
          enddo
        enddo
      enddo
        l = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = l + 1
          call MPI_IRECV(ata(1,l),ii_pe(mp,nproc)*jj*kt,MTYPEI2,
     &                   idproc(mp,nproc), 9941,
     &                   mpicomm, mpireqb(l), mpierr)
        enddo
        BARRIER
        l = 0
        do mp= mpe_1(nproc)+1,mpe_e(nproc)
          l = l + 1
          call MPI_WAIT(mpireqb(l), mpistat, mpierr)
        do k= 1,kt
          do j= 1,jj
            do i= 1,ii_pe(mp,nproc)
            outm(i0_pe(mp,nproc)+i,j,k) = ata(i+(j-1)*ii_pe(mp,nproc)
     &                                   +(k-1)*jj*ii_pe(mp,nproc),l)
            enddo
          enddo
       enddo
        enddo
      else  !mproc>1
      do k= 1,kt
        do j= 1,jj
          do i= 1,ii
            at(i+(j-1)*ii+(k-1)*ii*jj) = inm(i,j,k)
          enddo
        enddo
      enddo
        BARRIER
        call MPI_SEND(at,ii*jj*kt,MTYPEI2,
     &                idproc(mpe_1(nproc),nproc), 9941,
     &                mpicomm, mpierr)
      endif
#endif
      return
      end subroutine xcgetrowint2


