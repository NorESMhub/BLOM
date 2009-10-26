      program micom
c
c --- ---------------------------------------------
c --- miami isopycnic coordinate ocean model
c
c --- NERSC version 2.0, based on MICOM version 2.8
c --- Contact person: mats@nersc.no
c --- ---------------------------------------------
c
      use mod_xc
      use data_mct
c
      implicit none
c
#include "common_blocks.h"
#ifdef MPI
#  include "mpif.h"
#endif
c
      integer i
c
#ifdef MPI
      call mpi_init
      mpicom_mct=mpi_comm_world
#endif
c
      open (unit=nfu,file='ocn_in',status='old',action='read')
      read (nfu,'(i6)') i
      read (nfu,*)
      read (nfu,*)
      read (nfu,*)
      read (nfu,'(a)') runid_mct    ! Experiment name
      close (unit=nfu)
      if (i.eq.0) then
        runtype_mct='initial'
      else
        runtype_mct='continue'
      endif
c
c --- initialize the model
      call micom_init
c
c --- advance the model from time step nstep1 to nstep2
      micom_loop: do
        call micom_step
c       if (nstep.eq.nstep2) exit micom_loop
        if (nstep.eq.nstep1+2) exit micom_loop
      enddo micom_loop
c
c --- write check sum of layer thickness
      call chksummsk(dp(1-nbdy,1-nbdy,1+mod(nstep2,2)*kk),ip,1,'dp')
c
      if (mnproc.eq.1) then
        open (unit=nfu,file='run.status',status='unknown')
        write (nfu,*) 'success'
        close (unit=nfu)
      endif
      call xcstop('(normal)')
             stop '(normal)'
c
      end
