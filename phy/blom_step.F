! ------------------------------------------------------------------------------
! Copyright (C) 2008-2020 Mats Bentsen, Mehmet Ilicak
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

      subroutine blom_step
c
c --- ------------------------------------------------------------------
c --- integrate a model time step
c --- ------------------------------------------------------------------
c
      use mod_config, only: expcnf
      use mod_time, only: date, nday_of_year, nstep1,
     .                    nstep, nstep_in_day, delt1, step_time
      use mod_timing, only: total_time, total_xio_time,
     .                      auxil_total_time, getfrc_total_time,
     .                      tmsmt1_total_time, advdif_total_time,
     .                      sfcstr_total_time, momtum_total_time,
     .                      pgforc_total_time, barotp_total_time,
     .                      pbcor2_total_time, convec_total_time,
     .                      diapfl_total_time, thermf_total_time,
     .                      mxlayr_total_time, tmsmt2_total_time,
     .                      diaacc_total_time, io_total_time,
     .                      get_time
      use mod_xc, only: lp, mnproc, xctilr, xcsum
      use mod_swabs, only: updswa
      use mod_tmsmt, only: tmsmt1, tmsmt2
      use mod_eddtra, only: eddtra
      use mod_advect, only: advect
      use mod_pbcor, only: pbcor1, pbcor2
      use mod_pgforc, only: pgforc
      use mod_momtum, only: momtum
      use mod_mxlayr, only: mxlayr
      use mod_barotp, only: barotp
      use mod_cmnfld, only: cmnfld
      use mod_forcing, only: fwbbal
      use mod_budget, only: budget_sums, budget_output
      use mod_eddtra, only: eddtra
      use mod_momtum, only: momtum
      use mod_chkvar, only: chkvar
      use mod_dia
c
      implicit none
c
      real q
      integer i,m,n,mm,nn,k1m,k1n
c
      real total_step_time,
     .     auxil_time ,
     .     getfrc_time,
     .     tmsmt1_time,
     .     advdif_time,
     .     sfcstr_time,
     .     momtum_time,
     .     pgforc_time,
     .     barotp_time,
     .     pbcor2_time,
     .     convec_time,
     .     diapfl_time,
     .     thermf_time,
     .     mxlayr_time,
     .     tmsmt2_time,
     .     diaacc_time,
     .     io_time
c
c --- letter 'm' refers to mid-time level (example: dp(i,j,km) )
c --- letter 'n' refers to old and new time level
c
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
      mm=(m-1)*kk
      nn=(n-1)*kk
      k1m=1+mm
      k1n=1+nn
c
      call budget_sums(1,n,nn)
c
      call step_time
c
c --- ------------------------------------------------------------------
c --- Update some flux halos the first time step of a day to reproduce
c --- results after restart with arctic
c --- ------------------------------------------------------------------
c
      if (nreg.eq.2.and.mod(nstep,nstep_in_day).eq.1) then
        if (mnproc.eq.1) write (lp,*) 'blom_step: update flux halos'
        call xctilr(uflx(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_uv)
        call xctilr(utflx(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_uv)
        call xctilr(usflx(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_uv)
        call xctilr(vflx(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_vv)
        call xctilr(vtflx(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_vv)
        call xctilr(vsflx(1-nbdy,1-nbdy,k1n), 1,kk, 1,1, halo_vv)
      endif
c
      auxil_time=get_time()
c
c --- ------------------------------------------------------------------
c --- Get forcing
c --- ------------------------------------------------------------------
c
      call getfrc
c
c --- ------------------------------------------------------------------
c --- Update arrays related to shortwave radiation absorption.
c --- ------------------------------------------------------------------
c
      call updswa
c
      getfrc_time=get_time()
c
cdiag write (lp,*) 'tmsmt1...'
      call tmsmt1(m,n,mm,nn,k1m,k1n)
      tmsmt1_time=get_time()
c
cdiag write (lp,*) 'advdif...'
      call difest(m,n,mm,nn,k1m,k1n)
      call eddtra(m,n,mm,nn,k1m,k1n)
      call advect(m,n,mm,nn,k1m,k1n)
      call pbcor1(m,n,mm,nn,k1m,k1n)
      call diffus(m,n,mm,nn,k1m,k1n)
      advdif_time=get_time()
c
      call budget_sums(2,n,nn)
      auxil_time=auxil_time+get_time()
c
cdiag write (lp,*) 'sfcstr...'
      call sfcstr(m,n,mm,nn,k1m,k1n)
      sfcstr_time=get_time()
c
cdiag write (lp,*) 'pgforc...'
      call pgforc(m,n,mm,nn,k1m,k1n)
      pgforc_time=get_time()
c
cdiag write (lp,*) 'momtum...'
      call momtum(m,n,mm,nn,k1m,k1n)
      momtum_time=get_time()
c
cdiag write (lp,*) 'convec...'
      call convec(m,n,mm,nn,k1m,k1n)
      convec_time=get_time()
c
      call budget_sums(3,n,nn)
      auxil_time=auxil_time+get_time()
c
cdiag write (lp,*) 'diapfl...'
      call diapfl(m,n,mm,nn,k1m,k1n)
      diapfl_time=get_time()
c
      call budget_sums(4,n,nn)
      auxil_time=auxil_time+get_time()
c
cdiag write (lp,*) 'thermf...'
      call thermf(m,n,mm,nn,k1m,k1n)
      thermf_time=get_time()
c
cdiag write (lp,*) 'mxlayr...'
      call mxlayr(m,n,mm,nn,k1m,k1n)
      mxlayr_time=get_time()
c
#ifdef TRC
c --- update tracer due to non-passive processes
      call updtrc(m,n,mm,nn,k1m,k1n)
#endif
c
      call budget_sums(5,n,nn)
      auxil_time=auxil_time+get_time()
c
cdiag write (lp,*) 'barotp...'
      call barotp(m,n,mm,nn,k1m,k1n)
      barotp_time=get_time()
c
cdiag write (lp,*) 'pbcor2...'
      call pbcor2(m,n,mm,nn,k1m,k1n)
      pbcor2_time=get_time()
c
      call budget_sums(6,m,mm)
      auxil_time=auxil_time+get_time()
c
cdiag write (lp,*) 'tmsmt2...'
      call tmsmt2(m,n,mm,nn,k1m,k1n)
      tmsmt2_time=get_time()
c
      call cmnfld(m,n,mm,nn,k1m,k1n)
c
      call budget_sums(7,m,mm)
c
      call fwbbal(m,n,mm,nn,k1m,k1n)
c
      call budget_output(m)
c
      auxil_time=auxil_time+get_time()
c
c ----------------------------------------------------------------------
c
c --- output and diagnostic calculations
c
c ----------------------------------------------------------------------
c
      call chkvar(m,n,mm,nn,k1m,k1n)
c
      call diaacc(m,n,mm,nn,k1m,k1n)
      diaacc_time=get_time()
c
      if (mod(nstep,nstep_in_day).eq.0.and.nday_of_year.eq.1) then
c
c --- ------------------------------------------------------------------
c --- - output diagnosed heat and salt flux
c --- ------------------------------------------------------------------
c
        call wdiflx
c
      endif
c
c
c --- ------------------------------------------------------------------
c --- - output of BLOM diagnostics
c --- ------------------------------------------------------------------
c
      do i=1,nphy 
        if (((diagann_phy(i).and.nday_of_year.eq.1.or.diagmon_phy(i)
     .        .and.date%day.eq.1).and.mod(nstep,nstep_in_day).eq.0).or.
     .      .not.(diagann_phy(i).or.diagmon_phy(i)).and.
     .      mod(nstep+.5,diagfq_phy(i)).lt.1.) 
     .    call diaout(i,m,n,mm,nn,k1m,k1n)
      enddo
c
c --- update total time spent by various tasks
      auxil_total_time=auxil_total_time+auxil_time
      getfrc_total_time=getfrc_total_time+getfrc_time
      tmsmt1_total_time=tmsmt1_total_time+tmsmt1_time
      advdif_total_time=advdif_total_time+advdif_time
      sfcstr_total_time=sfcstr_total_time+sfcstr_time
      momtum_total_time=momtum_total_time+momtum_time
      pgforc_total_time=pgforc_total_time+pgforc_time
      barotp_total_time=barotp_total_time+barotp_time
      pbcor2_total_time=pbcor2_total_time+pbcor2_time
      convec_total_time=convec_total_time+convec_time
      diapfl_total_time=diapfl_total_time+diapfl_time
      thermf_total_time=thermf_total_time+thermf_time
      mxlayr_total_time=mxlayr_total_time+mxlayr_time
      tmsmt2_total_time=tmsmt2_total_time+tmsmt2_time
      diaacc_total_time=diaacc_total_time+diaacc_time
c
      if (((rstann.and.nday_of_year.eq.1.or.rstmon.and.date%day.eq.1)
     .     .and.mod(nstep,nstep_in_day).eq.0).or.
     .    .not.(rstann.or.rstmon).and.mod(nstep+.5,rstfrq).lt.1.) then
c
        if (expcnf.ne.'cesm') then
c
c --- ------------------------------------------------------------------
c --- --- output restart files
c --- ------------------------------------------------------------------
c
          call restart_wt
c
        endif
c
        io_time=get_time()
c
c --- ------------------------------------------------------------------
c --- - write timing diagnostics to standard out
c --- ------------------------------------------------------------------
c
        io_total_time=io_total_time+io_time
        total_step_time=auxil_time +getfrc_time+tmsmt1_time+advdif_time
     .                 +sfcstr_time+momtum_time+pgforc_time+barotp_time
     .                 +pbcor2_time+convec_time+diapfl_time+thermf_time
     .                 +mxlayr_time+tmsmt2_time+diaacc_time+io_time
        total_time=total_time+total_step_time
        total_xio_time=total_xio_time+total_step_time-io_time
c
        if (mnproc.eq.1) then
          write (lp,'(f12.4,a,i8)')
     .      total_step_time, '  sec for step ', nstep
          write (lp,'(f12.4,a,i8)')
     .      total_time/(nstep-nstep1),' Avg Time'
          write (lp,'(f12.4,a,i8)')
     .      total_xio_time/(nstep-nstep1),' Avg Time excluding IO'
          write (lp,'(f12.4,a,i8)')
     .      total_time,' Tot Time with contributions:'
          q=100./total_time
          write (lp,'(f12.4,a,i8)') auxil_total_time*q ,'% auxil '
          write (lp,'(f12.4,a,i8)') getfrc_total_time*q,'% getfrc'
          write (lp,'(f12.4,a,i8)') tmsmt1_total_time*q,'% tmsmt1'
          write (lp,'(f12.4,a,i8)') advdif_total_time*q,'% advdif'
          write (lp,'(f12.4,a,i8)') sfcstr_total_time*q,'% sfcstr'
          write (lp,'(f12.4,a,i8)') momtum_total_time*q,'% momtum'
          write (lp,'(f12.4,a,i8)') pgforc_total_time*q,'% pgforc'
          write (lp,'(f12.4,a,i8)') barotp_total_time*q,'% barotp'
          write (lp,'(f12.4,a,i8)') pbcor2_total_time*q,'% pbcor2'
          write (lp,'(f12.4,a,i8)') convec_total_time*q,'% convec'
          write (lp,'(f12.4,a,i8)') diapfl_total_time*q,'% diapfl'
          write (lp,'(f12.4,a,i8)') thermf_total_time*q,'% thermf'
          write (lp,'(f12.4,a,i8)') mxlayr_total_time*q,'% mxlayr'
          write (lp,'(f12.4,a,i8)') tmsmt2_total_time*q,'% tmsmt2'
          write (lp,'(f12.4,a,i8)') diaacc_total_time*q,'% diaacc'
          write (lp,'(f12.4,a,i8)') io_total_time*q    ,'% IO'
        endif
c
      else
c
c --- ------------------------------------------------------------------
c --- - write time spent for current time step
c --- ------------------------------------------------------------------
c
        io_time=get_time()
        io_total_time=io_total_time+io_time
        total_step_time=auxil_time +getfrc_time+tmsmt1_time+advdif_time
     .                 +sfcstr_time+momtum_time+pgforc_time+barotp_time
     .                 +pbcor2_time+convec_time+diapfl_time+thermf_time
     .                 +mxlayr_time+tmsmt2_time+diaacc_time+io_time
        total_time=total_time+total_step_time
        total_xio_time=total_xio_time+total_step_time-io_time
c
        if (mnproc.eq.1) then
          write (lp,'(f12.4,a,i8)') total_step_time, '  sec for step ',
     .                              nstep
        endif
c
      endif
c
      delt1=baclin+baclin
c
      return
      end
