! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, M. Bentsen
!
! This file is part of BLOM/iHAMOCC.
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
! along with BLOM. If not, see https://www.gnu.org/licenses/.


SUBROUTINE hamocc_step(m,n,mm,nn,k1m,k1n)
!
! --- ------------------------------------------------------------------
! --- perform one HAMOCC step
! --- ------------------------------------------------------------------
!
  use mod_xc,         only: idm,jdm,kdm,nbdy
  use mod_time,       only: date,nday_of_year,nstep,nstep_in_day
  use mod_grid,       only: plat
  use mod_state,      only: temp,saln
  use mod_forcing,    only: swa,slp,abswnd,atmco2,flxco2,flxdms,                &
       &                    atmbrf,flxbrf
  use mod_seaice,     only: ficem
  use mo_bgcmean,     only: nbgc,bgcwrt, diagfq_bgc,diagmon_bgc,                &
       &                    diagann_bgc
  use mo_intfcblom,   only: bgc_dx,bgc_dy,bgc_dp,bgc_rho,omask,                 &
       &                    blom2hamocc,hamocc2blom
  use mo_read_rivin,  only: rivflx
  use mo_read_fedep,  only: get_fedep
  use mo_read_ndep,   only: get_ndep
  use mo_read_oafx,   only: get_oafx
  use mo_read_pi_ph,  only: get_pi_ph,pi_ph
  use mo_control_bgc, only: with_dmsph

  implicit none

  integer, intent(in) :: m,n,mm,nn,k1m,k1n

  integer :: l,ldtday
  real    :: ndep(idm,jdm)
  real    :: dust(idm,jdm)
  real    :: oafx(idm,jdm)      

  call trc_limitc(nn)

  call blom2hamocc(m,n,mm,nn)

  ldtday = mod(nstep,nstep_in_day)

  do l=1,nbgc
     bgcwrt(l)=.false.
     if (((diagann_bgc(l).and.nday_of_year.eq.1.or.diagmon_bgc(l)               &
          &   .and.date%day.eq.1).and.mod(nstep,nstep_in_day).eq.0).or.         &
          &   .not.(diagann_bgc(l).or.diagmon_bgc(l)).and.                      &
          &   mod(nstep+.5,diagfq_bgc(l)).lt.1.)                                &
          &   bgcwrt(l)=.true.
  enddo

  call get_fedep(idm,jdm,date%month,dust)
  call get_ndep(idm,jdm,date%year,date%month,omask,ndep)
  call get_oafx(idm,jdm,date%year,date%month,omask,oafx)
  if(with_dmsph) call get_pi_ph(idm,jdm,date%month)

  call hamocc4bcm(idm,jdm,kdm,nbdy,                                             &
       &   date%year,date%month,date%day,ldtday,                                &
       &   bgc_dx,bgc_dy,bgc_dp,bgc_rho,plat,omask,                             &
       &   dust,rivflx,ndep,oafx,pi_ph,                                         &
       &   swa,ficem,slp,abswnd,                                                &
       &   temp(1-nbdy,1-nbdy,1+nn),saln(1-nbdy,1-nbdy,1+nn),                   &
       &   atmco2,flxco2,flxdms,atmbrf,flxbrf)

  !
  ! --- accumulate fields and write output
  !
  call accfields(idm,jdm,kdm,bgc_dx,bgc_dy,bgc_dp,omask)

  call hamocc2blom(m,n,mm,nn)

END SUBROUTINE hamocc_step
