! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger, I. Kriest,
!                     A. Moree, C. Heinze
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

module mo_vertical_fluxes

  implicit none
  private

  public :: sinking

contains

  subroutine get_ws(i,j,k,wpoc,wcal,wopal,wdust)
    !***********************************************************************************************
    ! Select sinking velocity for particulate tracers
    !
    use mo_kind,         only: rp
    use mo_biomod,       only: wmass
    use mo_param_bgc,    only: wmin,wmax,wlin,wcal_const,wdust_const,wopal_const,wpoc_const,dustsink
    use mo_vgrid,        only: ptiestu
    use mo_control_bgc,  only: use_AGG,use_WLIN,use_M4AGO
    use mo_ihamocc4m4ago,only: ws_agg

    ! Arguments
    integer, intent(in)    :: i,j,k                 ! indices for spatial domain to be used
    real(rp),intent(inout) :: wpoc,wcal,wopal,wdust ! sinking velocities of particulates

    if (use_AGG) then
       wpoc  = wmass(i,j,k)
       wcal  = wmass(i,j,k)
       wopal = wmass(i,j,k)
       wdust = dustsink
    else if (use_WLIN) then
       wpoc  = min(wmin+wlin*ptiestu(i,j,k), wmax)
       wcal  = wcal_const
       wopal = wopal_const
       wdust = wdust_const
    else if (use_M4AGO) then
       wpoc   = ws_agg(i,j,k)
       wcal   = ws_agg(i,j,k)
       wopal  = ws_agg(i,j,k)
       wdust  = ws_agg(i,j,k)
    else
       wpoc   = wpoc_const
       wcal   = wcal_const
       wopal  = wopal_const
       wdust  = wdust_const
    endif

  end subroutine get_ws


  subroutine sinking(kpie,kpje,kpke,pddpo,omask)
    !***********************************************************************************************
    ! Particulates sinking and sedimentation
    !
    use mo_kind,          only: rp
    use mo_carbch,        only: ocetra
    use mo_param_bgc,     only: rcar,wmin
    use mo_param1_bgc,    only: idet,icalc,iopal,ifdust,iphy,inos,iadust,idet13,idet14,icalc13,    &
                                icalc14,isco212,isco213,isco214,ialkali,inatalkali,inatsco212,     &
                                isilica,iphy13,iphy14,inatcalc
    use mo_control_bgc,   only: use_sedbypass,use_AGG,use_WLIN,use_M4AGO,use_natDIC,use_cisonew,   &
                                use_PBGC_OCNP_TIMESTEP
    use mo_sedmnt,        only: prcaca,produs,prorca,silpro,pror13,pror14,prca13,prca14,prnatcaca
    use mo_vgrid,         only: dp_min,dp_min_sink,k0100,k0500,k1000,k2000,k4000
    use mo_biomod,        only: bsiflx0100,bsiflx0500,bsiflx1000,bsiflx2000,bsiflx4000,bsiflx_bot, &
                                calflx0100,calflx0500,calflx1000,calflx2000,calflx4000,calflx_bot, &
                                carflx0100,carflx0500,carflx1000,carflx2000,carflx4000,carflx_bot, &
                                dustflx0100,dustflx0500,dustflx1000,dustflx2000,dustflx4000,       &
                                dustflx_bot,aggregate,dustagg
    use mo_ihamocc4m4ago, only: ws_agg

    ! Arguments
    integer, intent(in) :: kpie                         ! 1st dimension of model grid.
    integer, intent(in) :: kpje                         ! 2nd dimension of model grid.
    integer, intent(in) :: kpke                         ! 3rd (vertical) dimension of model grid.
    real(rp),intent(in) :: pddpo(kpie,kpje,kpke)        ! size of grid cell (3rd dimension) [m].
    real(rp),intent(in) :: omask(kpie,kpje)             ! land/ocean mask (1=ocean)

    ! Local parameters and variables
    integer, parameter :: nsinkmax = 12                 ! maximum number of sinking tracers
    integer  :: i,j,k
    integer  :: is,kdonor
    real(rp) :: dz
    real(rp) :: tco(nsinkmax),tcn(nsinkmax),q(nsinkmax)
    real(rp) :: wpoc, wcal, wopal, wdust
    real(rp) :: wpocd,wcald,wopald,wdustd
    ! use_AGG
    real(rp) :: wnos,wnosd,dagg
    ! sedbypass
    real(rp) :: florca,flcaca,flsil,flnatcaca
    real(rp) :: flor13,flor14,flca13,flca14

    ! Set output fields to zero
    carflx0100 (:,:) = 0._rp
    carflx0500 (:,:) = 0._rp
    carflx1000 (:,:) = 0._rp
    carflx2000 (:,:) = 0._rp
    carflx4000 (:,:) = 0._rp
    carflx_bot (:,:) = 0._rp
    bsiflx0100 (:,:) = 0._rp
    bsiflx0500 (:,:) = 0._rp
    bsiflx1000 (:,:) = 0._rp
    bsiflx2000 (:,:) = 0._rp
    bsiflx4000 (:,:) = 0._rp
    bsiflx_bot (:,:) = 0._rp
    calflx0100 (:,:) = 0._rp
    calflx0500 (:,:) = 0._rp
    calflx1000 (:,:) = 0._rp
    calflx2000 (:,:) = 0._rp
    calflx4000 (:,:) = 0._rp
    calflx_bot (:,:) = 0._rp
    dustflx0100(:,:) = 0._rp
    dustflx0500(:,:) = 0._rp
    dustflx1000(:,:) = 0._rp
    dustflx2000(:,:) = 0._rp
    dustflx4000(:,:) = 0._rp
    dustflx_bot(:,:) = 0._rp

    ! implicit method for sinking of particles:
    ! C(k,T+dt)=C(k,T) + (w*dt/ddpo(k))*(C(k-1,T+1)-C(k,T+1))
    ! -->
    ! C(k,T+dt)=(ddpo(k)*C(k,T)+w*dt*C(k-1,T+dt))/(ddpo(k)+w*dt)
    ! sedimentation=w*dt*C(ks,T+dt)
    !
    !$OMP PARALLEL DO PRIVATE(kdonor,wpoc,wpocd,wcal,wcald,wopal,wopald,wdust,wdustd,tco,tcn,q,wnos,wnosd,dagg,i,k) ORDERED
    do j = 1,kpje
      do i = 1,kpie

        tco(:) = 0.0_rp
        tcn(:) = 0.0_rp

        if(omask(i,j) > 0.5_rp) then

          kdonor = 1
          do k = 1,kpke
            !$OMP ORDERED
            ! Sum up total column inventory before sinking scheme
            if( pddpo(i,j,k) > dp_min ) then
              tco( 1) = tco( 1) + ocetra(i,j,k,idet  )*pddpo(i,j,k)
              tco( 2) = tco( 2) + ocetra(i,j,k,icalc )*pddpo(i,j,k)
              if (use_natDIC) then
                tco( 3) = tco( 3) + ocetra(i,j,k,inatcalc)*pddpo(i,j,k)
              endif
              tco( 4) = tco( 4) + ocetra(i,j,k,iopal )*pddpo(i,j,k)
              tco( 5) = tco( 5) + ocetra(i,j,k,ifdust)*pddpo(i,j,k)
              if (use_AGG) then
                tco( 6) = tco( 6) + ocetra(i,j,k,iphy  )*pddpo(i,j,k)
                tco( 7) = tco( 7) + ocetra(i,j,k,inos  )*pddpo(i,j,k)
                tco( 8) = tco( 8) + ocetra(i,j,k,iadust)*pddpo(i,j,k)
              endif
              if (use_cisonew) then
                tco( 9) = tco( 9) + ocetra(i,j,k,idet13 )*pddpo(i,j,k)
                tco(10) = tco(10) + ocetra(i,j,k,idet14 )*pddpo(i,j,k)
                tco(11) = tco(11) + ocetra(i,j,k,icalc13)*pddpo(i,j,k)
                tco(12) = tco(12) + ocetra(i,j,k,icalc14)*pddpo(i,j,k)
              endif
            endif

            if(pddpo(i,j,k) > dp_min_sink) then
              call get_ws(i,j,k,     wpoc, wcal, wopal, wdust)
              call get_ws(i,j,kdonor,wpocd,wcald,wopald,wdustd)
              if (use_AGG) then
                dagg   = dustagg(i,j,k)
              else
                dagg   = 0._rp
              endif

              if( k == 1 ) then
                wpocd  = 0.0_rp
                wcald  = 0.0_rp
                wopald = 0.0_rp
                wdustd = 0.0_rp
                if (use_AGG) then
                  wnosd  = 0.0_rp
                else if (use_WLIN) then
                  wpoc = wmin
                else if (use_M4AGO) then
                  wpoc = ws_agg(i,j,k)
                endif
              endif

              ocetra(i,j,k,idet)   = (ocetra(i,j,k,     idet)  *pddpo(i,j,k)                       &
                                   +  ocetra(i,j,kdonor,idet)  *wpocd) / (pddpo(i,j,k)+wpoc)
              ocetra(i,j,k,icalc)  = (ocetra(i,j,k,     icalc) *pddpo(i,j,k)                       &
                   &               +  ocetra(i,j,kdonor,icalc) *wcald) / (pddpo(i,j,k)+wcal)
              ocetra(i,j,k,iopal)  = (ocetra(i,j,k,     iopal) *pddpo(i,j,k)                       &
                   &               +  ocetra(i,j,kdonor,iopal) *wopald)/ (pddpo(i,j,k)+wopal)
              ocetra(i,j,k,ifdust) = (ocetra(i,j,k,     ifdust)*pddpo(i,j,k)                       &
                   &               +  ocetra(i,j,kdonor,ifdust)*wdustd)/ (pddpo(i,j,k)+wdust)      &
                                   -  dagg
              if (use_cisonew) then
                ocetra(i,j,k,idet13)  = (ocetra(i,j,k,     idet13) *pddpo(i,j,k)                   &
                     &                +  ocetra(i,j,kdonor,idet13) *wpocd) / (pddpo(i,j,k)+wpoc)
                ocetra(i,j,k,idet14)  = (ocetra(i,j,k,     idet14) *pddpo(i,j,k)                   &
                     &                +  ocetra(i,j,kdonor,idet14) *wpocd) / (pddpo(i,j,k)+wpoc)
                ocetra(i,j,k,icalc13) = (ocetra(i,j,k,     icalc13)*pddpo(i,j,k)                   &
                     &                +  ocetra(i,j,kdonor,icalc13)*wcald) / (pddpo(i,j,k)+wcal)
                ocetra(i,j,k,icalc14) = (ocetra(i,j,k,     icalc14)*pddpo(i,j,k)                   &
                     &                +  ocetra(i,j,kdonor,icalc14)*wcald) / (pddpo(i,j,k)+wcal)
              endif
              if (use_natDIC) then
                ocetra(i,j,k,inatcalc)= (ocetra(i,j,k,     inatcalc)*pddpo(i,j,k)                  &
                     &                +  ocetra(i,j,kdonor,inatcalc)*wcald) / (pddpo(i,j,k)+wcal)
              endif
              if (use_AGG) then
                ocetra(i,j,k,iphy)    = (ocetra(i,j,k,     iphy)*pddpo(i,j,k)                      &
                     &                +  ocetra(i,j,kdonor,iphy)*wpocd) / (pddpo(i,j,k)+wpoc)
                ocetra(i,j,k,inos)    = (ocetra(i,j,k,     inos)*pddpo(i,j,k)                      &
                     &                +  ocetra(i,j,kdonor,inos)*wnosd) / (pddpo(i,j,k)+wnos)      &
                                      -  aggregate(i,j,k)
                ocetra(i,j,k,iadust)  = (ocetra(i,j,k,    iadust)*pddpo(i,j,k)                     &
                     &                +  ocetra(i,j,kdonor,iadust)*wpocd) / (pddpo(i,j,k)+wpoc)    &
                                      +  dagg
              endif
              kdonor = k

            else if( pddpo(i,j,k) > dp_min ) then

              ocetra(i,j,k,idet)   = ocetra(i,j,kdonor,idet)
              ocetra(i,j,k,icalc)  = ocetra(i,j,kdonor,icalc)
              if (use_cisonew) then
                ocetra(i,j,k,idet13) = ocetra(i,j,kdonor,idet13)
                ocetra(i,j,k,idet14) = ocetra(i,j,kdonor,idet14)
                ocetra(i,j,k,icalc13) = ocetra(i,j,kdonor,icalc13)
                ocetra(i,j,k,icalc14) = ocetra(i,j,kdonor,icalc14)
              endif
              if (use_natDIC) then
                ocetra(i,j,k,inatcalc) = ocetra(i,j,kdonor,inatcalc)
              endif
              ocetra(i,j,k,iopal)  = ocetra(i,j,kdonor,iopal)
              ocetra(i,j,k,ifdust) = ocetra(i,j,kdonor,ifdust)
              if (use_AGG) then
                ocetra(i,j,k,iphy)   = ocetra(i,j,kdonor,iphy)
                ocetra(i,j,k,inos)   = ocetra(i,j,kdonor,inos)
                ocetra(i,j,k,iadust) = ocetra(i,j,kdonor,iadust)
              endif

            endif  ! pddpo > dp_min_sink

            ! Sum up total column inventory after sinking scheme
            ! flux to sediment added after kpke-loop
            if( pddpo(i,j,k) > dp_min ) then
              tcn( 1) = tcn( 1) + ocetra(i,j,k,idet  )*pddpo(i,j,k)
              tcn( 2) = tcn( 2) + ocetra(i,j,k,icalc )*pddpo(i,j,k)
              if (use_natDIC) then
                tcn( 3) = tcn( 3) + ocetra(i,j,k,inatcalc)*pddpo(i,j,k)
              endif
              tcn( 4) = tcn( 4) + ocetra(i,j,k,iopal )*pddpo(i,j,k)
              tcn( 5) = tcn( 5) + ocetra(i,j,k,ifdust)*pddpo(i,j,k)
              if (use_AGG) then
                tcn( 6) = tcn( 6) + ocetra(i,j,k,iphy  )*pddpo(i,j,k)
                tcn( 7) = tcn( 7) + ocetra(i,j,k,inos  )*pddpo(i,j,k)
                tcn( 8) = tcn( 8) + ocetra(i,j,k,iadust)*pddpo(i,j,k)
              endif
              if (use_cisonew) then
                tcn( 9) = tcn( 9) + ocetra(i,j,k,idet13 )*pddpo(i,j,k)
                tcn(10) = tcn(10) + ocetra(i,j,k,idet14 )*pddpo(i,j,k)
                tcn(11) = tcn(11) + ocetra(i,j,k,icalc13)*pddpo(i,j,k)
                tcn(12) = tcn(12) + ocetra(i,j,k,icalc14)*pddpo(i,j,k)
              endif
            endif
            !$OMP END ORDERED
          enddo  ! loop k=1,kpke


          ! Add fluxes to sediment to new total column inventory
          tcn( 1) = tcn( 1) + ocetra(i,j,kdonor,idet  )*wpoc
          tcn( 2) = tcn( 2) + ocetra(i,j,kdonor,icalc )*wcal
          if (use_natDIC) then
            tcn( 3) = tcn( 3) + ocetra(i,j,kdonor,inatcalc)*wcal
          endif
          tcn( 4) = tcn( 4) + ocetra(i,j,kdonor,iopal )*wopal
          tcn( 5) = tcn( 5) + ocetra(i,j,kdonor,ifdust)*wdust
          if (use_AGG) then
            tcn( 6) = tcn( 6) + ocetra(i,j,kdonor,iphy  )*wpoc
            tcn( 7) = tcn( 7) + ocetra(i,j,kdonor,inos  )*wnos
            tcn( 8) = tcn( 8) + ocetra(i,j,kdonor,iadust)*wpoc
          endif
          if (use_cisonew) then
            tcn( 9) = tcn( 9) + ocetra(i,j,kdonor,idet13 )*wpoc
            tcn(10) = tcn(10) + ocetra(i,j,kdonor,idet14 )*wpoc
            tcn(11) = tcn(11) + ocetra(i,j,kdonor,icalc13)*wcal
            tcn(12) = tcn(12) + ocetra(i,j,kdonor,icalc14)*wcal
          endif

          ! Do columnwise multiplicative mass conservation correction
          q(:) = 1.0_rp
          do is = 1,nsinkmax
            if( tco(is) > 1.e-12_rp .and. tcn(is) > 1.e-12_rp ) q(is) = tco(is)/tcn(is)
          enddo
          do k = 1,kpke
            if( pddpo(i,j,k) > dp_min ) then
              ocetra(i,j,k,idet  ) = ocetra(i,j,k,idet  )*q(1)
              ocetra(i,j,k,icalc ) = ocetra(i,j,k,icalc )*q(2)
              if (use_natDIC) then
                ocetra(i,j,k,inatcalc) = ocetra(i,j,k,inatcalc)*q(3)
              endif
              ocetra(i,j,k,iopal ) = ocetra(i,j,k,iopal )*q(4)
              ocetra(i,j,k,ifdust) = ocetra(i,j,k,ifdust)*q(5)
              if (use_AGG) then
                ocetra(i,j,k,iphy  ) = ocetra(i,j,k,iphy  )*q(6)
                ocetra(i,j,k,inos  ) = ocetra(i,j,k,inos  )*q(7)
                ocetra(i,j,k,iadust) = ocetra(i,j,k,iadust)*q(8)
              endif
              if (use_cisonew) then
                ocetra(i,j,k,idet13 ) = ocetra(i,j,k,idet13 )*q(9)
                ocetra(i,j,k,idet14 ) = ocetra(i,j,k,idet14 )*q(10)
                ocetra(i,j,k,icalc13) = ocetra(i,j,k,icalc13)*q(11)
                ocetra(i,j,k,icalc14) = ocetra(i,j,k,icalc14)*q(12)
              endif
            endif
          enddo

          ! Fluxes to sediment, layers thinner than dp_min_sink are ignored.
          ! Note that kdonor=kbo(i,j) by definition since kbo is the lowermost
          ! layer thicker than dp_min_sink.
          if (use_AGG) then
            prorca(i,j) = ocetra(i,j,kdonor,iphy  )*wpoc  + ocetra(i,j,kdonor,idet  )*wpoc
            prcaca(i,j) = ocetra(i,j,kdonor,icalc )*wcal
            silpro(i,j) = ocetra(i,j,kdonor,iopal )*wopal
            produs(i,j) = ocetra(i,j,kdonor,ifdust)*wdust + ocetra(i,j,kdonor,iadust)*wpoc

            if (use_cisonew) then
              pror13(i,j) = ocetra(i,j,kdonor,iphy13)*wpoc + ocetra(i,j,kdonor,idet13)*wpoc
              pror14(i,j) = ocetra(i,j,kdonor,iphy14)*wpoc + ocetra(i,j,kdonor,idet14)*wpoc
              prca13(i,j) = ocetra(i,j,kdonor,icalc13)*wcal
              prca14(i,j) = ocetra(i,j,kdonor,icalc14)*wcal
            endif
            if (use_natDIC) then
              prnatcaca(i,j) = ocetra(i,j,kdonor,inatcalc)*wcal
            endif
          else
            prorca(i,j) = ocetra(i,j,kdonor,idet  )*wpoc
            prcaca(i,j) = ocetra(i,j,kdonor,icalc )*wcal
            silpro(i,j) = ocetra(i,j,kdonor,iopal )*wopal
            produs(i,j) = ocetra(i,j,kdonor,ifdust)*wdust
            if (use_cisonew) then
              pror13(i,j) = ocetra(i,j,kdonor,idet13 )*wpoc
              prca13(i,j) = ocetra(i,j,kdonor,icalc13)*wcal
              pror14(i,j) = ocetra(i,j,kdonor,idet14 )*wpoc
              prca14(i,j) = ocetra(i,j,kdonor,icalc14)*wcal
            endif
            if (use_natDIC) then
              prnatcaca(i,j) = ocetra(i,j,kdonor,inatcalc)*wcal
            endif
          endif

        endif  ! omask > 0.5
      enddo    ! loop i=1,kpie
    enddo    ! loop j=1,kpje
    !$OMP END PARALLEL DO


    ! Calculate mass sinking flux for carbon, opal and calcium carbonate
    ! through the 100 m, 500 m, 1000 m, 2000 m, and 4000 m depth surfaces. These
    ! fluxes are intentionally calculated using values at the NEW timelevel
    ! to be fully consistent with the implicit sinking scheme

    !$OMP PARALLEL DO PRIVATE(i,k,wpoc,wcal,wopal)
    do j = 1,kpje
      do i = 1,kpie
        if(omask(i,j) > 0.5_rp) then

          ! 100 m
          k = k0100(i,j)
          if(k > 0) then
            call get_ws(i,j,k,wpoc,wcal,wopal,wdust)

            if (use_AGG) then
              carflx0100(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx0100(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx0100(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx0100(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx0100(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx0100(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 500 m
          k = k0500(i,j)
          if(k > 0) then
            call get_ws(i,j,k,wpoc,wcal,wopal,wdust)

            if (use_AGG) then
              carflx0500(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx0500(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx0500(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx0500(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx0500(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx0500(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 1000 m
          k = k1000(i,j)
          if(k > 0) then
            call get_ws(i,j,k,wpoc,wcal,wopal,wdust)

            if (use_AGG) then
              carflx1000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx1000(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx1000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx1000(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx1000(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx1000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 2000 m
          k = k2000(i,j)
          if(k > 0) then
            call get_ws(i,j,k,wpoc,wcal,wopal,wdust)

            if (use_AGG) then
              carflx2000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx2000(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx2000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx2000(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx2000(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx2000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! 4000 m
          k = k4000(i,j)
          if(k > 0) then
            call get_ws(i,j,k,wpoc,wcal,wopal,wdust)

            if (use_AGG) then
              carflx4000(i,j) = (ocetra(i,j,k,idet)+ocetra(i,j,k,iphy))*rcar*wpoc
              dustflx4000(i,j)= ocetra(i,j,k,ifdust)*wdust + ocetra(i,j,k,iadust)*wpoc
            else
              carflx4000(i,j) = ocetra(i,j,k,idet)*rcar*wpoc
              dustflx4000(i,j)= ocetra(i,j,k,ifdust)*wdust
            endif
            bsiflx4000(i,j) = ocetra(i,j,k,iopal)*wopal
            calflx4000(i,j) = ocetra(i,j,k,icalc)*wcal
          endif

          ! bottom fluxes
          carflx_bot(i,j) = prorca(i,j)*rcar
          bsiflx_bot(i,j) = silpro(i,j)
          calflx_bot(i,j) = prcaca(i,j)
          dustflx_bot(i,j)= produs(i,j)

        endif ! omask > 0.5
      enddo
    enddo
    !$OMP END PARALLEL DO

    if (use_sedbypass) then

      ! If sediment bypass is activated, fluxes to the sediment are distributed
      ! over the water column. Detritus is kept as detritus, while opal and CaCO3
      ! are remineralised instantanously

      !$OMP PARALLEL DO PRIVATE(dz,florca,flcaca,flnatcaca,flsil,flor13,flor14,flca13,flca14,i,k) ORDERED
      do j=1,kpje
        do i = 1,kpie
          if(omask(i,j) > 0.5_rp) then

            ! calculate depth of water column
            dz = 0.0_rp
            do k = 1,kpke
              !$OMP ORDERED
              if( pddpo(i,j,k) > dp_min ) dz = dz+pddpo(i,j,k)
              !$OMP END ORDERED
            enddo

            florca = prorca(i,j)/dz
            flcaca = prcaca(i,j)/dz
            flsil = silpro(i,j)/dz
            prorca(i,j) = 0._rp
            prcaca(i,j) = 0._rp
            silpro(i,j) = 0._rp
            if (use_cisonew) then
              flor13 = pror13(i,j)/dz
              flor14 = pror13(i,j)/dz
              flca13 = prca13(i,j)/dz
              flca14 = prca14(i,j)/dz
              pror13(i,j) = 0._rp
              pror14(i,j) = 0._rp
              prca13(i,j) = 0._rp
              prca14(i,j) = 0._rp
            endif
            if (use_natDIC) then
              flnatcaca = prnatcaca(i,j)/dz
            endif

            do k = 1,kpke
              if( pddpo(i,j,k) <= dp_min ) cycle

              ocetra(i,j,k,idet)    = ocetra(i,j,k,idet)+florca
              ocetra(i,j,k,ialkali) = ocetra(i,j,k,ialkali)+2._rp*flcaca
              ocetra(i,j,k,isco212) = ocetra(i,j,k,isco212)+flcaca
              ocetra(i,j,k,isilica) = ocetra(i,j,k,isilica)+flsil
              if (use_cisonew) then
                ocetra(i,j,k,idet13)  = ocetra(i,j,k,idet13)+flor13
                ocetra(i,j,k,idet14)  = ocetra(i,j,k,idet14)+flor14
                ocetra(i,j,k,isco213) = ocetra(i,j,k,isco213)+flca13
                ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)+flca14
              endif
              if (use_natDIC) then
                ocetra(i,j,k,inatalkali) = ocetra(i,j,k,inatalkali)+2._rp*flnatcaca
                ocetra(i,j,k,inatsco212) = ocetra(i,j,k,inatsco212)+flnatcaca
              endif
            enddo ! k=1,kpke

          endif ! omask > 0.5
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif ! use_sedbypass

  end subroutine sinking
end module mo_vertical_fluxes

