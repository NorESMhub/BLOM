! ------------------------------------------------------------------------------
! Copyright (C) 2002-2025 Mats Bentsen, Mehmet Ilicak
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

module mod_thermf_ben02

  use mod_constants, only: spcifh, t0deg, alpha0, epsilt, onem, &
                           g2kg, kg2g
  use mod_time,      only: nday_in_year, nday_of_year, nstep, &
                           nstep_in_day, baclin, &
                           xmi, l1mi, l2mi, l3mi, l4mi, l5mi
  use mod_xc
  use mod_vcoord,    only: vcoord_tag, vcoord_isopyc_bulkml
  use mod_grid,      only: scp2, plat, area
  use mod_state,     only: dp, temp, saln, p
  use mod_swtfrz,    only: swtfrz
  use mod_forcing,   only: sref, tflxap, sflxap, tflxdi, sflxdi, &
                           nflxdi, aptflx, apsflx, ditflx, disflx, &
                           sstclm, ricclm, sssclm, trxday, srxday, &
                           trxdpt, srxdpt, trxlim, srxlim, srxbal, &
                           swa, nsf, hmltfz, lip, sop, eva, rnf, rfi, &
                           fmltfz, sfl, ustarw, surflx, surrlx, &
                           sswflx, salflx, brnflx, salrlx, ustar, &
                           salt_corr, trc_corr, t_rs_nonloc, s_rs_nonloc
  use mod_swabs,     only: swal2, swfc2
  use mod_ben02,     only: tsi_tda, tml_tda, sml_tda, alb_tda, fice_tda, &
                           tsi, ntda, dfl, albw, alb, &
                           rnfins, rnfres, nrfets, rhowat
  use mod_thdysi,    only: tsrfm, ticem, albi_f, albi_m, albs_f, &
                           albs_m, rhoice, rhosnw, rkice, rksnw, fusi, &
                           fuss, fice_max, tice_m, tsnw_m, hice_nhmn, &
                           hice_shmn, sagets, sice, cwi, cuc
  use mod_seaice,    only: ficem, hicem, hsnwm, ustari, iagem
  use mod_utility,   only: util1, util2, util3
  use mod_checksum,  only: csdiag, chksum
  use mod_intp1d,    only: intp1d
  use mod_tracers,   only: ntr, itrtke, itrgls, trc, trflx
  use mod_diffusion, only: difdia
  use mod_tke,       only: gls_cmu0, Zos, gls_p, gls_m, gls_n, vonKar
  use mod_ifdefs,    only: use_TRC, use_TKE, use_GLS

  implicit none
  private

  public thermf_ben02

contains


  subroutine thermf_ben02(m,n,mm,nn,k1m,k1n)

    ! NERSC version of thermf.

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vrtsfl
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ntr) :: trflxc_aw
    integer :: i,j,k,l,m1,m2,m3,m4,m5,ntld,kn,kl,nt
    real :: dt,cpsw,rnf_fac,sag_fac,y
    real :: dpotl,hotl,totl,sotl,tice_f,hice_min,fice,hice,hsnw,tsrf
    real :: fice0,hice0,hsnw0,qsww,qnsw,tice,albi,tsmlt,albi_h,qswi,dh
    real :: qsnwf,fcond,qdamp,qsmlt,qo2i,qbot,swfac,dtml,q,volice,df,dvi
    real :: dvs,fwflx,sstc,rice,dpmxl,hmxl,tmxl,trxflx,pbot,dprsi,sssc
    real :: smxl,srxflx,sflxc,totsrp,totsrn,qp,qn,trflxc

    ! Due to conservation, the ratio of ice and snow density must be equal to
    ! the ratio of ice and snow heat of fusion
    if (abs(fuss/fusi-rhosnw/rhoice) > epsilt) then
      if (mnproc == 1) then
        write (lp,*) &
             'thermf: check consistency between snow/ice densities'
        write (lp,*) &
             'and heat of fusion!'
        stop
      end if
    end if

    ! Set various constants
    dt = baclin                         ! Time step
    cpsw = spcifh                       ! Specific heat of seawater
    rnf_fac = baclin/real(nrfets*86400) ! Runoff reservoar detrainment rate
    sag_fac = exp(-sagets*dt)           ! Snow aging rate

    ! Set parameters for time interpolation when applying diagnosed heat and
    ! salt relaxation fluxes
    y = (nday_of_year-1+mod(nstep,nstep_in_day)/real(nstep_in_day))*48. &
        /real(nday_in_year)
    m3 = int(y)+1
    y = y-real(m3-1)
    m1 = mod(m3+45,48)+1
    m2 = mod(m3+46,48)+1
    m4 = mod(m3   ,48)+1
    m5 = mod(m3+ 1,48)+1

    ! Time level for diagnosing heat and salt relaxation fluxes
    ntld = m3

    if (ditflx.or.disflx) nflxdi(ntld) = nflxdi(ntld)+1

    !$omp parallel do private( &
    !$omp l,i,dpotl,hotl,totl,sotl,tice_f,hice_min,fice,hice,hsnw,tsrf, &
    !$omp fice0,hice0,hsnw0,qsww,qnsw,tice,albi,tsmlt,albi_h,qswi,dh,qsnwf, &
    !$omp fcond,qdamp,qsmlt,qo2i,qbot,swfac,dtml,q,volice,df,dvi,dvs,fwflx, &
    !$omp sstc,rice,dpmxl,hmxl,tmxl,trxflx,pbot,dprsi,kn,kl,sssc,smxl, &
    !$omp srxflx, nt) &
    !$omp shared(xmi,y)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! Initialize variables describing the state of the ocean top layer,
          ! the mixed layer and ice/snow fraction
          dpotl = dp(i,j,k1n)
          hotl = dpotl/onem
          totl = temp(i,j,k1n)+t0deg
          sotl = saln(i,j,k1n)

          fice = ficem(i,j)
          hice = hicem(i,j)
          hsnw = hsnwm(i,j)
          tsrf = tsrfm(i,j)

          fice0 = fice
          hice0 = hice
          hsnw0 = hsnw

          ! Freezing point of sea water (in k)
          tice_f = swtfrz(p(i,j,1),sotl)+t0deg

          ! Minmimum ice thickness
          if (plat(i,j) > 0.) then
            hice_min = hice_nhmn
          else
            hice_min = hice_shmn
          end if

          if     (fice*hice < 1.e-5) then

            ! ------------------------------------------------------------------
            ! At most, a small amount of ice is present in the grid cell. Melt
            ! the remainders of ice and snow.
            ! ------------------------------------------------------------------

            hice = 0.
            hsnw = 0.
            fice = 0.

            ! Mean albedo of grid cell
            alb(i,j) = albw(i,j)

            ! Solar heat flux that enters the ocean
            qsww = swa(i,j)

            ! Non solar heat flux that enters the ocean
            qnsw = nsf(i,j)

            ! Set surface temperature and ice surface temperature
            tsrf = totl
            tice = totl

          else

            ! ------------------------------------------------------------------
            ! Do thermodynamics for an ice slab
            ! ------------------------------------------------------------------

            if (fice*hsnw > 1.e-3) then

              ! Set various variables in the case of a snow layer

              ! Albedo
              if (tsrf > tsnw_m-.1) then
                albi = albs_m
              else
                albi = albs_f
              end if

              ! Surface melting temperature
              tsmlt = tsnw_m

            else

              ! Set various variables in the case a thin or non existent snow
              ! layer

              ! Albedo
              albi_h = .065+.44*hice**.28
              if (tsrf > tice_m-.1) then
                albi = min(albi_m,albi_h)
              else
                albi = min(albi_f,albi_h)
              end if

              ! Surface melting temperature
              tsmlt = tice_m

            end if

            ! Mean albedo of the grid cell
            alb(i,j) = albi*fice+albw(i,j)*(1.-fice)

            ! Short wave radiation trough the ice covered fraction
            qswi = swa(i,j)*(1.-albi)/(1.-alb(i,j))

            ! Solar heat flux trough the open water fraction
            qsww = swa(i,j)*(1.-albw(i,j))/(1.-alb(i,j))

            ! Update snow thickness due to precipitation
            dh = sop(i,j)*dt/rhosnw
            hsnw = hsnw+dh

            ! Heat flux from snow to ice to balance the latent heat of snow fall
            qsnwf = dh*fuss/dt

            ! Conductive factor in snow and ice layer
            fcond = rkice*rksnw/(rksnw*hice+rkice*hsnw)

            ! Find the snow surface temperature and the non solar heat flux that
            ! enters the open water fraction
            if (abs(fcond-dfl(i,j)*(2.-fice)) < 1.e-3) then
              tsrf = tice_f+(qswi+nsf(i,j))/fcond
              qnsw = nsf(i,j)
              qdamp = 0.
            else
              tsrf = (qswi+nsf(i,j)-dfl(i,j)*(tsi(i,j)+(1.-fice)*totl) &
                   +fcond*tice_f) &
                   /(fcond-dfl(i,j)*(2.-fice))
              qnsw = nsf(i,j)+dfl(i,j)*fice*(totl-min(tsrf,tsmlt))
              qdamp = dfl(i,j)*(min(tsrf,tsmlt)-tsi(i,j))
            end if

            ! If the new surface temperature is above the snow melting
            ! temperature, determine the heat that goes to melting at the
            ! surface
            if (tsrf > tsmlt) then
              tsrf = tsmlt
              qsmlt = qswi+nsf(i,j) &
                   +dfl(i,j)*((1.-fice)*(tsrf-totl)+tsrf-tsi(i,j)) &
                   +fcond*(tice_f-tsrf)
            else
              qsmlt = 0.
            end if

            ! Set ice surface temperature
            tice = tice_f-fcond*(tice_f-tsrf)*hice/rkice

            ! Heat flux from ocean to ice (Maykut and McPhee 1995)
            qo2i = rhowat*cpsw*cwi*max(ustari(i,j),.2e-2)*min(tice_f-totl,0.) &
                 + cuc*max(tice_f-totl,0.)

            ! Heat budget at bottom of ice
            qbot = -fcond*(tice_f-tsrf)-qo2i-qdamp+qsnwf

            ! Update snow thickness due to melting
            dh = -qsmlt*dt/fuss
            if (hsnw+dh < 0.) then
              qsmlt = qsmlt-hsnw*fuss/dt
              hsnw = 0.
            else
              qsmlt = 0.
              hsnw = hsnw+dh
            end if

            ! Update ice thickness due to melting/freezing
            hice = max(0.,hice-(qbot+qsmlt)*dt/fusi)

            ! Convert snow to ice due to aging
            hice = hice+hsnw*(1.-sag_fac)*rhosnw/rhoice
            hsnw = hsnw*sag_fac

            ! Convert snow to ice if snow load is larger than the updrift of ice
            dh = (hsnw*rhosnw-hice*(rhowat-rhoice))/rhowat
            if (dh > 0.) then
              hice = hice+dh
              hsnw = hsnw-dh*rhoice/rhosnw
            end if

          end if

          ! --------------------------------------------------------------------
          ! Do thermodynamics for open water fraction of the grid cell
          ! --------------------------------------------------------------------

          ! Predict temperature change in mixed layer after a leapfrog time step
          ! due to heat fluxes
          swfac = 1.-swfc2(i,j)*exp(-hotl/swal2(i,j))
          dtml = (swfac*qsww+qnsw)*2.*dt/(cpsw*rhowat*hotl)

          if     (totl+dtml < tice_f) then

            ! Heat flux required to change the mixed layer temperature to the
            ! freezing point after a leapfrog time step
            q = .5*(tice_f-totl)*cpsw*rhowat*hotl/dt

            ! Ice volume that has to freeze to balance the heat budget
            volice = -(qsww+qnsw-q)*(1.-fice)*dt/fusi

            if (volice > epsilt) then

              ! New ice in the lead is formed with a specified thickness.
              ! Estimate the change in ice fraction
              df = volice/hice_min

              ! Redistribute ice and snow over an updated ice fraction
              hice = (hice*fice+volice)/min(fice_max,fice+df)
              hsnw = hsnw*fice/min(fice_max,fice+df)
              fice = min(fice_max,fice+df)

            end if

          else if (swfac*qsww+qnsw > 0.) then

            ! If the lead is warming, let the fraction  (1 - fice)  go to warm
            ! the lead, and the fraction  fice  to melt ice laterally
            fice = fice-(swfac*qsww+qnsw)*fice*dt &
                        /max(hice*fusi+hsnw*fuss,epsilt)
            if (fice < 0.) then
              fice = 0.
              hice = 0.
              hsnw = 0.
            end if

          end if

          ! --------------------------------------------------------------------
          ! Store the updated ice/snow state variables
          ! --------------------------------------------------------------------

          ficem(i,j) = fice
          hicem(i,j) = hice
          hsnwm(i,j) = hsnw
          tsrfm(i,j) = tsrf
          ticem(i,j) = tice

          ! --------------------------------------------------------------------
          ! Accumutate variables to produce averages in flux calculations
          ! --------------------------------------------------------------------

          alb_tda(i,j) = alb_tda(i,j)+alb(i,j)
          tml_tda(i,j) = tml_tda(i,j)+totl
          sml_tda(i,j) = sml_tda(i,j)+sotl
          fice_tda(i,j) = fice_tda(i,j)+fice
          tsi_tda(i,j) = tsi_tda(i,j)+tsrf

          ! --------------------------------------------------------------------
          ! Compute fluxes of heat and salt to the ocean
          ! --------------------------------------------------------------------

          ! Ice volume change
          dvi = hice*fice-hice0*fice0

          ! Snow volume change
          dvs = hsnw*fice-hsnw0*fice0

          ! Accumulate the runoff in a reservoar to delay the discharge into the
          ! ocean (by nrfets days approximately 1/e of runoff added will by
          ! discharged).
          rnfres(i,j) = rnfres(i,j)+rnfins(i,j)
          rnf(i,j) = rnfres(i,j)*rnf_fac
          rnfres(i,j) = rnfres(i,j)*(1.-rnf_fac)

          ! Fresh water flux due to melting/freezing [kg m-2 s-1] (positive
          ! downwards)
          fmltfz(i,j) = -(dvi*rhoice+dvs*rhosnw)/dt

          ! Fresh water flux [kg m-2 s-1] (positive downwards)
          fwflx = eva(i,j)+lip(i,j)+sop(i,j)+rnf(i,j)+rfi(i,j)+fmltfz(i,j)

          ! Salt flux [kg m-2 s-1] (positive downwards)
          sfl(i,j) = -sice*dvi*rhoice/dt*g2kg

          ! Salt flux due to brine rejection of freezing sea ice [kg m-2 s-1]
          ! (positive downwards)
          brnflx(i,j) = max(0.,-sotl*fmltfz(i,j)*g2kg+sfl(i,j))

          ! Virtual salt flux [kg m-2 s-1] (positive downwards)
          vrtsfl(i,j) = -sotl*fwflx*g2kg

          ! Store area weighted correction to the virtual salt flux. When the
          ! global average correction is applied, the virtual salt flux is
          ! globally consistent with a salt flux based on some reference
          ! salinity and any salinity limiting compensated for.
          util1(i,j) = -(sref*fwflx*g2kg &
                        +vrtsfl(i,j) &
                        +salt_corr(i,j)*g2kg/(2.*dt)) &
                        *scp2(i,j)

          ! Reset salt correction
          salt_corr(i,j) = 0.

          ! Heat flux due to melting/freezing [W m-2] (positive downwards)
          hmltfz(i,j) = (dvi*fusi+dvs*fuss)/dt

          ! Total heat flux in BLOM units [W m-2] (positive upwards)
          surflx(i,j) = -(swa(i,j)+nsf(i,j)+hmltfz(i,j))

          ! Short-wave heat flux in BLOM units [W m-2] (positive
          ! upwards)
          sswflx(i,j) = -qsww*(1.-fice0)

          if (use_TRC) then
            ! ------------------------------------------------------------------
            ! Tracer fluxes (positive downwards)
            ! ------------------------------------------------------------------

            do nt = 1,ntr
              if (use_TKE) then
                if (nt == itrtke) then
                  trflx(nt,i,j) = 0.
                  trflxc_aw(i,j,nt) = 0.
                  cycle
                end if
                if (use_GLS) then
                  if (nt == itrgls) then
                    trflx(nt,i,j) = -gls_n*difdia(i,j,1)*(gls_cmu0**gls_p) &
                         *(trc(i,j,k1n,itrtke)**gls_m) &
                         *(vonKar**gls_n)*Zos**(gls_n-1.)
                    trflxc_aw(i,j,nt) = 0.
                    cycle
                  end if
                else
                  if (nt == itrgls) then
                    trflx(nt,i,j) = 0.
                    trflxc_aw(i,j,nt) = 0.
                    cycle
                  end if
                end if
              end if
              trflx(nt,i,j) = -trc(i,j,k1n,nt)*fwflx
              trflxc_aw(i,j,nt) = -(trflx(nt,i,j) &
                                   +trc_corr(i,j,nt)/(2.*dt)) &
                                   *scp2(i,j)
              trc_corr(i,j,nt) = 0.
            end do
          end if

          ! --------------------------------------------------------------------
          ! Relaxation fluxes
          ! --------------------------------------------------------------------

          surrlx(i,j) = 0.

          ! If  trxday>0 , apply relaxation towards observed sst
          if (trxday > epsilt) then
            sstc = intp1d(sstclm(i,j,l1mi),sstclm(i,j,l2mi), &
                          sstclm(i,j,l3mi),sstclm(i,j,l4mi), &
                          sstclm(i,j,l5mi),xmi)
            rice = intp1d(ricclm(i,j,l1mi),ricclm(i,j,l2mi), &
                          ricclm(i,j,l3mi),ricclm(i,j,l4mi), &
                          ricclm(i,j,l5mi),xmi)
            sstc = (1.-rice)*max(sstc,tice_f)+rice*tice_f
            if (vcoord_tag == vcoord_isopyc_bulkml) then
              dpmxl = dp(i,j,1+nn)+dp(i,j,2+nn)
              hmxl = dpmxl/onem
              tmxl = (temp(i,j,1+nn)*dp(i,j,1+nn) &
                     +temp(i,j,2+nn)*dp(i,j,2+nn))/dpmxl+t0deg
              trxflx = spcifh*min(hmxl,trxdpt)/(trxday*86400.) &
                       *min(trxlim,max(-trxlim,sstc-tmxl))/alpha0
            else
              pbot = p(i,j,1)
              do k = 1,kk
                kn = k+nn
                pbot = pbot+dp(i,j,kn)
              end do
              dprsi = 1./min(trxdpt*onem,pbot-p(i,j,1))
              t_rs_nonloc(i,j,1) = 1.
              tmxl = 0.
              do k = 1,kk
                kn = k+nn
                t_rs_nonloc(i,j,k+1) = t_rs_nonloc(i,j,k)-dp(i,j,kn)*dprsi
                if (t_rs_nonloc(i,j,k+1) < 0.) then
                  tmxl = tmxl+temp(i,j,kn)*t_rs_nonloc(i,j,k)+t0deg
                  exit
                else
                  tmxl = tmxl+temp(i,j,kn)*(t_rs_nonloc(i,j,k  ) &
                                           -t_rs_nonloc(i,j,k+1))
                end if
              end do
              do kl = k,kk
                t_rs_nonloc(i,j,kl+1) = 0.
              end do
              trxflx = spcifh*trxdpt/(trxday*86400.) &
                       *min(trxlim,max(-trxlim,sstc-tmxl))/alpha0
            end if
            surrlx(i,j) = -trxflx
          else
            trxflx = 0.
          end if

          ! If aptflx=.true., apply diagnosed relaxation flux
          if (aptflx) then
            surrlx(i,j) = surrlx(i,j) &
                        - intp1d(tflxap(i,j,m1),tflxap(i,j,m2),tflxap(i,j,m3), &
                                 tflxap(i,j,m4),tflxap(i,j,m5),y)
          end if

          ! If ditflx=.true., diagnose relaxation flux by accumulating the
          ! relaxation flux
          if (ditflx) then
            tflxdi(i,j,ntld) = tflxdi(i,j,ntld)+trxflx
          end if

          salrlx(i,j) = 0.

          ! If  srxday>0 , apply relaxation towards observed sss
          if (srxday > epsilt) then
            sssc = intp1d(sssclm(i,j,l1mi),sssclm(i,j,l2mi), &
                          sssclm(i,j,l3mi),sssclm(i,j,l4mi), &
                          sssclm(i,j,l5mi),xmi)
            if (vcoord_tag == vcoord_isopyc_bulkml) then
              dpmxl = dp(i,j,1+nn)+dp(i,j,2+nn)
              hmxl = dpmxl/onem
              smxl = (saln(i,j,1+nn)*dp(i,j,1+nn) &
                     +saln(i,j,2+nn)*dp(i,j,2+nn))/dpmxl
              srxflx = min(hmxl,srxdpt)/(srxday*86400.) &
                       *min(srxlim,max(-srxlim,sssc-smxl))/alpha0
            else
              pbot = p(i,j,1)
              do k = 1,kk
                kn = k+nn
                pbot = pbot+dp(i,j,kn)
              end do
              dprsi = 1./min(srxdpt*onem,pbot-p(i,j,1))
              s_rs_nonloc(i,j,1) = 1.
              smxl = 0.
              do k = 1,kk
                kn = k+nn
                s_rs_nonloc(i,j,k+1) = s_rs_nonloc(i,j,k)-dp(i,j,kn)*dprsi
                if (s_rs_nonloc(i,j,k+1) < 0.) then
                  smxl = smxl+saln(i,j,kn)*s_rs_nonloc(i,j,k)
                  exit
                else
                  smxl = smxl+saln(i,j,kn)*(s_rs_nonloc(i,j,k  ) &
                                           -s_rs_nonloc(i,j,k+1))
                end if
              end do
              do kl = k,kk
                s_rs_nonloc(i,j,kl+1) = 0.
              end do
              srxflx = srxdpt/(srxday*86400.) &
                       *min(srxlim,max(-srxlim,sssc-smxl))/alpha0
            end if
            salrlx(i,j) = -srxflx
            util2(i,j) = max(0.,salrlx(i,j))*scp2(i,j)
            util3(i,j) = min(0.,salrlx(i,j))*scp2(i,j)
          else
            srxflx = 0.
          end if

          ! If apsflx=.true., apply diagnosed relaxation flux
          if (apsflx) then
            salrlx(i,j) = salrlx(i,j) &
                        - intp1d(sflxap(i,j,m1),sflxap(i,j,m2),sflxap(i,j,m3), &
                                 sflxap(i,j,m4),sflxap(i,j,m5),y)
          end if

          ! If disflx=.true., diagnose relaxation flux by accumulating the
          ! relaxation flux
          if (disflx) then
            sflxdi(i,j,ntld) = sflxdi(i,j,ntld)+srxflx
          end if

          ! --------------------------------------------------------------------
          ! Update age of ice
          ! --------------------------------------------------------------------

          if (fice*hice < 1.e-5) then
            iagem(i,j) = 0.
          else
            iagem(i,j) = (iagem(i,j)+dt/86400.)*(1.-max(0.,dvi)/(fice*hice))
          end if

          ! --------------------------------------------------------------------
          ! Friction velocity (m/s)
          ! --------------------------------------------------------------------

          ustar(i,j) = (min(ustari(i,j),.8e-2)*fice0+ustarw(i,j)*(1.-fice0))

        end do
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------------
    ! Apply the virtual salt flux correction and the compute the total
    ! salt flux by combining the virtual and true salt flux. Also
    ! convert salt fluxes used later to unit [g m-2 s-1] and positive
    ! upwards.
    ! ------------------------------------------------------------------

    call xcsum(sflxc,util1,ips)
    sflxc = sflxc/area
    if (mnproc == 1) then
      write (lp,'(a,e15.7)') ' thermf: sflxc', sflxc
    end if

    !$omp parallel do private(l,i)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          salflx(i,j) = -(vrtsfl(i,j)+sflxc+sfl(i,j))*kg2g
          brnflx(i,j) = -brnflx(i,j)*kg2g
        end do
      end do
    end do
    !$omp end parallel do

    ! If  srxday>0  and  srxbal=.true. , balance the sss relaxation flux
    ! so the net input of salt in grid cells connected to the world
    ! ocean is zero
    if (srxday > epsilt.and.srxbal) then
      call xcsum(totsrp,util2,ipwocn)
      call xcsum(totsrn,util3,ipwocn)
      if (abs(totsrp-totsrn) > 0.) then
        qp = 2.*totsrn/(totsrn-totsrp)
        qn = 2.*totsrp/(totsrp-totsrn)
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              if (ipwocn(i,j) == 1) then
                salrlx(i,j) = qp*max(0.,salrlx(i,j)) + qn*min(0.,salrlx(i,j))
              end if
            end do
          end do
        end do
        !$omp end parallel do
      end if
    end if

    if (use_TRC) then
      do nt = 1,ntr
        if (use_TKE) then
          if (nt == itrtke.or.nt == itrgls) cycle
        end if
        call xcsum(trflxc,trflxc_aw(1-nbdy,1-nbdy,nt),ips)
        trflxc = trflxc/area
        !$omp parallel do private(l,i)
        do j = 1,jj
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              trflx(nt,i,j) = -(trflx(nt,i,j)+trflxc)
            end do
          end do
        end do
        !$omp end parallel do
      end do
    end if

    ! --------------------------------------------------------------------------
    ! Number of accumulated fields for flux calculations
    ! --------------------------------------------------------------------------

    ntda = ntda+1

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'thermf_ben02:'
      end if
      call chksum(alb     , 1, halo_ps, 'alb'     )
      call chksum(ficem   , 1, halo_ps, 'ficem'   )
      call chksum(hicem   , 1, halo_ps, 'hicem'   )
      call chksum(hsnwm   , 1, halo_ps, 'hsnwm'   )
      call chksum(tsrfm   , 1, halo_ps, 'tsrfm'   )
      call chksum(ticem   , 1, halo_ps, 'ticem'   )
      call chksum(alb_tda , 1, halo_ps, 'alb_tda' )
      call chksum(tml_tda , 1, halo_ps, 'tml_tda' )
      call chksum(sml_tda , 1, halo_ps, 'sml_tda' )
      call chksum(fice_tda, 1, halo_ps, 'fice_tda')
      call chksum(tsi_tda , 1, halo_ps, 'tsi_tda' )
      call chksum(rnfres  , 1, halo_ps, 'rnfres'  )
      call chksum(surflx  , 1, halo_ps, 'surflx'  )
      call chksum(sswflx  , 1, halo_ps, 'sswflx'  )
      call chksum(salflx  , 1, halo_ps, 'salflx'  )
      call chksum(brnflx  , 1, halo_ps, 'brnflx'  )
      call chksum(surrlx  , 1, halo_ps, 'surrlx'  )
      call chksum(salrlx  , 1, halo_ps, 'salrlx'  )
      call chksum(iagem   , 1, halo_ps, 'iagem'   )
      call chksum(ustar   , 1, halo_ps, 'ustar'   )
      if (vcoord_tag /= vcoord_isopyc_bulkml) then
        call chksum(t_rs_nonloc, kk+1, halo_ps, 't_rs_nonloc')
        call chksum(s_rs_nonloc, kk+1, halo_ps, 's_rs_nonloc')
      end if
    end if

  end subroutine thermf_ben02

end module mod_thermf_ben02
