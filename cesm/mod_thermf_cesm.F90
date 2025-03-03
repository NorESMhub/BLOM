! ------------------------------------------------------------------------------
! Copyright (C) 2008-2025 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_thermf_cesm

  use mod_constants, only: grav, spcifh, t0deg, alpha0, epsilt, onem, &
                           g2kg, kg2g
  use mod_time,      only: nstep, nstep_in_day, nday_in_year, &
                           nday_of_year, baclin, &
                           xmi, l1mi, l2mi, l3mi, l4mi, l5mi
  use mod_xc
  use mod_vcoord,    only: vcoord_tag, vcoord_isopyc_bulkml
  use mod_grid,      only: scp2, area
  use mod_state,     only: dp, temp, saln, p
  use mod_swtfrz,    only: swtfrz
  use mod_forcing,   only: sref, tflxap, sflxap, tflxdi, sflxdi, &
                           nflxdi, aptflx, apsflx, ditflx, disflx, &
                           sstclm, ricclm, sssclm, trxday, srxday, &
                           trxdpt, srxdpt, trxlim, srxlim, srxbal, &
                           swa, nsf, hmltfz, lip, sop, eva, rnf, rfi, &
                           fmltfz, sfl, ustarw, surflx, surrlx, &
                           sswflx, salflx, brnflx, salrlx, ustar, &
                           salt_corr, trc_corr, t_rs_nonloc, s_rs_nonloc, &
                           sss_stream, sst_stream, ice_stream, &
                           use_stream_relaxation
  use mod_cesm,      only: hmlt, frzpot, mltpot
  use mod_utility,   only: util1, util2, util3
  use mod_checksum,  only: csdiag, chksummsk
  use mod_tracers,   only: ntr, itrtke, itrgls, trc, trflx
  use mod_diffusion, only: difdia
  use mod_tke,       only: gls_cmu0, Zos, gls_p, gls_m, gls_n, vonKar
  use mod_intp1d,    only: intp1d
  use mod_ifdefs,    only: use_TRC, use_TKE, use_GLS

  implicit none
  private

  public :: thermf_cesm

contains

  subroutine thermf_cesm(m,n,mm,nn,k1m,k1n)

    ! NERSC version of thermf. To be used when coupled to CESM

    ! Arguments
    integer, intent(in) :: m,n,mm,nn,k1m,k1n

    ! Local variables
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tfrz, tfrzm, vrtsfl
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ntr) :: trflxc_aw
    integer :: i,j,k,l,m1,m2,m3,m4,m5,ntld,kn,kl,nt
    real :: y,dpotl,hotl,totl,sotl,tice_f,fwflx,sstc,rice,dpmxl,hmxl
    real :: tmxl,trxflx,pbot,dprsi,sssc,smxl,srxflx,sflxc
    real :: totsrp,totsrn,qp,qn,trflxc

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

    ! Compute freezing temperatures of sea water
    tfrz(:,:) = swtfrz(p(:,:,1),saln(:,:,k1n))
    tfrzm(:,:) = swtfrz(p(:,:,1),.5*(saln(:,:,k1m)+saln(:,:,k1n)))

    if (ditflx.or.disflx) nflxdi(ntld) = nflxdi(ntld)+1

    !$omp parallel do private( &
    !$omp l,i,dpotl,hotl,totl,sotl,tice_f,fwflx,sstc,rice,dpmxl,hmxl,tmxl, &
    !$omp trxflx,pbot,dprsi,kn,kl,sssc,smxl,srxflx,nt)
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))

          ! --------------------------------------------------------------------
          ! Set some quantities
          ! --------------------------------------------------------------------

          ! ocean top layer quantities
          dpotl = dp(i,j,k1n)
          hotl = dpotl/onem
          totl = temp(i,j,k1n)+t0deg
          sotl = saln(i,j,k1n)

          tice_f = tfrz(i,j)+t0deg

          ! --------------------------------------------------------------------
          ! Fresh water and salt fluxes
          ! --------------------------------------------------------------------

          ! Fresh water flux [kg m-2 s-1] (positive downwards)
          fwflx = eva(i,j)+lip(i,j)+sop(i,j)+rnf(i,j)+rfi(i,j)+fmltfz(i,j)

          ! Salt flux due to brine rejection of freezing sea ice [kg m-2 m-1]
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
                        +salt_corr(i,j)*g2kg/(2.*baclin)) &
                        *scp2(i,j)

          ! Reset salt correction
          salt_corr(i,j) = 0.

          ! --------------------------------------------------------------------
          ! Heat fluxes
          ! --------------------------------------------------------------------

          ! Freezing/melting potential [J m-2]. A positive flux means the ocean
          ! surface has a temperature below freezing temperature and must be
          ! heated. Note the freezing potential is multiplied by 1/2 due to the
          ! leap-frog time stepping. The melting potential uses time averaged
          ! quantities since it is not accumulated.
          frzpot(i,j) = max(0.,tice_f-totl)*spcifh*dpotl/(2.*grav)
          mltpot(i,j) = min(0.,tfrzm(i,j)-.5*(temp(i,j,k1m)+temp(i,j,k1n))) &
                        *spcifh*.5*(dp(i,j,k1m)+dp(i,j,k1n))/grav

          ! Heat flux due to melting/freezing [W m-2] (positive downwards)
          hmltfz(i,j) = hmlt(i,j)+frzpot(i,j)/baclin

          ! Total heat flux in BLOM units [W cm-2] (positive upwards)
          surflx(i,j) = -(swa(i,j)+nsf(i,j)+hmltfz(i,j))

          ! Short-wave heat flux in BLOM units [W cm-2] (positive
          ! upwards)
          sswflx(i,j) = -swa(i,j)

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
                                   +trc_corr(i,j,nt)/(2.*baclin)) &
                                   *scp2(i,j)
              trc_corr(i,j,nt) = 0.
            end do
          end if

          ! --------------------------------------------------------------------
          ! Relaxation fluxes
          ! --------------------------------------------------------------------

          surrlx(i,j) = 0.

          ! If  trxday>0 , apply relaxation towards observed sst
          if (trxday > epsilt ) then
            if (use_stream_relaxation) then
               sstc = sst_stream(i,j)
               rice = ice_stream(i,j)
               sstc = (1.-rice)*max(sstc,tice_f) + rice*tice_f
            else
               sstc = intp1d(sstclm(i,j,l1mi),sstclm(i,j,l2mi), &
                             sstclm(i,j,l3mi),sstclm(i,j,l4mi), &
                             sstclm(i,j,l5mi),xmi)
               rice = intp1d(ricclm(i,j,l1mi),ricclm(i,j,l2mi), &
                             ricclm(i,j,l3mi),ricclm(i,j,l4mi), &
                             ricclm(i,j,l5mi),xmi)
               sstc = (1.-rice)*max(sstc,tice_f)+rice*tice_f
            end if
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
            if (use_stream_relaxation) then
               sssc = sss_stream(i,j)
            else
               sssc = intp1d(sssclm(i,j,l1mi),sssclm(i,j,l2mi), &
                             sssclm(i,j,l3mi),sssclm(i,j,l4mi),sssclm(i,j,l5mi),xmi)
            end if
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
          ! Friction velocity (cm/s)
          ! --------------------------------------------------------------------

          ustar(i,j) = ustarw(i,j)

        end do
      end do
    end do
    !$omp end parallel do

    ! ------------------------------------------------------------------
    ! Apply the virtual salt flux correction and the compute the total
    ! salt flux by combining the virtual and true salt flux. Also
    ! convert salt fluxes used later to unit [10e-3 g cm-2 s-1] and
    ! positive upwards.
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
          salflx(i,j) = -(vrtsfl(i,j)+sflxc+sfl(i,j)) &
                         *(kg2g)
          brnflx(i,j) = -brnflx(i,j)*(kg2g)
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

    if (csdiag) then
      if (mnproc == 1) then
        write (lp,*) 'thermf_cesm:'
      end if
      call chksummsk(surflx,ip,1,'surflx')
      call chksummsk(sswflx,ip,1,'sswflx')
      call chksummsk(salflx,ip,1,'salflx')
      call chksummsk(brnflx,ip,1,'brnflx')
      call chksummsk(surrlx,ip,1,'surrlx')
      call chksummsk(salrlx,ip,1,'salrlx')
      call chksummsk(hmltfz,ip,1,'hmltfz')
      call chksummsk(ustar,ip,1,'ustar')
      call chksummsk(frzpot,ip,1,'frzpot')
      call chksummsk(mltpot,ip,1,'mltpot')
      if (vcoord_tag /= vcoord_isopyc_bulkml) then
        call chksummsk(t_rs_nonloc, ip, kk+1, 't_rs_nonloc')
        call chksummsk(s_rs_nonloc, ip, kk+1, 's_rs_nonloc')
      end if
    end if

  end subroutine thermf_cesm

end module mod_thermf_cesm
