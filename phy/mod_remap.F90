! ------------------------------------------------------------------------------
! Copyright (C) 2006-2024 Mats Bentsen, Mehmet Ilicak, Mariana Vertenstein
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

module mod_remap

  ! ------------------------------------------------------------------
  ! This module contains variables and procedures related to advection
  ! of layer pressure thickness and tracers by incremental remapping.
  ! ------------------------------------------------------------------

  ! NOTE: natr was not included below - so the ifdef for ATRC was never tested
  ! since it would not compile

  use mod_types,     only: r8
  use mod_xc
  use mod_tracers,   only: ntr, itrtke, itrgls, natr, trc ! TRC
  use mod_ifdefs,    only: use_TRC, use_ATRC, use_TKE, use_TKEADV

  implicit none
  private

  ! Parameters:
  real(r8), parameter :: &
       dpeps = 1.e-12_r8           ! Small layer pressure thickness (equivalent
                                   ! to approximately 10-16 m) [g cm-1 s-2].
  real(r8), parameter :: &
       treps = 1.e-14_r8           ! Small tracer concentration. (for use_TRC and use_ATRC)

  public :: remap_eitvel, remap_eitflx

  private :: triint, penint

contains

  ! ------------------------------------------------------------------
  ! Private procedures.
  ! ------------------------------------------------------------------

  subroutine triint(ac,x1,y1,x2,y2,x3,y3,&
                    a,ax,ay,axx,ayy,axy, &
                    axxx,ayyy,axxy,axyy)

    ! Arguments
    real, intent(in)  :: ac,x1,y1,x2,y2,x3,y3
    real, intent(out) :: a,ax,ay,axx,ayy,axy
    real, intent(out) :: axxx,ayyy,axxy,axyy

    ! Local variables
    real :: r1_3,r1_6,r1_12
    parameter (r1_3=1./3.,r1_6=1./6.,r1_12 = 1./12.)
    real :: r1_10,r1_30
    parameter (r1_10=.1,r1_30 = 1./30.)
    real :: xx,yy,xy1,xy2,xy3,xy

    xx = x1*x2+x2*x3+x1*x3
    yy = y1*y2+y2*y3+y1*y3
    xy1 = x1*y1
    xy2 = x2*y2
    xy3 = x3*y3
    xy = xy1+xy2+xy3

    a = .5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))*ac

    ax = r1_3*(x1+x2+x3)
    ay = r1_3*(y1+y2+y3)
    axx = r1_6*(9.*ax*ax-xx)
    ayy = r1_6*(9.*ay*ay-yy)
    axy = r1_12*(9.*ax*ay+xy)
    if (use_TRC .and. use_ATRC) then
      axxx = r1_10*((18.*axx-3.*xx)*ax+x1*x2*x3)
      ayyy = r1_10*((18.*ayy-3.*yy)*ay+y1*y2*y3)
      axxy = r1_30*(18.*axx*ay+3.*ax*xy+x1*xy1+x2*xy2+x3*xy3)
      axyy = r1_30*(18.*ayy*ax+3.*ay*xy+y1*xy1+y2*xy2+y3*xy3)
    end if

    ax = ax*a
    ay = ay*a
    axx = axx*a
    ayy = ayy*a
    axy = axy*a
    if (use_TRC .and. use_ATRC) then
      axxx = axxx*a
      ayyy = ayyy*a
      axxy = axxy*a
      axyy = axyy*a
    end if

  end subroutine triint

  subroutine penint(ac,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5, &
                    a,ax,ay,axx,ayy,axy, &
                    axxx,ayyy,axxy,axyy)

    ! Arguments
    real, intent(in)  :: ac,x1,y1,x2,y2,x3,y3,x4,y4,x5,y5
    real, intent(out) :: a,ax,ay,axx,ayy,axy
    real, intent(out) :: axxx,ayyy,axxy,axyy

    ! Local variables
    real, parameter :: r1_3 = 1./3.
    real, parameter :: r1_6 = 1./6.
    real, parameter :: r1_12 = 1./12.
    real, parameter :: r1_10 = .1
    real, parameter :: r1_30 = 1./30.
    real :: xx123,yy123,xx135,yy135,xx345,yy345,xy1,xy2,xy3,xy4,xy5
    real :: xy123,xy135,xy345,a123,a135,a345,ax123,ax135,ax345
    real :: ay123,ay135,ay345,axx123,axx135,axx345,ayy123,ayy135,ayy345
    real :: axy123,axy135,axy345
    real :: axxx123,axxx135,axxx345,ayyy123,ayyy135,ayyy345
    real :: axxy123,axxy135,axxy345,axyy123,axyy135,axyy345

    xx123 = x1*x2+x2*x3+x1*x3
    yy123 = y1*y2+y2*y3+y1*y3
    xx135 = x1*x3+x3*x5+x1*x5
    yy135 = y1*y3+y3*y5+y1*y5
    xx345 = x3*x4+x4*x5+x3*x5
    yy345 = y3*y4+y4*y5+y3*y5

    xy1 = x1*y1
    xy2 = x2*y2
    xy3 = x3*y3
    xy4 = x4*y4
    xy5 = x5*y5

    xy123 = xy1+xy2+xy3
    xy135 = xy1+xy3+xy5
    xy345 = xy3+xy4+xy5

    a123 = .5*((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))*ac
    a135 = .5*((x3-x1)*(y5-y1)-(y3-y1)*(x5-x1))*ac
    a345 = .5*((x4-x3)*(y5-y3)-(y4-y3)*(x5-x3))*ac

    ax123 = r1_3*(x1+x2+x3)
    ay123 = r1_3*(y1+y2+y3)
    ax135 = r1_3*(x1+x3+x5)
    ay135 = r1_3*(y1+y3+y5)
    ax345 = r1_3*(x3+x4+x5)
    ay345 = r1_3*(y3+y4+y5)

    axx123 = r1_6*(9.*ax123*ax123-xx123)
    ayy123 = r1_6*(9.*ay123*ay123-yy123)
    axy123 = r1_12*(9.*ax123*ay123+xy123)
    axx135 = r1_6*(9.*ax135*ax135-xx135)
    ayy135 = r1_6*(9.*ay135*ay135-yy135)
    axy135 = r1_12*(9.*ax135*ay135+xy135)
    axx345 = r1_6*(9.*ax345*ax345-xx345)
    ayy345 = r1_6*(9.*ay345*ay345-yy345)
    axy345 = r1_12*(9.*ax345*ay345+xy345)

    if (use_TRC .and. use_ATRC) then
      axxx123 = r1_10*((18.*axx123-3.*xx123)*ax123+x1*x2*x3)
      ayyy123 = r1_10*((18.*ayy123-3.*yy123)*ay123+y1*y2*y3)
      axxy123 = r1_30*(18.*axx123*ay123+3.*ax123*xy123 &
               +x1*xy1+x2*xy2+x3*xy3)
      axyy123 = r1_30*(18.*ayy123*ax123+3.*ay123*xy123 &
               +y1*xy1+y2*xy2+y3*xy3)
      axxx135 = r1_10*((18.*axx135-3.*xx135)*ax135+x1*x3*x5)
      ayyy135 = r1_10*((18.*ayy135-3.*yy135)*ay135+y1*y3*y5)
      axxy135 = r1_30*(18.*axx135*ay135+3.*ax135*xy135 &
               +x1*xy1+x3*xy3+x5*xy5)
      axyy135 = r1_30*(18.*ayy135*ax135+3.*ay135*xy135 &
               +y1*xy1+y3*xy3+y5*xy5)
      axxx345 = r1_10*((18.*axx345-3.*xx345)*ax345+x3*x4*x5)
      ayyy345 = r1_10*((18.*ayy345-3.*yy345)*ay345+y3*y4*y5)
      axxy345 = r1_30*(18.*axx345*ay345+3.*ax345*xy345 &
               +x3*xy3+x4*xy4+x5*xy5)
      axyy345 = r1_30*(18.*ayy345*ax345+3.*ay345*xy345 &
               +y3*xy3+y4*xy4+y5*xy5)
    end if

    a = a123+a135+a345
    ax = ax123*a123+ax135*a135+ax345*a345
    ay = ay123*a123+ay135*a135+ay345*a345
    axx = axx123*a123+axx135*a135+axx345*a345
    ayy = ayy123*a123+ayy135*a135+ayy345*a345
    axy = axy123*a123+axy135*a135+axy345*a345

    if (use_TRC .and. use_ATRC) then
      axxx = axxx123*a123+axxx135*a135+axxx345*a345
      ayyy = ayyy123*a123+ayyy135*a135+ayyy345*a345
      axxy = axxy123*a123+axxy135*a135+axxy345*a345
      axyy = axyy123*a123+axyy135*a135+axyy345*a345
    end if

  end subroutine penint

  !---------------------------------------------------------------
  ! Public procedures.
  !---------------------------------------------------------------

  subroutine remap_eitvel(scuy,scvx,scp2i,scp2,pbmin,pbu,pbv,plo, &
                          u,v,dt,mrg,dp,temp,saln,uflx,vflx, &
                          utflx,vtflx,usflx,vsflx,k)

    !---------------------------------------------------------------
    ! Advection of layer pressure thickness and tracers by incremental
    ! remapping.
    !---------------------------------------------------------------

    ! Argument variables:
    !   scuy   - length of cell boundary with u-point as midpoint.
    !   scvx   - length of cell boundary with v-point as midpoint.
    !   scp2i  - inverse of grid cell area.
    !   scp2   - grid cell area.
    !   pbmin  - minimum bottom pressure of a grid cell and its
    !            neighbors.
    !   pbu    - bottom pressure at u-point, defined as
    !            min(pb(i-1,j),pb(i,j)).
    !   pbv    - bottom pressure at v-point, defined as
    !            min(pb(i,j-1),pb(i,j)).
    !   plo    - lower interface pressure of layer pressure thickness.
    !   u      - u-component of velocity.
    !   v      - v-component of velocity.
    !   dt     - time step.
    !   temp   - temperature.
    !   saln   - salinity.
    !   uflx   - u-component of mass flux.
    !   vflx   - v-component of mass flux.
    !   utflx  - u-component of heat flux.
    !   vtflx  - v-component of heat flux.
    !   usflx  - u-component of salt flux.
    !   vsflx  - v-component of salt flux.
    !   mrg    - margin of halo that must be valid upon return.
    !   k      - layer index

    ! Arguments
    integer, intent(in) :: mrg
    integer, intent(in) :: k
    real,    intent(in) :: dt
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: scuy
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: scvx
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: scp2i
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: scp2
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbmin
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbu
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbv
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: plo
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: u
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: v
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dp
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: temp
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: saln
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: uflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: utflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vtflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: usflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vsflx

    ! Local variables.
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         pup,dx,dy,xd,yd,tx,ty,td,sx,sy,sd,cu,cv,cuc,cvc,fdu,fdv,ftu,ftv, &
         fsu,fsv
    real :: dxi,dyi,dpw,dpe,dps,dpn,dpsw,dpse,dpc,dpnw,dpne, &
         dgmx,dfmx,dfmn,q,q1,q2,q3,q4,tgmx,tgmn,tfmx,tfmn, &
         sgmx,sgmn,sfmx,sfmn,xm,ym,xc0,xc1,yc0,yc1,x2,y2,x4,y4, &
         a,ax,ay,axx,ayy,axy,dl,fd,qx,qy
    integer :: i,j,l,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne,nw
    real, dimension(ntr-natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         trx,try,trd,ftru,ftrv
    real, dimension(natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         ag,agx,agy,agd,fagu,fagv
    real :: xdt,ydt,axxx,ayyy,axxy,axyy,qxx,qyy,qxy,fdt
    integer :: nt,nat

    !---------------------------------------------------------------
    ! General information:
    !   Logical arrangment of variables is as follows: Layer pressure
    !   thickness dp(i,j), as an example of a scalar variable, is the
    !   mean layer pressure thickness of grid cell (i,j). Velocity
    !   component u(i,j) is located at the midpoint of the cell boundary
    !   separating grid cells (i-1,j) and (i,j). Velocity component
    !   v(i,j) is located at the midpoint of the cell boundary
    !   separating grid cells (i,j-1) and (i,j). A corner variable with
    !   index (i,j) is located at the common grid cell corner of grid
    !   cells (i-1,j-1), (i,j-1), (i-1,j), and (i,j).
    !
    !   The divergence of the velocity field is defined as follows:
    !     (u(i+1,j)*scuy(i+1,j)-u(i,j)*scuy(i,j)
    !     +v(i,j+1)*scvy(i,j+1)-v(i,j)*scvx(i,j))*scp2i(i,j)
    !   By construction, the "fluxing areas" used in obtaining fluxes
    !   trough cell boundaries containing u(i,j) and v(i,j), are equal
    !   to u(i,j)*scuy(i,j)*dt and v(i,j)*scvx(i,j)*dt, respectively.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! Add small number to density field and initialize some variables.
    !---------------------------------------------------------------

    do j = 1-mrg-2,jj+mrg+2
      do l = 1,isp(j)
        do i = max(1-mrg-2,ifp(j,l)),min(ii+mrg+2,ilp(j,l))
          dp(i,j) = max(0.,dp(i,j))+dpeps
          pup(i,j) = plo(i,j)-dp(i,j)
        end do
      end do
      do i = 1-mrg-1,ii+mrg+1
        fdu(i,j) = 0.
        fdv(i,j) = 0.
        ftu(i,j) = 0.
        ftv(i,j) = 0.
        fsu(i,j) = 0.
        fsv(i,j) = 0.
        if (use_TRC) then
          if (use_ATRC) then
            do nt = 1,ntr-natr
              if (use_TKE .and. .not. use_TKEADV) then
                if (nt == itrtke.or.nt == itrgls) cycle
              end if
              ftru(nt,i,j) = 0.
              ftrv(nt,i,j) = 0.
            end do
            do nt = 1,natr
              fagu(nt,i,j) = 0.
              fagv(nt,i,j) = 0.
            end do
          else
            do nt = 1,ntr
              if (use_TKE .and. .not. use_TKEADV) then
                if (nt == itrtke.or.nt == itrgls) cycle
              end if
              ftru(nt,i,j) = 0.
              ftrv(nt,i,j) = 0.
            end do
          end if
          cu(i,j) = 0.
          cv(i,j) = 0.
        end if
      end do
    end do

    if (use_TRC .and. use_ATRC) then
      do nt = 1,natr
        nat = ntr-natr+nt
        do j = 1-mrg-2,jj+mrg+2
          do l = 1,isp(j)
            do i = max(1-mrg-2,ifp(j,l)),min(ii+mrg+2,ilp(j,l))
              trc(i,j,k,nt) = max(0.,trc(i,j,k,nt))+treps
              ag(nt,i,j) = trc(i,j,k,nat)/trc(i,j,k,nt)
            end do
          end do
        end do
      end do
    end if

    !---------------------------------------------------------------
    ! Compute limited gradients, center of mass coordinates, and
    ! non-dimensional velocities.
    !---------------------------------------------------------------

    do j = 1-mrg-1,jj+mrg+1

      do l = 1,isp(j)
        do i = max(1-mrg-1,ifp(j,l)),min(ii+mrg+1,ilp(j,l))

          ! Define indices for grid cell neighbors, ensuring that only wet
          ! points are used.
          iw = i-iu(i  ,j)
          ie = i+iu(i+1,j)
          js = j-iv(i,j  )
          jn = j+iv(i,j+1)
          isw = i*(1-ip(iw,js))+iw*ip(iw,js)
          jsw = j*(1-ip(iw,js))+js*ip(iw,js)
          ise = i*(1-ip(ie,js))+ie*ip(ie,js)
          jse = j*(1-ip(ie,js))+js*ip(ie,js)
          inw = i*(1-ip(iw,jn))+iw*ip(iw,jn)
          jnw = j*(1-ip(iw,jn))+jn*ip(iw,jn)
          ine = i*(1-ip(ie,jn))+ie*ip(ie,jn)
          jne = j*(1-ip(ie,jn))+jn*ip(ie,jn)

          dxi = 1./max(1,ie-iw)
          dyi = 1./max(1,jn-js)

          ! Compute limited gradient for layer pressure thickness and
          ! center of mass coordinate.
          dpsw = max(dpeps,min(pbmin(i,j)-pup(isw,jsw),dp(isw,jsw)))
          dps  = max(dpeps,min(pbmin(i,j)-pup(i  ,js ),dp(i  ,js )))
          dpse = max(dpeps,min(pbmin(i,j)-pup(ise,jse),dp(ise,jse)))
          dpw  = max(dpeps,min(pbmin(i,j)-pup(iw ,j  ),dp(iw ,j  )))
          dpc  = max(dpeps,min(pbmin(i,j)-pup(i  ,j  ),dp(i  ,j  )))
          dpe  = max(dpeps,min(pbmin(i,j)-pup(ie ,j  ),dp(ie ,j  )))
          dpnw = max(dpeps,min(pbmin(i,j)-pup(inw,jnw),dp(inw,jnw)))
          dpn  = max(dpeps,min(pbmin(i,j)-pup(i  ,jn ),dp(i  ,jn )))
          dpne = max(dpeps,min(pbmin(i,j)-pup(ine,jne),dp(ine,jne)))
          dx(i,j) = (dpe-dpw)*dxi
          dy(i,j) = (dpn-dps)*dyi
          dgmx = .5*(abs(dx(i,j))+abs(dy(i,j)))
          dfmx = max(0.,max(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
          dfmn = min(0.,min(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
          if (dfmx > 0..and.dfmn < 0.) then
            q = min(dfmx/max(dfmx,dgmx),dfmn/min(dfmn,-dgmx))
            dx(i,j) = dx(i,j)*q
            dy(i,j) = dy(i,j)*q
            xd(i,j) = dx(i,j)/(12.*dp(i,j))
            yd(i,j) = dy(i,j)/(12.*dp(i,j))
          else
            dx(i,j) = 0.
            dy(i,j) = 0.
            xd(i,j) = 0.
            yd(i,j) = 0.
          end if

          ! Compute limited gradients for temperature, salinity, and
          ! density
          tx(i,j) = (temp(ie,j)-temp(iw,j))*dxi
          ty(i,j) = (temp(i,jn)-temp(i,js))*dyi
          q1 = tx(i,j)*(-.5-xd(i,j))
          q2 = tx(i,j)*( .5-xd(i,j))
          q3 = ty(i,j)*(-.5-yd(i,j))
          q4 = ty(i,j)*( .5-yd(i,j))
          tgmx = max(q1,q2)+max(q3,q4)
          tgmn = min(q1,q2)+min(q3,q4)
          tfmx = max(0.,max(temp(isw,jsw),temp(i  ,js ), &
               temp(ise,jse),temp(iw ,j  ), &
               temp(ie ,j  ),temp(inw,jnw), &
               temp(i  ,jn ),temp(ine,jne)) &
               -temp(i,j))
          tfmn = min(0.,min(temp(isw,jsw),temp(i  ,js ), &
               temp(ise,jse),temp(iw ,j  ), &
               temp(ie ,j  ),temp(inw,jnw), &
               temp(i  ,jn ),temp(ine,jne)) &
               -temp(i,j))
          if (tfmx > 0..and.tfmn < 0.) then
            q = min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
            tx(i,j) = tx(i,j)*q
            ty(i,j) = ty(i,j)*q
            td(i,j) = temp(i,j)-tx(i,j)*xd(i,j)-ty(i,j)*yd(i,j)
          else
            tx(i,j) = 0.
            ty(i,j) = 0.
            td(i,j) = temp(i,j)
          end if

          sx(i,j) = (saln(ie,j)-saln(iw,j))*dxi
          sy(i,j) = (saln(i,jn)-saln(i,js))*dyi
          q1 = sx(i,j)*(-.5-xd(i,j))
          q2 = sx(i,j)*( .5-xd(i,j))
          q3 = sy(i,j)*(-.5-yd(i,j))
          q4 = sy(i,j)*( .5-yd(i,j))
          sgmx = max(q1,q2)+max(q3,q4)
          sgmn = min(q1,q2)+min(q3,q4)
          sfmx = max(0.,max(saln(isw,jsw),saln(i  ,js ), &
               saln(ise,jse),saln(iw ,j  ), &
               saln(ie ,j  ),saln(inw,jnw), &
               saln(i  ,jn ),saln(ine,jne)) &
               -saln(i,j))
          sfmn = min(0.,min(saln(isw,jsw),saln(i  ,js ), &
               saln(ise,jse),saln(iw ,j  ), &
               saln(ie ,j  ),saln(inw,jnw), &
               saln(i  ,jn ),saln(ine,jne)) &
               -saln(i,j))
          if (sfmx > 0..and.sfmn < 0.) then
            q = min(sfmx/max(sfmx,sgmx),sfmn/min(sfmn,sgmn))
            sx(i,j) = sx(i,j)*q
            sy(i,j) = sy(i,j)*q
            sd(i,j) = saln(i,j)-sx(i,j)*xd(i,j)-sy(i,j)*yd(i,j)
          else
            sx(i,j) = 0.
            sy(i,j) = 0.
            sd(i,j) = saln(i,j)
          end if
          if (use_TRC) then
            if (use_ATRC) then
              ! Compute limited gradient for tracers.
              do nt = 1,ntr-natr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trx(nt,i,j) = (trc(ie,j,k,nt)-trc(iw,j,k,nt))*dxi
                try(nt,i,j) = (trc(i,jn,k,nt)-trc(i,js,k,nt))*dyi
                q1 = trx(nt,i,j)*(-.5-xd(i,j))
                q2 = trx(nt,i,j)*( .5-xd(i,j))
                q3 = try(nt,i,j)*(-.5-yd(i,j))
                q4 = try(nt,i,j)*( .5-yd(i,j))
                tgmx = max(q1,q2)+max(q3,q4)
                tgmn = min(q1,q2)+min(q3,q4)
                tfmx = max(0.,max(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt), &
                     trc(ise,jse,k,nt),trc(iw ,j  ,k,nt), &
                     trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt), &
                     trc(i  ,jn ,k,nt),trc(ine,jne,k,nt)) &
                     -trc(i,j,k,nt))
                tfmn = min(0.,min(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt), &
                     trc(ise,jse,k,nt),trc(iw ,j  ,k,nt), &
                     trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt), &
                     trc(i  ,jn ,k,nt),trc(ine,jne,k,nt)) &
                     -trc(i,j,k,nt))
                if (tfmx > 0..and.tfmn < 0.) then
                  q = min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                  trx(nt,i,j) = trx(nt,i,j)*q
                  try(nt,i,j) = try(nt,i,j)*q
                  trd(nt,i,j) = trc(i,j,k,nt) &
                       -trx(nt,i,j)*xd(i,j)-try(nt,i,j)*yd(i,j)
                else
                  trx(nt,i,j) = 0.
                  try(nt,i,j) = 0.
                  trd(nt,i,j) = trc(i,j,k,nt)
                end if
              end do

              ! Compute limited gradient for age tracers.
              do nt = 1,natr
                nat = ntr-natr+nt
                agx(nt,i,j) = (ag(nt,ie,j)-ag(nt,iw,j))*dxi
                agy(nt,i,j) = (ag(nt,i,jn)-ag(nt,i,js))*dyi
                q = 1./(12.*trc(i,j,k,nt))
                xdt = (12.*xd(i,j)*trd(nt,i,j)+trx(nt,i,j))*q
                ydt = (12.*yd(i,j)*trd(nt,i,j)+try(nt,i,j))*q
                q1 = agx(nt,i,j)*(-.5-xdt)
                q2 = agx(nt,i,j)*( .5-xdt)
                q3 = agy(nt,i,j)*(-.5-ydt)
                q4 = agy(nt,i,j)*( .5-ydt)
                tgmx = max(q1,q2)+max(q3,q4)
                tgmn = min(q1,q2)+min(q3,q4)
                tfmx = max(0.,max(ag(nt,isw,jsw),ag(nt,i  ,js ), &
                                  ag(nt,ise,jse),ag(nt,iw ,j  ), &
                                  ag(nt,ie ,j  ),ag(nt,inw,jnw), &
                                  ag(nt,i  ,jn ),ag(nt,ine,jne)) &
                              -ag(nt,i,j))
                tfmn = min(0.,min(ag(nt,isw,jsw),ag(nt,i  ,js ), &
                                  ag(nt,ise,jse),ag(nt,iw ,j  ), &
                                  ag(nt,ie ,j  ),ag(nt,inw,jnw), &
                                  ag(nt,i  ,jn ),ag(nt,ine,jne)) &
                             -ag(nt,i,j))
                if (tfmx > 0..and.tfmn < 0.) then
                  q = min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                  agx(nt,i,j) = agx(nt,i,j)*q
                  agy(nt,i,j) = agy(nt,i,j)*q
                  agd(nt,i,j) = ag(nt,i,j)-agx(nt,i,j)*xdt-agy(nt,i,j)*ydt
                else
                  agx(nt,i,j) = 0.
                  agy(nt,i,j) = 0.
                  agd(nt,i,j) = ag(nt,i,j)
                end if
              end do
            else

              ! Compute limited gradient for tracers.
              do nt = 1,ntr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trx(nt,i,j) = (trc(ie,j,k,nt)-trc(iw,j,k,nt))*dxi
                try(nt,i,j) = (trc(i,jn,k,nt)-trc(i,js,k,nt))*dyi
                q1 = trx(nt,i,j)*(-.5-xd(i,j))
                q2 = trx(nt,i,j)*( .5-xd(i,j))
                q3 = try(nt,i,j)*(-.5-yd(i,j))
                q4 = try(nt,i,j)*( .5-yd(i,j))
                tgmx = max(q1,q2)+max(q3,q4)
                tgmn = min(q1,q2)+min(q3,q4)
                tfmx = max(0.,max(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt), &
                                  trc(ise,jse,k,nt),trc(iw ,j  ,k,nt), &
                                  trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt), &
                                  trc(i  ,jn ,k,nt),trc(ine,jne,k,nt)) &
                             -trc(i,j,k,nt))
                tfmn = min(0.,min(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt), &
                                  trc(ise,jse,k,nt),trc(iw ,j  ,k,nt), &
                                  trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt), &
                                  trc(i  ,jn ,k,nt),trc(ine,jne,k,nt)) &
                             -trc(i,j,k,nt))
                if (tfmx > 0..and.tfmn < 0.) then
                  q = min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                  trx(nt,i,j) = trx(nt,i,j)*q
                  try(nt,i,j) = try(nt,i,j)*q
                  trd(nt,i,j) = trc(i,j,k,nt)&
                               -trx(nt,i,j)*xd(i,j)-try(nt,i,j)*yd(i,j)
                else
                  trx(nt,i,j) = 0.
                  try(nt,i,j) = 0.
                  trd(nt,i,j) = trc(i,j,k,nt)
                end if
              end do
            end if
          end if

        end do
      end do
    end do

    ! Compute non-dimensional velocities.

    do j = 1-mrg-1,jj+mrg+1
      do l = 1,isu(j)
        do i = max(1-mrg,ifu(j,l)),min(ii+mrg+1,ilu(j,l))
          if (u(i,j) > 0.) then
            cu(i,j) = u(i,j)*dt*scuy(i,j)*scp2i(i-1,j)
          else
            cu(i,j) = u(i,j)*dt*scuy(i,j)*scp2i(i  ,j)
          end if
        end do
      end do
    end do

    do j = 1-mrg,jj+mrg+1
      do l = 1,isv(j)
        do i = max(1-mrg-1,ifv(j,l)),min(ii+mrg+1,ilv(j,l))
          if (v(i,j) > 0.) then
            cv(i,j) = v(i,j)*dt*scvx(i,j)*scp2i(i,j-1)
          else
            cv(i,j) = v(i,j)*dt*scvx(i,j)*scp2i(i,j  )
          end if
        end do
      end do
    end do

    !---------------------------------------------------------------
    ! Compute corner velocities. The velocity components are computed as
    ! the harmonic mean of the nearest C-grid velocity components with
    ! the following exeptions: The corner velocity component is set to
    ! zero if the nearest C-grid components have different sign, or one
    ! or tree of the neighboring grid cells are wet, or two neighbors
    ! are wet and are arranged diagonally. This construction of corner
    ! velocities will ensure that the entire fluxing area is located
    ! upwind of the cell boundary.
    !---------------------------------------------------------------

    do j = 1-mrg,jj+mrg+1
      do i = 1-mrg,ii+mrg+1
        nw = ip(i-1,j-1)+ip(i,j-1)+ip(i-1,j)+ip(i,j)
        if     (nw == 4) then
          if (cu(i,j-1)*cu(i,j) <= 0.) then
            cuc(i,j) = 0.
          else
            cuc(i,j) = 2.*cu(i,j-1)*cu(i,j)/(cu(i,j-1)+cu(i,j))
          end if
          if (cv(i-1,j)*cv(i,j) <= 0.) then
            cvc(i,j) = 0.
          else
            cvc(i,j) = 2.*cv(i-1,j)*cv(i,j)/(cv(i-1,j)+cv(i,j))
          end if
        else if (nw == 2) then
          if     (ip(i-1,j-1)+ip(i  ,j-1) == 2) then
            cuc(i,j) = cu(i,j-1)
            cvc(i,j) = 0.
          else if (ip(i-1,j  )+ip(i  ,j  ) == 2) then
            cuc(i,j) = cu(i,j  )
            cvc(i,j) = 0.
          else if (ip(i-1,j-1)+ip(i-1,j  ) == 2) then
            cuc(i,j) = 0.
            cvc(i,j) = cv(i-1,j)
          else if (ip(i  ,j-1)+ip(i  ,j  ) == 2) then
            cuc(i,j) = 0.
            cvc(i,j) = cv(i  ,j)
          else
            cuc(i,j) = 0.
            cvc(i,j) = 0.
          end if
        else
          cuc(i,j) = 0.
          cvc(i,j) = 0.
        end if
      end do
    end do

    !---------------------------------------------------------------
    ! Compute cell boundary fluxes.
    !---------------------------------------------------------------

    ! - u-components of fluxes.

    do j = 1-mrg,jj+mrg

      do l = 1,isu(j)
        do i = max(1-mrg,ifu(j,l)),min(ii+mrg+1,ilu(j,l))

          ! Assuming coordinate [0,0] at the u-point, the non-dimensional
          ! fluxing area is defined as the area of a polygon with vertices
          ! [0,1/2], [-cuc(i,j+1),-cvc(i,j+1)+1/2], [xm,ym],
          ! [-cuc(i,j),-cvc(i,j)-1/2], and [0,-1/2]. The vertex [xm,ym] is
          ! defined so that the polygon area is equal to cu(i,j).

          ym = -.5*(cvc(i,j)+cvc(i,j+1))
          xm = ((ym+.5)*cuc(i,j)-(ym-.5)*cuc(i,j+1)-2.*cu(i,j)) &
               /(1.+cvc(i,j)-cvc(i,j+1))

          if (cu(i,j) > 0.) then

            if (cvc(i,j) > 0.) then

              ! Add contributions from grid cell (i-1,j-1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [xc1+1/2,1/2], [-cuc(i,j)+1/2,-cvc(i,j)+1/2], and
              ! [1/2,1/2].

              xc0 = (xm*cvc(i,j)-cuc(i,j)*(ym+.5))/(cvc(i,j)+ym+.5)
              xc1 = xc0*scp2(i-1,j)*scp2i(i-1,j-1)
              x4 = xc0+.5
              y4 = -.5
              call triint(scp2(i-1,j-1), &
                   xc1+.5,.5,-cuc(i,j)+.5,-cvc(i,j)+.5,.5,.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)
              dl = min(dp(i-1,j-1),max(0.,pbu(i,j)-pup(i-1,j-1)))
              fd = a*dl+ax*dx(i-1,j-1)+ay*dy(i-1,j-1)
              fdu(i,j) = fdu(i,j)+fd
              qx = ax*dl+axx*dx(i-1,j-1)+axy*dy(i-1,j-1)
              qy = ay*dl+axy*dx(i-1,j-1)+ayy*dy(i-1,j-1)
              ftu(i,j) = ftu(i,j)+fd*td(i-1,j-1) &
                        +qx*tx(i-1,j-1)+qy*ty(i-1,j-1)
              fsu(i,j) = fsu(i,j)+fd*sd(i-1,j-1) &
                        +qx*sx(i-1,j-1)+qy*sy(i-1,j-1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i-1,j-1)+axxy*dy(i-1,j-1)
                  qyy = ayy*dl+axyy*dx(i-1,j-1)+ayyy*dy(i-1,j-1)
                  qxy = axy*dl+axxy*dx(i-1,j-1)+axyy*dy(i-1,j-1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i-1,j-1) &
                         +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                    ftru(nt,i,j) = ftru(nt,i,j)+fdt
                    fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i-1,j-1) &
                                  +(qx *trd(nt,i-1,j-1) &
                                   +qxx*trx(nt,i-1,j-1) &
                                   +qxy*try(nt,i-1,j-1))*agx(nt,i-1,j-1) &
                                  +(qy *trd(nt,i-1,j-1) &
                                   +qxy*trx(nt,i-1,j-1) &
                                   +qyy*try(nt,i-1,j-1))*agy(nt,i-1,j-1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j-1) &
                                  +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j-1) &
                                  +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                  end do
                end if
              end if
            else
              x4 = -cuc(i,j)+.5
              y4 = -cvc(i,j)-.5
            end if

            if (cvc(i,j+1) < 0.) then

              ! Add contributions from grid cell (i-1,j+1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [xc1+1/2,-1/2], [1/2,-1/2], and
              ! [-cuc(i,j+1)+1/2,-cvc(i,j+1)-1/2].

              xc0 = (xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
              xc1 = xc0*scp2(i-1,j)*scp2i(i-1,j+1)
              x2 = xc0+.5
              y2 = .5
              call triint(scp2(i-1,j+1), &
                   xc1+.5,-.5,.5,-.5,-cuc(i,j+1)+.5,-cvc(i,j+1)-.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)
              dl = min(dp(i-1,j+1),max(0.,pbu(i,j)-pup(i-1,j+1)))
              fd = a*dl+ax*dx(i-1,j+1)+ay*dy(i-1,j+1)
              fdu(i,j) = fdu(i,j)+fd
              qx = ax*dl+axx*dx(i-1,j+1)+axy*dy(i-1,j+1)
              qy = ay*dl+axy*dx(i-1,j+1)+ayy*dy(i-1,j+1)
              ftu(i,j) = ftu(i,j)+fd*td(i-1,j+1) &
                        +qx*tx(i-1,j+1)+qy*ty(i-1,j+1)
              fsu(i,j) = fsu(i,j)+fd*sd(i-1,j+1) &
                        +qx*sx(i-1,j+1)+qy*sy(i-1,j+1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i-1,j+1)+axxy*dy(i-1,j+1)
                  qyy = ayy*dl+axyy*dx(i-1,j+1)+ayyy*dy(i-1,j+1)
                  qxy = axy*dl+axxy*dx(i-1,j+1)+axyy*dy(i-1,j+1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i-1,j+1) &
                         +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                    ftru(nt,i,j) = ftru(nt,i,j)+fdt
                    fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i-1,j+1) &
                                  +(qx *trd(nt,i-1,j+1) &
                                   +qxx*trx(nt,i-1,j+1) &
                                   +qxy*try(nt,i-1,j+1))*agx(nt,i-1,j+1) &
                                  +(qy *trd(nt,i-1,j+1) &
                                   +qxy*trx(nt,i-1,j+1) &
                                   +qyy*try(nt,i-1,j+1))*agy(nt,i-1,j+1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j+1) &
                                  +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j+1) &
                                  +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                  end do
                end if
              end if
            else
              x2 = -cuc(i,j+1)+.5
              y2 = -cvc(i,j+1)+.5
            end if

            !-- Add contributions from grid cell (i-1,j). Assuming
            !-- coordinate [0,0] at the cell center, the contributions are
            !-- flux integrals over the pentagon with vertices [1/2,1/2],
            !-- [x2,y2], [xm+1/2,ym], [x4,y4], and [1/2,-1/2].

            call penint(scp2(i-1,j), &
                 .5,.5,x2,y2,xm+.5,ym,x4,y4,.5,-.5, &
                 a,ax,ay,axx,ayy,axy, &
                 axxx,ayyy,axxy,axyy)

            dl = min(dp(i-1,j),max(0.,pbu(i,j)-pup(i-1,j)))
            fd = a*dl+ax*dx(i-1,j)+ay*dy(i-1,j)
            fdu(i,j) = fdu(i,j)+fd
            qx = ax*dl+axx*dx(i-1,j)+axy*dy(i-1,j)
            qy = ay*dl+axy*dx(i-1,j)+ayy*dy(i-1,j)
            ftu(i,j) = ftu(i,j)+fd*td(i-1,j)+qx*tx(i-1,j)+qy*ty(i-1,j)
            fsu(i,j) = fsu(i,j)+fd*sd(i-1,j)+qx*sx(i-1,j)+qy*sy(i-1,j)
            if (use_TRC) then
              if (use_ATRC) then
                qxx = axx*dl+axxx*dx(i-1,j)+axxy*dy(i-1,j)
                qyy = ayy*dl+axyy*dx(i-1,j)+ayyy*dy(i-1,j)
                qxy = axy*dl+axxy*dx(i-1,j)+axyy*dy(i-1,j)
                do nt = 1,natr
                  fdt = fd*trd(nt,i-1,j) &
                       +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                  ftru(nt,i,j) = ftru(nt,i,j)+fdt
                  fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i-1,j) &
                                +(qx *trd(nt,i-1,j) &
                                 +qxx*trx(nt,i-1,j) &
                                 +qxy*try(nt,i-1,j))*agx(nt,i-1,j) &
                                +(qy *trd(nt,i-1,j) &
                                 +qxy*trx(nt,i-1,j) &
                                 +qyy*try(nt,i-1,j))*agy(nt,i-1,j)
                end do
                do nt = natr+1,ntr-natr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j) &
                                +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                end do
              else
                do nt = 1,ntr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j) &
                                +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                end do
              end if
            end if

          else

            if (cvc(i,j) > 0.) then

              ! Add contributions from grid cell (i,j-1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [xc1-1/2,1/2], [-cuc(i,j)-1/2,-cvc(i,j)+1/2], and
              ! [-1/2,1/2].

              xc0 = (xm*cvc(i,j)-cuc(i,j)*(ym+.5))/(cvc(i,j)+ym+.5)
              xc1 = xc0*scp2(i,j)*scp2i(i,j-1)
              x4 = xc0-.5
              y4 = -.5

              call triint(scp2(i,j-1), &
                   xc1-.5,.5,-cuc(i,j)-.5,-cvc(i,j)+.5,-.5,.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i,j-1),max(0.,pbu(i,j)-pup(i,j-1)))
              fd = a*dl+ax*dx(i,j-1)+ay*dy(i,j-1)
              fdu(i,j) = fdu(i,j)+fd
              qx = ax*dl+axx*dx(i,j-1)+axy*dy(i,j-1)
              qy = ay*dl+axy*dx(i,j-1)+ayy*dy(i,j-1)
              ftu(i,j) = ftu(i,j)+fd*td(i,j-1) &
                        +qx*tx(i,j-1)+qy*ty(i,j-1)
              fsu(i,j) = fsu(i,j)+fd*sd(i,j-1) &
                        +qx*sx(i,j-1)+qy*sy(i,j-1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i,j-1)+axxy*dy(i,j-1)
                  qyy = ayy*dl+axyy*dx(i,j-1)+ayyy*dy(i,j-1)
                  qxy = axy*dl+axxy*dx(i,j-1)+axyy*dy(i,j-1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i,j-1) &
                         +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                    ftru(nt,i,j) = ftru(nt,i,j)+fdt
                    fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i,j-1) &
                                  +(qx *trd(nt,i,j-1) &
                                   +qxx*trx(nt,i,j-1) &
                                   +qxy*try(nt,i,j-1))*agx(nt,i,j-1) &
                                  +(qy *trd(nt,i,j-1) &
                                   +qxy*trx(nt,i,j-1) &
                                   +qyy*try(nt,i,j-1))*agy(nt,i,j-1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j-1) &
                                  +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j-1) &
                                  +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                  end do
                end if
              end if
            else
              x4 = -cuc(i,j)-.5
              y4 = -cvc(i,j)-.5
            end if

            if (cvc(i,j+1) < 0.) then

              ! Add contributions from grid cell (i,j+1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [xc1-1/2,-1/2], [-1/2,-1/2], and
              ! [-cuc(i,j+1)-1/2,-cvc(i,j+1)-1/2].

              xc0 = (xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
              xc1 = xc0*scp2(i,j)*scp2i(i,j+1)
              x2 = xc0-.5
              y2 = .5

              call triint(scp2(i,j+1), &
                   xc1-.5,-.5,-.5,-.5,-cuc(i,j+1)-.5,-cvc(i,j+1)-.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i,j+1),max(0.,pbu(i,j)-pup(i,j+1)))
              fd = a*dl+ax*dx(i,j+1)+ay*dy(i,j+1)
              fdu(i,j) = fdu(i,j)+fd
              qx = ax*dl+axx*dx(i,j+1)+axy*dy(i,j+1)
              qy = ay*dl+axy*dx(i,j+1)+ayy*dy(i,j+1)
              ftu(i,j) = ftu(i,j)+fd*td(i,j+1) &
                   +qx*tx(i,j+1)+qy*ty(i,j+1)
              fsu(i,j) = fsu(i,j)+fd*sd(i,j+1) &
                   +qx*sx(i,j+1)+qy*sy(i,j+1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i,j+1)+axxy*dy(i,j+1)
                  qyy = ayy*dl+axyy*dx(i,j+1)+ayyy*dy(i,j+1)
                  qxy = axy*dl+axxy*dx(i,j+1)+axyy*dy(i,j+1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i,j+1) &
                         +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                    ftru(nt,i,j) = ftru(nt,i,j)+fdt
                    fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i,j+1) &
                                  +(qx *trd(nt,i,j+1) &
                                   +qxx*trx(nt,i,j+1) &
                                   +qxy*try(nt,i,j+1))*agx(nt,i,j+1) &
                                  +(qy *trd(nt,i,j+1) &
                                   +qxy*trx(nt,i,j+1) &
                                   +qyy*try(nt,i,j+1))*agy(nt,i,j+1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j+1) &
                                  +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j+1) &
                                  +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                  end do
                end if
              end if
            else
              x2 = -cuc(i,j+1)-.5
              y2 = -cvc(i,j+1)+.5
            end if

            !-- Add contributions from grid cell (i,j). Assuming
            !-- coordinate [0,0] at the cell center, the contributions are
            !-- flux integrals over the pentagon with vertices [-1/2,1/2],
            !-- [x2,y2], [xm-1/2,ym], [x4,y4], and [-1/2,-1/2].

            call penint(scp2(i,j), &
                 -.5,.5,x2,y2,xm-.5,ym,x4,y4,-.5,-.5, &
                 a,ax,ay,axx,ayy,axy, &
                 axxx,ayyy,axxy,axyy)

            dl = min(dp(i,j),max(0.,pbu(i,j)-pup(i,j)))
            fd = a*dl+ax*dx(i,j)+ay*dy(i,j)
            fdu(i,j) = fdu(i,j)+fd
            qx = ax*dl+axx*dx(i,j)+axy*dy(i,j)
            qy = ay*dl+axy*dx(i,j)+ayy*dy(i,j)
            ftu(i,j) = ftu(i,j)+fd*td(i,j)+qx*tx(i,j)+qy*ty(i,j)
            fsu(i,j) = fsu(i,j)+fd*sd(i,j)+qx*sx(i,j)+qy*sy(i,j)

            if (use_TRC) then
              if (use_ATRC) then
                qxx = axx*dl+axxx*dx(i,j)+axxy*dy(i,j)
                qyy = ayy*dl+axyy*dx(i,j)+ayyy*dy(i,j)
                qxy = axy*dl+axxy*dx(i,j)+axyy*dy(i,j)
                do nt = 1,natr
                  fdt = fd*trd(nt,i,j) &
                       +qx*trx(nt,i,j)+qy*try(nt,i,j)
                  ftru(nt,i,j) = ftru(nt,i,j)+fdt
                  fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i,j) &
                                +(qx *trd(nt,i,j) &
                                 +qxx*trx(nt,i,j) &
                                 +qxy*try(nt,i,j))*agx(nt,i,j) &
                                +(qy *trd(nt,i,j) &
                                 +qxy*trx(nt,i,j) &
                                 +qyy*try(nt,i,j))*agy(nt,i,j)
                end do
                do nt = natr+1,ntr-natr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j) &
                       +qx*trx(nt,i,j)+qy*try(nt,i,j)
                end do
              else
                do nt = 1,ntr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j) &
                                +qx*trx(nt,i,j)+qy*try(nt,i,j)
                end do
              end if
            end if

          end if

          ! u-component of mass, heat and salt flux.
          uflx(i,j) = uflx(i,j)+fdu(i,j)
          utflx(i,j) = utflx(i,j)+ftu(i,j)
          usflx(i,j) = usflx(i,j)+fsu(i,j)

        end do
      end do

    end do

    ! v-components of fluxes.

    do j = 1-mrg,jj+mrg+1

      do l = 1,isv(j)
        do i = max(1-mrg,ifv(j,l)),min(ii+mrg,ilv(j,l))

          ! Assuming coordinate [0,0] at the v-point, the non-dimensional
          ! fluxing area is defined as the area of a polygon with vertices
          ! [-1/2,0], [-cuc(i,j)-1/2,-cvc(i,j)], [xm,ym],
          ! [-cuc(i+1,j)+1/2,-cvc(i+1,j)], and [1/2,0]. The vertex [xm,ym]
          ! is defined so that the polygon area is equal to cv(i,j).

          xm = -.5*(cuc(i,j)+cuc(i+1,j))
          ym = ((xm+.5)*cvc(i,j)-(xm-.5)*cvc(i+1,j)-2.*cv(i,j)) &
               /(1.+cuc(i,j)-cuc(i+1,j))

          if (cv(i,j) > 0) then

            if (cuc(i,j) > 0.) then

              ! Add contributions from grid cell (i-1,j-1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [1/2,yc1+1/2], [1/2,1/2], and
              ! [-cuc(i,j)+1/2,-cvc(i,j)+1/2].

              yc0 = (ym*cuc(i,j)-cvc(i,j)*(xm+.5))/(cuc(i,j)+xm+.5)
              yc1 = yc0*scp2(i,j-1)*scp2i(i-1,j-1)
              x2 = -.5
              y2 = yc0+.5

              call triint(scp2(i-1,j-1), &
                   .5,yc1+.5,.5,.5,-cuc(i,j)+.5,-cvc(i,j)+.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i-1,j-1),max(0.,pbv(i,j)-pup(i-1,j-1)))
              fd = a*dl+ax*dx(i-1,j-1)+ay*dy(i-1,j-1)
              fdv(i,j) = fdv(i,j)+fd
              qx = ax*dl+axx*dx(i-1,j-1)+axy*dy(i-1,j-1)
              qy = ay*dl+axy*dx(i-1,j-1)+ayy*dy(i-1,j-1)
              ftv(i,j) = ftv(i,j)+fd*td(i-1,j-1) &
                        +qx*tx(i-1,j-1)+qy*ty(i-1,j-1)
              fsv(i,j) = fsv(i,j)+fd*sd(i-1,j-1) &
                        +qx*sx(i-1,j-1)+qy*sy(i-1,j-1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i-1,j-1)+axxy*dy(i-1,j-1)
                  qyy = ayy*dl+axyy*dx(i-1,j-1)+ayyy*dy(i-1,j-1)
                  qxy = axy*dl+axxy*dx(i-1,j-1)+axyy*dy(i-1,j-1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i-1,j-1) &
                         +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                    fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i-1,j-1) &
                                  +(qx *trd(nt,i-1,j-1) &
                                   +qxx*trx(nt,i-1,j-1) &
                                   +qxy*try(nt,i-1,j-1))*agx(nt,i-1,j-1) &
                                  +(qy *trd(nt,i-1,j-1) &
                                   +qxy*trx(nt,i-1,j-1) &
                                   +qyy*try(nt,i-1,j-1))*agy(nt,i-1,j-1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i-1,j-1) &
                                  +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i-1,j-1) &
                                  +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                  end do
                end if
              end if
            else
              x2 = -cuc(i,j)-.5
              y2 = -cvc(i,j)+.5
            end if

            if (cuc(i+1,j) < 0.) then

              ! Add contributions from grid cell (i+1,j-1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [-1/2,yc1+1/2], [-cuc(i+1,j)-1/2,-cvc(i+1,j)+1/2], and
              ! [-1/2,1/2].

              yc0 = (ym*cuc(i+1,j)-cvc(i+1,j)*(xm-.5))/(cuc(i+1,j)+xm-.5)
              yc1 = yc0*scp2(i,j-1)*scp2i(i+1,j-1)
              x4 = .5
              y4 = yc0+.5

              call triint(scp2(i+1,j-1), &
                   -.5,yc1+.5,-cuc(i+1,j)-.5,-cvc(i+1,j)+.5,-.5,.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i+1,j-1),max(0.,pbv(i,j)-pup(i+1,j-1)))
              fd = a*dl+ax*dx(i+1,j-1)+ay*dy(i+1,j-1)
              fdv(i,j) = fdv(i,j)+fd
              qx = ax*dl+axx*dx(i+1,j-1)+axy*dy(i+1,j-1)
              qy = ay*dl+axy*dx(i+1,j-1)+ayy*dy(i+1,j-1)
              ftv(i,j) = ftv(i,j)+fd*td(i+1,j-1) &
                        +qx*tx(i+1,j-1)+qy*ty(i+1,j-1)
              fsv(i,j) = fsv(i,j)+fd*sd(i+1,j-1) &
                        +qx*sx(i+1,j-1)+qy*sy(i+1,j-1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i+1,j-1)+axxy*dy(i+1,j-1)
                  qyy = ayy*dl+axyy*dx(i+1,j-1)+ayyy*dy(i+1,j-1)
                  qxy = axy*dl+axxy*dx(i+1,j-1)+axyy*dy(i+1,j-1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i+1,j-1) &
                         +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                    fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i+1,j-1) &
                                  +(qx *trd(nt,i+1,j-1) &
                                   +qxx*trx(nt,i+1,j-1) &
                                   +qxy*try(nt,i+1,j-1))*agx(nt,i+1,j-1) &
                                  +(qy *trd(nt,i+1,j-1) &
                                   +qxy*trx(nt,i+1,j-1) &
                                   +qyy*try(nt,i+1,j-1))*agy(nt,i+1,j-1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i+1,j-1) &
                                  +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i+1,j-1) &
                                  +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
                  end do
                end if
              end if
            else
              x4 = -cuc(i+1,j)+.5
              y4 = -cvc(i+1,j)+.5
            end if

            !-- Add contributions from grid cell (i,j-1). Assuming
            !-- coordinate [0,0] at the cell center, the contributions are
            !-- flux integrals over the pentagon with vertices [-1/2,1/2],
            !-- [x2,y2], [xm,ym+1/2], [x4,y4], and [1/2,1/2].

            call penint(scp2(i,j-1), &
                 -.5,.5,x2,y2,xm,ym+.5,x4,y4,.5,.5, &
                 a,ax,ay,axx,ayy,axy, &
                 axxx,ayyy,axxy,axyy)

            dl = min(dp(i,j-1),max(0.,pbv(i,j)-pup(i,j-1)))
            fd = a*dl+ax*dx(i,j-1)+ay*dy(i,j-1)
            fdv(i,j) = fdv(i,j)+fd
            qx = ax*dl+axx*dx(i,j-1)+axy*dy(i,j-1)
            qy = ay*dl+axy*dx(i,j-1)+ayy*dy(i,j-1)
            ftv(i,j) = ftv(i,j)+fd*td(i,j-1)+qx*tx(i,j-1)+qy*ty(i,j-1)
            fsv(i,j) = fsv(i,j)+fd*sd(i,j-1)+qx*sx(i,j-1)+qy*sy(i,j-1)
            if (use_TRC) then
              if (use_ATRC) then
                qxx = axx*dl+axxx*dx(i,j-1)+axxy*dy(i,j-1)
                qyy = ayy*dl+axyy*dx(i,j-1)+ayyy*dy(i,j-1)
                qxy = axy*dl+axxy*dx(i,j-1)+axyy*dy(i,j-1)
                do nt = 1,natr
                  fdt = fd*trd(nt,i,j-1)+qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                  fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i,j-1) &
                                +(qx *trd(nt,i,j-1) &
                                 +qxx*trx(nt,i,j-1) &
                                 +qxy*try(nt,i,j-1))*agx(nt,i,j-1) &
                                +(qy *trd(nt,i,j-1) &
                                 +qxy*trx(nt,i,j-1) &
                                 +qyy*try(nt,i,j-1))*agy(nt,i,j-1)
                end do
                do nt = natr+1,ntr-natr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i,j-1) &
                                +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                end do
              else
                do nt = 1,ntr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i,j-1) &
                                +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                end do
              end if
            end if
          else

            if (cuc(i,j) > 0.) then

              ! Add contributions from grid cell (i-1,j). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [1/2,yc1-1/2], [1/2,-1/2], and
              ! [-cuc(i,j)+1/2,-cvc(i,j)-1/2].

              yc0 = (ym*cuc(i,j)-cvc(i,j)*(xm+.5))/(cuc(i,j)+xm+.5)
              yc1 = yc0*scp2(i,j)*scp2i(i-1,j)
              x2 = -.5
              y2 = yc0-.5

              call triint(scp2(i-1,j), &
                   .5,yc1-.5,.5,-.5,-cuc(i,j)+.5,-cvc(i,j)-.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i-1,j),max(0.,pbv(i,j)-pup(i-1,j)))
              fd = a*dl+ax*dx(i-1,j)+ay*dy(i-1,j)
              fdv(i,j) = fdv(i,j)+fd
              qx = ax*dl+axx*dx(i-1,j)+axy*dy(i-1,j)
              qy = ay*dl+axy*dx(i-1,j)+ayy*dy(i-1,j)
              ftv(i,j) = ftv(i,j)+fd*td(i-1,j) &
                        +qx*tx(i-1,j)+qy*ty(i-1,j)
              fsv(i,j) = fsv(i,j)+fd*sd(i-1,j) &
                        +qx*sx(i-1,j)+qy*sy(i-1,j)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i-1,j)+axxy*dy(i-1,j)
                  qyy = ayy*dl+axyy*dx(i-1,j)+ayyy*dy(i-1,j)
                  qxy = axy*dl+axxy*dx(i-1,j)+axyy*dy(i-1,j)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i-1,j)+qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                    fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i-1,j) &
                                  +(qx *trd(nt,i-1,j) &
                                   +qxx*trx(nt,i-1,j) &
                                   +qxy*try(nt,i-1,j))*agx(nt,i-1,j) &
                                  +(qy *trd(nt,i-1,j) &
                                   +qxy*trx(nt,i-1,j) &
                                   +qyy*try(nt,i-1,j))*agy(nt,i-1,j)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke .or. nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i-1,j) &
                                  +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i-1,j) &
                                  +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                  end do
                end if
              end if
            else
              x2 = -cuc(i,j)-.5
              y2 = -cvc(i,j)-.5
            end if

            if (cuc(i+1,j) < 0.) then

              ! Add contributions from grid cell (i+1,j). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [-1/2,yc1-1/2], [-cuc(i+1,j)-1/2,-cvc(i+1,j)-1/2], and
              ! [-1/2,-1/2].

              yc0 = (ym*cuc(i+1,j)-cvc(i+1,j)*(xm-.5))/(cuc(i+1,j)+xm-.5)
              yc1 = yc0*scp2(i,j)*scp2i(i+1,j)
              x4 = .5
              y4 = yc0-.5

              call triint(scp2(i+1,j), &
                   -.5,yc1-.5,-cuc(i+1,j)-.5,-cvc(i+1,j)-.5,-.5,-.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i+1,j),max(0.,pbv(i,j)-pup(i+1,j)))
              fd = a*dl+ax*dx(i+1,j)+ay*dy(i+1,j)
              fdv(i,j) = fdv(i,j)+fd
              qx = ax*dl+axx*dx(i+1,j)+axy*dy(i+1,j)
              qy = ay*dl+axy*dx(i+1,j)+ayy*dy(i+1,j)
              ftv(i,j) = ftv(i,j)+fd*td(i+1,j) &
                        +qx*tx(i+1,j)+qy*ty(i+1,j)
              fsv(i,j) = fsv(i,j)+fd*sd(i+1,j) &
                        +qx*sx(i+1,j)+qy*sy(i+1,j)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i+1,j)+axxy*dy(i+1,j)
                  qyy = ayy*dl+axyy*dx(i+1,j)+ayyy*dy(i+1,j)
                  qxy = axy*dl+axxy*dx(i+1,j)+axyy*dy(i+1,j)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i+1,j) &
                         +qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                    fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i+1,j) &
                                  +(qx *trd(nt,i+1,j) &
                                   +qxx*trx(nt,i+1,j) &
                                   +qxy*try(nt,i+1,j))*agx(nt,i+1,j) &
                                  +(qy *trd(nt,i+1,j) &
                                   +qxy*trx(nt,i+1,j) &
                                   +qyy*try(nt,i+1,j))*agy(nt,i+1,j)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i+1,j) &
                                  +qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i+1,j) &
                         +qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
                  end do
                end if
              end if
            else
              x4 = -cuc(i+1,j)+.5
              y4 = -cvc(i+1,j)-.5
            end if

            !-- Add contributions from grid cell (i,j). Assuming
            !-- coordinate [0,0] at the cell center, the contributions are
            !-- flux integrals over the pentagon with vertices [-1/2,-1/2],
            !-- [x2,y2], [xm,ym-1/2], [x4,y4], and [1/2,-1/2].


            call penint(scp2(i,j), &
                 -.5,-.5,x2,y2,xm,ym-.5,x4,y4,.5,-.5, &
                 a,ax,ay,axx,ayy,axy, &
                 axxx,ayyy,axxy,axyy)

            dl = min(dp(i,j),max(0.,pbv(i,j)-pup(i,j)))
            fd = a*dl+ax*dx(i,j)+ay*dy(i,j)
            fdv(i,j) = fdv(i,j)+fd
            qx = ax*dl+axx*dx(i,j)+axy*dy(i,j)
            qy = ay*dl+axy*dx(i,j)+ayy*dy(i,j)
            ftv(i,j) = ftv(i,j)+fd*td(i,j) &
                 +qx*tx(i,j)+qy*ty(i,j)
            fsv(i,j) = fsv(i,j)+fd*sd(i,j) &
                 +qx*sx(i,j)+qy*sy(i,j)
            if (use_TRC) then
              if (use_ATRC) then
                qxx = axx*dl+axxx*dx(i,j)+axxy*dy(i,j)
                qyy = ayy*dl+axyy*dx(i,j)+ayyy*dy(i,j)
                qxy = axy*dl+axxy*dx(i,j)+axyy*dy(i,j)
                do nt = 1,natr
                  fdt = fd*trd(nt,i,j) &
                       +qx*trx(nt,i,j)+qy*try(nt,i,j)
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                  fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i,j) &
                                +(qx *trd(nt,i,j) &
                                 +qxx*trx(nt,i,j) &
                                 +qxy*try(nt,i,j))*agx(nt,i,j) &
                                +(qy *trd(nt,i,j) &
                                 +qxy*trx(nt,i,j) &
                                 +qyy*try(nt,i,j))*agy(nt,i,j)
                end do
                do nt = natr+1,ntr-natr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i,j) &
                       +qx*trx(nt,i,j)+qy*try(nt,i,j)
                end do
              else
                do nt = 1,ntr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i,j) &
                                +qx*trx(nt,i,j)+qy*try(nt,i,j)
                end do
              end if
            end if

          end if

          ! v-component of mass, heat and salt flux.
          vflx(i,j) = fdv(i,j)
          vtflx(i,j) = ftv(i,j)
          vsflx(i,j) = fsv(i,j)

        end do
      end do

    end do

    !---------------------------------------------------------------
    ! Update fields
    !---------------------------------------------------------------

    do j = 1-mrg,jj+mrg
      do l = 1,isp(j)
        do i = max(1-mrg,ifp(j,l)),min(ii+mrg,ilp(j,l))
          q = dp(i,j)
          dp(i,j) = q-(fdu(i+1,j)-fdu(i,j) &
               +fdv(i,j+1)-fdv(i,j))*scp2i(i,j)
          temp(i,j) = (q*temp(i,j) &
                      -(ftu(i+1,j)-ftu(i,j) &
                       +ftv(i,j+1)-ftv(i,j))*scp2i(i,j)) &
                     /dp(i,j)
          saln(i,j) = (q*saln(i,j) &
                      -(fsu(i+1,j)-fsu(i,j) &
                       +fsv(i,j+1)-fsv(i,j))*scp2i(i,j)) &
                     /dp(i,j)
          if (use_TRC) then
            if (use_ATRC) then
              do nt = 1,natr
                nat = ntr-natr+nt
                trc(i,j,k,nt)= &
                     max(0.,(q*trc(i,j,k,nt) &
                     -(ftru(nt,i+1,j)-ftru(nt,i,j) &
                      +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j)) &
                     /dp(i,j)-treps)
                trc(i,j,k,nat) = (q*trc(i,j,k,nat) &
                                 -(fagu(nt,i+1,j)-fagu(nt,i,j) &
                                  +fagv(nt,i,j+1)-fagv(nt,i,j))*scp2i(i,j)) &
                                /dp(i,j)
              end do
              do nt = natr+1,ntr-natr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trc(i,j,k,nt) = (q*trc(i,j,k,nt) &
                                -(ftru(nt,i+1,j)-ftru(nt,i,j) &
                                 +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j)) &
                               /dp(i,j)
              end do
            else
              do nt = 1,ntr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trc(i,j,k,nt) = (q*trc(i,j,k,nt) &
                                -(ftru(nt,i+1,j)-ftru(nt,i,j) &
                                 +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j)) &
                               /dp(i,j)
              end do
            end if
          end if
          dp(i,j) = max(0.,dp(i,j)-dpeps)
        end do
      end do
    end do

  end subroutine remap_eitvel

  !---------------------------------------------------------------

  subroutine remap_eitflx(scuy,scvx,scp2i,scp2,pbmin,pbu,pbv,plo, &
                          u,v,umfl,vmfl,dt,mrg,dp,temp,saln, &
                          uflx,vflx,utflx,vtflx,usflx,vsflx,k)

    !---------------------------------------------------------------
    ! Advection of layer pressure thickness and tracers by incremental
    ! remapping.
    !---------------------------------------------------------------

    ! Argument variables:
    !   scuy   - length of cell boundary with u-point as midpoint.
    !   scvx   - length of cell boundary with v-point as midpoint.
    !   scp2i  - inverse of grid cell area.
    !   scp2   - grid cell area.
    !   pbmin  - minimum bottom pressure of a grid cell and its
    !            neighbors.
    !   pbu    - bottom pressure at u-point, defined as
    !            min(pb(i-1,j),pb(i,j)).
    !   pbv    - bottom pressure at v-point, defined as
    !            min(pb(i,j-1),pb(i,j)).
    !   plo    - lower interface pressure of layer pressure thickness.
    !   u      - u-component of velocity.
    !   v      - v-component of velocity.
    !   umfl   - u-component of mass flux to be applied in the advection.
    !   vmfl   - v-component of mass flux to be applied in the advection.
    !   dt     - time step.
    !   dp     - layer pressure thickness.
    !   temp   - temperature.
    !   saln   - salinity.
    !   uflx   - u-component of total mass flux applied.
    !   vflx   - v-component of total mass flux applied.
    !   utflx  - u-component of heat flux.
    !   vtflx  - v-component of heat flux.
    !   usflx  - u-component of salt flux.
    !   vsflx  - v-component of salt flux.
    !   mrg    - margin of halo that must be valid upon return.
    !   k      - layer index

    ! Arguments
    integer, intent(in) :: k
    integer, intent(in) :: mrg
    real,    intent(in) :: dt
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: scuy
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: scvx
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: scp2i
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: scp2
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbmin
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbu
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: pbv
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: plo
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: u
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: v
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: umfl
    real,    intent(in) ,   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vmfl
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: dp
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: temp
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: saln
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: uflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: utflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vtflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: usflx
    real,    intent(out),   dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: vsflx

    ! Local variables.
    real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
         pup,dx,dy,xd,yd,tx,ty,td,sx,sy,sd,cu,cv,cuc,cvc,fdu,fdv,ftu,ftv, &
         fsu,fsv
    real :: dxi,dyi,dpw,dpe,dps,dpn,dpsw,dpse,dpc,dpnw,dpne, &
         dgmx,dfmx,dfmn,q,q1,q2,q3,q4,tgmx,tgmn,tfmx,tfmn, &
         sgmx,sgmn,sfmx,sfmn,dlm,dlp,aa,mflpos,mflneg,ca, &
         xm,ym,xc0,xc1,yc0,yc1,x2,y2,x4,y4, &
         a,ax,ay,axx,ayy,axy,dl,fd,qx,qy
    integer :: i,j,l,iw,ie,js,jn,isw,jsw,ise,jse,inw,jnw,ine,jne,nw
    real, dimension(ntr-natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: trx
    real, dimension(ntr-natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: try
    real, dimension(ntr-natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: trd
    real, dimension(ntr-natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ftru
    real, dimension(ntr-natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ftrv
    real, dimension(natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ag
    real, dimension(natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: agx
    real, dimension(natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: agy
    real, dimension(natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: agd
    real, dimension(natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fagu
    real, dimension(natr,1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: fagv
    real :: xdt,ydt,axxx,ayyy,axxy,axyy,qxx,qyy,qxy,fdt
    integer :: nt,nat

    !---------------------------------------------------------------
    ! General information:
    !   Logical arrangment of variables is as follows: Layer pressure
    !   thickness dp(i,j), as an example of a scalar variable, is the
    !   mean layer pressure thickness of grid cell (i,j). Velocity
    !   component u(i,j) is located at the midpoint of the cell boundary
    !   separating grid cells (i-1,j) and (i,j). Velocity component
    !   v(i,j) is located at the midpoint of the cell boundary
    !   separating grid cells (i,j-1) and (i,j). A corner variable with
    !   index (i,j) is located at the common grid cell corner of grid
    !   cells (i-1,j-1), (i,j-1), (i-1,j), and (i,j).
    !
    !   The divergence of the velocity field is defined as follows:
    !     (u(i+1,j)*scuy(i+1,j)-u(i,j)*scuy(i,j)
    !     +v(i,j+1)*scvy(i,j+1)-v(i,j)*scvx(i,j))*scp2i(i,j)
    !   By construction, the "fluxing areas" used in obtaining fluxes
    !   trough cell boundaries containing u(i,j) and v(i,j), are equal
    !   to u(i,j)*scuy(i,j)*dt and v(i,j)*scvx(i,j)*dt, respectively.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! Add small number to density field and initialize some variables.
    !---------------------------------------------------------------

    do j = 1-mrg-2,jj+mrg+2
      do l = 1,isp(j)
        do i = max(1-mrg-2,ifp(j,l)),min(ii+mrg+2,ilp(j,l))
          dp(i,j) = max(0.,dp(i,j))+dpeps
          pup(i,j) = plo(i,j)-dp(i,j)
        end do
      end do
      do i = 1-mrg-1,ii+mrg+1
        fdu(i,j) = 0.
        fdv(i,j) = 0.
        ftu(i,j) = 0.
        ftv(i,j) = 0.
        fsu(i,j) = 0.
        fsv(i,j) = 0.
        if (use_TRC) then
          if (use_ATRC) then
            do nt = 1,ntr-natr
              if (use_TKE .and. .not. use_TKEADV) then
                if (nt == itrtke.or.nt == itrgls) cycle
              end if
              ftru(nt,i,j) = 0.
              ftrv(nt,i,j) = 0.
            end do
            do nt = 1,natr
              fagu(nt,i,j) = 0.
              fagv(nt,i,j) = 0.
            end do
          else
            do nt = 1,ntr
              if (use_TKE .and. .not. use_TKEADV) then
                if (nt == itrtke.or.nt == itrgls) cycle
              end if
              ftru(nt,i,j) = 0.
              ftrv(nt,i,j) = 0.
            end do
          end if
        end if
        cu(i,j) = 0.
        cv(i,j) = 0.
      end do
    end do
    if (use_TRC .and. use_ATRC) then
      do nt = 1,natr
        nat = ntr-natr+nt
        do j = 1-mrg-2,jj+mrg+2
          do l = 1,isp(j)
            do i = max(1-mrg-2,ifp(j,l)),min(ii+mrg+2,ilp(j,l))
              trc(i,j,k,nt) = max(0.,trc(i,j,k,nt))+treps
              ag(nt,i,j) = trc(i,j,k,nat)/trc(i,j,k,nt)
            end do
          end do
        end do
      end do
    end if

    !---------------------------------------------------------------
    ! Compute limited gradients, center of mass coordinates, and
    ! non-dimensional velocities.
    !---------------------------------------------------------------

    do j = 1-mrg-1,jj+mrg+1

      do l = 1,isp(j)
        do i = max(1-mrg-1,ifp(j,l)),min(ii+mrg+1,ilp(j,l))

          ! Define indices for grid cell neighbors, ensuring that only wet
          ! points are used.
          iw = i-iu(i  ,j)
          ie = i+iu(i+1,j)
          js = j-iv(i,j  )
          jn = j+iv(i,j+1)
          isw = i*(1-ip(iw,js))+iw*ip(iw,js)
          jsw = j*(1-ip(iw,js))+js*ip(iw,js)
          ise = i*(1-ip(ie,js))+ie*ip(ie,js)
          jse = j*(1-ip(ie,js))+js*ip(ie,js)
          inw = i*(1-ip(iw,jn))+iw*ip(iw,jn)
          jnw = j*(1-ip(iw,jn))+jn*ip(iw,jn)
          ine = i*(1-ip(ie,jn))+ie*ip(ie,jn)
          jne = j*(1-ip(ie,jn))+jn*ip(ie,jn)

          dxi = 1./max(1,ie-iw)
          dyi = 1./max(1,jn-js)

          ! Compute limited gradient for layer pressure thickness and
          ! center of mass coordinate.
          dpsw = max(dpeps,min(pbmin(i,j)-pup(isw,jsw),dp(isw,jsw)))
          dps  = max(dpeps,min(pbmin(i,j)-pup(i  ,js ),dp(i  ,js )))
          dpse = max(dpeps,min(pbmin(i,j)-pup(ise,jse),dp(ise,jse)))
          dpw  = max(dpeps,min(pbmin(i,j)-pup(iw ,j  ),dp(iw ,j  )))
          dpc  = max(dpeps,min(pbmin(i,j)-pup(i  ,j  ),dp(i  ,j  )))
          dpe  = max(dpeps,min(pbmin(i,j)-pup(ie ,j  ),dp(ie ,j  )))
          dpnw = max(dpeps,min(pbmin(i,j)-pup(inw,jnw),dp(inw,jnw)))
          dpn  = max(dpeps,min(pbmin(i,j)-pup(i  ,jn ),dp(i  ,jn )))
          dpne = max(dpeps,min(pbmin(i,j)-pup(ine,jne),dp(ine,jne)))
          dx(i,j) = (dpe-dpw)*dxi
          dy(i,j) = (dpn-dps)*dyi
          dgmx = .5*(abs(dx(i,j))+abs(dy(i,j)))
          dfmx = max(0.,max(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
          dfmn = min(0.,min(dpsw,dps,dpse,dpw,dpe,dpnw,dpn,dpne)-dpc)
          if (dfmx > 0..and.dfmn < 0.) then
            q = min(dfmx/max(dfmx,dgmx),dfmn/min(dfmn,-dgmx))
            dx(i,j) = dx(i,j)*q
            dy(i,j) = dy(i,j)*q
            xd(i,j) = dx(i,j)/(12.*dp(i,j))
            yd(i,j) = dy(i,j)/(12.*dp(i,j))
          else
            dx(i,j) = 0.
            dy(i,j) = 0.
            xd(i,j) = 0.
            yd(i,j) = 0.
          end if

          ! Compute limited gradients for temperature, salinity, and
          ! density
          tx(i,j) = (temp(ie,j)-temp(iw,j))*dxi
          ty(i,j) = (temp(i,jn)-temp(i,js))*dyi
          q1 = tx(i,j)*(-.5-xd(i,j))
          q2 = tx(i,j)*( .5-xd(i,j))
          q3 = ty(i,j)*(-.5-yd(i,j))
          q4 = ty(i,j)*( .5-yd(i,j))
          tgmx = max(q1,q2)+max(q3,q4)
          tgmn = min(q1,q2)+min(q3,q4)
          tfmx = max(0.,max(temp(isw,jsw),temp(i  ,js ), &
                            temp(ise,jse),temp(iw ,j  ), &
                            temp(ie ,j  ),temp(inw,jnw), &
                            temp(i  ,jn ),temp(ine,jne)) &
                         -temp(i,j))
          tfmn = min(0.,min(temp(isw,jsw),temp(i  ,js ), &
                            temp(ise,jse),temp(iw ,j  ), &
                            temp(ie ,j  ),temp(inw,jnw), &
                            temp(i  ,jn ),temp(ine,jne)) &
                         -temp(i,j))
          if (tfmx > 0..and.tfmn < 0.) then
            q = min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
            tx(i,j) = tx(i,j)*q
            ty(i,j) = ty(i,j)*q
            td(i,j) = temp(i,j)-tx(i,j)*xd(i,j)-ty(i,j)*yd(i,j)
          else
            tx(i,j) = 0.
            ty(i,j) = 0.
            td(i,j) = temp(i,j)
          end if

          sx(i,j) = (saln(ie,j)-saln(iw,j))*dxi
          sy(i,j) = (saln(i,jn)-saln(i,js))*dyi
          q1 = sx(i,j)*(-.5-xd(i,j))
          q2 = sx(i,j)*( .5-xd(i,j))
          q3 = sy(i,j)*(-.5-yd(i,j))
          q4 = sy(i,j)*( .5-yd(i,j))
          sgmx = max(q1,q2)+max(q3,q4)
          sgmn = min(q1,q2)+min(q3,q4)
          sfmx = max(0.,max(saln(isw,jsw),saln(i  ,js ), &
                            saln(ise,jse),saln(iw ,j  ), &
                            saln(ie ,j  ),saln(inw,jnw), &
                            saln(i  ,jn ),saln(ine,jne)) &
                         -saln(i,j))
          sfmn = min(0.,min(saln(isw,jsw),saln(i  ,js ), &
                            saln(ise,jse),saln(iw ,j  ), &
                            saln(ie ,j  ),saln(inw,jnw), &
                            saln(i  ,jn ),saln(ine,jne)) &
                         -saln(i,j))
          if (sfmx > 0..and.sfmn < 0.) then
            q = min(sfmx/max(sfmx,sgmx),sfmn/min(sfmn,sgmn))
            sx(i,j) = sx(i,j)*q
            sy(i,j) = sy(i,j)*q
            sd(i,j) = saln(i,j)-sx(i,j)*xd(i,j)-sy(i,j)*yd(i,j)
          else
            sx(i,j) = 0.
            sy(i,j) = 0.
            sd(i,j) = saln(i,j)
          end if
          if (use_TRC) then
            if (use_ATRC) then

              ! Compute limited gradient for tracers.
              do nt = 1,ntr-natr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trx(nt,i,j) = (trc(ie,j,k,nt)-trc(iw,j,k,nt))*dxi
                try(nt,i,j) = (trc(i,jn,k,nt)-trc(i,js,k,nt))*dyi
                q1 = trx(nt,i,j)*(-.5-xd(i,j))
                q2 = trx(nt,i,j)*( .5-xd(i,j))
                q3 = try(nt,i,j)*(-.5-yd(i,j))
                q4 = try(nt,i,j)*( .5-yd(i,j))
                tgmx = max(q1,q2)+max(q3,q4)
                tgmn = min(q1,q2)+min(q3,q4)
                tfmx = max(0.,max(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt), &
                                  trc(ise,jse,k,nt),trc(iw ,j  ,k,nt), &
                                  trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt), &
                                  trc(i  ,jn ,k,nt),trc(ine,jne,k,nt)) &
                               -trc(i,j,k,nt))
                tfmn = min(0.,min(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt), &
                                  trc(ise,jse,k,nt),trc(iw ,j  ,k,nt), &
                                  trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt), &
                                  trc(i  ,jn ,k,nt),trc(ine,jne,k,nt)) &
                               -trc(i,j,k,nt))
                if (tfmx > 0..and.tfmn < 0.) then
                  q = min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                  trx(nt,i,j) = trx(nt,i,j)*q
                  try(nt,i,j) = try(nt,i,j)*q
                  trd(nt,i,j) = trc(i,j,k,nt) &
                               -trx(nt,i,j)*xd(i,j)-try(nt,i,j)*yd(i,j)
                else
                  trx(nt,i,j) = 0.
                  try(nt,i,j) = 0.
                  trd(nt,i,j) = trc(i,j,k,nt)
                end if
              end do

              ! Compute limited gradient for age tracers.
              do nt = 1,natr
                nat = ntr-natr+nt
                agx(nt,i,j) = (ag(nt,ie,j)-ag(nt,iw,j))*dxi
                agy(nt,i,j) = (ag(nt,i,jn)-ag(nt,i,js))*dyi
                q = 1./(12.*trc(i,j,k,nt))
                xdt = (12.*xd(i,j)*trd(nt,i,j)+trx(nt,i,j))*q
                ydt = (12.*yd(i,j)*trd(nt,i,j)+try(nt,i,j))*q
                q1 = agx(nt,i,j)*(-.5-xdt)
                q2 = agx(nt,i,j)*( .5-xdt)
                q3 = agy(nt,i,j)*(-.5-ydt)
                q4 = agy(nt,i,j)*( .5-ydt)
                tgmx = max(q1,q2)+max(q3,q4)
                tgmn = min(q1,q2)+min(q3,q4)
                tfmx = max(0.,max(ag(nt,isw,jsw),ag(nt,i  ,js ), &
                                  ag(nt,ise,jse),ag(nt,iw ,j  ), &
                                  ag(nt,ie ,j  ),ag(nt,inw,jnw), &
                                  ag(nt,i  ,jn ),ag(nt,ine,jne)) &
                               -ag(nt,i,j))
                tfmn = min(0.,min(ag(nt,isw,jsw),ag(nt,i  ,js ), &
                                  ag(nt,ise,jse),ag(nt,iw ,j  ), &
                                  ag(nt,ie ,j  ),ag(nt,inw,jnw), &
                                  ag(nt,i  ,jn ),ag(nt,ine,jne)) &
                               -ag(nt,i,j))
                if (tfmx > 0..and.tfmn < 0.) then
                  q = min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                  agx(nt,i,j) = agx(nt,i,j)*q
                  agy(nt,i,j) = agy(nt,i,j)*q
                  agd(nt,i,j) = ag(nt,i,j)-agx(nt,i,j)*xdt-agy(nt,i,j)*ydt
                else
                  agx(nt,i,j) = 0.
                  agy(nt,i,j) = 0.
                  agd(nt,i,j) = ag(nt,i,j)
                end if
              end do
            else

              ! Compute limited gradient for tracers.
              do nt = 1,ntr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trx(nt,i,j) = (trc(ie,j,k,nt)-trc(iw,j,k,nt))*dxi
                try(nt,i,j) = (trc(i,jn,k,nt)-trc(i,js,k,nt))*dyi
                q1 = trx(nt,i,j)*(-.5-xd(i,j))
                q2 = trx(nt,i,j)*( .5-xd(i,j))
                q3 = try(nt,i,j)*(-.5-yd(i,j))
                q4 = try(nt,i,j)*( .5-yd(i,j))
                tgmx = max(q1,q2)+max(q3,q4)
                tgmn = min(q1,q2)+min(q3,q4)
                tfmx = max(0.,max(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt), &
                                  trc(ise,jse,k,nt),trc(iw ,j  ,k,nt), &
                                  trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt), &
                                  trc(i  ,jn ,k,nt),trc(ine,jne,k,nt)) &
                               -trc(i,j,k,nt))
                tfmn = min(0.,min(trc(isw,jsw,k,nt),trc(i  ,js ,k,nt), &
                                  trc(ise,jse,k,nt),trc(iw ,j  ,k,nt), &
                                  trc(ie ,j  ,k,nt),trc(inw,jnw,k,nt), &
                                  trc(i  ,jn ,k,nt),trc(ine,jne,k,nt)) &
                               -trc(i,j,k,nt))
                if (tfmx > 0..and.tfmn < 0.) then
                  q = min(tfmx/max(tfmx,tgmx),tfmn/min(tfmn,tgmn))
                  trx(nt,i,j) = trx(nt,i,j)*q
                  try(nt,i,j) = try(nt,i,j)*q
                  trd(nt,i,j) = trc(i,j,k,nt) &
                               -trx(nt,i,j)*xd(i,j)-try(nt,i,j)*yd(i,j)
                else
                  trx(nt,i,j) = 0.
                  try(nt,i,j) = 0.
                  trd(nt,i,j) = trc(i,j,k,nt)
                end if
              end do
            end if
          end if

        end do
      end do

    end do

    ! Compute non-dimensional velocities.

    do j = 1-mrg-1,jj+mrg+1
      do l = 1,isu(j)
        do i = max(1-mrg,ifu(j,l)),min(ii+mrg+1,ilu(j,l))
          dlm = min(dp(i-1,j),max(0.,pbu(i,j)-pup(i-1,j)))
          dlp = min(dp(i  ,j),max(0.,pbu(i,j)-pup(i  ,j)))
          aa = u(i,j)*dt*scuy(i,j)
          if (dlm < dpeps) then
            mflpos = 0.
          else
            mflpos= umfl(i,j)+(dlm+.5*dx(i-1,j))*aa
          end if
          if (dlp < dpeps) then
            mflneg = 0.
          else
            mflneg = -umfl(i,j)-(dlp-.5*dx(i  ,j))*aa
          end if
          if (mflpos > 0..and.mflpos > mflneg) then
            ca = aa*scp2i(i-1,j)
            q = dlm+.5*dx(i-1,j)*(1.-ca)
            if (abs(2.*dx(i-1,j)*umfl(i,j)*scp2i(i-1,j)) < 1.e-8*q*q) then
              cu(i,j) = ca+umfl(i,j)*scp2i(i-1,j)/q
            else
              !diag
              if (q*q-2.*dx(i-1,j)*umfl(i,j)*scp2i(i-1,j) < 0.) then
                write (lp,*) 'remap: u pos error: ', &
                     q*q-2.*dx(i-1,j)*umfl(i,j)*scp2i(i-1,j), &
                     'dlm:',dlm,'dx',dx(i-1,j),'ca',ca,'umfl',umfl(i,j), &
                     'scp2i',scp2i(i-1,j)
                call flush(lp)
              end if
              !diag
              cu(i,j) = ca &
                   +(q-sqrt(max(0.,q*q-2.*dx(i-1,j)*umfl(i,j)*scp2i(i-1,j)))) &
                   /dx(i-1,j)
            end if
          else if (mflneg > 0.) then
            ca = aa*scp2i(i,j)
            q = dlp-.5*dx(i,j)*(1.+ca)
            if (abs(2.*dx(i,j)*umfl(i,j)*scp2i(i,j)) < 1.e-8*q*q) then
              cu(i,j) = ca+umfl(i,j)*scp2i(i,j)/q
            else
              !diag
              if (q*q-2.*dx(i,j)*umfl(i,j)*scp2i(i,j) < 0.) then
                write (lp,*) 'remap: u neg error: ', &
                     q*q-2.*dx(i,j)*umfl(i,j)*scp2i(i,j), &
                     'dlp:',dlp,'dx',dx(i,j),'ca',ca,'umfl',umfl(i,j), &
                     'scp2i',scp2i(i,j)
                call flush(lp)
              end if
              !diag
              cu(i,j) = ca &
                   +(q-sqrt(max(0.,q*q-2.*dx(i,j)*umfl(i,j)*scp2i(i,j))))/dx(i,j)
            end if
          else
            cu(i,j) = 0.
          end if
          cu(i,j) = max(-.25,min(.25,cu(i,j)))
        end do
      end do
    end do

    do j = 1-mrg,jj+mrg+1
      do l = 1,isv(j)
        do i = max(1-mrg-1,ifv(j,l)),min(ii+mrg+1,ilv(j,l))
          dlm = min(dp(i,j-1),max(0.,pbv(i,j)-pup(i,j-1)))
          dlp = min(dp(i,j  ),max(0.,pbv(i,j)-pup(i,j  )))
          aa = v(i,j)*dt*scvx(i,j)
          if (dlm < dpeps) then
            mflpos = 0.
          else
            mflpos= vmfl(i,j)+(dlm+.5*dy(i,j-1))*aa
          end if
          if (dlp < dpeps) then
            mflneg = 0.
          else
            mflneg = -vmfl(i,j)-(dlp-.5*dy(i,j  ))*aa
          end if
          if (mflpos > 0..and.mflpos > mflneg) then
            ca = aa*scp2i(i,j-1)
            q = dlm+.5*dy(i,j-1)*(1.-ca)
            if (abs(2.*dy(i,j-1)*vmfl(i,j)*scp2i(i,j-1)) < 1.e-8*q*q) then
              cv(i,j) = ca+vmfl(i,j)*scp2i(i,j-1)/q
            else
              !diag
              if (q*q-2.*dy(i,j-1)*vmfl(i,j)*scp2i(i,j-1) < 0.) then
                write (lp,*) 'remap: v pos error: ', &
                     q*q-2.*dy(i,j-1)*vmfl(i,j)*scp2i(i,j-1), &
                     'dlm:',dlm,'dy',dy(i,j-1),'ca',ca,'vmfl',vmfl(i,j), &
                     'scp2i',scp2i(i,j-1)
                call flush(lp)
              end if
              !diag
              cv(i,j) = ca &
                   +(q-sqrt(max(0.,q*q-2.*dy(i,j-1)*vmfl(i,j)*scp2i(i,j-1)))) &
                   /dy(i,j-1)
            end if
          else if (mflneg > 0.) then
            ca = aa*scp2i(i,j  )
            q = dlp-.5*dy(i,j  )*(1.+ca)
            if (abs(2.*dy(i,j  )*vmfl(i,j)*scp2i(i,j  )) < 1.e-8*q*q) then
              cv(i,j) = ca+vmfl(i,j)*scp2i(i,j  )/q
            else
              !diag
              if (q*q-2.*dy(i,j  )*vmfl(i,j)*scp2i(i,j  ) < 0.) then
                write (lp,*) 'remap: v neg error: ', &
                     q*q-2.*dy(i,j  )*vmfl(i,j)*scp2i(i,j  ), &
                     'dlp:',dlp,'dy',dy(i,j  ),'ca',ca,'vmfl',vmfl(i,j), &
                     'scp2i',scp2i(i,j  )
                call flush(lp)
              end if
              !diag
              cv(i,j) = ca &
                   +(q-sqrt(max(0.,q*q-2.*dy(i,j  )*vmfl(i,j)*scp2i(i,j  )))) &
                   /dy(i,j  )
            end if
          else
            cv(i,j) = 0.
          end if
          cv(i,j) = max(-.25,min(.25,cv(i,j)))
        end do
      end do
    end do

    !---------------------------------------------------------------
    ! Compute corner velocities. The velocity components are computed as
    ! the harmonic mean of the nearest C-grid velocity components with
    ! the following exeptions: The corner velocity component is set to
    ! zero if the nearest C-grid components have different sign, or one
    ! or tree of the neighboring grid cells are wet, or two neighbors
    ! are wet and are arranged diagonally. This construction of corner
    ! velocities will ensure that the entire fluxing area is located
    ! upwind of the cell boundary.
    !---------------------------------------------------------------

    do j = 1-mrg,jj+mrg+1
      do i = 1-mrg,ii+mrg+1
        nw = ip(i-1,j-1)+ip(i,j-1)+ip(i-1,j)+ip(i,j)
        if     (nw == 4) then
          if (cu(i,j-1)*cu(i,j) <= 0.) then
            cuc(i,j) = 0.
          else
            cuc(i,j) = 2.*cu(i,j-1)*cu(i,j)/(cu(i,j-1)+cu(i,j))
          end if
          if (cv(i-1,j)*cv(i,j) <= 0.) then
            cvc(i,j) = 0.
          else
            cvc(i,j) = 2.*cv(i-1,j)*cv(i,j)/(cv(i-1,j)+cv(i,j))
          end if
        else if (nw == 2) then
          if     (ip(i-1,j-1)+ip(i,j-1) == 2) then
            cuc(i,j) = cu(i,j-1)
            cvc(i,j) = 0.
          else if (ip(i-1,j  )+ip(i,j  ) == 2) then
            cuc(i,j) = cu(i,j  )
            cvc(i,j) = 0.
          else if (ip(i-1,j-1)+ip(i-1,j  ) == 2) then
            cuc(i,j) = 0.
            cvc(i,j) = cv(i-1,j)
          else if (ip(i,j-1)+ip(i,j  ) == 2) then
            cuc(i,j) = 0.
            cvc(i,j) = cv(i,j)
          else
            cuc(i,j) = 0.
            cvc(i,j) = 0.
          end if
        else
          cuc(i,j) = 0.
          cvc(i,j) = 0.
        end if
      end do
    end do

    !---------------------------------------------------------------
    ! Compute cell boundary fluxes.
    !---------------------------------------------------------------

    ! u-components of fluxes.

    do j = 1-mrg,jj+mrg

      do l = 1,isu(j)
        do i = max(1-mrg,ifu(j,l)),min(ii+mrg+1,ilu(j,l))

          ! Assuming coordinate [0,0] at the u-point, the non-dimensional
          ! fluxing area is defined as the area of a polygon with vertices
          ! [0,1/2], [-cuc(i,j+1),-cvc(i,j+1)+1/2], [xm,ym],
          ! [-cuc(i,j),-cvc(i,j)-1/2], and [0,-1/2]. The vertex [xm,ym] is
          ! defined so that the polygon area is equal to cu(i,j).

          ym = -.5*(cvc(i,j)+cvc(i,j+1))
          xm = ((ym+.5)*cuc(i,j)-(ym-.5)*cuc(i,j+1)-2.*cu(i,j)) &
               /(1.+cvc(i,j)-cvc(i,j+1))

          if (cu(i,j) > 0.) then

            if (cvc(i,j) > 0.) then

              ! Add contributions from grid cell (i-1,j-1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [xc1+1/2,1/2], [-cuc(i,j)+1/2,-cvc(i,j)+1/2], and
              ! [1/2,1/2].

              xc0 = (xm*cvc(i,j)-cuc(i,j)*(ym+.5))/(cvc(i,j)+ym+.5)
              xc1 = xc0*scp2(i-1,j)*scp2i(i-1,j-1)
              x4 = xc0+.5
              y4 = -.5

              call triint(scp2(i-1,j-1), &
                   xc1+.5,.5,-cuc(i,j)+.5,-cvc(i,j)+.5,.5,.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i-1,j-1),max(0.,pbu(i,j)-pup(i-1,j-1)))
              fd = a*dl+ax*dx(i-1,j-1)+ay*dy(i-1,j-1)
              fdu(i,j) = fdu(i,j)+fd
              qx = ax*dl+axx*dx(i-1,j-1)+axy*dy(i-1,j-1)
              qy = ay*dl+axy*dx(i-1,j-1)+ayy*dy(i-1,j-1)
              ftu(i,j) = ftu(i,j)+fd*td(i-1,j-1) &
                   +qx*tx(i-1,j-1)+qy*ty(i-1,j-1)
              fsu(i,j) = fsu(i,j)+fd*sd(i-1,j-1) &
                   +qx*sx(i-1,j-1)+qy*sy(i-1,j-1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i-1,j-1)+axxy*dy(i-1,j-1)
                  qyy = ayy*dl+axyy*dx(i-1,j-1)+ayyy*dy(i-1,j-1)
                  qxy = axy*dl+axxy*dx(i-1,j-1)+axyy*dy(i-1,j-1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i-1,j-1) &
                         +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                    ftru(nt,i,j) = ftru(nt,i,j)+fdt
                    fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i-1,j-1) &
                         +(qx *trd(nt,i-1,j-1) &
                         +qxx*trx(nt,i-1,j-1) &
                         +qxy*try(nt,i-1,j-1))*agx(nt,i-1,j-1) &
                         +(qy*trd(nt,i-1,j-1) &
                         +qxy*trx(nt,i-1,j-1) &
                         +qyy*try(nt,i-1,j-1))*agy(nt,i-1,j-1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j-1) &
                         +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j-1) &
                         +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                  end do
                end if
              end if
            else
              x4 = -cuc(i,j)+.5
              y4 = -cvc(i,j)-.5
            end if

            if (cvc(i,j+1) < 0.) then

              ! Add contributions from grid cell (i-1,j+1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [xc1+1/2,-1/2], [1/2,-1/2], and
              ! [-cuc(i,j+1)+1/2,-cvc(i,j+1)-1/2].

              xc0 = (xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
              xc1 = xc0*scp2(i-1,j)*scp2i(i-1,j+1)
              x2 = xc0+.5
              y2 = .5

              call triint(scp2(i-1,j+1), &
                   xc1+.5,-.5,.5,-.5,-cuc(i,j+1)+.5,-cvc(i,j+1)-.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i-1,j+1),max(0.,pbu(i,j)-pup(i-1,j+1)))
              fd = a*dl+ax*dx(i-1,j+1)+ay*dy(i-1,j+1)
              fdu(i,j) = fdu(i,j)+fd
              qx = ax*dl+axx*dx(i-1,j+1)+axy*dy(i-1,j+1)
              qy = ay*dl+axy*dx(i-1,j+1)+ayy*dy(i-1,j+1)
              ftu(i,j) = ftu(i,j)+fd*td(i-1,j+1) &
                   +qx*tx(i-1,j+1)+qy*ty(i-1,j+1)
              fsu(i,j) = fsu(i,j)+fd*sd(i-1,j+1) &
                   +qx*sx(i-1,j+1)+qy*sy(i-1,j+1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i-1,j+1)+axxy*dy(i-1,j+1)
                  qyy = ayy*dl+axyy*dx(i-1,j+1)+ayyy*dy(i-1,j+1)
                  qxy = axy*dl+axxy*dx(i-1,j+1)+axyy*dy(i-1,j+1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i-1,j+1) &
                         +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                    ftru(nt,i,j) = ftru(nt,i,j)+fdt
                    fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i-1,j+1) &
                         +(qx *trd(nt,i-1,j+1) &
                         +qxx*trx(nt,i-1,j+1) &
                         +qxy*try(nt,i-1,j+1))*agx(nt,i-1,j+1) &
                         +(qy *trd(nt,i-1,j+1) &
                         +qxy*trx(nt,i-1,j+1) &
                         +qyy*try(nt,i-1,j+1))*agy(nt,i-1,j+1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j+1) &
                         +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j+1) &
                         +qx*trx(nt,i-1,j+1)+qy*try(nt,i-1,j+1)
                  end do
                end if
              end if
            else
              x2 = -cuc(i,j+1)+.5
              y2 = -cvc(i,j+1)+.5
            end if

            !-- Add contributions from grid cell (i-1,j). Assuming
            !-- coordinate [0,0] at the cell center, the contributions are
            !-- flux integrals over the pentagon with vertices [1/2,1/2],
            !-- [x2,y2], [xm+1/2,ym], [x4,y4], and [1/2,-1/2].

            call penint(scp2(i-1,j), &
                 .5,.5,x2,y2,xm+.5,ym,x4,y4,.5,-.5, &
                 a,ax,ay,axx,ayy,axy, &
                 axxx,ayyy,axxy,axyy)

            dl = min(dp(i-1,j),max(0.,pbu(i,j)-pup(i-1,j)))
            fd = a*dl+ax*dx(i-1,j)+ay*dy(i-1,j)
            fdu(i,j) = fdu(i,j)+fd
            qx = ax*dl+axx*dx(i-1,j)+axy*dy(i-1,j)
            qy = ay*dl+axy*dx(i-1,j)+ayy*dy(i-1,j)
            ftu(i,j) = ftu(i,j)+fd*td(i-1,j) &
                 +qx*tx(i-1,j)+qy*ty(i-1,j)
            fsu(i,j) = fsu(i,j)+fd*sd(i-1,j) &
                 +qx*sx(i-1,j)+qy*sy(i-1,j)
            if (use_TRC) then
              if (use_ATRC) then
                qxx = axx*dl+axxx*dx(i-1,j)+axxy*dy(i-1,j)
                qyy = ayy*dl+axyy*dx(i-1,j)+ayyy*dy(i-1,j)
                qxy = axy*dl+axxy*dx(i-1,j)+axyy*dy(i-1,j)
                do nt = 1,natr
                  fdt = fd*trd(nt,i-1,j) &
                       +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                  ftru(nt,i,j) = ftru(nt,i,j)+fdt
                  fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i-1,j) &
                       +(qx *trd(nt,i-1,j) &
                       +qxx*trx(nt,i-1,j) &
                       +qxy*try(nt,i-1,j))*agx(nt,i-1,j) &
                       +(qy *trd(nt,i-1,j) &
                       +qxy*trx(nt,i-1,j) &
                       +qyy*try(nt,i-1,j))*agy(nt,i-1,j)
                end do
                do nt = natr+1,ntr-natr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j) &
                       +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                end do
              else
                do nt = 1,ntr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i-1,j) &
                       +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                end do
              end if
            end if

          else

            if (cvc(i,j) > 0.) then

              ! Add contributions from grid cell (i,j-1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [xc1-1/2,1/2], [-cuc(i,j)-1/2,-cvc(i,j)+1/2], and
              ! [-1/2,1/2].

              xc0 = (xm*cvc(i,j)-cuc(i,j)*(ym+.5))/(cvc(i,j)+ym+.5)
              xc1 = xc0*scp2(i,j)*scp2i(i,j-1)
              x4 = xc0-.5
              y4 = -.5

              call triint(scp2(i,j-1), &
                   xc1-.5,.5,-cuc(i,j)-.5,-cvc(i,j)+.5,-.5,.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i,j-1),max(0.,pbu(i,j)-pup(i,j-1)))
              fd = a*dl+ax*dx(i,j-1)+ay*dy(i,j-1)
              fdu(i,j) = fdu(i,j)+fd
              qx = ax*dl+axx*dx(i,j-1)+axy*dy(i,j-1)
              qy = ay*dl+axy*dx(i,j-1)+ayy*dy(i,j-1)
              ftu(i,j) = ftu(i,j)+fd*td(i,j-1)+qx*tx(i,j-1)+qy*ty(i,j-1)
              fsu(i,j) = fsu(i,j)+fd*sd(i,j-1)+qx*sx(i,j-1)+qy*sy(i,j-1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i,j-1)+axxy*dy(i,j-1)
                  qyy = ayy*dl+axyy*dx(i,j-1)+ayyy*dy(i,j-1)
                  qxy = axy*dl+axxy*dx(i,j-1)+axyy*dy(i,j-1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i,j-1) &
                         +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                    ftru(nt,i,j) = ftru(nt,i,j)+fdt
                    fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i,j-1) &
                                  +(qx *trd(nt,i,j-1) &
                                   +qxx*trx(nt,i,j-1) &
                                   +qxy*try(nt,i,j-1))*agx(nt,i,j-1) &
                                  +(qy *trd(nt,i,j-1) &
                                   +qxy*trx(nt,i,j-1) &
                                   +qyy*try(nt,i,j-1))*agy(nt,i,j-1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j-1) &
                         +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j-1) &
                         +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                  end do
                end if
              end if
            else
              x4 = -cuc(i,j)-.5
              y4 = -cvc(i,j)-.5
            end if

            if (cvc(i,j+1) < 0.) then

              ! Add contributions from grid cell (i,j+1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [xc1-1/2,-1/2], [-1/2,-1/2], and
              ! [-cuc(i,j+1)-1/2,-cvc(i,j+1)-1/2].

              xc0 = (xm*cvc(i,j+1)-cuc(i,j+1)*(ym-.5))/(cvc(i,j+1)+ym-.5)
              xc1 = xc0*scp2(i,j)*scp2i(i,j+1)
              x2 = xc0-.5
              y2 = .5

              call triint(scp2(i,j+1), &
                   xc1-.5,-.5,-.5,-.5,-cuc(i,j+1)-.5,-cvc(i,j+1)-.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i,j+1),max(0.,pbu(i,j)-pup(i,j+1)))
              fd = a*dl+ax*dx(i,j+1)+ay*dy(i,j+1)
              fdu(i,j) = fdu(i,j)+fd
              qx = ax*dl+axx*dx(i,j+1)+axy*dy(i,j+1)
              qy = ay*dl+axy*dx(i,j+1)+ayy*dy(i,j+1)
              ftu(i,j) = ftu(i,j)+fd*td(i,j+1)+qx*tx(i,j+1)+qy*ty(i,j+1)
              fsu(i,j) = fsu(i,j)+fd*sd(i,j+1)+qx*sx(i,j+1)+qy*sy(i,j+1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i,j+1)+axxy*dy(i,j+1)
                  qyy = ayy*dl+axyy*dx(i,j+1)+ayyy*dy(i,j+1)
                  qxy = axy*dl+axxy*dx(i,j+1)+axyy*dy(i,j+1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i,j+1)+qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                    ftru(nt,i,j) = ftru(nt,i,j)+fdt
                    fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i,j+1) &
                                  +(qx *trd(nt,i,j+1) &
                                   +qxx*trx(nt,i,j+1) &
                                   +qxy*try(nt,i,j+1))*agx(nt,i,j+1) &
                                  +(qy *trd(nt,i,j+1) &
                                   +qxy*trx(nt,i,j+1) &
                                   +qyy*try(nt,i,j+1))*agy(nt,i,j+1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j+1) &
                                  +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j+1) &
                                  +qx*trx(nt,i,j+1)+qy*try(nt,i,j+1)
                  end do
                end if
              end if
            else
              x2 = -cuc(i,j+1)-.5
              y2 = -cvc(i,j+1)+.5
            end if

            !-- Add contributions from grid cell (i,j). Assuming
            !-- coordinate [0,0] at the cell center, the contributions are
            !-- flux integrals over the pentagon with vertices [-1/2,1/2],
            !-- [x2,y2], [xm-1/2,ym], [x4,y4], and [-1/2,-1/2].

            call penint(scp2(i,j), &
                 -.5,.5,x2,y2,xm-.5,ym,x4,y4,-.5,-.5, &
                 a,ax,ay,axx,ayy,axy, &
                 axxx,ayyy,axxy,axyy)

            dl = min(dp(i,j),max(0.,pbu(i,j)-pup(i,j)))
            fd = a*dl+ax*dx(i,j)+ay*dy(i,j)
            fdu(i,j) = fdu(i,j)+fd
            qx = ax*dl+axx*dx(i,j)+axy*dy(i,j)
            qy = ay*dl+axy*dx(i,j)+ayy*dy(i,j)
            ftu(i,j) = ftu(i,j)+fd*td(i,j)+qx*tx(i,j)+qy*ty(i,j)
            fsu(i,j) = fsu(i,j)+fd*sd(i,j)+qx*sx(i,j)+qy*sy(i,j)
            if (use_TRC) then
              if (use_ATRC) then
                qxx = axx*dl+axxx*dx(i,j)+axxy*dy(i,j)
                qyy = ayy*dl+axyy*dx(i,j)+ayyy*dy(i,j)
                qxy = axy*dl+axxy*dx(i,j)+axyy*dy(i,j)
                do nt = 1,natr
                  fdt = fd*trd(nt,i,j) &
                       +qx*trx(nt,i,j)+qy*try(nt,i,j)
                  ftru(nt,i,j) = ftru(nt,i,j)+fdt
                  fagu(nt,i,j) = fagu(nt,i,j)+fdt*agd(nt,i,j) &
                                +(qx *trd(nt,i,j) &
                                 +qxx*trx(nt,i,j) &
                                 +qxy*try(nt,i,j))*agx(nt,i,j) &
                                +(qy *trd(nt,i,j) &
                                 +qxy*trx(nt,i,j) &
                                 +qyy*try(nt,i,j))*agy(nt,i,j)
                end do
                do nt = natr+1,ntr-natr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j) &
                                +qx*trx(nt,i,j)+qy*try(nt,i,j)
                end do
              else
                do nt = 1,ntr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftru(nt,i,j) = ftru(nt,i,j)+fd*trd(nt,i,j) &
                                +qx*trx(nt,i,j)+qy*try(nt,i,j)
                end do
              end if
            end if

          end if

          ! u-component of mass, heat and salt flux.
           uflx(i,j) = uflx(i,j)+fdu(i,j)
          utflx(i,j) = utflx(i,j)+ftu(i,j)
          usflx(i,j) = usflx(i,j)+fsu(i,j)

        end do
      end do

    end do

    ! v-components of fluxes.

    do j = 1-mrg,jj+mrg+1

      do l = 1,isv(j)
        do i = max(1-mrg,ifv(j,l)),min(ii+mrg,ilv(j,l))

          ! Assuming coordinate [0,0] at the v-point, the non-dimensional
          ! fluxing area is defined as the area of a polygon with vertices
          ! [-1/2,0], [-cuc(i,j)-1/2,-cvc(i,j)], [xm,ym],
          ! [-cuc(i+1,j)+1/2,-cvc(i+1,j)], and [1/2,0]. The vertex [xm,ym]
          ! is defined so that the polygon area is equal to cv(i,j).

          xm = -.5*(cuc(i,j)+cuc(i+1,j))
          ym = ((xm+.5)*cvc(i,j)-(xm-.5)*cvc(i+1,j)-2.*cv(i,j)) &
               /(1.+cuc(i,j)-cuc(i+1,j))

          if (cv(i,j) > 0) then

            if (cuc(i,j) > 0.) then

              ! Add contributions from grid cell (i-1,j-1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [1/2,yc1+1/2], [1/2,1/2], and
              ! [-cuc(i,j)+1/2,-cvc(i,j)+1/2].

              yc0 = (ym*cuc(i,j)-cvc(i,j)*(xm+.5))/(cuc(i,j)+xm+.5)
              yc1 = yc0*scp2(i,j-1)*scp2i(i-1,j-1)
              x2 = -.5
              y2 = yc0+.5

              call triint(scp2(i-1,j-1), &
                   .5,yc1+.5,.5,.5,-cuc(i,j)+.5,-cvc(i,j)+.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i-1,j-1),max(0.,pbv(i,j)-pup(i-1,j-1)))
              fd = a*dl+ax*dx(i-1,j-1)+ay*dy(i-1,j-1)
              fdv(i,j) = fdv(i,j)+fd
              qx = ax*dl+axx*dx(i-1,j-1)+axy*dy(i-1,j-1)
              qy = ay*dl+axy*dx(i-1,j-1)+ayy*dy(i-1,j-1)
              ftv(i,j) = ftv(i,j)+fd*td(i-1,j-1) &
                   +qx*tx(i-1,j-1)+qy*ty(i-1,j-1)
              fsv(i,j) = fsv(i,j)+fd*sd(i-1,j-1) &
                   +qx*sx(i-1,j-1)+qy*sy(i-1,j-1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i-1,j-1)+axxy*dy(i-1,j-1)
                  qyy = ayy*dl+axyy*dx(i-1,j-1)+ayyy*dy(i-1,j-1)
                  qxy = axy*dl+axxy*dx(i-1,j-1)+axyy*dy(i-1,j-1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i-1,j-1)&
                         +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                    fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i-1,j-1) &
                                  +(qx *trd(nt,i-1,j-1) &
                                   +qxx*trx(nt,i-1,j-1) &
                                   +qxy*try(nt,i-1,j-1))*agx(nt,i-1,j-1) &
                                  +(qy *trd(nt,i-1,j-1) &
                                   +qxy*trx(nt,i-1,j-1) &
                                   +qyy*try(nt,i-1,j-1))*agy(nt,i-1,j-1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i-1,j-1) &
                                  +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i-1,j-1) &
                                  +qx*trx(nt,i-1,j-1)+qy*try(nt,i-1,j-1)
                  end do
                end if
              end if
            else
              x2 = -cuc(i,j)-.5
              y2 = -cvc(i,j)+.5
            end if

            if (cuc(i+1,j) < 0.) then

              ! Add contributions from grid cell (i+1,j-1). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [-1/2,yc1+1/2], [-cuc(i+1,j)-1/2,-cvc(i+1,j)+1/2], and
              ! [-1/2,1/2].

              yc0 = (ym*cuc(i+1,j)-cvc(i+1,j)*(xm-.5))/(cuc(i+1,j)+xm-.5)
              yc1 = yc0*scp2(i,j-1)*scp2i(i+1,j-1)
              x4 = .5
              y4 = yc0+.5

              call triint(scp2(i+1,j-1), &
                   -.5,yc1+.5,-cuc(i+1,j)-.5,-cvc(i+1,j)+.5,-.5,.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i+1,j-1),max(0.,pbv(i,j)-pup(i+1,j-1)))
              fd = a*dl+ax*dx(i+1,j-1)+ay*dy(i+1,j-1)
              fdv(i,j) = fdv(i,j)+fd
              qx = ax*dl+axx*dx(i+1,j-1)+axy*dy(i+1,j-1)
              qy = ay*dl+axy*dx(i+1,j-1)+ayy*dy(i+1,j-1)
              ftv(i,j) = ftv(i,j)+fd*td(i+1,j-1) &
                        +qx*tx(i+1,j-1)+qy*ty(i+1,j-1)
              fsv(i,j) = fsv(i,j)+fd*sd(i+1,j-1) &
                        +qx*sx(i+1,j-1)+qy*sy(i+1,j-1)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i+1,j-1)+axxy*dy(i+1,j-1)
                  qyy = ayy*dl+axyy*dx(i+1,j-1)+ayyy*dy(i+1,j-1)
                  qxy = axy*dl+axxy*dx(i+1,j-1)+axyy*dy(i+1,j-1)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i+1,j-1) &
                         +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                    fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i+1,j-1) &
                                  +(qx *trd(nt,i+1,j-1) &
                                   +qxx*trx(nt,i+1,j-1) &
                                   +qxy*try(nt,i+1,j-1))*agx(nt,i+1,j-1) &
                                  +(qy *trd(nt,i+1,j-1) &
                                   +qxy*trx(nt,i+1,j-1) &
                                   +qyy*try(nt,i+1,j-1))*agy(nt,i+1,j-1)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i+1,j-1) &
                                  +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i+1,j-1) &
                                  +qx*trx(nt,i+1,j-1)+qy*try(nt,i+1,j-1)
                  end do
                end if
              end if
            else
              x4 = -cuc(i+1,j)+.5
              y4 = -cvc(i+1,j)+.5
            end if

            !-- Add contributions from grid cell (i,j-1). Assuming
            !-- coordinate [0,0] at the cell center, the contributions are
            !-- flux integrals over the pentagon with vertices [-1/2,1/2],
            !-- [x2,y2], [xm,ym+1/2], [x4,y4], and [1/2,1/2].

            call penint(scp2(i,j-1), &
                 -.5,.5,x2,y2,xm,ym+.5,x4,y4,.5,.5, &
                 a,ax,ay,axx,ayy,axy, &
                 axxx,ayyy,axxy,axyy)

            dl = min(dp(i,j-1),max(0.,pbv(i,j)-pup(i,j-1)))
            fd = a*dl+ax*dx(i,j-1)+ay*dy(i,j-1)
            fdv(i,j) = fdv(i,j)+fd
            qx = ax*dl+axx*dx(i,j-1)+axy*dy(i,j-1)
            qy = ay*dl+axy*dx(i,j-1)+ayy*dy(i,j-1)
            ftv(i,j) = ftv(i,j)+fd*td(i,j-1)+qx*tx(i,j-1)+qy*ty(i,j-1)
            fsv(i,j) = fsv(i,j)+fd*sd(i,j-1)+qx*sx(i,j-1)+qy*sy(i,j-1)
            if (use_TRC) then
              if (use_ATRC) then
                qxx = axx*dl+axxx*dx(i,j-1)+axxy*dy(i,j-1)
                qyy = ayy*dl+axyy*dx(i,j-1)+ayyy*dy(i,j-1)
                qxy = axy*dl+axxy*dx(i,j-1)+axyy*dy(i,j-1)
                do nt = 1,natr
                  fdt = fd*trd(nt,i,j-1)+qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                  fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i,j-1) &
                                +(qx *trd(nt,i,j-1) &
                                 +qxx*trx(nt,i,j-1) &
                                 +qxy*try(nt,i,j-1))*agx(nt,i,j-1) &
                                +(qy *trd(nt,i,j-1) &
                                 +qxy*trx(nt,i,j-1) &
                                 +qyy*try(nt,i,j-1))*agy(nt,i,j-1)
                end do
                do nt = natr+1,ntr-natr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i,j-1) &
                                +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                end do
              else
                do nt = 1,ntr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i,j-1) &
                                +qx*trx(nt,i,j-1)+qy*try(nt,i,j-1)
                end do
              end if
            end if

          else

            if (cuc(i,j) > 0.) then

              ! Add contributions from grid cell (i-1,j). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [1/2,yc1-1/2], [1/2,-1/2], and
              ! [-cuc(i,j)+1/2,-cvc(i,j)-1/2].

              yc0 = (ym*cuc(i,j)-cvc(i,j)*(xm+.5))/(cuc(i,j)+xm+.5)
              yc1 = yc0*scp2(i,j)*scp2i(i-1,j)
              x2 = -.5
              y2 = yc0-.5

              call triint(scp2(i-1,j), &
                   .5,yc1-.5,.5,-.5,-cuc(i,j)+.5,-cvc(i,j)-.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i-1,j),max(0.,pbv(i,j)-pup(i-1,j)))
              fd = a*dl+ax*dx(i-1,j)+ay*dy(i-1,j)
              fdv(i,j) = fdv(i,j)+fd
              qx = ax*dl+axx*dx(i-1,j)+axy*dy(i-1,j)
              qy = ay*dl+axy*dx(i-1,j)+ayy*dy(i-1,j)
              ftv(i,j) = ftv(i,j)+fd*td(i-1,j) &
                        +qx*tx(i-1,j)+qy*ty(i-1,j)
              fsv(i,j) = fsv(i,j)+fd*sd(i-1,j) &
                        +qx*sx(i-1,j)+qy*sy(i-1,j)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i-1,j)+axxy*dy(i-1,j)
                  qyy = ayy*dl+axyy*dx(i-1,j)+ayyy*dy(i-1,j)
                  qxy = axy*dl+axxy*dx(i-1,j)+axyy*dy(i-1,j)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i-1,j)+qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                    fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i-1,j) &
                                  +(qx *trd(nt,i-1,j) &
                                   +qxx*trx(nt,i-1,j) &
                                   +qxy*try(nt,i-1,j))*agx(nt,i-1,j) &
                                  +(qy *trd(nt,i-1,j) &
                                   +qxy*trx(nt,i-1,j) &
                                   +qyy*try(nt,i-1,j))*agy(nt,i-1,j)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i-1,j) &
                                  +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i-1,j) &
                                  +qx*trx(nt,i-1,j)+qy*try(nt,i-1,j)
                  end do
                end if
              end if
            else
              x2 = -cuc(i,j)-.5
              y2 = -cvc(i,j)-.5
            end if

            if (cuc(i+1,j) < 0.) then

              ! Add contributions from grid cell (i+1,j). Assuming
              ! coordinate [0,0] at the cell center, the contributions are
              ! flux integrals over the triangle with vertices
              ! [-1/2,yc1-1/2], [-cuc(i+1,j)-1/2,-cvc(i+1,j)-1/2], and
              ! [-1/2,-1/2].

              yc0 = (ym*cuc(i+1,j)-cvc(i+1,j)*(xm-.5))/(cuc(i+1,j)+xm-.5)
              yc1 = yc0*scp2(i,j)*scp2i(i+1,j)
              x4 = .5
              y4 = yc0-.5

              call triint(scp2(i+1,j), &
                   -.5,yc1-.5,-cuc(i+1,j)-.5,-cvc(i+1,j)-.5,-.5,-.5, &
                   a,ax,ay,axx,ayy,axy, &
                   axxx,ayyy,axxy,axyy)

              dl = min(dp(i+1,j),max(0.,pbv(i,j)-pup(i+1,j)))
              fd = a*dl+ax*dx(i+1,j)+ay*dy(i+1,j)
              fdv(i,j) = fdv(i,j)+fd
              qx = ax*dl+axx*dx(i+1,j)+axy*dy(i+1,j)
              qy = ay*dl+axy*dx(i+1,j)+ayy*dy(i+1,j)
              ftv(i,j) = ftv(i,j)+fd*td(i+1,j) &
                        +qx*tx(i+1,j)+qy*ty(i+1,j)
              fsv(i,j) = fsv(i,j)+fd*sd(i+1,j) &
                        +qx*sx(i+1,j)+qy*sy(i+1,j)
              if (use_TRC) then
                if (use_ATRC) then
                  qxx = axx*dl+axxx*dx(i+1,j)+axxy*dy(i+1,j)
                  qyy = ayy*dl+axyy*dx(i+1,j)+ayyy*dy(i+1,j)
                  qxy = axy*dl+axxy*dx(i+1,j)+axyy*dy(i+1,j)
                  do nt = 1,natr
                    fdt = fd*trd(nt,i+1,j)+qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                    fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i+1,j) &
                                  +(qx *trd(nt,i+1,j) &
                                   +qxx*trx(nt,i+1,j) &
                                   +qxy*try(nt,i+1,j))*agx(nt,i+1,j) &
                                  +(qy *trd(nt,i+1,j) &
                                   +qxy*trx(nt,i+1,j) &
                                   +qyy*try(nt,i+1,j))*agy(nt,i+1,j)
                  end do
                  do nt = natr+1,ntr-natr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i+1,j) &
                                  +qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
                  end do
                else
                  do nt = 1,ntr
                    if (use_TKE .and. .not. use_TKEADV) then
                      if (nt == itrtke.or.nt == itrgls) cycle
                    end if
                    ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i+1,j) &
                                  +qx*trx(nt,i+1,j)+qy*try(nt,i+1,j)
                  end do
                end if
              end if
            else
              x4 = -cuc(i+1,j)+.5
              y4 = -cvc(i+1,j)-.5
            end if

            !-- Add contributions from grid cell (i,j). Assuming
            !-- coordinate [0,0] at the cell center, the contributions are
            !-- flux integrals over the pentagon with vertices [-1/2,-1/2],
            !-- [x2,y2], [xm,ym-1/2], [x4,y4], and [1/2,-1/2].

            call penint(scp2(i,j), &
                 -.5,-.5,x2,y2,xm,ym-.5,x4,y4,.5,-.5, &
                 a,ax,ay,axx,ayy,axy, &
                 axxx,ayyy,axxy,axyy)

            dl = min(dp(i,j),max(0.,pbv(i,j)-pup(i,j)))
            fd = a*dl+ax*dx(i,j)+ay*dy(i,j)
            fdv(i,j) = fdv(i,j)+fd
            qx = ax*dl+axx*dx(i,j)+axy*dy(i,j)
            qy = ay*dl+axy*dx(i,j)+ayy*dy(i,j)
            ftv(i,j) = ftv(i,j)+fd*td(i,j) &
                    +qx*tx(i,j)+qy*ty(i,j)
            fsv(i,j) = fsv(i,j)+fd*sd(i,j) &
                    +qx*sx(i,j)+qy*sy(i,j)
            if (use_TRC) then
              if (use_ATRC) then
                qxx = axx*dl+axxx*dx(i,j)+axxy*dy(i,j)
                qyy = ayy*dl+axyy*dx(i,j)+ayyy*dy(i,j)
                qxy = axy*dl+axxy*dx(i,j)+axyy*dy(i,j)
                do nt = 1,natr
                  fdt = fd*trd(nt,i,j)+qx*trx(nt,i,j)+qy*try(nt,i,j)
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fdt
                  fagv(nt,i,j) = fagv(nt,i,j)+fdt*agd(nt,i,j) &
                                +(qx *trd(nt,i,j) &
                                 +qxx*trx(nt,i,j) &
                                 +qxy*try(nt,i,j))*agx(nt,i,j) &
                                +(qy *trd(nt,i,j) &
                                 +qxy*trx(nt,i,j) &
                                 +qyy*try(nt,i,j))*agy(nt,i,j)
                end do
                do nt = natr+1,ntr-natr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i,j) &
                                +qx*trx(nt,i,j)+qy*try(nt,i,j)
                end do
              else
                do nt = 1,ntr
                  if (use_TKE .and. .not. use_TKEADV) then
                    if (nt == itrtke.or.nt == itrgls) cycle
                  end if
                  ftrv(nt,i,j) = ftrv(nt,i,j)+fd*trd(nt,i,j) &
                                +qx*trx(nt,i,j)+qy*try(nt,i,j)
                end do
              end if
            end if

          end if

          ! v-component of mass, heat and salt flux.
          vflx(i,j) = fdv(i,j)
          vtflx(i,j) = ftv(i,j)
          vsflx(i,j) = fsv(i,j)

        end do
      end do

    end do

    !---------------------------------------------------------------
    ! Update fields
    !---------------------------------------------------------------

    do j = 1-mrg,jj+mrg
      do l = 1,isp(j)
        do i = max(1-mrg,ifp(j,l)),min(ii+mrg,ilp(j,l))
          q = dp(i,j)
          dp(i,j) = q-(fdu(i+1,j)-fdu(i,j) &
                      +fdv(i,j+1)-fdv(i,j))*scp2i(i,j)
          temp(i,j) = (q*temp(i,j) &
               -(ftu(i+1,j)-ftu(i,j) &
                +ftv(i,j+1)-ftv(i,j))*scp2i(i,j)) &
               /dp(i,j)
          saln(i,j) = (q*saln(i,j) &
               -(fsu(i+1,j)-fsu(i,j) &
                +fsv(i,j+1)-fsv(i,j))*scp2i(i,j)) &
               /dp(i,j)
          if (use_TRC) then
            if (use_ATRC) then
              do nt = 1,natr
                nat = ntr-natr+nt
                trc(i,j,k,nt)= &
                     max(0.,(q*trc(i,j,k,nt) &
                     -(ftru(nt,i+1,j)-ftru(nt,i,j) &
                      +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j)) &
                     /dp(i,j)-treps)
                trc(i,j,k,nat) = (q*trc(i,j,k,nat) &
                     -(fagu(nt,i+1,j)-fagu(nt,i,j) &
                      +fagv(nt,i,j+1)-fagv(nt,i,j))*scp2i(i,j)) &
                     /dp(i,j)
              end do
              do nt = natr+1,ntr-natr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trc(i,j,k,nt) = (q*trc(i,j,k,nt) &
                     -(ftru(nt,i+1,j)-ftru(nt,i,j) &
                      +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j)) &
                     /dp(i,j)
              end do
            else
              do nt = 1,ntr
                if (use_TKE .and. .not. use_TKEADV) then
                  if (nt == itrtke.or.nt == itrgls) cycle
                end if
                trc(i,j,k,nt) = (q*trc(i,j,k,nt) &
                     -(ftru(nt,i+1,j)-ftru(nt,i,j) &
                      +ftrv(nt,i,j+1)-ftrv(nt,i,j))*scp2i(i,j)) &
                     /dp(i,j)
              end do
            end if
          end if
          dp(i,j) = max(0.,dp(i,j)-dpeps)
        end do
      end do
    end do

  end subroutine remap_eitflx

end module mod_remap
