! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
! Copyright (C) 2003  I. Kriest
! Copyright (C) 2020  K. Assmann, J. Tjiputra, J. Schwinger
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


module mo_sedshi

  implicit none
  private

  public :: sedshi

contains

  subroutine sedshi(kpie,kpje,omask,kplyear)

    !***********************************************************************************************
    ! Sediment shifting
    !
    ! Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    ! Modified:
    ! S.Legutke,        *MPI-MaD, HH*    10.04.01
    !  - rename ssssil(i,j,k)=sedlay(i,j,k,issssil) etc.
    ! I. Kriest         *MPI-Met, HH*,   27.05.03
    !***********************************************************************************************

    use mo_kind,        only: rp
    use mo_sedmnt,      only: burial,calfa,clafa,oplfa,orgfa,porsol,sedlay,seddw,solfu
    use mo_param_bgc,   only: rcar,sec_per_year,sec_per_day
    use mo_param1_bgc,  only: isssc12,issssil,issso12,issster,ks,nsedtra,isssc13,isssc14,          &
                              issso13,issso14,issso12_age,nsedtra_woage
    use mo_carbch,      only: sedfluxb
    use mo_control_bgc, only: use_cisonew,use_sediment_quality,dtbgc,                              &
                            & do_sedspinup,sedspin_yr_s,sedspin_yr_e,sedspin_ncyc,ldyn_sed_age

    ! Arguments
    integer, intent(in) :: kpie
    integer, intent(in) :: kpje
    real(rp),intent(in) :: omask(kpie,kpje)
    integer, intent(in) :: kplyear                                         ! current year.

    ! Local variables
    integer  :: i,j,k,l,iv
    real(rp) :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
    real(rp) :: wsed(kpie,kpje), fulsed(kpie,kpje)
    real(rp) :: sedlo,uebers,seddef,spresent,buried
    real(rp) :: refill,frac
    real(rp) :: eps=epsilon(1._rp)
    real(rp) :: acc_time=0._rp

    sedfluxb(:,:,:) = 0._rp

    if(do_sedspinup .and. kplyear>=sedspin_yr_s .and. kplyear<=sedspin_yr_e) then
      ! accumulated time spent due to sediment acceleration
      acc_time = sec_per_day*sedspin_ncyc/sec_per_year ! *dtbgc/dtbgc
    endif

    ! DOWNWARD SHIFTING
    ! shift solid sediment sediment downwards, if layer is full, i.e., if
    ! the volume filled by the four constituents poc, opal, caco3, clay
    ! is more than porsol*seddw
    ! the outflow of layer i is given by sedlay(i)*porsol(i)*seddw(i), it is
    ! distributed in the layer below over a volume of porsol(i+1)*seddw(i+1)

    do k=1,ks-1

      !$OMP PARALLEL DO PRIVATE(i,sedlo)
      do j=1,kpje
        do i=1,kpie
          if(omask(i,j) > 0.5_rp) then
            !ka          if(bolay(i,j).gt.0._rp) then
            sedlo  = rcar*orgfa*sedlay(i,j,k,issso12) &
                 & +      calfa*sedlay(i,j,k,isssc12) &
                 & +      oplfa*sedlay(i,j,k,issssil) &
                 & +      clafa*sedlay(i,j,k,issster)
            ! "full sediment has sedlo=1
            wsed(i,j)=max(0._rp,(sedlo-1._rp)/(abs(sedlo)+1.e-10_rp))
          endif
        enddo !end i-loop
      enddo !end j-loop
      !$OMP END PARALLEL DO

      ! filling downward  (accumulation)
      do iv=1,nsedtra_woage
        !$OMP PARALLEL DO PRIVATE(i,uebers)
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j) > 0.5_rp) then
              !ka          if(bolay(i,j).gt.0._rp) then
              uebers=wsed(i,j)*sedlay(i,j,k,iv)
              if (use_sediment_quality .and. iv == issso12 .and. ldyn_sed_age) then
                sedlay(i,j,k+1,issso12_age) = ( uebers                                             &
                 & *(seddw(k)*porsol(i,j,k))/(seddw(k+1)*porsol(i,j,k+1))*sedlay(i,j,k,issso12_age)&
                 &     + sedlay(i,j,k,issso12)*sedlay(i,j,k,issso12_age))                          &
                 & / (uebers*(seddw(k)*porsol(i,j,k))/(seddw(k+1)*porsol(i,j,k+1))                 &
                       + sedlay(i,j,k,issso12)+eps)
              endif
              sedlay(i,j,k  ,iv)=sedlay(i,j,k  ,iv)-uebers
              sedlay(i,j,k+1,iv)=sedlay(i,j,k+1,iv)+uebers   &
                               *(seddw(k)*porsol(i,j,k))/(seddw(k+1)*porsol(i,j,k+1))
            endif
          enddo !end i-loop
        enddo !end j-loop
        !$OMP END PARALLEL DO
      enddo !end iv-loop

    enddo !end k-loop

    ! store amount lost from last sediment layer - this is a kind of
    ! permanent burial in deep consolidated layer, and this stuff is
    ! effectively lost from the whole ocean+sediment(+atmosphere) system.
    ! Would have to be supplied by river runoff or simple addition e.g.
    ! to surface layers in the long range. Can be supplied again if a
    ! sediment column has a deficiency in volume.

    !$OMP PARALLEL DO PRIVATE(i,sedlo)
    do j=1,kpje
      do i=1,kpie
        if(omask(i,j) > 0.5_rp) then
          !ka          if(bolay(i,j).gt.0._rp) then
          sedlo  = rcar*orgfa*sedlay(i,j,ks,issso12)  &
               & +      calfa*sedlay(i,j,ks,isssc12)  &
               & +      oplfa*sedlay(i,j,ks,issssil)  &
               & +      clafa*sedlay(i,j,ks,issster)
          wsed(i,j)=max(0._rp,(sedlo-1._rp)/(abs(sedlo)+1.e-10_rp))
        endif
      enddo !end i-loop
    enddo !end j-loop
    !$OMP END PARALLEL DO

    do iv=1,nsedtra_woage
      !$OMP PARALLEL DO PRIVATE(i,uebers)
      do j=1,kpje
        do i=1,kpie
          if(omask(i,j) > 0.5_rp) then
            !ka          if(bolay(i,j).gt.0._rp) then
            uebers=wsed(i,j)*sedlay(i,j,ks,iv)
            if (use_sediment_quality .and. iv == issso12  .and. ldyn_sed_age) then
              burial(i,j,issso12_age) = (uebers*seddw(ks)*porsol(i,j,ks)*sedlay(i,j,ks,issso12_age)&
                                      &   + burial(i,j,issso12)*burial(i,j,issso12_age))           &
                                      & /(uebers*seddw(ks)*porsol(i,j,ks) + burial(i,j,issso12)+eps)
            endif
            sedlay(i,j,ks ,iv)=sedlay(i,j,ks ,iv)-uebers
            burial(i,j,iv)=burial(i,j,iv)+uebers*seddw(ks)*porsol(i,j,ks)
            sedfluxb(i,j,iv) = uebers*seddw(ks)*porsol(i,j,ks)
          endif
        enddo !end i-loop
      enddo !end j-loop
      !$OMP END PARALLEL DO
    enddo !end iv-loop

    ! now the loading nowhere excceds 1

    ! digging from below in case of erosion
    ! UPWARD SHIFTING
    ! shift solid sediment sediment upwards, if total sediment volume is less
    ! than required, i.e., if the volume filled by the four constituents
    ! poc, opal, caco3, claycik (integrated over total sediment column)
    ! is less than porsol*seddw (integrated over total sediment column)
    ! first, the last box is filled from below with total required volume;
    ! then, successively, the following layers are filled upwards.
    ! if there is not enough solid matter to fill the column, add clay.

    !$OMP PARALLEL DO PRIVATE(i)
    do j=1,kpje
      do i=1,kpie
        fulsed(i,j)=0._rp
      enddo !end i-loop
    enddo !end j-loop
    !$OMP END PARALLEL DO

    ! determine how the total sediment column is filled
    do k=1,ks
      !$OMP PARALLEL DO PRIVATE(i,sedlo)
      do j=1,kpje
        do i=1,kpie
          if(omask(i,j) > 0.5_rp) then
            !ka        if(bolay(i,j).gt.0._rp) then
            sedlo  = rcar*orgfa*sedlay(i,j,k,issso12)  &
                 & +      calfa*sedlay(i,j,k,isssc12)  &
                 & +      oplfa*sedlay(i,j,k,issssil)  &
                 & +      clafa*sedlay(i,j,k,issster)
            fulsed(i,j)=fulsed(i,j)+porsol(i,j,k)*seddw(k)*sedlo
          endif
        enddo !end i-loop
      enddo !end j-loop
      !$OMP END PARALLEL DO
    enddo !end k-loop

    ! shift the sediment deficiency from the deepest (burial)
    ! layer into layer ks

    !$OMP PARALLEL DO PRIVATE(i,seddef,spresent,buried,refill,frac)
    do j=1,kpje
      do i=1,kpie
        if(omask(i,j) > 0.5_rp) then
          !ka      if(bolay(i,j).gt.0._rp) then

          ! deficiency to fully loaded sediment packed in sedlay(i,j,ks)
          ! this is the volume required from the buried layer

          seddef=solfu(i,j)-fulsed(i,j)

          ! total volume of solid constituents in buried layer
          spresent = rcar*orgfa*burial(i,j,issso12)  &
               &   +      calfa*burial(i,j,isssc12)  &
               &   +      oplfa*burial(i,j,issssil)  &
               &   +      clafa*burial(i,j,issster)

          ! determine whether an additional amount of clay is needed in the burial
          ! layer to fill the whole sediment; I assume that there is an infinite
          ! supply of clay from below
          burial(i,j,issster) = burial(i,j,issster) + max(0._rp,seddef-spresent)/clafa

          ! determine new volume of buried layer
          buried = rcar*orgfa*burial(i,j,issso12)  &
               & +      calfa*burial(i,j,isssc12)  &
               & +      oplfa*burial(i,j,issssil)  &
               & +      clafa*burial(i,j,issster)

          ! fill the last active layer
          refill=seddef/(buried+1.e-10_rp)
          frac = porsol(i,j,ks)*seddw(ks)

          if (use_sediment_quality  .and. ldyn_sed_age) then
            ! Update burial POC age [yrs] - NOTE that sedshi is called once per day!
            burial(i,j,issso12_age)    = burial(i,j,issso12_age) + sec_per_day/sec_per_year + acc_time
            sedlay(i,j,ks,issso12_age) = (refill*burial(i,j,issso12)/frac * burial(i,j,issso12_age)&
                                       &    + sedlay(i,j,ks,issso12)*sedlay(i,j,ks,issso12_age))   &
                                       & /(refill*burial(i,j,issso12)/frac + sedlay(i,j,ks,issso12)+eps)
          endif

          sedlay(i,j,ks,issso12)=sedlay(i,j,ks,issso12)+refill*burial(i,j,issso12)/frac
          sedlay(i,j,ks,isssc12)=sedlay(i,j,ks,isssc12)+refill*burial(i,j,isssc12)/frac
          sedlay(i,j,ks,issssil)=sedlay(i,j,ks,issssil)+refill*burial(i,j,issssil)/frac
          sedlay(i,j,ks,issster)=sedlay(i,j,ks,issster)+refill*burial(i,j,issster)/frac

          if (use_cisonew) then
            sedlay(i,j,ks,issso13)=sedlay(i,j,ks,issso13)+refill*burial(i,j,issso13)/frac
            sedlay(i,j,ks,isssc13)=sedlay(i,j,ks,isssc13)+refill*burial(i,j,isssc13)/frac
            sedlay(i,j,ks,issso14)=sedlay(i,j,ks,issso14)+refill*burial(i,j,issso14)/frac
            sedlay(i,j,ks,isssc14)=sedlay(i,j,ks,isssc14)+refill*burial(i,j,isssc14)/frac
          endif

          ! account for refluxes to get net-burial fluxes for output:
          sedfluxb(i,j,issso12) = sedfluxb(i,j,issso12) - refill*burial(i,j,issso12)
          sedfluxb(i,j,isssc12) = sedfluxb(i,j,isssc12) - refill*burial(i,j,isssc12)
          sedfluxb(i,j,issssil) = sedfluxb(i,j,issssil) - refill*burial(i,j,issssil)
          sedfluxb(i,j,issster) = sedfluxb(i,j,issster) - refill*burial(i,j,issster)

          ! account for losses in buried sediment
          burial(i,j,issso12) = burial(i,j,issso12)-refill*burial(i,j,issso12)
          burial(i,j,isssc12) = burial(i,j,isssc12)-refill*burial(i,j,isssc12)
          burial(i,j,issssil) = burial(i,j,issssil)-refill*burial(i,j,issssil)
          burial(i,j,issster) = burial(i,j,issster)-refill*burial(i,j,issster)
          if (use_cisonew) then
            burial(i,j,issso13) = burial(i,j,issso13)-refill*burial(i,j,issso13)
            burial(i,j,isssc13) = burial(i,j,isssc13)-refill*burial(i,j,isssc13)
            burial(i,j,issso14) = burial(i,j,issso14)-refill*burial(i,j,issso14)
            burial(i,j,isssc14) = burial(i,j,isssc14)-refill*burial(i,j,isssc14)
          endif
        endif
      enddo !end i-loop
    enddo !end j-loop
    !$OMP END PARALLEL DO

    !     redistribute overload of layer ks
    do  k=ks,2,-1
      !$OMP PARALLEL DO PRIVATE(i,sedlo)
      do j=1,kpje
        do i=1,kpie
          if(omask(i,j) > 0.5_rp) then
            !ka        if(bolay(i,j).gt.0._rp) then
            sedlo  = rcar*orgfa*sedlay(i,j,k,issso12)  &
                 & +      calfa*sedlay(i,j,k,isssc12)  &
                 & +      oplfa*sedlay(i,j,k,issssil)  &
                 & +      clafa*sedlay(i,j,k,issster)
            wsed(i,j)=max(0._rp,(sedlo-1._rp)/(abs(sedlo)+1.e-10_rp))
          endif
        enddo !end i-loop
      enddo !end j-loop
      !$OMP END PARALLEL DO

      do iv=1,nsedtra_woage
        !$OMP PARALLEL DO PRIVATE(i,uebers,frac)
        do j=1,kpje
          do i=1,kpie
            if(omask(i,j) > 0.5_rp) then
              !ka        if(bolay(i,j).gt.0._rp) then
              uebers=sedlay(i,j,k,iv)*wsed(i,j)
              frac=porsol(i,j,k)*seddw(k)/(porsol(i,j,k-1)*seddw(k-1))
              if (use_sediment_quality .and. iv == issso12 .and. ldyn_sed_age) then
                sedlay(i,j,k-1,issso12_age) = (uebers*frac*sedlay(i,j,k,issso12_age)               &
                                            &+ sedlay(i,j,k-1,issso12)*sedlay(i,j,k-1,issso12_age))&
                                            & / (uebers*frac + sedlay(i,j,k-1,issso12)+eps)
              endif
              sedlay(i,j,k,iv)=sedlay(i,j,k,iv)-uebers
              sedlay(i,j,k-1,iv)=sedlay(i,j,k-1,iv)+uebers*frac
            endif
          enddo !end i-loop
        enddo !end j-loop
        !$OMP END PARALLEL DO
      enddo !end iv-loop

    enddo  !end k-loop

  end subroutine sedshi

end module mo_sedshi
