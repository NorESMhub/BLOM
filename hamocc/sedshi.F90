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


      SUBROUTINE SEDSHI(kpie,kpje,omask)
!**********************************************************************
!
!**** *SEDSHI* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - rename ssssil(i,j,k)=sedlay(i,j,k,issssil) etc.
!     I. Kriest         *MPI-Met, HH*,   27.05.03
!     - change specific weights for opal, CaCO3, POC
!     - include upward transport 
!     Purpose
!     -------
!     .
!
!     Method
!     -------
!     .
!
!**   Interface.
!     ----------
!
!     *CALL*       *SEDSHI*
!
!     Externals
!     ---------
!     none.
!
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod
      use mo_param1_bgc 

      implicit none

      INTEGER :: kpie,kpje,i,j,k,l,iv
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje),omask(kpie,kpje)
      REAL :: wsed(kpie,kpje), fulsed(kpie,kpje)
      REAL :: sedlo,uebers,seddef,spresent,buried
      REAL :: refill,frac

! DOWNWARD SHIFTING	 
! shift solid sediment sediment downwards, if layer is full, i.e., if 
! the volume filled by the four constituents poc, opal, caco3, clay
! is more than porsol*seddw
! the outflow of layer i is given by sedlay(i)*porsol(i)*seddw(i), it is
! distributed in the layer below over a volume of porsol(i+1)*seddw(i+1)

      do k=1,ks-1

!$OMP PARALLEL DO PRIVATE(sedlo) 
        do j=1,kpje
        do i=1,kpie
          if(omask(i,j).gt.0.5) then
!ka          if(bolay(i,j).gt.0.) then
             sedlo = orgfa*rcar*sedlay(i,j,k,issso12)                  &
     &              +calfa*sedlay(i,j,k,isssc12)                       &
     &              +oplfa*sedlay(i,j,k,issssil)                       &
     &              +clafa*sedlay(i,j,k,issster)
! "full sediment has sedlo=1
             wsed(i,j)=max(0.,(sedlo-1.)/(sedlo+1.e-10))
          endif
       enddo !end i-loop
       enddo !end j-loop
!$OMP END PARALLEL DO

! filling downward  (accumulation)
        do iv=1,nsedtra
!$OMP PARALLEL DO PRIVATE(uebers) 
        do j=1,kpje
        do i=1,kpie
          if(omask(i,j).gt.0.5) then
!ka          if(bolay(i,j).gt.0.) then
            uebers=wsed(i,j)*sedlay(i,j,k,iv)
            sedlay(i,j,k  ,iv)=sedlay(i,j,k  ,iv)-uebers
            sedlay(i,j,k+1,iv)=sedlay(i,j,k+1,iv)+uebers               &
     &        *(seddw(k)*porsol(k))/(seddw(k+1)*porsol(k+1))
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

!$OMP PARALLEL DO PRIVATE(sedlo)  
       do j=1,kpje
       do i=1,kpie
          if(omask(i,j).gt.0.5) then
!ka          if(bolay(i,j).gt.0.) then
             sedlo = orgfa*rcar*sedlay(i,j,ks,issso12)                 &
     &              +calfa*sedlay(i,j,ks,isssc12)                      &
     &              +oplfa*sedlay(i,j,ks,issssil)                      &
     &              +clafa*sedlay(i,j,ks,issster)
             wsed(i,j)=max(0.,(sedlo-1.)/(sedlo+1.e-10))
          endif
      enddo !end i-loop
      enddo !end j-loop
!$OMP END PARALLEL DO

      do iv=1,nsedtra
!$OMP PARALLEL DO PRIVATE(uebers)  
      do j=1,kpje
      do i=1,kpie
          if(omask(i,j).gt.0.5) then
!ka          if(bolay(i,j).gt.0.) then
            uebers=wsed(i,j)*sedlay(i,j,k,iv)
            sedlay(i,j,ks ,iv)=sedlay(i,j,ks ,iv)-uebers
            burial(i,j,iv)=burial(i,j,iv)+uebers*seddw(k)*porsol(k)
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

!$OMP PARALLEL DO 
      do j=1,kpje
      do i=1,kpie
        fulsed(i,j)=0.
      enddo !end i-loop
      enddo !end j-loop
!$OMP END PARALLEL DO
      
! determine how the total sediment column is filled 
      do k=1,ks
!$OMP PARALLEL DO PRIVATE(sedlo) 
      do j=1,kpje
      do i=1,kpie
        if(omask(i,j).gt.0.5) then
!ka        if(bolay(i,j).gt.0.) then
          sedlo=orgfa*rcar*sedlay(i,j,k,issso12)                       &
     &         +calfa*sedlay(i,j,k,isssc12)                            &
     &         +oplfa*sedlay(i,j,k,issssil)                            &
     &         +clafa*sedlay(i,j,k,issster)
          fulsed(i,j)=fulsed(i,j)+porsol(k)*seddw(k)*sedlo
        endif
      enddo !end i-loop
      enddo !end j-loop
!$OMP END PARALLEL DO
      enddo !end k-loop

! shift the sediment deficiency from the deepest (burial) 
! layer into layer ks
!$OMP PARALLEL DO                                          &
!$OMP&PRIVATE(seddef,spresent,buried,refill,frac) 
      do j=1,kpje
      do i=1,kpie
      if(omask(i,j).gt.0.5) then
!ka      if(bolay(i,j).gt.0.) then

! deficiency to fully loaded sediment packed in sedlay(i,j,ks)
! this is the volume required from the buried layer

        seddef=solfu-fulsed(i,j)

! total volume of solid constituents in buried layer
        spresent=orgfa*rcar*burial(i,j,issso12)                        &
     &         +calfa*burial(i,j,isssc12)                              &
     &         +oplfa*burial(i,j,issssil)                              &
     &         +clafa*burial(i,j,issster)

! determine whether an additional amount of clay is needed in the burial
! layer to fill the whole sediment; I assume that there is an infinite
! supply of clay from below
        burial(i,j,issster) = burial(i,j,issster)                      &
     &                      + MAX(0.,seddef-spresent)/clafa

! determine new volume of buried layer
        buried=orgfa*rcar*burial(i,j,issso12)                          &
     &        +calfa*burial(i,j,isssc12)                               &
     &        +oplfa*burial(i,j,issssil)                               &
     &        +clafa*burial(i,j,issster)

! fill the last active layer
        refill=seddef/(buried+1.e-10) 
        frac = porsol(ks)*seddw(ks) !changed k to ks, ik
        
        sedlay(i,j,ks,issso12)=sedlay(i,j,ks,issso12)                  &
     &                        +refill*burial(i,j,issso12)/frac
        sedlay(i,j,ks,isssc12)=sedlay(i,j,ks,isssc12)                  &
     &                        +refill*burial(i,j,isssc12)/frac
        sedlay(i,j,ks,issssil)=sedlay(i,j,ks,issssil)                  &
     &                        +refill*burial(i,j,issssil)/frac
        sedlay(i,j,ks,issster)=sedlay(i,j,ks,issster)                  &
     &                        +refill*burial(i,j,issster)/frac

! account for losses in buried sediment
        burial(i,j,issso12) = burial(i,j,issso12)                      &
     &                      - refill*burial(i,j,issso12)
        burial(i,j,isssc12) = burial(i,j,isssc12)                      &
     &                      - refill*burial(i,j,isssc12)
        burial(i,j,issssil) = burial(i,j,issssil)                      &
     &                      - refill*burial(i,j,issssil)
        burial(i,j,issster) = burial(i,j,issster)                      &
     &                      - refill*burial(i,j,issster)

      endif
      enddo !end i-loop
      enddo !end j-loop
!$OMP END PARALLEL DO

!     redistribute overload of layer ks
      do  k=ks,2,-1
!$OMP PARALLEL DO PRIVATE(sedlo) 
      do j=1,kpje
      do i=1,kpie
        if(omask(i,j).gt.0.5) then
!ka        if(bolay(i,j).gt.0.) then
          sedlo=orgfa*rcar*sedlay(i,j,k,issso12)                       &
     &         +calfa*sedlay(i,j,k,isssc12)                            &
     &         +oplfa*sedlay(i,j,k,issssil)                            &
     &         +clafa*sedlay(i,j,k,issster)
          wsed(i,j)=max(0.,(sedlo-1.)/(sedlo+1.e-10))
        endif
      enddo !end i-loop
      enddo !end j-loop
!$OMP END PARALLEL DO

      do iv=1,4
!$OMP PARALLEL DO PRIVATE(uebers,frac) 
      do j=1,kpje
      do i=1,kpie
        if(omask(i,j).gt.0.5) then
!ka        if(bolay(i,j).gt.0.) then
          uebers=sedlay(i,j,k,iv)*wsed(i,j)
          frac=porsol(k)*seddw(k)/(porsol(k-1)*seddw(k-1))
          sedlay(i,j,k,iv)=sedlay(i,j,k,iv)-uebers
          sedlay(i,j,k-1,iv)=sedlay(i,j,k-1,iv)+uebers*frac
#ifdef cisonew
          if(iv.eq.issso12)then
            sedlay(i,j,k,issso13)  =sedlay(i,j,k,issso13)-uebers
            sedlay(i,j,k-1,issso13)=sedlay(i,j,k-1,issso13)+uebers*frac
            sedlay(i,j,k,issso14)  =sedlay(i,j,k,issso14)-uebers
            sedlay(i,j,k-1,issso14)=sedlay(i,j,k-1,issso14)+uebers*frac
          endif
          if(iv.eq.isssc12)then
            sedlay(i,j,k,isssc13)  =sedlay(i,j,k,isssc13)-uebers
            sedlay(i,j,k-1,isssc13)=sedlay(i,j,k-1,isssc13)+uebers*frac
            sedlay(i,j,k,isssc14)  =sedlay(i,j,k,isssc14)-uebers
            sedlay(i,j,k-1,isssc14)=sedlay(i,j,k-1,isssc14)+uebers*frac
          endif
#endif
        endif
      enddo !end i-loop
      enddo !end j-loop
!$OMP END PARALLEL DO
      enddo !end iv-loop

      enddo  !end k-loop


      RETURN
      END
