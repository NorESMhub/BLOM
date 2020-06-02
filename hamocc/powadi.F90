! Copyright (C) 2001  Ernst Maier-Reimer, S. Legutke
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


      SUBROUTINE POWADI(j,kpie,kpje,solrat,sedb1,sediso,bolven,omask)
!**********************************************************************
!
!**** *POWADI* - vertical diffusion with simultaneous dissolution.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Purpose
!     -------
!     .
!
!     Method
!     -------
!     implicit discretisation.
!
!**   Interface.
!     ----------
!
!     *CALL*       *POWADI(j,solrat,sedb1,sediso)*
!
!     Input  solrat : dissolution rate
!     =====       j : zonal grid index
!             sedb1 : tracer at entry
!
!     Output: sediso: diffused tracer at exit
!     ======
!
!     *PARAMETER*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
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
      USE mo_control_bgc
      use mo_vgrid, only: bolay

      implicit none

      REAL :: sedb1(kpie,0:ks),sediso(kpie,0:ks)
      REAL :: solrat(kpie,ks)
      REAL :: bolven(kpie)
      REAL :: TREDSY(kpie,0:ks,3)
      REAL :: omask(kpie,kpje)
      INTEGER :: kpie,kpje,i,j,k,l
      REAL :: asu,alo
      
      DO 1321 k=1,ks
         asu=sedict*seddzi(k)*porwah(k)
         alo=0.
         IF(k.LT.ks)alo=sedict*seddzi(k+1)*porwah(k+1)
         DO 1321 i=1,kpie
            tredsy(i,k,1) = -asu
            tredsy(i,k,3) = -alo 
            tredsy(i,k,2) =  seddw(k)*porwat(k) - tredsy(i,k,1)       &
     &          - tredsy(i,k,3) + solrat(i,k)*porwat(k)*seddw(k)
 1321    CONTINUE

         k=0
         asu=0.
         alo=sedict*seddzi(1)*porwah(1)
         DO 1421 i=1,kpie
          IF(omask(i,j).GT.0.5) THEN
              tredsy(i,k,1) = -asu
              tredsy(i,k,3) = -alo 
              tredsy(i,k,2) = bolven(i)*bolay(i,j)                   &
     &                        -tredsy(i,k,1)-tredsy(i,k,3)
          ELSE
              tredsy(i,k,1) = 0
              tredsy(i,k,3) = 0
              tredsy(i,k,2) = 0
          ENDIF
	  
1421  CONTINUE

       DO 132 K=1,KS
          DO 133 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN
                tredsy(i,k-1,1) = tredsy(i,k,1) / tredsy(i,k-1,2)
                tredsy(i,k,2)   = tredsy(i,k,2)                       &
     &              - tredsy(i,k-1,3) * tredsy(i,k,1) / tredsy(i,k-1,2)
             ENDIF
 133      CONTINUE
 132   CONTINUE
 
       DO 135 k=1,ks
       DO 135 i=1,kpie
          sedb1(i,k) = sedb1(i,k) - tredsy(i,k-1,1) * sedb1(i,k-1)
135    CONTINUE

       k=ks
       DO 136 i=1,kpie
          IF(omask(i,j).GT.0.5) sediso(i,k) = sedb1(i,k) / tredsy(i,k,2)
136    CONTINUE

       DO 137 k=1,ks
          l=ks-k
          DO 137 i=1,kpie
         IF(omask(i,j).GT.0.5) sediso(i,l) =                           &
     &           ( sedb1(i,l) - tredsy(i,l,3) * sediso(i,l+1) )       &
     &           / tredsy(i,l,2)
 137      CONTINUE


      RETURN
      END
