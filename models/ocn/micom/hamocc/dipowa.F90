      SUBROUTINE DIPOWA(kpie,kpje,kpke,pdlxp,pdlyp,omask)

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/dipowa.f90,v $\\
!$Revision: 1.2.20.1.2.1 $\\
!$Date: 2006/04/03 11:27:49 $\\
!
!****************************************************************
!
!**** *DIPOWA* - 'diffusion of pore water'
!      vertical diffusion of sediment pore water tracers
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - all npowtra-1 properties are diffused in 1 go.
!     js: not mass conserving check c13/powtra/ocetra
!
!     Purpose
!     -------
!     calculate vertical diffusion of sediment pore water properties
!     and diffusive flux through the ocean/sediment interface.
!     integration.
!
!     Method
!     -------
!     implicit formulation;
!     constant diffusion coefficient : 1.e-9 set in BODENSED.
!     diffusion coefficient : zcoefsu/zcoeflo for upper/lower
!     sediment layer boundary.
!
!**   Interface.
!     ----------
!
!     *CALL*       *DIPOWA*
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

      implicit none

      INTEGER :: i,j,k,l,iv
      INTEGER :: kpie,kpje,kpke
      integer :: iv_oc                                ! index of ocetra in powtra loop

      REAL :: sedb1(kpie,0:ks,npowtra)                ! ????
      REAL :: zcoefsu(0:ks),zcoeflo(0:ks)             ! diffusion coefficients (upper/lower)
      REAL :: TREDSY(kpie,0:kpke,3)                   ! redsy for 'reduced system'?
      
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje),omask(kpie,kpje)
      REAL :: aprior                                  ! start value of oceanic tracer in bottom layer
      
!ik accelerated sediment
!ik needed for boundary layer ventilation in fast sediment routine      

      REAL :: bolven(kpie)                ! bottom layer ventilation rate

      zcoefsu( 0)=0.0
      DO  k=1,ks   
         zcoefsu(k  )=-sedict*seddzi(k)*porwah(k)    ! sediment diffusion coefficient * 1/dz * fraction of pore water at half depths
         zcoeflo(k-1)=-sedict*seddzi(k)*porwah(k)    ! why the same ?
      ENDDO
      zcoeflo(ks)=0.0                                ! diffusion coefficient for bottom sediment layer

!$OMP PARALLEL DO                            &                
!$OMP&PRIVATE(bolven,tredsy,sedb1,aprior,iv_oc)
      DO 11000 j=1,kpje

! calculate bottom ventilation rate for scaling of sediment-water exchange
      do 1170 i=1,kpie      
        bolven(i) = 1.
1170  continue


      k=0
      DO 1421 i=1,kpie
         tredsy(i,k,1) = zcoefsu(k)
         tredsy(i,k,3) = zcoeflo(k)
         tredsy(i,k,2) = bolven(i)*bolay(i,j) - tredsy(i,k,1) - tredsy(i,k,3) ! dz(kbo) - diff upper - diff lower
1421  CONTINUE

      k=0
      DO 1422 iv=1,npowtra      ! loop over pore water tracers
      iv_oc=iv
#ifdef __c_isotopes
      if(iv.eq.ipowc13) iv_oc=isco213
      if(iv.eq.ipowc14) iv_oc=isco214
#endif
      DO 1422 i=1,kpie
         sedb1(i,k,iv)=0.
         IF(omask(i,j).GT.0.5)                                     &
!         IF(bolay(i,j).GT.0.)                                     &
     &      sedb1(i,k,iv)=ocetra(i,j,kbo(i,j),iv_oc)*bolay(i,j)*bolven(i)     ! tracer_concentration(kbo) * dz(kbo)
 1422 CONTINUE

      DO 1321 k=1,ks
         DO 1321 i=1,kpie     
            tredsy(i,k,1) = zcoefsu(k)
            tredsy(i,k,3) = zcoeflo(k)
            tredsy(i,k,2) = seddw(k)*porwat(k) -tredsy(i,k,1) -tredsy(i,k,3)      
 1321    CONTINUE

      DO 1322 iv=1,npowtra
      DO 1322 k=1,ks
      DO 1322 i=1,kpie     
         sedb1(i,k,iv)=powtra(i,j,k,iv)*porwat(k)*seddw(k) ! tracer_concentration(k[1:ks]) * porewater fraction(k) * dz(k)
 1322 CONTINUE

         DO 132 k=1,ks
            DO 133 i=1,kpie
               IF(omask(i,j).GT.0.5) THEN     
!               IF(bolay(i,j).GT.0.) THEN
                  tredsy(i,k-1,1) = tredsy(i,k,1) / tredsy(i,k-1,2)    ! this overwrites tredsy(k=0) for k=1
!                                   diff upper    / conc (k-1)
                  tredsy(i,k,2)   = tredsy(i,k,2)                      &
     &                             -tredsy(i,k-1,3)*tredsy(i,k,1) /tredsy(i,k-1,2)
!                    concentration -diff lower     * diff upper   /conc(k-1)
               ENDIF
 133        CONTINUE
 132     CONTINUE
 
! diffusion from above
      DO 135 iv=1,npowtra
         DO 135 k=1,ks
         DO 135 i=1,kpie
            sedb1(i,k,iv) = sedb1(i,k,iv)                        &
     &          - tredsy(i,k-1,1) * sedb1(i,k-1,iv)       
 135     CONTINUE

! sediment bottom layer
         k=ks
      DO 136 iv=1,npowtra
         DO 136 i=1,kpie
         IF(omask(i,j).GT.0.5) THEN     
!         IF(bolay(i,j).GT.0.) THEN
               powtra(i,j,k,iv) = sedb1(i,k,iv) / tredsy(i,k,2)
            ENDIF
 136     CONTINUE

!    call maschk(kpie,kpje,kpke,23)
! sediment column
      DO 137 iv=1,npowtra
         DO 137 k=1,ks-1
            l=ks-k
            DO 137 i=1,kpie
                IF(omask(i,j).GT.0.5) THEN     
!               IF(bolay(i,j).GT.0.) THEN
                  powtra(i,j,l,iv) = ( sedb1(i,l,iv)            &
     &                - tredsy(i,l,3) * powtra(i,j,l+1,iv) )    &
     &                / tredsy(i,l,2)
               ENDIF
 137    CONTINUE
!    call maschk(kpie,kpje,kpke,24)

! sediment ocean interface
      DO 139 iv=1,npowtra! caution - the following assumes same indecees for ocetra and powtra test npowa_base 071106
                            ! check mo_param1_bgc.f90 for consistency
      iv_oc=iv
#ifdef __c_isotopes
      if(iv.eq.ipowc13) iv_oc=isco213
      if(iv.eq.ipowc14) iv_oc=isco214
#endif
        DO 139 i=1,kpie
          l=0
           IF(omask(i,j).GT.0.5) THEN     
!          IF(bolay(i,j).GT.0.) THEN

           aprior = ocetra(i,j,kbo(i,j),iv_oc)
           ocetra(i,j,kbo(i,j),iv_oc) =                                  &
     &         ( sedb1(i,l,iv) - tredsy(i,l,3) * powtra(i,j,l+1,iv) ) &
     &         / tredsy(i,l,2) 

           sedfluxo(i,j,iv)=sedfluxo(i,j,iv)                          &    !used in inventory_bgc/maschk (diagnostics) 
     &                     +ocetra(i,j,kbo(i,j),iv)-aprior

           ENDIF
 139  CONTINUE
!     call maschk(kpie,kpje,kpke,25)

   
11000 CONTINUE    ! j loop  

      RETURN
     END
