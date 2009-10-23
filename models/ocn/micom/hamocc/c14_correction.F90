!     c14 correction for biological effects to save spin-up time
!     insert in aufr_bgc.f90 ONCE! 
!     check if all necessary parameters are known at time of call to aufr
!     satoxy is computed in chemcon (called before aufr, seems ok.)

      subroutine c14_correction(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,psao,ptho,omask)

#ifdef __c_isotopes
      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc

      INTEGER  kpie,kpje,kpke
      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: ptho (kpie,kpje,kpke)
      REAL :: psao (kpie,kpje,kpke)
      REAL :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      real  aq,cq,sq                        ! global average alkalinity, 12C, salinity
      real  volcell
      real  totvol                          ! total volume of wet cells
      real  an (kpie,kpje,kpke)             ! salinity-normalised alkalinity
      real  cn (kpie,kpje,kpke)             ! salinity-normalised 12C
      real  orgrem(kpie,kpje,kpke)          ! remineralized organic carbon
      real  caldis(kpie,kpje,kpke)          ! remineralized carbonate shells
      real  omask(kpie,kpje)          ! remineralized carbonate shells
      totvol =0.

      do k=1,kpke
      do j=1,kpje
      do i=1,kpie

       if (omask(i,j).gt.0.5 .and. ptho(i,j,k).lt.10.) then               ! wet cells with T < 10 degC

       volcell = pdlxp(i,j) * pdlyp(i,j) * pddpo(i,j,k)

       aq = aq + ocetra(i,j,k,ialkali) * volcell     ! sum up alkalinity 
       cq = cq + ocetra(i,j,k,isco212) * volcell     !        DIC
       sq = sq + psao(i,j,k)           * volcell     !        salinity

       totvol = totvol + volcell

       endif

      enddo
      enddo
      enddo

! global averaging
      aq = aq / totvol
      cq = cq / totvol
      sq = sq / totvol

      do k=1,kpke
      do j=1,kpje
      do i=1,kpie

       if (omask(i,j).gt.0.5) then

! normalisation
       an(i,j,k) = ocetra(i,j,k,ialkali) * sq / psao(i,j,k)     ! alkalinity
       cn(i,j,k) = ocetra(i,j,k,isco212) * sq / psao(i,j,k)     ! 12C

       endif

      enddo
      enddo
      enddo

! add biogenic 14C
      do k=1,kpke
      do j=1,kpje
      do i=1,kpie
       if (omask(i,j).gt.0.5) then

!                              delta C           delta Alk
       orgrem(i,j,k) = (2*(cn(i,j,k) - cq) - (an(i,j,k) - aq)) / 1.7
! 2nd approach for orgrem (AOU-based, add in full)
       orgrem(i,j,k) = (satoxy(i,j,k) - ocetra(i,j,k,ioxygen)) * 122./172.

       caldis(i,j,k) = orgrem(i,j,k) - (cn(i,j,k) - cq)

       ocetra(i,j,k,isco214) = ocetra(i,j,k,isco214)               &
     &                       + 0.9 * orgrem(i,j,k)                 &
!    &                       + orgrem(i,j,k)                       & ! for AOU-based orgrem
     &                       + 0.95* caldis(i,j,k)                   ! comment out for AOU based orgrem (?)

       endif

      enddo
      enddo
      enddo

! end 14C correction
#endif /*__c_isotopes*/

      return
      end





