c***********************************************************************
      subroutine relcon (lrhs)
c***********************************************************************
c
c     this routine calculates those parts of the momentum
c     equations that do not depend on the u and v values at the
c     previous time step, and the coefficients for the u and v
c     values at the grid point (the diagonal terms).
c     this subroutine sets up the equations to be solved by  relax.
c
c-----------------------------------------------------------------------
      implicit none
c-----------------------------------------------------------------------
c Minimum ice mass for calculation of momentum balance
c-----------------------------------------------------------------------
      integer lrhs
c
#include "param_icedyn.h"
#include "common_icedyn.h"
c
      real umax
      parameter ( umax=1.0)                        ! [m/s]
c
c-----------------------------------------------------------------------
c Number of patches in the overrelaxation routine
c-----------------------------------------------------------------------
      integer NP
      parameter ( NP=1 )
c-----------------------------------------------------------------------
c Patch dimension
c-----------------------------------------------------------------------
      integer PD
      parameter ( PD=(MP+NP-1)/NP )
c-----------------------------------------------------------------------
c Actual number of patches used
c-----------------------------------------------------------------------
      integer NPP
      parameter ( NPP=(MP+PD-1)/PD)
c-----------------------------------------------------------------------
c Local variables
c-----------------------------------------------------------------------
      real tmpu, tmpv, tmp1, tmp2, tauax, tauay,
     .     px, py, wta,verr
      integer i,j,jm1,jp1,k,ib,lr,ln,n,pdd,jj
      real UU(L,-1:PD,NPP),VV(L,-1:PD,NPP)

c-----------------------------------------------------------------------
c Start
c-----------------------------------------------------------------------

c$OMP PARALLEL DO PRIVATE(jm1)
      do j=0,M
        jm1=mod(j-1+MP,MP)
        do i=1,L
c-----------------------------------------------------------------------
c Calculate ice mass per area at vector grid points
c-----------------------------------------------------------------------
            amas(i,j) = rhoice*0.25*( h(i-1,j  )+h(i,j  )
     &                               +h(i-1,jm1)+h(i,jm1))
c-----------------------------------------------------------------------
c Very small ice masses are neglected and treated as zero
c-----------------------------------------------------------------------
          if (amas(i,j).lt.amasmin) then
            amas (i,j)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c-----------------------------------------------------------------------
c Calculate those parts of the diagonal terms that depend on the viscosities
c-----------------------------------------------------------------------
      call bcoef(lrhs)
c---------------------------------------------------------------------
c Main loop for all inner grid points
c---------------------------------------------------------------------
c$OMP PARALLEL DO PRIVATE(jm1,px,py,tmp1,tauax,tauay)
      do 600 j=0,m
         jm1=mod(j-1+mp,mp)
         do 600 i=1,l
c---------------------------------------------------------------------
c Pressure term d/dx (- p/2 )  to u equation
c---------------------------------------------------------------------
          px = -0.5*((p(i,j)+p(i,jm1))/(pn(i,j)+pn(i,jm1))
     2       -(p(i-1,j)+p(i-1,jm1))/(pn(i-1,j)+pn(i-1,jm1)))/dx
c---------------------------------------------------------------------
c Pressure term d/dy ( -p/2 )  to v equation
c---------------------------------------------------------------------
          py = -0.5*((p(i,j)+p(i-1,j))/(pm(i,j)+pm(i-1,j))
     2       -(p(i,jm1)+p(i-1,jm1))/(pm(i,jm1)+pm(i-1,jm1)))/dy
c---------------------------------------------------------------------
c Divide fluxes by the area at the grid point to get proper divergences
c---------------------------------------------------------------------
          bu(i,j)=bu(i,j)*pmpnv(i,j)
          bv(i,j)=bv(i,j)*pmpnv(i,j)
          px=px*pmpnv(i,j)
          py=py*pmpnv(i,j)
c---------------------------------------------------------------------
c Add pressure term to F
c---------------------------------------------------------------------
          fx(i,j)=px
          fy(i,j)=py
c---------------------------------------------------------------------
c Local rate of change  -m (du/dt)
c
c This term is usually small. M.Harder
c---------------------------------------------------------------------
c-----------------------------------------------------------------------
c Oceanic drag from rotated geostrophic velocities.
c This will also have contributions to the asymetric term and
c to the fx and fy force terms.
c
c Note: term is - cd *rho*abs(uwat-u)*(uwat-u), two - signs cancel
c       for the diagonal term.
c
c Attention: if daily forcing is applied, the velocity from the previous
c time step may not be a good guess for the new velocity.
c Markus Harder, August 1994
c-----------------------------------------------------------------------
          tmp1 = (uwat(i,j)-u(i,j,lrhs))**2
     1          +(vwat(i,j)-v(i,j,lrhs))**2
          tmp1 = max(tmp1,4e-4)      ! prevent very small drag
          tmp1 = rhowat*cdwat*sqrt(tmp1)
          bu(i,j) = bu(i,j) + tmp1*coswat
          bv(i,j) = bv(i,j) + tmp1*coswat
          asy(i,j) = sinwat(i,j)*tmp1
          fx(i,j) = fx(i,j) + tmp1*(coswat*uwat(i,j)
     1                             -sinwat(i,j)*vwat(i,j))
          fy(i,j) = fy(i,j) + tmp1*(sinwat(i,j)*uwat(i,j)
     1                             +coswat*vwat(i,j))
c---------------------------------------------------------------------
c Coriolis term
c---------------------------------------------------------------------
         asy(i,j) = asy(i,j)+amas(i,j)*f(i,j)
c---------------------------------------------------------------------
c Contribution to the x and y forces due to surface tilt
c---------------------------------------------------------------------
        fx(i,j)=fx(i,j)-amas(i,j)*g*gradhx(i,j)
        fy(i,j)=fy(i,j)-amas(i,j)*g*gradhy(i,j)
c---------------------------------------------------------------------
c Stress on the top of the ice to the x and y forces
c---------------------------------------------------------------------
         tauax=coswin*tauxw(i,j)-sinwin(i,j)*tauyw(i,j)
         tauay=sinwin(i,j)*tauxw(i,j)+coswin*tauyw(i,j)
         fx(i,j)=fx(i,j)+tauax
         fy(i,j)=fy(i,j)+tauay
c---------------------------------------------------------------------
c Contribution from the rate of change that comes from the old time step.
c Note: when on the rhs, this becomes a positive term.
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c Get the reciprocals of bu and bv
C---------------------------------------------------------------------
          bu(i,j) = 1./bu(i,j)
          bv(i,j) = 1./bv(i,j)
  600 continue
c$OMP END PARALLEL DO
      return
C***********************************************************************
      entry relax (LR, LN)
C***********************************************************************
c   Note that RELAX is called twice for each time step IIC.
c   The first  call is indicated by LR=3.
c   The second call is indicated by LN=3.
C     THE VELOCITY ARRAYS WILL BE USED AS FOLLOWS:
C     U OR V (I,J,LN) = VELOCITIES AT NEXT TIME STEP
C     U OR V (I,J,LR) = VELOCITIES USED FOR CALCULATING THE RHS OF THE
C     MOMENTUM EQUATIONS
c
c Extensive documentation done by Markus Harder, April 1993
c
c Major modifications done by Markus Harder, April 1994
c      RELAX runs faster now:
c      The whole program is accelerated by a factor of about 1.2
c
c Further modifications by Markus Harder, September 1994
c
C     THIS SUBROUTINE CALCULATES VELOCITIES AT THE NEXT TIME STEP BY
C     FIRST ESTIMATING THE VELOCITIES USING SEQUENTIAL OVERRELAXATION
C     WITH CHEYBESHEV ACCELERATION (SEE HOCKNEY AND JESSHOPE (1981),
C     PG 334-341).  IN THIS CASE WE ARE UPDATING A CHECKERBOARD ALTERNATING
C     COLORS, IE FIRST RED, THEN BLACK.
c-----------------------------------------------------------------------
C MAKE A FIRST GUESS THAT THE NEW VELOCITIES EQUAL THE old ADVECTIVE
C     VELOCITIES
c
c Attention: if daily forcing is applied, the velocity from the previous
c time step may not be a good guess for the new velocity.
c Markus Harder, August 1994
c-----------------------------------------------------------------------

c$OMP PARALLEL DO
        do j=0,M
          do i=1,L
            u(i,j,ln)=u(i,j,lr)
            v(i,j,ln)=v(i,j,lr)
          enddo
        enddo
c$OMP END PARALLEL DO

c-----------------------------------------------------------------------
C     CALCULATE THE DENOMINATOR OF THE SOLUTION EQUATION
C     MAKE SURE THAT WE DO NOT GET A GARBAGE NUMBER WHERE WE DON'T
C     REALLY NEED TO DO CALCULATION BY USING THE V-MASK
c     -> set DEN = 1.0 where vm=0
c-----------------------------------------------------------------------

c$OMP PARALLEL DO
      do 20 j=0,m
        do 20 i=1,l
   20   den(i,j) = 1./(1.+(vm(i,j)*asy(i,j))**2*
     1                         bu(i,j)*bv(i,j))
c$OMP END PARALLEL DO
c-----------------------------------------------------------------------
C     NOW START THE OVERRELAXATION PROCEDURE
C     WE ARE SOLVING THE EQUATIONS:
C     U - (ASY * BU) * V = (RU + FX) * BU
C     V + (ASY * BV) * U = (RV + FY) * BV
C     WHERE ASY IS THE AXYMETRICAL TERM (CORIOLIS + WATER DRAG),
C           BU AND BV ARE THE RECIPROCALS OF THE DIAGONAL TERMS,
C           FX AND FY ARE THOSE PARTS OF THE MOMENTUM EQUATION THAT
C           DO NOT DEPENDED ON THE VELOCITIES AT THE NEW TIME STEP,
C           AND RU AND RV WILL DEPEND ON THE VALUES AT THE OTHER
C           GRID POINTS, BUT NOT DIRECTLY ON THE LOCAL GRID POINT.
C
c-----------------------------------------------------------------------
            mrelax=      0            ! iteration step counter
      verr=      1000.            ! large dummy value for first step
c-----------------------------------------------------------------------
c Main loop of RELAX
c-----------------------------------------------------------------------
  120      continue


c-----------------------------------------------------------------------
c Stop if we have not found a solution after MMAX iterations
c-----------------------------------------------------------------------
      if (mrelax .gt. mmax) goto 400
c-----------------------------------------------------------------------
c Do not use overrelaxation for first and for late iteration steps
c-----------------------------------------------------------------------
      if ((verr.lt.vrmax*vrwt).or.(mrelax.gt.mwt).or.(mrelax.le.1)) then
         wta = 1.0
      else
         wta = wt
      end if
c=======================================================================
c The checkerboard pattern
c ========================
c
c The velocity fields u(i,j), v(i,j) and similar fields with vector
c components (i=1..L,j=1..M) may be represented as a checkerboard with
c 2 alternating colors, e.g. black (b) and red (r).
c
c The lower left corner (1,1) is defined to be black.
c If the sum (i+j) of the indices of a vector field, e.g. u(i,j),
c is even, the corresponding field is black.
c Otherwise, i.e. if (i+j) is odd, the corresponding field is red.
c
c Thus the black fields are
c   (1,1), (1,3), (1,5), ..., (2,2), (2,4), (2,6), ...
c and the red fields are
c   (1,2), (1,4), (1,6), ..., (2,1), (2,3), (2,5), ...
c
c   j            . . .                    . . .
c
c     rrr bbb rrr bbb rrr bbb rrr bbb
c   4 rrr bbb rrr bbb rrr bbb rrr bbb     . . .
c     rrr bbb rrr bbb rrr bbb rrr bbb
c
c     bbb rrr bbb rrr bbb rrr bbb rrr
c   3 bbb rrr bbb rrr bbb rrr bbb rrr     . . .
c     bbb rrr bbb rrr bbb rrr bbb rrr
c
c     rrr bbb rrr bbb rrr bbb rrr bbb
c   2 rrr bbb rrr bbb rrr bbb rrr bbb     . . .
c     rrr bbb rrr bbb rrr bbb rrr bbb
c
c     bbb rrr bbb rrr bbb rrr bbb rrr
c   1 bbb rrr bbb rrr bbb rrr bbb rrr     . . .
c     bbb rrr bbb rrr bbb rrr bbb rrr
c      1   2   3   4   5   6   7   8      . . .     i
c=======================================================================
c-----------------------------------------------------------------------
C     DO RELAXATION ALTERNATING COLORS FOR CHECKERBOARD PATTERN
c
c K=0      black fields
c K=1      red fields
c-----------------------------------------------------------------------
c Apply checkerboard pattern:
c K=0: black fields
c   J even      I=2,4,6,...,48            (for LM=48)
c   J odd      I=3,5,7,...,47            (for LM=48)
c
c K=1: red fields
c   J odd      I=2,4,6,...,48            (for LM=48)
c   J even      I=3,5,7,...,47            (for LM=48)
c-----------------------------------------------------------------------

      verr=0.
      DO 171 K = 0,1
c$OMP PARALLEL DO PRIVATE(PDD,JM1,J,JP1,IB,TMP1,TMP2,TMPU,TMPV)
c$OMP+  REDUCTION(MAX:VERR)
      DO 172 N=1,NPP
      PDD=MOD(MIN(N*PD-1,M),PD)

c-----------------------------------------------------------------------
c Copy velocity into patch arrays
c-----------------------------------------------------------------------
      JM1=MOD((N-1)*PD-1+MP,MP)
      DO I=1,L
        UU(I,-1,N)=U(I,JM1,LN)
        VV(I,-1,N)=V(I,JM1,LN)
      ENDDO
      DO JJ=0,PDD
        J=(N-1)*PD+JJ
        DO I=1,L
          UU(I,JJ,N)=U(I,J,LN)
          VV(I,JJ,N)=V(I,J,LN)
        ENDDO
      ENDDO
      JP1=MOD((N-1)*PD+PDD+1,MP)
      DO I=1,L
        UU(I,PDD+1,N)=U(I,JP1,LN)
        VV(I,PDD+1,N)=V(I,JP1,LN)
      ENDDO

      DO 173 JJ=0,PDD
        J=(N-1)*PD+JJ
        JM1=MOD(J-1+MP,MP)
        JP1=MOD(J+1,MP)


c-----------------------------------------------------------------------
c Calculate index IB to start I-iteration with
c   K=0:      IB=2 for even J            J=2,4,6,...12 (for MM=13)
c            IB=3 for odd  J            J=3,5,7,...13 (for MM=13)
c
c   K=1:      IB=3 for even J
c            IB=2 for odd  J
c-----------------------------------------------------------------------
        IB=MOD(J+K+2,2)+2
        do 173 I=IB,LM,2
c-----------------------------------------------------------------------
C     ADD D/DX ( ZETA * (DU/DX + DV/DY) )  TO U EQUATION
c-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C     CALCULATE D/DX(ZETA * DU/DX)
C-----------------------------------------------------------------------
           RU(I,J) = (ZETA(I,J)  +ZETA(I,JM1))  *UU(I+1,JJ,N)*sx2p(i,j)
     &             +(ZETA(I-1,J)+ZETA(I-1,JM1))*UU(I-1,JJ,N)*sx2m(i,j)


C-----------------------------------------------------------------------
C     CALCULATE D/DX( ZETA * DV/DY)
C     NOTE THIS WILL INCLUDE THE TERM AT THE GRID POINT I,J, SINCE
C     FOR ALL CASES THIS WILL BE USED FOR THE OTHER EQUATION, IE
C     IF CALLED WITH U IT WILL BE FOR THE V EQUATION.
C-----------------------------------------------------------------------
          RU(I,J) = RU(I,J) + SXY*(ZETA(I,J)*
     1     ((VV(I,JJ+1,N)  +VV(I+1,JJ+1,N))
     2     -(VV(I,JJ,N)    +VV(I+1,JJ,N)))
     3                     +ZETA(I,JM1)*
     4     ((VV(I,JJ,N)    +VV(I+1,JJ,N))
     5     -(VV(I,JJ-1,N)  +VV(I+1,JJ-1,N)))
     6                     -ZETA(I-1,J)*
     1     ((VV(I-1,JJ+1,N)+VV(I,JJ+1,N))
     2     -(VV(I-1,JJ,N)  +VV(I,JJ,N)))
     3                     -ZETA(I-1,JM1)*
     4     ((VV(I-1,JJ,N)  +VV(I,JJ,N))
     5     -(VV(I-1,JJ-1,N)+VV(I,JJ-1,N))))

c-----------------------------------------------------------------------
C     ADD D/DY ( ZETA * (DU/DX + DV/DY) )  TO V EQUATION
c-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     CALCULATE D/DY( ZETA * DV/DY )
C-----------------------------------------------------------------------
          RV(I,J) = (ZETA(I,J)  +ZETA(I-1,J))  *VV(I,JJ+1,N)*sy2p(i,j)
     &             +(ZETA(I,JM1)+ZETA(I-1,JM1))*VV(I,JJ-1,N)*sy2m(i,j)
C-----------------------------------------------------------------------
C     CALCULATE D/DY( ZETA * DV/DX)
C-----------------------------------------------------------------------
          RV(I,J) = RV(I,J) + SXY*(ZETA(I,J)*
     1     ((UU(I+1,JJ,N)  +UU(I+1,JJ+1,N))
     2     -(UU(I,JJ,N)    +UU(I,JJ+1,N)))
     3                     +ZETA(I-1,J)*
     4     ((UU(I,JJ,N)    +UU(I,JJ+1,N))
     5     -(UU(I-1,JJ,N)  +UU(I-1,JJ+1,N)))
     6                     -ZETA(I,JM1)*
     1     ((UU(I+1,JJ-1,N)+UU(I+1,JJ,N))
     2     -(UU(I,JJ-1,N)  +UU(I,JJ,N)))
     3                     -ZETA(I-1,JM1)*
     4     ((UU(I,JJ-1,N)  +UU(I,JJ,N))
     5     -(UU(I-1,JJ-1,N)+UU(I-1,JJ,N))))
c-----------------------------------------------------------------------
C     ADD D/DX ( ETA * (DU/DX - DV/DY) ) TO U EQUATION
c-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     CALCULATE D/DX(ETA * DU/DX)
C-----------------------------------------------------------------------
      RU(I,J)=RU(I,J)
     &            +(ETA(I,J)  +ETA(I,JM1))  *UU(I+1,JJ,N)*sx2p(i,j)
     &            +(ETA(I-1,J)+ETA(I-1,JM1))*UU(I-1,JJ,N)*sx2m(i,j)

C-----------------------------------------------------------------------
C     CALCULATE D/DX( ETA * DV/DY)
C     NOTE THIS WILL INCLUDE THE TERM AT THE GRID POINT I,J, SINCE
C     FOR ALL CASES THIS WILL BE USED FOR THE OTHER EQUATION, IE
C     IF CALLED WITH U IT WILL BE FOR THE V EQUATION.
C-----------------------------------------------------------------------
          RU(I,J) = RU(I,J) - SXY*(ETA(I,J)*      ! MINUS SIGN !
     1     ((VV(I,JJ+1,N)  +VV(I+1,JJ+1,N))
     2     -(VV(I,JJ,N)    +VV(I+1,JJ,N)))
     3                     +ETA(I,JM1)*
     4     ((VV(I,JJ,N)    +VV(I+1,JJ,N))
     5     -(VV(I,JJ-1,N)  +VV(I+1,JJ-1,N)))
     6                     -ETA(I-1,J)*
     1     ((VV(I-1,JJ+1,N)+VV(I,JJ+1,N))
     2     -(VV(I-1,JJ,N)  +VV(I,JJ,N)))
     3                     -ETA(I-1,JM1)*
     4     ((VV(I-1,JJ,N)  +VV(I,JJ,N))
     5     -(VV(I-1,JJ-1,N)+VV(I,JJ-1,N))))

c-----------------------------------------------------------------------
C     ADD D/DY ( ETA * (DV/DY - DU/DX) ) TO V EQUATION
c-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     CALCULATE D/DY( ETA * DV/DY )
C-----------------------------------------------------------------------
      RV(I,J)=RV(I,J)
     &            +(ETA(I,J)  +ETA(I-1,J))  *VV(I,JJ+1,N)*sy2p(i,j)
     &            +(ETA(I,JM1)+ETA(I-1,JM1))*VV(I,JJ-1,N)*sy2m(i,j)
C-----------------------------------------------------------------------
C     CALCULATE D/DY( ETA * DU/DX)
C-----------------------------------------------------------------------
          RV(I,J) = RV(I,J) - SXY*(ETA(I,J)*      ! MINUS SIGN !
     1     ((UU(I+1,JJ,N)  +UU(I+1,JJ+1,N))
     2     -(UU(I,JJ,N)    +UU(I,JJ+1,N)))
     3                     +ETA(I-1,J)*
     4     ((UU(I,JJ,N)    +UU(I,JJ+1,N))
     5     -(UU(I-1,JJ,N)  +UU(I-1,JJ+1,N)))
     6                     -ETA(I,JM1)*
     1     ((UU(I+1,JJ-1,N)+UU(I+1,JJ,N))
     2     -(UU(I,JJ-1,N)  +UU(I,JJ,N)))
     3                     -ETA(I-1,JM1)*
     4     ((UU(I,JJ-1,N)  +UU(I,JJ,N))
     5     -(UU(I-1,JJ-1,N)+UU(I-1,JJ,N))))
c-----------------------------------------------------------------------
C     ADD D/DY ( ETA * (DU/DY + DV/DX) ) TO U EQUATION
c-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     CALCULATE D/DY( ETA * DU/DY )
C-----------------------------------------------------------------------
      RU(I,J)=RU(I,J)
     &            +(ETA(I,J)  +ETA(I-1,J))  *UU(I,JJ+1,N)*sy2p(i,j)
     &            +(ETA(I,JM1)+ETA(I-1,JM1))*UU(I,JJ-1,N)*sy2m(i,j)

C-----------------------------------------------------------------------
C     CALCULATE D/DY( ETA * DV/DX)
C-----------------------------------------------------------------------
          RU(I,J) = RU(I,J) + SXY*(ETA(I,J)*
     1     ((VV(I+1,JJ,N)  +VV(I+1,JJ+1,N))
     2     -(VV(I,JJ,N)    +VV(I,JJ+1,N)))
     3                     +ETA(I-1,J)*
     4     ((VV(I,JJ,N)    +VV(I,JJ+1,N))
     5     -(VV(I-1,JJ,N)  +VV(I-1,JJ+1,N)))
     6                     -ETA(I,JM1)*
     1     ((VV(I+1,JJ-1,N)+VV(I+1,JJ,N))
     2     -(VV(I,JJ-1,N)  +VV(I,JJ,N)))
     3                     -ETA(I-1,JM1)*
     4     ((VV(I,JJ-1,N)  +VV(I,JJ,N))
     5     -(VV(I-1,JJ-1,N)+VV(I-1,JJ,N))))

c-----------------------------------------------------------------------
C     ADD D/DX ( ETA * (DU/DY + DV/DX) ) TO V EQUATION
c-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     CALCULATE D/DX(ETA * DV/DX)
C-----------------------------------------------------------------------
      RV(I,J)=RV(I,J)
     &            +(ETA(I,J)  +ETA(I,JM1))  *VV(I+1,JJ,N)*sx2p(i,j)
     &            +(ETA(I-1,J)+ETA(I-1,JM1))*VV(I-1,JJ,N)*sx2m(i,j)
C-----------------------------------------------------------------------
C     CALCULATE D/DX( ETA * DU/DY)
C     NOTE THIS WILL INCLUDE THE TERM AT THE GRID POINT I,J, SINCE
C     FOR ALL CASES THIS WILL BE USED FOR THE OTHER EQUATION, IE
C     IF CALLED WITH U IT WILL BE FOR THE V EQUATION.
C-----------------------------------------------------------------------
          RV(I,J) = RV(I,J) + SXY*(ETA(I,J)*
     1     ((UU(I,JJ+1,N)  +UU(I+1,JJ+1,N))
     2     -(UU(I,JJ,N)    +UU(I+1,JJ,N)))
     3                     +ETA(I,JM1)*
     4     ((UU(I,JJ,N)    +UU(I+1,JJ,N))
     5     -(UU(I,JJ-1,N)  +UU(I+1,JJ-1,N)))
     6                     -ETA(I-1,J)*
     1     ((UU(I-1,JJ+1,N)+UU(I,JJ+1,N))
     2     -(UU(I-1,JJ,N)  +UU(I,JJ,N)))
     3                     -ETA(I-1,JM1)*
     4     ((UU(I-1,JJ,N)  +UU(I,JJ,N))
     5     -(UU(I-1,JJ-1,N)+UU(I,JJ-1,N))))
c-----------------------------------------------------------------------
C     ADD IN THE ADVECTION TERMS AND DIVIDE BY THE AREA AROUND
C     THE GRID POINT TO GET FLUX DIVERGENCES
cmh      CALL MADV(LR,LN,K)      ! old code
c
c Advection of momentum is neglected. M. Harder, 25-JAN-95
c-----------------------------------------------------------------------
C     ADD IN - M * U * D/DX (U)
cmh          RU(I,J) = RU(I,J)-AMAS(I,J)*U(I,J,LR)*
cmh     1     (UU(I+1,JJ,N)/(PN(I,J)+PN(I,J-1))
cmh     2     -UU(I-1,JJ,N)/(PN(I-1,J)+PN(I-1,J-1)))/DX
C     ADD IN - M * V * D/DY (U)
cmh     2                     -AMAS(I,J)*V(I,J,LR)*
cmh     3     (UU(I,JJ+1,N)/(PM(I,J)+PM(I-1,J))
cmh     4     -UU(I,JJ-1,N)/(PM(I,J-1)+PM(I-1,J-1)))/DY
C
C     ADD IN - M * U D/DX (V)
cmh          RV(I,J) = RV(I,J)-AMAS(I,J)*U(I,J,LR)*
cmh     1     (VV(I+1,JJ,N)/(PN(I,J)+PN(I,J-1))
cmh     2     -VV(I-1,JJ,N)/(PN(I-1,J)+PN(I-1,J-1)))/DX
C     ADD IN - M * V D/DY (V)
cmh     2                     -AMAS(I,J)*V(I,J,LR)*
cmh     3     (VV(I,JJ+1,N)/(PM(I,J)+PM(I-1,J))
cmh     4     -VV(I,JJ-1,N)/(PM(I,J-1)+PM(I-1,J-1)))/DY

          RU(I,J) = RU(I,J)*pmpnv(I,J)
          RV(I,J) = RV(I,J)*pmpnv(I,J)

c-----------------------------------------------------------------------
C     ADD IN THE OTHER PARTS OF THE MOMENTUM EQUATION AND DIVIDE BY
C     THE DIAGONAL TERM
c-----------------------------------------------------------------------
         TMP1 = VM(I,J)*(RU(I,J)+FX(I,J))*BU(I,J)
         TMP2 = VM(I,J)*(RV(I,J)+FY(I,J))*BV(I,J)
c-----------------------------------------------------------------------
C      NOW SOLVE FOR THE NEXT GUESS OF THE VELOCITIES
c-----------------------------------------------------------------------
         TMPU = (TMP1+ASY(I,J)*BU(I,J)*TMP2)*
     1                 DEN(I,J)
         TMPV = (TMP2-ASY(I,J)*BV(I,J)*TMP1)*
     1                 DEN(I,J)
c-----------------------------------------------------------------------
c Calculate the difference to the old velocities,
c and multiply by overrelaxation factor WTA for a faster approach
c to the solution.
c-----------------------------------------------------------------------
         TMPU = WTA*(TMPU-UU(I,JJ,N))
         TMPV = WTA*(TMPV-VV(I,JJ,N))
c-----------------------------------------------------------------------
c Keep track of the maximum error of actual solution -> VERR
c The goal is to reach an equilibrium state where the modifications
c TMPU,TMPV to the velocities U and V are very small.
c-----------------------------------------------------------------------
         VERR=max(abs(TMPU),VERR)
         VERR=max(abs(TMPV),VERR)
c-----------------------------------------------------------------------
C     UPDATE THIS HALF OF THE VELOCITIES
C     (INCLUDING HALF THE OUTER BOUNDARY)
c-----------------------------------------------------------------------
         UU(I,JJ,N) =  UU(I,JJ,N)+TMPU
         VV(I,JJ,N) =  VV(I,JJ,N)+TMPV
C-----------------------------------------------------------------------
C     SET VELOCITIES TO 0 USING VELOCITY MASK
c
cmh ( replaces call of BCSV )
C-----------------------------------------------------------------------
         UU(I,JJ,N) = VM(I,J)*UU(I,JJ,N)
         VV(I,JJ,N) = VM(I,J)*VV(I,JJ,N)
c-----------------------------------------------------------------------
  173      continue            ! I,J loop
  172      continue            ! N loop
c$OMP END PARALLEL DO

c$OMP PARALLEL DO PRIVATE(PDD,J)
      DO 174 N=1,NPP
      PDD=MOD(MIN(N*PD-1,M),PD)

c-----------------------------------------------------------------------
c Copy velocity from patch arrays
c-----------------------------------------------------------------------
      DO JJ=0,PDD
        J=(N-1)*PD+JJ
        DO I=1,L
          U(I,J,LN)=UU(I,JJ,N)
          V(I,J,LN)=VV(I,JJ,N)
        ENDDO
      ENDDO

  174      continue
c$OMP END PARALLEL DO

  171      continue            ! K   loop
c-----------------------------------------------------------------------
c Increment iteration counter
c-----------------------------------------------------------------------
      mrelax=mrelax+1
c-----------------------------------------------------------------------
c We have found a new solution if the error is small enough.
c Otherwise, we repeat the iteration loop.
c-----------------------------------------------------------------------
      if (verr.gt.vrmax) goto 120
c-----------------------------------------------------------------------
c We have found a new solution for the coefficients that have been
c set up by subroutine relcon. Some of these coefficients are
c NOT really constants but depend on the velocity field
c (especially for the quadratic water drag term),
c These coefficients are correct for our first guess of the new velocity
c field but are often inaccurate for the new velocity field.
c To fix this problem, an iteration with several calls of RELAX could be used.
c
c The direct effect of improperly set coefficients is that very large
c velocities may occur.
c Physically, these large velocities make no sense. Numerically, they
c mean trouble and at least a slower approach to an accurate solution.
c
c So these large velocities should be cut off. However, this is an
c ARTIFICIAL modification of the momentum equation. The final solution
c for the velocities should not contain any extremly large speeds.
c Here, it is decided to count the number of values that are cut off,
c and to print them to a statistics file. If the final solution requires
c any cut-offs, it should be seriously taken as a sign that the momentum
c equation has not been solved correctly.
c
c      Markus Harder, September 1994
c-----------------------------------------------------------------------

      numax=0
      numaxi=0

c$OMP PARALLEL DO REDUCTION(+:numax,numaxi)
      do j=0,M
        do i=1,L
          if ((abs(u(i,j,ln)).gt.umax).or.(abs(v(i,j,ln)).gt.umax)
     &             .and.(vm(i,j).gt.0)) then
            numax=numax+1
            if (amas(i,j).gt.0.) numaxi=numaxi+1
          endif
        enddo
      enddo
c$OMP END PARALLEL DO

c$OMP PARALLEL DO
      do j=0,M
        do i=1,L
          u(i,j,ln)=max(u(i,j,ln),-umax)
          u(i,j,ln)=min(u(i,j,ln), umax)
          v(i,j,ln)=max(v(i,j,ln),-umax)
          v(i,j,ln)=min(v(i,j,ln), umax)
        enddo
      enddo
c$OMP END PARALLEL DO
c-----------------------------------------------------------------------
c Finished
c-----------------------------------------------------------------------
        nosolution=0
      return
c-----------------------------------------------------------------------
c We have not found a solution after MMAX iterations
c
c This does not necessarily mean that no solution exists.
c Rather, it indicates in many cases that the velocity field
c of the solution is strongly different to the velocity field
c that we have used to set up the equation system of RELAX.
c
c One way to handle this problem may be to take the velocities
c predicted by RELAX, despite they are not accurate solutions,
c as a better guess for the new velocities. The coefficients of RELAX
c should be recalculated (by RELCON) with this new guess,
c then RELAX called again.
c
c This method will work quite well if the problem behaves well but
c requires just some more iteration steps due to rapidly varying velocities
c (especially when the forcing functions vary strongly from one time step
c to the next).
c This method will not work if the problem is ill-posed, e.g. if numerical
c instabilities occur due to a bad-behaving rheology.
c
c In this routine, we do not decide which action to take if no solution
c is found. We set a flag indicating the problem and leave the action
c to the calling routine.
c-----------------------------------------------------------------------
  400      continue
        nosolution=1
      return
      end            ! subroutine relax
