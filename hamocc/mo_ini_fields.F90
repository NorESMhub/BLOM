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

module mo_ini_fields

  implicit none
  private

  public :: ini_fields_ocean
  public :: ini_fields_atm

contains

  subroutine ini_fields_atm(kpie,kpje)
    !***********************************************************************************************
    ! Initialize bgc variables for atmosphere (used if not coupled to an atmosphere model).
    !
    !  Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !  Modified
    !  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
    !  -split the original BELEG_BGC in two parts, BELEG_PARM (NOW MO_PARAM_BGC) and BELEG_VARS
    !***********************************************************************************************

    use mo_control_bgc, only: use_natDIC,use_cisonew,use_BROMO
    use mo_param1_bgc,  only: iatmco2,iatmo2,iatmn2,iatmnco2,iatmc13,iatmc14,iatmbromo
    use mo_param_bgc,   only: atm_o2,atm_n2,atm_co2_nat,atm_c13,atm_c14,c14fac,atm_bromo
    use mo_carbch,      only: atm,atm_co2

    ! Initialise atmosphere fields. We use a 2D representation of atmospheric
    ! fields for simplicity, even for cases where actually only a scalar value
    ! is used. The overhead of this is small. If an atm-field is present in
    ! restart file (if BOXATM is activated), this will be overwritten later.

    ! Arguments
    integer, intent(in) :: kpie,kpje

    ! local variables
    integer             :: i,j

    do j=1,kpje
      do i=1,kpie
        atm(i,j,iatmco2)  = atm_co2
        atm(i,j,iatmo2)   = atm_o2
        atm(i,j,iatmn2)   = atm_n2
        if (use_natDIC) then
          atm(i,j,iatmnco2) = atm_co2_nat
        endif
        if (use_cisonew) then
          atm(i,j,iatmc13)  = atm_c13
          atm(i,j,iatmc14)  = atm_c14/c14fac
        endif
        if (use_BROMO) then
          atm(i,j,iatmbromo)= atm_bromo
        endif
      enddo
    enddo
  end subroutine ini_fields_atm


  subroutine ini_fields_ocean(kpaufr,kpie,kpje,kpke,kbnd,pddpo,prho,omask,pglon,pglat)

    !***********************************************************************************************
    ! Initialize bgc variables.
    !
    !  Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
    !
    !  Modified
    !  J.Schwinger,        *NORCE Climate, Bergen*    2020-05-19
    !  -split the original BELEG_BGC in two parts, BELEG_PARM (NOW MO_PARAM_BGC) and BELEG_VARS
    !***********************************************************************************************

    use mo_carbch,      only: co2star,co3,hi,ocetra
    use mo_param_bgc,   only: fesoly,cellmass,fractdim,bifr13_ini,bifr14_ini,c14fac,re1312,re14to
    use mo_biomod,      only: abs_oce
    use mo_control_bgc, only: rmasks,use_FB_BGC_OCE,use_cisonew,use_AGG,use_CFC,use_natDIC,        &
                              use_BROMO, use_sedbypass, use_dom
    use mo_param1_bgc,  only: ialkali,ian2o,iano3,icalc,idet,idicsat,idms,idoc,ifdust,igasnit,     &
                              iiron,iopal,ioxygen,iphosph,iphy,iprefalk,iprefdic,iprefo2,iprefpo4, &
                              isco212,isilica,izoo,iadust,inos,ibromo,icfc11,icfc12,isf6,          &
                              icalc13,icalc14,idet13,idet14,idoc13,idoc14,iphy13,iphy14,           &
                              isco213,isco214,izoo13,izoo14,safediv,inatcalc,                      &
                              idocsl,idocsr,idocr,iprefdoc,iprefdocsl,iprefdocsr,iprefdocr,        &
                              ipowaal,ipowaic,ipowaox,ipowaph,ipowasi,ipown2,ipowno3,isssc12,      &
                              issso12,issssil,issster,ks,nsedtra,ipowc13,ipowc13,issso13,issso13,  &
                              isssc13,ipowc14,isssc14,issso14
    use mo_vgrid,       only: kmle,kbo
    use mo_carbch,      only: nathi,natco3
    use mo_sedmnt,      only: sedhpl,burial,powtra,sedlay
    use mo_profile_gd,  only: profile_gd

    ! Arguments
    integer, intent(in) :: kpaufr                                   ! 1/0 flag, 1 indicating a restart run
    integer, intent(in) :: kpie                                     ! 1st dimension of model grid.
    integer, intent(in) :: kpje                                     ! 2nd dimension of model grid.
    integer, intent(in) :: kpke                                     ! 3rd (vertical) dimension of model grid.
    integer, intent(in) :: kbnd                                     ! nb of halo grid points
    real,    intent(in) :: pddpo(kpie,kpje,kpke)                    ! size of grid cell (3rd dimension) [m].
    real,    intent(in) :: prho (kpie,kpje,kpke)                    ! density [g/cm^3].
    real,    intent(in) :: omask(kpie,kpje)                         ! ocean mask.
    real,    intent(in) :: pglon(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd) ! longitude of grid cell [deg].
    real,    intent(in) :: pglat(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd) ! latitude  of grid cell [deg].

    ! local variables
    integer :: i,j,k,l
    real    :: rco213,rco214,beta13,beta14 ! cisonew
    real    :: snow ! AGG

    if (use_FB_BGC_OCE) then
      do k=1,kpke
        do j=1,kpje
          do i=1,kpie
            abs_oce(i,j,k)=1.
          enddo
        enddo
      enddo
    endif
    !
    ! Initialisation of ocean tracers and sediment
    !
    ! Initialise ocean tracers with WOA and GLODAP data. This is done even in case
    ! of a restart since some tracers (e.g. C-isotopes) might not be in the restart
    ! file and aufr.f90 instead expects an initialised field.
    call profile_gd(kpie,kpje,kpke,kbnd,pglon,pglat,omask)

    ! If this is a restart run initialisation is done in aufr.F90
    if (kpaufr==1) RETURN  !jt

    do k=1,kpke
      do j=1,kpje
        do i=1,kpie
          if (omask(i,j) > 0.5 ) then
            ! convert WOA tracers kmol/m^3 -> mol/kg; GLODAP dic and alk
            ! are already in mol/kg. We need these units here, since after
            ! initialisation the tracer field is passed to the ocean model
            ! first where units are mol/kg.
            ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph)/prho(i,j,k)
            ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen)/prho(i,j,k)
            ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)  /prho(i,j,k)
            ocetra(i,j,k,isilica) = ocetra(i,j,k,isilica)/prho(i,j,k)
            if (use_cisonew) then
              ! d13C based on Eide data is read in above (profile_gd)
              ! Convert to 13C using model initial (ie GLODAP) total C
              ! If restarting, this is redone with model total C from restart in aufr_bgc.F90
              beta13=ocetra(i,j,k,isco213)/1000.+1.
              ocetra(i,j,k,isco213) = ocetra(i,j,k,isco212)*beta13*re1312/(1.+beta13*re1312)

              ! 14C is read in as small delta14C (calculated from R. Key, 2003 and Eide et al. 2017)
              ! Convert to 14C using model total C, and normalize by c14fac to prevent numerical errors
              beta14=ocetra(i,j,k,isco214)/1000.+1.
              ocetra(i,j,k,isco214) = ocetra(i,j,k,isco212)*beta14*re14to/c14fac
            endif
          endif
        enddo
      enddo
    enddo

    ! Initialise remaining ocean tracers
    do k=1,kpke
      do j=1,kpje
        do i=1,kpie
          if (omask(i,j) > 0.5) then
            ocetra(i,j,k,igasnit)=1.e-10
            ocetra(i,j,k,idoc)   =1.e-8
            ocetra(i,j,k,iphy)   =1.e-8
            ocetra(i,j,k,izoo)   =1.e-8
            ocetra(i,j,k,idet)   =1.e-8
            ocetra(i,j,k,icalc)  =0.
            ocetra(i,j,k,iopal)  =1.e-8
            ocetra(i,j,k,ian2o)  =0.
            ocetra(i,j,k,idms)   =0.
            ocetra(i,j,k,ifdust) =0.
            ocetra(i,j,k,iiron)  =fesoly
            ocetra(i,j,k,iprefo2)=0.
            ocetra(i,j,k,iprefpo4)=0.
            ocetra(i,j,k,iprefalk)=0.
            ocetra(i,j,k,iprefdic)=0.
            ocetra(i,j,k,idicsat)=1.e-8
            hi(i,j,k)            =1.e-8
            co3(i,j,k)           =0.
            co2star(i,j,k)       =20.e-6
            if (use_AGG) then
              ! calculate initial numbers from mass, to start with appropriate size distribution
              snow = (ocetra(i,j,k,iphy)+ocetra(i,j,k,idet))*1.e+6
              ocetra(i,j,k,inos)   = snow / cellmass / (FractDim+1.)
              ocetra(i,j,k,iadust) =0.
            endif
            if (use_CFC) then
              ocetra(i,j,k,icfc11)   =0.
              ocetra(i,j,k,icfc12)   =0.
              ocetra(i,j,k,isf6)     =0.
            endif
            if (use_natDIC) then
              nathi(i,j,k)           =1.e-8
              natco3(i,j,k)          =0.
              ocetra(i,j,k,inatcalc) =0.
            endif
            if (use_cisonew) then
              rco213=ocetra(i,j,k,isco213)/(ocetra(i,j,k,isco212)+safediv)
              rco214=ocetra(i,j,k,isco214)/(ocetra(i,j,k,isco212)+safediv)
              ocetra(i,j,k,iphy13) =ocetra(i,j,k,iphy)*rco213*bifr13_ini
              ocetra(i,j,k,iphy14) =ocetra(i,j,k,iphy)*rco214*bifr14_ini
              ocetra(i,j,k,izoo13) =ocetra(i,j,k,izoo)*rco213*bifr13_ini
              ocetra(i,j,k,izoo14) =ocetra(i,j,k,izoo)*rco214*bifr14_ini
              ocetra(i,j,k,idoc13) =ocetra(i,j,k,idoc)*rco213*bifr13_ini
              ocetra(i,j,k,idoc14) =ocetra(i,j,k,idoc)*rco214*bifr14_ini
              ocetra(i,j,k,idet13) =ocetra(i,j,k,idet)*rco213*bifr13_ini
              ocetra(i,j,k,idet14) =ocetra(i,j,k,idet)*rco214*bifr14_ini
              ocetra(i,j,k,icalc13)=ocetra(i,j,k,icalc)*rco213
              ocetra(i,j,k,icalc14)=ocetra(i,j,k,icalc)*rco214
            endif
            if (use_BROMO) then
              ! Initialise to 0,01 pmol L-1 (Stemmler et al., 2015) => mol/kg
              ocetra(i,j,k,ibromo)= 1.e-14/prho(i,j,k)
            endif
            if (use_dom) then
              ocetra(i,j,k,idocsl)   =1.e-8
              ocetra(i,j,k,idocsr)   =1.e-8
              ocetra(i,j,k,idocr )   =1.e-8
              ocetra(i,j,k,iprefdoc) =0.
              ocetra(i,j,k,iprefdocsl) =0.
              ocetra(i,j,k,iprefdocsr) =0.
              ocetra(i,j,k,iprefdocr) =0.
            endif
          endif ! omask > 0.5
        enddo
      enddo
    enddo

    ! Initialise preformed tracers in the mixed layer; note that the
    ! whole field has been initialised to zero above
    do j=1,kpje
      do i=1,kpie
        if (omask(i,j) > 0.5) then
          ocetra(i,j,1:kmle(i,j),iprefo2)  = ocetra(i,j,1:kmle(i,j),ioxygen)
          ocetra(i,j,1:kmle(i,j),iprefpo4) = ocetra(i,j,1:kmle(i,j),iphosph)
          ocetra(i,j,1:kmle(i,j),iprefalk) = ocetra(i,j,1:kmle(i,j),ialkali)
          ocetra(i,j,1:kmle(i,j),iprefdic) = ocetra(i,j,1:kmle(i,j),isco212)
          if (use_dom) then
            ocetra(i,j,1:kmle(i,j),iprefdoc) = ocetra(i,j,1:kmle(i,j),idoc)
            ocetra(i,j,1:kmle(i,j),iprefdocsl) = ocetra(i,j,1:kmle(i,j),idocsl)
            ocetra(i,j,1:kmle(i,j),iprefdocsr) = ocetra(i,j,1:kmle(i,j),idocsr)
            ocetra(i,j,1:kmle(i,j),iprefdocr) = ocetra(i,j,1:kmle(i,j),idocr)
          endif
        endif
      enddo
    enddo


    ! Initial values for sediment
    if (.not. use_sedbypass) then
      do  k=1,ks
        do  j=1,kpje
          do  i=1,kpie
            if (omask(i,j) > 0.5) then
              powtra(i,j,k,ipowaic)=ocetra(i,j,kbo(i,j),isco212)
              powtra(i,j,k,ipowaal)=ocetra(i,j,kbo(i,j),ialkali)
              powtra(i,j,k,ipowaph)=ocetra(i,j,kbo(i,j),iphosph)
              powtra(i,j,k,ipowaox)=ocetra(i,j,kbo(i,j),ioxygen)
              powtra(i,j,k,ipown2) =0.
              powtra(i,j,k,ipowno3)=ocetra(i,j,kbo(i,j),iano3)
              powtra(i,j,k,ipowasi)=ocetra(i,j,kbo(i,j),isilica)
              sedlay(i,j,k,issso12)=1.e-8
              sedlay(i,j,k,isssc12)=1.e-8
              sedlay(i,j,k,issster)=30.
              sedlay(i,j,k,issssil)=1.e-8
              sedhpl(i,j,k)        =hi(i,j,kbo(i,j))
              if (use_cisonew) then
                rco213=ocetra(i,j,kbo(i,j),isco213)/(ocetra(i,j,kbo(i,j),isco212)+safediv)
                rco214=ocetra(i,j,kbo(i,j),isco214)/(ocetra(i,j,kbo(i,j),isco212)+safediv)
                powtra(i,j,k,ipowc13)=powtra(i,j,k,ipowaic)*rco213*bifr13_ini
                powtra(i,j,k,ipowc14)=powtra(i,j,k,ipowaic)*rco214*bifr14_ini
                sedlay(i,j,k,issso13)=sedlay(i,j,k,issso12)*rco213*bifr13_ini
                sedlay(i,j,k,issso14)=sedlay(i,j,k,issso12)*rco214*bifr14_ini
                sedlay(i,j,k,isssc13)=sedlay(i,j,k,isssc12)*rco213
                sedlay(i,j,k,isssc14)=sedlay(i,j,k,isssc12)*rco214
              endif
            else
              powtra(i,j,k,ipowno3)=rmasks
              powtra(i,j,k,ipown2) =rmasks
              powtra(i,j,k,ipowaic)=rmasks
              powtra(i,j,k,ipowaal)=rmasks
              powtra(i,j,k,ipowaph)=rmasks
              powtra(i,j,k,ipowaox)=rmasks
              powtra(i,j,k,ipowasi)=rmasks
              sedlay(i,j,k,issso12)=rmasks
              sedlay(i,j,k,isssc12)=rmasks
              sedlay(i,j,k,issssil)=rmasks
              sedlay(i,j,k,issster)=rmasks
              sedlay(i,j,k,issssil)=rmasks
              sedhpl(i,j,k)        =rmasks
              if (use_cisonew) then
                powtra(i,j,k,ipowc13)=rmasks
                powtra(i,j,k,ipowc14)=rmasks
                sedlay(i,j,k,issso13)=rmasks
                sedlay(i,j,k,issso14)=rmasks
                sedlay(i,j,k,isssc13)=rmasks
                sedlay(i,j,k,isssc14)=rmasks
              endif
            endif
          enddo
        enddo
      enddo

      ! last and final sediment layer
      do  l=1,nsedtra
        do  j=1,kpje
          do  i=1,kpie
            burial(i,j,l)=0.
          enddo
        enddo
      enddo
    endif

  end subroutine ini_fields_ocean

end module mo_ini_fields
