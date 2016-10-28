c --- ------------------------------------------------------------------
c --- common blocks related to the model grids geographical coordinates
c --- ------------------------------------------------------------------
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,4) ::
     .  qclat,qclon,   ! grid coordinates q-cell corners
     .  pclat,pclon,   ! grid coordinates p-cell corners
     .  uclat,uclon,   ! grid coordinates u-cell corners
     .  vclat,vclon    ! grid coordinates v-cell corners
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  qlat,qlon,     ! grid coordinates q-points
     .  plat,plon,     ! grid coordinates p-points
     .  ulat,ulon,     ! grid coordinates u-points
     .  vlat,vlon,     ! grid coordinates v-points
     .  angle,         ! local angle of i-direction and eastward
                       ! direction at p-points
     .  cosang,sinang  ! cosine and sine of local angle of i-direction
                       ! and eastward direction at p-points

      integer ::
     .  nwp            ! number of wet grid cells
c
      common /geo/ qclat,qclon,pclat,pclon,uclat,uclon,vclat,vclon,
     .             qlat,qlon,plat,plon,ulat,ulon,vlat,vlon,angle,
     .             cosang,sinang,nwp
