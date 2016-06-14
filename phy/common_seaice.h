c --- ------------------------------------------------------------------
c --- common blocks for sea ice variables
c --- ------------------------------------------------------------------
c   
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  ficem,hicem,hsnwm,ustari,tauxice,tauyice,uicem,vicem,iagem
c 
      common /sivar/ ficem,hicem,hsnwm,ustari,tauxice,tauyice,uicem,
     .               vicem,iagem
