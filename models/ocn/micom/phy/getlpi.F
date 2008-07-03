      real function getlpi(temp,saln,phiu,phil,pu)
c
c --- get lower pressure interface of a layer knowing the temperature,
c --- salinity of the layer and the geopotential at upper and lower
c --- interface
c
      implicit none
c
      real temp,saln,phiu,phil,pu
c
      real pl,q,dphi,alpu,alpl
c
      real rho
      external rho
c
c --- first guess on pressure interface
      pl=pu-rho(pu,temp,saln)*(phil-phiu)
c
c --- improve the accuracy of the pressure interface by an
c --- iterative procedure
      q=1.
      do while (abs(q).gt.1.e-4)
        call delphi(pu,pl,temp,saln,dphi,alpu,alpl)
        q=(phil-phiu-dphi)/alpl
        pl=pl-q
      enddo
c
      getlpi=pl
c
      return
      end
