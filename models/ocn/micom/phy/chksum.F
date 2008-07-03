      subroutine chksum(a,kcsd,text)
c
c --- compute checksum of model field
c
      use mod_xc
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kcsd) :: a
      integer kcsd
      character text*(*)
c
      real, dimension(itdm,jtdm,kcsd) :: aa
      integer kcs
c
      integer crcfast
      external crcfast
c
      do kcs=1,kcsd
        call xcaget(aa(1,1,kcs),a(1-nbdy,1-nbdy,kcs),1)
      enddo
c
      if (mnproc.eq.1) then
        write (lp,'(3a,z8.8)') ' chksum: ',text,': 0x',
     .    crcfast(aa,itdm*jtdm*kcsd*8)
      endif
c
      return
      end
