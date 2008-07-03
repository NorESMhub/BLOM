      subroutine chksummsk(a,msk,kcsd,text)
c
c --- compute checksum of model field after multiplying with mask
c
      use mod_xc
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kcsd) :: a
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: msk
      integer kcsd
      character text*(*)
c
      real, dimension(itdm,jtdm,kcsd) :: aa
      real, dimension(itdm,jtdm) :: rrmsk
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: rmsk
      integer ics,jcs,kcs
c
      integer crcfast
      external crcfast
c
      do jcs=1,jj
        do ics=1,ii
          rmsk(ics,jcs)=msk(ics,jcs)
        enddo
      enddo
c
      do kcs=1,kcsd
        call xcaget(aa(1,1,kcs),a(1-nbdy,1-nbdy,kcs),1)
      enddo
      call xcaget(rrmsk,rmsk,1)
c
      do kcs=1,kcsd
        do jcs=1,jtdm
          do ics=1,itdm
            if (rrmsk(ics,jcs).lt..5) then
              aa(ics,jcs,kcs)=0.
            endif
          enddo
        enddo
      enddo
c
      if (mnproc.eq.1) then
        write (lp,'(3a,z8.8)') ' chksum: ',text,': 0x',
     .    crcfast(aa,itdm*jtdm*kcsd*8)
      endif
c
      return
      end
