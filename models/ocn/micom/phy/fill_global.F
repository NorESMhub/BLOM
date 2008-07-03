      subroutine fill_global(missing_value,fill_value,itype,field)
c
c --- ------------------------------------------------------------------
c --- Fill missing values by extrapolating values from neighbouring
c --- points. A value = fill_value will be modified.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      real missing_value,fill_value
      integer itype
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: field
c
      integer niter_max
      parameter (niter_max=100)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp
      real sum,done
      integer i,j,mbdy,niter,num,il,jl
c
      mbdy=1
      done=0.
      niter=0
c
      do while (done.eq.0.)
c
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(field, 1,1, nbdy,nbdy, itype)
        else
          mbdy=mbdy-1
        endif
c
        done=1.
        niter=niter+1
        if (niter.eq.niter_max) then
          if (mnproc.eq.1) then
            write(lp,*) 'fill_global: filling failed!'
          endif
          return
        endif
c
c$OMP PARALLEL DO PRIVATE(sum,num,jl,il) REDUCTION(min:done)
        do j=1-mbdy+1,jj+mbdy-1
          do i=1-mbdy+1,ii+mbdy-1
            if (abs(field(i,j)).eq.abs(fill_value)) then
              done=min(0.,done)
              sum=0.
              num=0
              do jl=j-1,j+1
                do il=i-1,i+1
                  if (abs(field(il,jl)).ne.abs(fill_value).and.
     .                abs(field(il,jl)).ne.abs(missing_value)) then
                    sum=sum+field(il,jl)
                    num=num+1
                  endif
                enddo
              enddo
              if (num.gt.0) then
                tmp(i,j)=sum/real(num)
              else
                tmp(i,j)=field(i,j)
              endif
            else
              tmp(i,j)=field(i,j)
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
        do j=1-mbdy+1,jj+mbdy-1
          do i=1-mbdy+1,ii+mbdy-1
            field(i,j)=tmp(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call xcminr(done)
c
      enddo
c
      return
      end
