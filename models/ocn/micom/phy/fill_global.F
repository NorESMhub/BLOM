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
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: tmp
      real sum,done,changed
      integer i,j,mbdy,num,il,jl
c
      mbdy=1
      done=0.
      changed=1.
c
      do while (changed.eq.1.)
c
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(field, 1,1, nbdy,nbdy, itype)
        else
          mbdy=mbdy-1
        endif
c
        done=1.
        changed=0.
c
c$OMP PARALLEL DO PRIVATE(sum,num,jl,il)
c$OMP+REDUCTION(min:done) REDUCTION(max:changed)
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
                changed=max(1.,changed)
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
        call xcmaxr(changed)
c
      enddo
c
      if (done.eq.0.) then
        if (mnproc.eq.1) then
          write(lp,*) 'fill_global: filling failed!'
        endif
        return
      endif
c
      return
      end
