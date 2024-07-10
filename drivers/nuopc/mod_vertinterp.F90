module mod_vertinterp

  use mod_types,     only: r8
  use mod_constants, only: epsil
  use mod_xc,        only: idm, jdm, jj, kk, nbdy, ifp, ilp, isp
  use mod_grid,      only: depths
  use mod_state,     only: dp
  use mod_time,      only: nstep

  implicit none
  public

contains

  subroutine vertinterp_weights(k, nlevs, vlevs, vlevs_bnds, ind1, ind2, weights)

    ! Arguments
    integer , intent(in) :: k
    integer , intent(in) :: nlevs
    real(r8), intent(in) :: vlevs(nlevs)
    real(r8), intent(in) :: vlevs_bnds(2,nlevs)
    integer , intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)     :: ind1
    integer , intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)     :: ind2
    real    , intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlevs) :: weights

    ! Local variables
    integer  :: m,mm
    integer  :: d,i,j,l,kl,km,kml,k1m
    real(r8) :: r,dzeps,dpeps
    real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk)    , save :: ztop
    real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk)    , save :: zbot
    real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlevs) , save :: dlevp
    real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)       , save :: pbath
    logical  :: iniflg = .true.

    ! Define thresholds
    dzeps = 1e1*epsilp
    dpeps = 1e5*epsilp

    ! Sort out stuff related to time stepping
    m = mod(nstep,2) + 1
    mm = (m-1)*kk
    km  = k+mm
    k1m = 1+mm

    ! Adjust bounds of levitus levels according to model bathymetry
    if (iniflg) then
      do j = 1,jj+1
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            pbath(i,j) = depths(i,j)
          end do
        end do
      end do
      do j = 1,jj+1
        do d = 1,nlevs
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              dlevp(i,j,d) = max(dzeps,min(pbath(i,j),vlevs_bnds(2,d))-vlevs_bnds(1,d))
            end do
          end do
        end do
      end do

      iniflg = .false.
    end if

    ! Compute top and bottom depths of density layers
    if (k == 1) then
      do j = 1,jj
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            zbot(i,j,1) = dp(i,j,k1m)
          end do
        end do
        do kl = 2,kk
          kml = kl+mm
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              zbot(i,j,kl) = zbot(i,j,kl-1) + dp(i,j,kml)
            end do
          end do
        end do
        do l = 1,isp(j)
          do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
            zbot(i,j,1) = zbot(i,j,1)*pbath(i,j)/zbot(i,j,kk)
            ztop(i,j,1) = 0.
            ind1(i,j) = 1
          end do
        end do
        do kl = 2,kk
          do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
              zbot(i,j,kl) = zbot(i,j,kl)*pbath(i,j)/zbot(i,j,kk)
              ztop(i,j,kl) = zbot(i,j,kl-1)
            end do
          end do
        end do
      end do
    end if

    ! Compute interpolation weights
    do j = 1,jj
      do l = 1,isp(j)
        do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
          ind2(i,j) = 0
          if (dp(i,j,km) > dpeps) then
            do d = ind1(i,j),nlevs
              if (vlevs_bnds(2,d) <= ztop(i,j,k)) then
                ind1(i,j) = d+1
                cycle
              else if (vlevs_bnds(1,d) >= zbot(i,j,k)) then
                exit
              end if
              ind2(i,j) = d
              weights(i,j,d) = ( min(zbot(i,j,k),vlevs_bnds(2,d))   &
                                -max(ztop(i,j,k),vlevs_bnds(1,d)) ) &
                                /dlevp(i,j,d)
            end do
          end if
        end do
      end do
    end do

  end subroutine vertinterp_weights

end module mod_vertinterp
