module mod_vertinterp

   use mod_types,     only: r8
   use mod_constants, only: epsilp
   use mod_xc,        only: idm, jdm, ii, jj, kk, nbdy, ifp, ilp, isp
   use mod_grid,      only: depths
   use mod_state,     only: dp
   use mod_time,      only: nstep, baclin

   implicit none
   private

   public :: vertinterp_wghts
   public :: vertinterp_accum

   real(r8), allocatable :: dlevp(:,:,:)
   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) :: ztop
   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kk) :: zbot
   real(r8), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)    :: pbath

contains

   ! ---------------------------------------------------------------------------
   subroutine vertinterp_wghts(k, nlevs, vlevs, vlevs_bnds, ind1, ind2, wghts)

      ! Arguments
      integer , intent(in) :: k
      integer , intent(in) :: nlevs
      real(r8), intent(in) :: vlevs(nlevs)
      real(r8), intent(in) :: vlevs_bnds(2,nlevs)
      integer , intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)     :: ind1
      integer , intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)     :: ind2
      real    , intent(out), dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlevs) :: wghts

      ! Local variables
      integer  :: m,mm
      integer  :: d,i,j,l,kl,km,kml,k1m
      real(r8) :: r,dzeps,dpeps
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

         allocate (dlevp(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlevs))
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

      ! Compute interpolation wghts
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
                     wghts(i,j,d) = ( min(zbot(i,j,k),vlevs_bnds(2,d))   &
                          -max(ztop(i,j,k),vlevs_bnds(1,d)) ) &
                          /dlevp(i,j,d)
                  end do
               end if
            end do
         end do
      end do

   end subroutine vertinterp_wghts

   ! ---------------------------------------------------------------------------
   subroutine vertinterp_accum(nlev_in, nlev_out, k, ind1, ind2, wghts, fld_in, fld_out)

      !---------------------------------------------------------------
      ! Description: accumulate layer fields mapped to levels
      !
      ! Arguments:
      !   int  nlev_in  (in)  : number of levels of input field
      !   int  nlev_out (in)  : number of levels of output field
      !   int  k        (in)  : layer index of input fld
      !   int  ind1     (in)  : index field for first accumulated level
      !   int  ind2     (in)  : index field for last accumulated level
      !   real wghts    (in)  : weights used for accumulation
      !   real fld_in   (in)  : input data used for accumulation
      !   real fld_out  (out) : output mapped and accumulated data
      !---------------------------------------------------------------

      integer  , intent(in) :: nlev_in
      integer  , intent(in) :: nlev_out
      integer  , intent(in) :: k
      integer  , intent(in)  , dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)          :: ind1
      integer  , intent(in)  , dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy)          :: ind2
      real(r8) , intent(in)  , dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlev_out) :: wghts
      real(r8) , intent(in)  , dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlev_in)  :: fld_in
      real(r8) , intent(out) , dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlev_out) :: fld_out

      ! Local variables
      integer :: d,i,j,l

      !$omp parallel do private(l,i,d)
      do j = 1,jj
         do l = 1,isp(j)
            do i = max(1,ifp(j,l)),min(ii,ilp(j,l))
               do d = ind1(i,j),ind2(i,j)
                  fld_out(i,j,d) = fld_out(i,j,d) + fld_in(i,j,k)*wghts(i,j,d)*baclin
               end do
            end do
         end do
      end do
      !$omp end parallel do

   end subroutine vertinterp_accum

end module mod_vertinterp
