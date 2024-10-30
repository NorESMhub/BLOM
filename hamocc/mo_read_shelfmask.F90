module mo_read_shelfmask

  implicit none
  private

  public :: read_shelfmask

  character(len=512),               public :: shelfsea_maskfile =''
  logical, allocatable,  protected, public :: shelfmask(:,:)

contains

  subroutine read_shelfmask(kpie,kpje,kbnd,pbath,omask)

    use mod_xc,         only: mnproc,xchalt
    use mod_dia,        only: iotype
    use mod_nctools,    only: ncfopn,ncread,ncfcls
    use mo_control_bgc, only: use_shelfsea_res_time,io_stdo_bgc
    use mo_param_bgc,   only: shelfbreak_depth

    !Arguments
    integer, intent(in) :: kpie,kpje,kbnd
    real,    intent(in) :: pbath(1-kbnd:kpie+kbnd,1-kbnd:kpje+kbnd) ! bathymetry fields - water depth [m]
    real,    intent(in) :: omask(kpie,kpje)                         ! land/ocean mask.

    ! Local variables
    logical :: file_exists=.false.
    integer :: i,j,errstat,dummymask(2)
    real    :: mask(kpie,kpje)

    mask = 0.

    if (.not.use_shelfsea_res_time) then
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,*) '******************************************************'
        write(io_stdo_bgc,*) 'read_shelfmask: shelf sea age tracer is not activated.'
      endif
      return
    endif

    ! Check if shelfsea mask file exists. If not, abort.
    inquire(file=shelfsea_maskfile,exist=file_exists)
    if (.not. file_exists .and. mnproc.eq.1) then
      write(io_stdo_bgc,*) ''
      write(io_stdo_bgc,*) '**********************************************************'
      write(io_stdo_bgc,*) 'read_shelfmask: Cannot find shelf sea region mask file... '
      write(io_stdo_bgc,*) '... using internal bathymetry data to reconstruct the mask'
    endif

    ! Allocate field to hold shelfsea mask
    if(mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variable shelfmask ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
    endif
    allocate(shelfmask(kpie,kpje),stat=errstat)
    allocate(mask(kpie,kpje),stat=errstat)
    if(errstat.ne.0) stop 'not enough memory shelfmask'
    shelfmask(:,:) = .false.

    if (file_exists) then
      ! read shelf sea mask from file
      if (mnproc.eq.1) then
        write(io_stdo_bgc,*) ''
        write(io_stdo_bgc,'(a)') 'read_shelfmask: read mask from ',trim(shelfsea_maskfile)
      endif
      call ncfopn(trim(shelfsea_maskfile),'r',' ',1,iotype)
      call ncread('shelfmask',mask,dummymask,0,0.)
      call ncfcls
    else
      ! reconstruct shelf sea mask from internal bathymetry
      !$OMP DO PARALLEL PRIVATE (i,j)
      do j=1,kpje
        do i=1,kpie
          if((omask(i,j) > 0.5) .and. (pbath(i,j) <= shelfbreak_depth)) then
            mask(i,j) = 1.
          endif
        enddo
      enddo
      !$OMP END PARALLEL DO
    endif

    !$OMP DO PARALLEL PRIVATE (i,j)
    do j = 1,kpje
      do i = 1,kpie
        if (1 == nint(mask(i,j))) then
          shelfmask(i,j) = .true.
        else
          shelfmask(i,j) = .false.
        endif
      enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine read_shelfmask


end module mo_read_shelfmask
