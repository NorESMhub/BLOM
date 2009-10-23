MODULE mo_kind

  ! L. Kornblueh, MPI, August 2001, added working precision and comments 

  IMPLICIT NONE

  ! Number model from which the SELECTED_*_KIND are requested:
  !
  !                   4 byte REAL      8 byte REAL
  !          CRAY:        -            precision =   13
  !                                    exponent  = 2465
  !          IEEE:    precision =  6   precision =   15  
  !                   exponent  = 37   exponent  =  307 
  !
  ! Most likely this are the only possible models.

  ! Floating point section 

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,37)  
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307)

  INTEGER, PARAMETER :: wp = dp   ! working precision

  ! Integer section

  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)

END MODULE mo_kind
