module kinds

!******************************************************************************
!
! Fortran module for setting data types.
! Author - Jacob A. Harvey
! Copyright - CRUNCH Lab Feb 2012 
!
!*******************************************************************************

  implicit none

  integer, parameter :: dp = selected_real_kind(14,300)
  integer, parameter :: ip = selected_int_kind(12)

! define vectors

  type vector
       real :: x , y , z
  end type vector

End Module kinds
