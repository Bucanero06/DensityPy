!>  Insert Doxygen comments here
!!  
!! 
program TestAngularMomentum

  use, intrinsic :: ISO_FORTRAN_ENV

  use ModuleMainInterface
  use ModuleString
  use ModuleAngularMomentum

  implicit none

  integer :: lmax


  call GetRunTimeParameters( lmax )
  
!!$  call TestClebschGordanCoefInteger    ( lmax )
!!$  call TestClebschGordanCoefHalfInteger( lmax )
!!$  call TestSixJSymbolsHalfInteger( lmax )


!!$  call CheckClebschGordanContractions_1st1B( lmax )
!!$  call CheckClebschGordanContractions_3rd1B( lmax )
!!$  call CheckClebschGordanContractions_1st2B( lmax )
!!$  call CheckClebschGordanContractions_2nd2B( lmax )
!!$  call CheckClebschGordanContractions_3rd2B( lmax )
!!$  call CheckClebschGordanContractions_4th2B( lmax )
!!$  call CheckClebschGordanContractions_6th2B( lmax )

!!$  call Check6jSpecialFormulas( lmax )

!!$  call CheckClebschGordanFormulas()
  
end program TestAngularMomentum
