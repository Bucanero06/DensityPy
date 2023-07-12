!! CONFIDENTIAL
!! Copyright (C) 2020 Luca Argenti, PhD - All Rights Reserved
!! email: luca.argenti@gmail.com
!! email: luca.argenti@ucf.edu
!! Luca Argenti is Associate Professor of Physics, Optics and Photonics
!! at the Department of Physics and the College of Optics
!! of the University of Central Florida
!! 4111 Libra Drive
!! Orlando, Florida, USA
!!
Module ModuleAstraCredit

  use, intrinsic :: ISO_FORTRAN_ENV
  implicit none
  private

  public :: PrintAstraName, PrintAstraCredit
  
contains

  subroutine PrintASTRAName()
    use ModuleAstraEnv
    call Execute_Command_Line("cat "//AstraEnv%GetLogoFile())
  end subroutine PrintASTRAName

  subroutine PrintAstraCredit()
    write(OUTPUT_UNIT,"(a)") "  A Close-Coupling program for single and double molecular ionization"
    write(OUTPUT_UNIT,"(a)") "  by Luca Argenti, Ph.D, University of Central Florida, FL, USA"
    write(OUTPUT_UNIT,"(a)") "  with the collaboration of"
    write(OUTPUT_UNIT,"(a)") "   - Jeppe Olsen, Ph.D, University of Aahrus, Denmark, EU"
    write(OUTPUT_UNIT,"(a)") "   - Barry Schneider, Ph.D, NIST, Gaithersburg, MD, USA"
    write(OUTPUT_UNIT,"(a)") "   - Juan Martin Randazzo, Ph.D, University of Central Florida, FL, USA"
    write(OUTPUT_UNIT,"(a)") "   - Ruben A. Fernandez Carbon,  University of Central Florida, FL, USA"
    write(OUTPUT_UNIT,"(a)") "   - ..."
    write(OUTPUT_UNIT,"(a)") "  with the support of the DOE Career Starting Grant DE-SC0020311"
    write(OUTPUT_UNIT,"(a)") "  'New correlated numerical methods for attosecond molecular "
    write(OUTPUT_UNIT,"(a)") "   single and double ionization'"
    write(OUTPUT_UNIT,"(a)")
  end subroutine PrintAstraCredit

end Module ModuleAstraCredit
