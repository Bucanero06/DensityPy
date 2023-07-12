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
!> \file
!!
!! Sets several constants that are used throughout the programs.

module ModuleConstants
  !> Sets the kind of the real or complex parameters to be used.
  integer        , parameter :: DOUBLE= kind(1d0)
  !> Greek \f$\pi\f$.
  real   (DOUBLE), parameter :: PI    = 3.1415926535897932d0
  !> Inverse of fine-structure constant.
  real   (DOUBLE), parameter :: CAU   = 137.03599970d0
  !> Fine-structure constant.
  real   (DOUBLE), parameter :: ALPHA = 0.0072973525664d0
  !> Zero complex number.
  complex(DOUBLE), parameter :: Z0    = (0.d0,0.d0)
  !> Real unit in complex number format.
  complex(DOUBLE), parameter :: Z1    = (1.d0,0.d0)
  !> Imaginary unit in complex number format.
  complex(DOUBLE), parameter :: Zi    = (0.d0,1.d0)
  !> Natural logarithm of 2
  real   (kind(1d0)), parameter  :: LN2 = 0.6931471805599453d0
  !> Greek \f$1/\sqrt(2\pi)\f$.
  real   (DOUBLE), parameter :: ONE_OVER_SQRT_2PI = 1.d0/sqrt(2.d0*3.1415926535897932d0)
  !> Bohr Radius in meters
  real   (DOUBLE), parameter :: BOHR_RADIUS = 5.2917721067d-11
  !> Boltzmann Constant in eV/K
  real   (DOUBLE), parameter :: BOLTZMANN_CONSTANT_eVK = 8.617333262d-5


  ! Conversion factors.

  !> Time conversion factor from a.u. to SI (seconds).
  real(kind(1d0)), parameter :: Energy_au_to_eV=27.21138602d0
  !> Time conversion factor from a.u. to SI (seconds).
  real(kind(1d0)), parameter :: TimeConvAUtoSI=2.418884326505d-17
  !> Time conversion factor from femtoseconds to a.u.
  real(kind(1d0)), parameter :: TimeConvFStoAU=1.d-15/TimeConvAUtoSI
  !> Energy conversion factor from a.u. to SI (Joules).
  real(kind(1d0)), parameter :: EnergyConvAUtoSI=4.35974417d-18
  !> Length conversion factor from a.u. to SI (meter).
  real(kind(1d0)), parameter :: LengthConvAUtoSI=5.291772192d-11
  !> Intensity conversion factor from PW/cm² to a.u.
  real(kind(1d0)), parameter :: IntensityConvPWCM2toAU=1d15/EnergyConvAUtoSI*TimeConvAUtoSI*(LengthConvAUtoSI)**2/1d-4
  !> Cross Section conversion factor from au (Bohr Radii square) to Mega Barn (1Mb=10^{-22}m^2)
  real(kind(1d0)), parameter :: RBOHR2_AU_to_MEGABARN = BOHR_RADIUS**2 / 1.d-22
  !> Boltzmann Constant in au/K
  real   (DOUBLE), parameter :: BOLTZMANN_CONSTANT_auK = BOLTZMANN_CONSTANT_eVK / Energy_au_to_eV


end module ModuleConstants
