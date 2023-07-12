! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)                               
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
MODULE phys_const
! definitions of physical constants
 use precisn

 implicit none
 
 !complex unit
 complex(kind=cfp), parameter :: imu=(0.0_cfp,1.0_cfp)

 !the number of digits in pi,etc. corresponds to quad precision, so these numbers can be used to define their quad precision equivalents
 real(kind=cfp), parameter :: pi=3.14159265358979323846264338327950288_cfp
 real(kind=cfp), parameter :: twopi=6.28318530717958647692528676655900577_cfp
 real(kind=cfp), parameter :: fourpi=12.5663706143591729538505735331180115_cfp
 real(kind=cfp), parameter :: inv_fourpi=1.0_cfp/fourpi
 real(kind=cfp), parameter :: pi5o2= 17.4934183276248628462628216798715538_cfp !pi^{5/2}
 real(kind=cfp), parameter :: pi54 = pi**1.25_cfp !pi^{5/4}
 real(kind=cfp), parameter :: rtpi = sqrt(pi)
 real(kind=cfp), parameter :: rtpi_half = rtpi*0.5_cfp

 real(kind=ep1), parameter :: qtwopi=6.28318530717958647692528676655900577_ep1

 !again, quad prec values
 real(kind=cfp), parameter :: roneh = 0.707106781186547524400844362104849039_cfp !sqrt(0.5)
 real(kind=cfp), parameter :: rtwo = 1.41421356237309504880168872420969808_cfp !sqrt(2)

 real(kind=cfp), parameter :: to_au=1.88972599_cfp !Conversion constant from angstroms to a.u. This is the value which Molpro uses when writing the Molden file.

 real(kind=cfp), parameter :: to_ev=27.21138469_cfp !Conversion constant from H to eV.

END MODULE phys_const
