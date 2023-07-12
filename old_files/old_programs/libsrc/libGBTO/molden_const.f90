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
module molden_const
use precisn

   !> Conversion constant from angstroms to a.u. This is the value which Molpro uses when writing the molden file.
   !> \todo move this constant to phys_const module?
   real(kind=cfp), parameter :: angstrom_to_au=1.889726131421507_cfp !1.88972599_cfp

   !> Threshold for distance of the CMS from the origin which will trigger recentering of the molecule to the actual center of mass in molden_mod. 
   real(kind=cfp), parameter :: cms_thrs = 10e-10_cfp

   !> Maximum GTO L supported by the Molden format.
   integer, parameter :: max_molden_l = 4

   !> Supported (by the Molden format) angular types of GTO basis functions.
   character(len=1), parameter :: gto_typ(0:max_molden_l) = (/'s','p','d','f','g'/)

   !> Length of the longest character string from the array ang_fact below.
   integer, parameter, private :: len_ang_fact = 4

   !> Definitions of the angular parts of the cartesian gaussians. The order of the angular factors is given by the Molden format.
   character(len=len_ang_fact), protected :: ang_fact(0:max_molden_l,1:15)

   data ang_fact(0,1)      /'s'/
   data ang_fact(1:1,1:3)  /'x','y','z'/
   data ang_fact(2:2,1:6)  /'xx','yy','zz','xy','xz','yz'/
   data ang_fact(3:3,1:10) /'xxx','yyy','zzz','xyy','xxy','xxz','xzz','yzz','yyz','xyz'/
   data ang_fact(4:4,1:15) /'xxxx','yyyy','zzzz','xxxy','xxxz','yyyx','yyyz','zzzx','zzzy',&
                            'xxyy','xxzz','yyzz','xxyz','yyxz','zzxy'/

   !> How many cartesian functions correspond to the cartesian shell with angular momentum L, i.e. the number of nonzero elements in each row of ang_fact.
   !> These numbers can be obtained using the nshell function from the cgto_hgp module. However, we prefer to spell the numbers out explicitly since Molden format can handle only a limited number of L
   !> values. These are given by the array gto_typ above: the gto_typ array and cart_shell must always be compatible.
   integer, parameter :: cart_shell(0:max_molden_l) = (/1,3,6,10,15/)

   !> Header defining the Molden file.
   character(len=*), parameter :: header_molden = '[Molden Format]'
   !> Header for the Nuclear data section.
   character(len=*), parameter :: header_atoms = '[Atoms]'
   !> Units switch.
   character(len=*), parameter :: str_angs = 'Angs'
   !> Units switch.
   character(len=*), parameter :: str_au = 'AU'
   !> Header for the GTO basis set section.
   character(len=*), parameter :: header_gto = '[GTO]'
   !> Flag for mixed sp basis set types.
   character(len=*), parameter :: str_sp = 'sp'
   !> Header for the Molecular orbital section
   character(len=*), parameter :: header_mo = '[MO]'
   !> String defining the line containing the orbital energy.
   character(len=*), parameter :: str_ene = 'Ene='
   !> String defining the line containing the orbital symmetry.
   character(len=*), parameter :: str_sym = 'Sym='
   !> String defining the line containing the particle's spin in the orbital.
   character(len=*), parameter :: str_spin = 'Spin='
   !> String defining the line containing the orbital occupation number.
   character(len=*), parameter :: str_occ = 'Occup='
   !> Flag for spherical GTO basis functions of the D-type.
   character(len=*), parameter :: str_5d = '[5D]'
   !> Flag for spherical GTO basis functions of the D-type and cartesian functions of the F-type.
   character(len=*), parameter :: str_5d10f = '[5D10F]'
   !> Flag for spherical GTO basis functions of the F-type.
   character(len=*), parameter :: str_7f = '[7F]'
   !> Flag for spherical GTO basis functions of the D and F-types.
   character(len=*), parameter :: str_5d7f = '[5D7F]'
   !> Flag for spherical GTO basis functions of the G-type.
   character(len=*), parameter :: str_9g = '[9G]'
   !> Flag specifying spin up in a given orbital.
   character(len=*), parameter :: str_alpha = 'Alpha'
   !> Flag specifying spin down in a given orbital.
   character(len=*), parameter :: str_beta = 'Beta'

   !> Line containing the orbital energy energy.
   character(len=*), parameter :: form_ene = '(1x,"' // str_ene // '",3x,f10.4)'
   !> Line containing the orbital symmetry.
   character(len=*), parameter :: form_sym = '(1x,"' // str_sym // '",3x,i4,".",i1)'
   !> Line containing the spin of the particle in the orbital.
   character(len=*), parameter :: form_spin = '(1x,"' // str_spin // '",1x,a)'
   !> Line containing the orbital occupation number.
   character(len=*), parameter :: form_occ = '(1x,"' // str_occ // '",3x,f8.6)'
   !> Line containing the information on the shell of GTOs: contracted GTO angular type, no. of primitives, 1.00.
   character(len=*), parameter :: form_gto_head = '(1x,a1,3x,i2,1x,f5.2)'
   !> Line containing the primitive GTO exponent and contraction coefficient.
   character(len=*), parameter :: form_gto_prim = '(d18.10,d18.10)'
   !> Line containing the nuclear data: element name, number, charge, center.
   character(len=*), parameter :: form_atom = '(a2,1x,i4,1x,i4,1x,3f20.10)' !'(a2,1x,i5,1x,i2,3(1x,f12.6))' !the second format is what the Molden program uses
   !> Line in the GTO basis section containing the nucleus index.
   character(len=*), parameter :: form_atom_id = '(1x,i2,1x,i1)'

   !> The cartesian and the corresponding spherical orbital coefficients for l .le. 1 must be the same (up to a rounding error) since the cartesian GTOs and the spherical GTOs for these shells are the same.
   !> This constant defines the precision required on the s,p orbital coefficients when the conversion between the cartesian/spherical representation is performed. It is used to perform sanity checks that
   !> the transformation was done correctly.
   real(kind=cfp), parameter :: s_p_conv_tol = 10e-13_cfp

   public process_ang_fact

contains

   !> This routine takes a single character string of the type ang_fact and returns the powers of the x, y, z factors in it.
   subroutine process_ang_fact(ang_fact,i,j,k)
      implicit none
      character(len=len_ang_fact), intent(in) :: ang_fact
      integer, intent(out) :: i, j, k

      integer :: it, real_len

         real_len = len_trim(ang_fact) !real length of the string in ang_fact: after removal of all unnecessary blanks

         i = 0 !x exponent
         j = 0 !y exponent
         k = 0 !z exponent
         do it=1,real_len
            if (ang_fact(it:it) .eq. 'x') i = i + 1
            if (ang_fact(it:it) .eq. 'y') j = j + 1
            if (ang_fact(it:it) .eq. 'z') k = k + 1
         enddo

   end subroutine process_ang_fact

end module molden_const
