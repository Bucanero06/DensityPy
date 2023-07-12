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
MODULE wigner_cf
use precisn

PRIVATE X_3j, Y_3j, Z_3j
PRIVATE X_6j, Y_6j, Z_6j

CONTAINS

FUNCTION C_G_cf(j1,m1,j2,m2,J,M)
!calculates C-G coefficient <j1,m1,j2,m2|JM> using recursion, i.e. it gives stable values even for large quantum numbers. Please note that half integers
!are not permitted yet, but it is easy to fix.
IMPLICIT NONE
INTEGER, INTENT(IN) :: j1, j2, J, m1, m2, M !assumed to be positive WHOLE integers
REAL(kind=cfp) :: C_G_cf

REAL(kind=cfp) :: w
INTEGER :: m_inp, min_j, j_d1, j_d2, j_d3, j_min, j_max
REAL(kind=cfp), ALLOCATABLE :: f(:)

   m_inp = -M !this is essential to get the correct C-G coefficient
   
   IF (m1 + m2 .ne. M .or. J < abs(j1-j2) .or. J > j1+j2 .or. abs(m1) > j1 .or. abs(m2) > j2 .or. abs(M) > J) THEN 
      C_G_cf = 0.0_cfp
      RETURN
   ELSE IF ((m1 == M .and. j1 == J .and. m2 == 0 .and. j2 == 0) .or. &
            (m2 == M .and. j2 == J .and. m1 == 0 .and. j1 == 0)) THEN !trivial case
      C_G_cf = 1.0_cfp
      RETURN
   ENDIF
   
   !try to minimize the computing time, by choosing the j, which gives the smallest number of allowed values and then permute the columns of W accordingly.
   j_d1 = j2+J - max(abs(j2-J),abs(m1))
   j_d2 = j1+J - max(abs(j1-J),abs(m2))
   j_d3 = j2+j1 - max(abs(j2-j1),abs(m_inp))
   
   min_j = min(j_d1,j_d2,j_d3)
   
   IF (min_j == j_d1) THEN

      CALL Wigner_3j(f,j2,J,m1,m2,m_inp,j_min,j_max)
      w=f(j1)
      C_G_cf = (-1)**(-j1+j2-M)*sqrt(2.0_cfp*J+1)*w
   
   ELSEIF (min_j == j_d2) THEN !account for odd permutation (we permute the first and the second column)

      CALL Wigner_3j(f,j1,J,m2,m1,m_inp,j_min,j_max)
      w=f(j2)
      C_G_cf = (-1)**(-j1+j2-M+j1+j2+J)*sqrt(2.0_cfp*J+1)*w  
   
   ELSE !min_j == j_d3; account for odd permutation (we permute the first and the third column)

      CALL Wigner_3j(f,j2,j1,m_inp,m2,m1,j_min,j_max)
      w=f(J)
      C_G_cf = (-1)**(-j1+j2-M+j1+j2+J)*sqrt(2.0_cfp*J+1)*w
     
   ENDIF

   DEALLOCATE(f)

END FUNCTION C_G_cf

!> This routine will be replaced by wigner_3j once I've figured out where the problem with my implementation is.
subroutine Wigner3j(w3j, jmin, jmax, j2, j3, m1, m2, m3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   This subroutine will calculate the Wigner 3j symbols
!
!      j  j2 j3
!      m1 m2 m3
!
!   for all allowable values of j. The returned values in the array j are 
!   calculated only for the limits
!
!      jmin = max(|j2-j3|, |m1|)
!      jmax = j2 + j3
!
!   To be non-zero, m1 + m2 + m3 = 0. In addition, it is assumed that all j and m are 
!   integers. Returned values have a relative error less than ~1.d-8 when j2 and j3 
!   are less than 103 (see below). In practice, this routine is probably usable up to 165.
!
!   This routine is based upon the stable non-linear recurence relations of Luscombe and 
!   Luban (1998) for the "non classical" regions near jmin and jmax. For the classical 
!   region, the standard three term recursion relationship is used (Schulten and Gordon 1975). 
!   Note that this three term recursion can be unstable and can also lead to overflows. Thus 
!   the values are rescaled by a factor "scalef" whenever the absolute value of the 3j coefficient 
!   becomes greater than unity. Also, the direction of the iteration starts from low values of j
!   to high values, but when abs(w3j(j+2)/w3j(j)) is less than one, the iteration will restart 
!   from high to low values. More efficient algorithms might be found for specific cases 
!   (for instance, when all m's are zero).
!
!   Verification: 
!
!   The results have been verified against this routine run in quadruple precision.
!   For 1.e7 acceptable random values of j2, j3, m2, and m3 between -200 and 200, the relative error
!   was calculated only for those 3j coefficients that had an absolute value greater than 
!   1.d-17 (values smaller than this are for all practical purposed zero, and can be heavily 
!   affected by machine roundoff errors or underflow). 853 combinations of parameters were found
!   to have relative errors greater than 1.d-8. Here I list the minimum value of max(j2,j3) for
!   different ranges of error, as well as the number of times this occured
!   
!   1.d-7 < error  <=1.d-8 = 103   # = 483
!   1.d-6 < error <= 1.d-7 =  116   # = 240
!   1.d-5 < error <= 1.d-6 =  165   # = 93
!   1.d-4 < error <= 1.d-5 = 167   # = 36
!
!   Many times (maybe always), the large relative errors occur when the 3j coefficient 
!   changes sign and is close to zero. (I.e., adjacent values are about 10.e7 times greater 
!   in magnitude.) Thus, if one does not need to know highly accurate values of the 3j coefficients
!   when they are almost zero (i.e., ~1.d-10) then this routine is probably usable up to about 160.
!
!   These results have also been verified for parameter values less than 100 using a code
!   based on the algorith of de Blanc (1987), which was originally coded by Olav van Genabeek, 
!   and modified by M. Fang (note that this code was run in quadruple precision, and
!   only calculates one coefficient for each call. I also have no idea if this code
!   was verified.) Maximum relative errors in this case were less than 1.d-8 for a large number
!   of values (again, only 3j coefficients greater than 1.d-17 were considered here).
!   
!   The biggest improvement that could be made in this routine is to determine when one should
!   stop iterating in the forward direction, and start iterating from high to low values. 
!
!   Calling parameters
!      IN   
!         j2, j3, m1, m2, m3    Integer values.
!      OUT   
!         w3j         Array of length jmax - jmin + 1.
!         jmin, jmax      Minimum and maximum values
!                  out output array.
!   Dependencies: None
!   
!   Written by Mark Wieczorek August (2004)
!
!   August 2009: Based on the suggestions of Roelof Rietbroek, the calculation of RS has been slightly
!   modified so that division by zero will not cause a run time crash (this behavior depends on how the 
!   compiler treats IEEE floating point exceptions). These values were never used in the original code 
!   when this did occur.
!
!   Copyright (c) 2005-2009, Mark A. Wieczorek
!   All rights reserved.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   use utils, only: xermsg
   implicit none
   integer, intent(in) ::   j2, j3, m1, m2, m3
   integer, intent(out) ::   jmin, jmax
   real(kind=cfp), intent(out) ::   w3j(:)
   real(kind=cfp) ::      wnmid, wpmid, scalef, denom, rs(0:j2+j3+1), wl(0:j2+j3+1), wu(0:j2+j3+1)
   real(kind=cfp) ::      xjmin, yjmin, yjmax, zjmax, xj, zj
   integer ::       j, jnum, jp, jn, k, flag1, flag2, jmid
   
   
   if (size(w3j) < j2+j3+1) then
      print*, "W3J must be dimensioned (J2+J3+1) where J2 and J3 are ", j2, j3
      print*, "Input array is dimensioned ", size(w3j)
      call xermsg('wigner_cf','Wigner3j','The input array w3j is too small; see output for details.',1,1)
   endif
   
   w3j = 0.0_cfp
   
   flag1 = 0
   flag2 = 0
   
   scalef = 1.0e3_cfp
   
   jmin = max(abs(j2-j3), abs(m1))
   jmax = j2 + j3
   jnum = jmax - jmin + 1

   if (abs(m2) > j2 .or. abs(m3) > j3) then
      return
   elseif (m1 + m2 + m3 /= 0) then
      return
   elseif (jmax < jmin) then
      return
   endif
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !   Only one term is present
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if (jnum == 1) then
      w3j(1) = 1.0_cfp / sqrt(2.0_cfp*jmin+1.0_cfp)
      if ( (w3j(1) < 0.0_cfp .and. (-1)**(j2-j3+m2+m3) > 0) .or. (w3j(1) > 0.0_cfp .and. (-1)**(j2-j3+m2+m3) < 0) ) w3j(1) = -w3j(1)
      return   
   endif
      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !    Calculate lower non-classical values for [jmin, jn]. If the second term
   !   can not be calculated because the recursion relationsips give rise to a
   !   1/0, then set flag1 to 1.  If all m's are zero, then this is not a problem 
   !   as all odd terms must be zero.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   rs = 0.0_cfp
   wl = 0.0_cfp
   
   xjmin = x_3j_cf(jmin)
   yjmin = y_3j_cf(jmin)
   
   if (m1 == 0 .and. m2 == 0 .and. m3 == 0) then      ! All m's are zero
   
      wl(jindex(jmin)) = 1.0_cfp
      wl(jindex(jmin+1)) = 0.0_cfp
      jn = jmin+1
      
   elseif (yjmin == 0.0_cfp) then            ! The second terms is either zero
   
      if (xjmin == 0.0_cfp) then         ! or undefined
         flag1 = 1
         jn = jmin
      else
         wl(jindex(jmin)) = 1.0_cfp
         wl(jindex(jmin+1)) = 0.0_cfp
         jn = jmin+1
      endif
      
   elseif ( xjmin * yjmin >= 0.0_cfp) then         ! The second term is outside of the non-classical region 
      wl(jindex(jmin)) = 1.0_cfp
      wl(jindex(jmin+1)) = -yjmin / xjmin
      jn = jmin+1
      
   else                     ! Calculate terms in the non-classical region
   
      rs(jindex(jmin)) = -xjmin / yjmin
      
      jn = jmax
      do j=jmin + 1, jmax-1, 1
         denom =  y_3j_cf(j) + z_3j_cf(j)*rs(jindex(j-1))
         xj = x_3j_cf(j)
         if (abs(xj) > abs(denom) .or. xj * denom >= 0.0_cfp .or. denom == 0.0_cfp) then
            jn = j-1
            exit
         else
            rs(jindex(j)) = -xj / denom
         endif
            
      enddo
      
      wl(jindex(jn)) = 1.0_cfp
      
      do k=1, jn - jmin, 1
         wl(jindex(jn-k)) = wl(jindex(jn-k+1)) * rs(jindex(jn-k))
      enddo
      
      if (jn == jmin) then               ! Calculate at least two terms so that
         wl(jindex(jmin+1)) = -yjmin / xjmin      ! these can be used in three term
         jn = jmin+1               ! recursion
         
      endif

   endif
   
   if (jn == jmax) then               ! All terms are calculated
   
      w3j(1:jnum) = wl(1:jnum)
      call normw3j
      call fixsign
      
      return

   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !    Calculate upper non-classical values for [jp, jmax].
   !   If the second last term can not be calculated because the
   !   recursion relations give a 1/0, then set flag2 to 1. 
   !   (Note, I don't think that this ever happens).
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   wu = 0.0_cfp
   
   yjmax = y_3j_cf(jmax)
   zjmax = z_3j_cf(jmax)
   
   if (m1 == 0 .and. m2 == 0 .and. m3 == 0) then
   
      wu(jindex(jmax)) = 1.0_cfp
      wu(jindex(jmax-1)) = 0.0_cfp
      jp = jmax-1
      
   elseif (yjmax == 0.0_cfp) then
   
      if (zjmax == 0.0_cfp) then
         flag2 = 1
         jp = jmax
      else
         wu(jindex(jmax)) = 1.0_cfp
         wu(jindex(jmax-1)) = - yjmax / zjmax
         jp = jmax-1
      endif
      
   elseif (yjmax * zjmax >= 0.0_cfp) then
   
      wu(jindex(jmax)) = 1.0_cfp
      wu(jindex(jmax-1)) = - yjmax / zjmax
      jp = jmax-1

   else
      rs(jindex(jmax)) = -zjmax / yjmax

      jp = jmin
      do j=jmax-1, jn, -1
         denom = y_3j_cf(j) + x_3j_cf(j)*rs(jindex(j+1))
         zj = z_3j_cf(j)
         if (abs(zj) > abs(denom) .or. zj * denom >= 0.0_cfp .or. denom == 0.0_cfp) then
            jp = j+1
            exit
         else
            rs(jindex(j)) = -zj / denom
         endif
      enddo   
      
      wu(jindex(jp)) = 1.0_cfp
      
      do k=1, jmax - jp, 1
         wu(jindex(jp+k)) = wu(jindex(jp+k-1))*rs(jindex(jp+k))
      enddo
      
      if (jp == jmax) then
         wu(jindex(jmax-1)) = - yjmax / zjmax
         jp = jmax-1
      endif
      
   endif
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !    Calculate classical terms for [jn+1, jp-1] using standard three
   !    term rercusion relationship. Start from both jn and jp and stop at the
   !    midpoint. If flag1 is set, then perform the recursion solely from high to
   !    low values. If flag2 is set, then perform the recursion solely from low to high.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   if (flag1 == 0) then
   
      jmid = (jn + jp)/2
      
      do j=jn, jmid - 1, 1         
         wl(jindex(j+1)) = - (z_3j_cf(j)*wl(jindex(j-1)) +y_3j_cf(j)*wl(jindex(j))) / x_3j_cf(j)
         
         if (abs(wl(jindex(j+1))) > 1.0_cfp) then            ! watch out for overflows.
            wl(jindex(jmin):jindex(j+1)) = wl(jindex(jmin):jindex(j+1)) / scalef
         endif
         
         ! if values are decreasing then stop upward iteration and start with the downward iteration.
         if (abs(wl(jindex(j+1)) / wl(jindex(j-1))) < 1.0_cfp .and. wl(jindex(j+1)) /= 0.0_cfp) then
            jmid = j+1
            exit
         endif
      enddo
      
      wnmid = wl(jindex(jmid))
    
      ! Make sure that the stopping midpoint value is not a zero, or close to it! 
      if (wl(jindex(jmid-1)) /= 0.0_cfp) then
         if (abs(wnmid/wl(jindex(jmid-1))) < 1.0e-6_cfp) then 
            wnmid = wl(jindex(jmid-1))      
            jmid = jmid - 1               
         endif
      endif
      
      
      do j=jp, jmid+1, -1
         wu(jindex(j-1)) = - (x_3j_cf(j)*wu(jindex(j+1)) + y_3j_cf(j)*wu(jindex(j)) ) / z_3j_cf(j)
         if (abs(wu(jindex(j-1))) > 1.0_cfp) then
            wu(jindex(j-1):jindex(jmax)) = wu(jindex(j-1):jindex(jmax)) / scalef
         endif
   
      enddo
      
      wpmid = wu(jindex(jmid))
      
      ! rescale two sequences to common midpoint
      
      if (jmid == jmax) then
         w3j(1:jnum) = wl(1:jnum)
      elseif (jmid == jmin) then
         w3j(1:jnum) = wu(1:jnum)
      else
         w3j(1:jindex(jmid)) = wl(1:jindex(jmid)) * wpmid / wnmid 
         w3j(jindex(jmid+1):jindex(jmax)) = wu(jindex(jmid+1):jindex(jmax))
      endif
      
   elseif (flag1 == 1 .and. flag2 == 0) then   ! iterature in downward direction only
      
      do j=jp, jmin+1, -1
         wu(jindex(j-1)) = - (x_3j_cf(j)*wu(jindex(j+1)) + y_3j_cf(j)*wu(jindex(j)) ) / z_3j_cf(j)
         if (abs(wu(jindex(j-1))) > 1) then
            wu(jindex(j-1):jindex(jmax)) = wu(jindex(j-1):jindex(jmax)) / scalef
         endif
      enddo
      
      w3j(1:jnum) = wu(1:jnum)
      
   elseif (flag2 == 1 .and. flag1 == 0) then   ! iterature in upward direction only
      
      do j=jn, jp-1, 1
         wl(jindex(j+1)) = - (z_3j_cf(j)*wl(jindex(j-1)) +y_3j_cf(j)*wl(jindex(j))) / x_3j_cf(j)
         if (abs(wl(jindex(j+1))) > 1) then
            wl(jindex(jmin):jindex(j+1)) = wl(jindex(jmin):jindex(j+1))/ scalef
         endif
      enddo
      
      w3j(1:jnum) = wl(1:jnum)
      
   elseif (flag1 == 1 .and. flag2 == 1) then

      call xermsg('wigner_cf','Wigner3j','Can not calculate function for input values, both flag1 and flag 2 are set.',2,1)

   endif

   
   call normw3j
   call fixsign
      
   
   contains
   
      integer function jindex(j)
         integer :: j
         jindex = j-jmin+1
      end function jindex
   
      real(kind=cfp) function a_3j_cf(j)
         integer :: j, i
         i = ((j)**2 - (j2-j3)**2) * ((j2+j3+1)**2 - (j)**2) * ((j)**2-(m1)**2)
         a_3j_cf = sqrt(real(i,kind=cfp))
      end function a_3j_cf
      
      real(kind=cfp) function y_3j_cf(j)
         integer :: j
         y_3j_cf = real(-(2*j+1) * ( (m1) * ((j2)*(j2+1) - (j3)*(j3+1) ) - (m3-m2)*(j)*(j+1) ),kind=cfp)
      end function y_3j_cf

      real(kind=cfp) function x_3j_cf(j)   
         integer :: j
         x_3j_cf = real(j,kind=cfp) * a_3j_cf(j+1)
      end function x_3j_cf
      
      real(kind=cfp) function z_3j_cf(j)
         integer :: j
         z_3j_cf = real(j+1,kind=cfp)*a_3j_cf(j)
      end function z_3j_cf
      
      subroutine normw3j
         real(kind=cfp):: norm
         integer j
         
         norm = 0.0_cfp
         do j = jmin, jmax
            norm = norm + real(2*j+1,kind=cfp) * w3j(jindex(j))**2
         enddo
         
         w3j(1:jnum) = w3j(1:jnum) / sqrt(norm)
         
      end subroutine normw3j
      
      subroutine fixsign
      
         if ( (w3j(jindex(jmax)) < 0.0_cfp .and. (-1)**(j2-j3+m2+m3) > 0) .or. &
              (w3j(jindex(jmax)) > 0.0_cfp .and. (-1)**(j2-j3+m2+m3) < 0) ) then
            w3j(1:jnum) = -w3j(1:jnum)
         endif
         
      end subroutine fixsign
      
end subroutine Wigner3j

SUBROUTINE Wigner_3j(f,j2,j3,m1,m2,m3,j_min,j_max)
use utils, only: xermsg
!calculates (recursively) the Wigner 3j symbol using the method of J. H. Luscombe and M. Luban: Phys. Rev. E, vol. 57, no. 6, 7274, 1998.
!the array f(:) is allocated like this f(j_min:j_max), where j_min and j_max are the limits of allowed j1 for which the 3j symbol is calculated.
!WARNING numerical instability may arise in case l1 > l2 > l3 and l1 large. I.e. it is best to call this routine in combination of l1,l2,l3 such that l1 is the smallest.
IMPLICIT NONE
INTEGER, INTENT(IN) :: j2, j3, m1, m2, m3 !assumed to be positive integers
INTEGER, INTENT(OUT) :: j_min, j_max
REAL(kind=cfp), ALLOCATABLE, INTENT(OUT) :: f(:)

INTEGER :: j, p, k, n_m, n_p
REAL(kind=cfp) :: norm, tmp, sp1, sm1, s_tmp, denom, nom
REAL(kind=cfp), ALLOCATABLE :: r(:), s(:), psi_m(:), psi_p(:)
LOGICAL :: from_b, from_a, phase
INTEGER :: err

   call xermsg ('wigner_cf', 'Wigner_3j', 'Use Wigner3j instead. This has some phase problems for very large L > 14.', 1, 1)

   j_max = j2+j3
   j_min = max(abs(j2-j3),abs(m1))
!   PRINT *,'j_min,j_max',j_min,j_max
!   WRITE(*,'("3j input",5i4)') j2,j3,m1,m2,m3

   IF (j_min > j_max) STOP "Error in Wigner_3j; Unphysical values of j_min and j_max."

   ALLOCATE(f(j_min:j_max),stat=err)
   if (err .ne. 0) call xermsg ('wigner_cf', 'Wigner_3j', 'Memory allocation failed; see Wigner_3j', err, 1)
   f(:) = 1.0_cfp
   
   !selection rules below. Please note that we are not checking if j1 + j2 + j3 is integer as required, i.e. this is an assuption about sanity of the user!
   IF (m1 + m2 + m3 .ne. 0 .or. abs(m2) > j2 .or. abs(m3) > j3) THEN
      f(:) = 0.0_cfp
      RETURN
   ENDIF

   IF (j_max == j_min) THEN !only one allowed value
      f(:) = 1.0_cfp/sqrt(2.0_cfp*j_max+1.0_cfp)
      IF ((-1)**(j2-j3+m2+m3) < 0) f(:)=-f(:)
      RETURN
   ENDIF
   
   ALLOCATE(r(j_min:j_max),s(j_min:j_max),psi_m(j_min:j_max),psi_p(j_min:j_max),stat=err)
   if (err .ne. 0) call xermsg ('wigner_cf', 'Wigner_3j', 'Memory allocation 2 failed; see Wigner_3j', err, 1)
   r(:)=0; s(:)=0

   from_b = .false.; from_a = .false.; phase = .false.
   
   !backward propagation in the 'outer' nonclassical region
   !if m1=0=m2=m3 then w(j) = 0 for j+j2+j3 = odd!!
   !if denom = 0 at the start then check if the backward three term recurrence can be started from here. If not then start it from below. 
   IF (m1 .eq. 0 .and. m2 .eq.0 .and. m3 .eq. 0) THEN

      n_p = j_max-1

      !if m1=0=m2=m3 then w(j) = 0 for j+j2+j3 = odd. For j = j_max and j_max = j2+j3 this means that j+j2+j3 is even, i.e. psi_p(j_max) .ne. 0
      phase = .true.
      psi_p(j_max) = 1
      psi_p(j_max-1) = 0
 
   ELSEIF (Y_3j(j_max, j2, j3, m1, m2, m3) .eq. 0.0_cfp) THEN !the two-term recurrence relation cannot be used (Y_3j is the denominator of the recurrence).

      IF (Z_3j(j_max, j2, j3, m1, m2, m3) .eq. 0.0_cfp) THEN !the three term backward recurrence relation cannot be used from this point
         n_p = j_max
         !set some logical variable saying that the three-term recurrence must be iterated from below:
         from_b = .true.
      ELSE
         phase = .true.
         n_p = j_max -1
         psi_p(j_max) = 1
         psi_p(n_p) = 0 !(= -Y/Z, but Y=0)
      ENDIF

   ELSE !no problem - we can iterate the two term recurrence relation until we hit either the classical region or one of the r(j) is 0.

      !PRINT *,'backward',j_max-1,j_min+1

      r(j_max) = -Z_3j(j_max, j2, j3, m1, m2, m3)/Y_3j(j_max, j2, j3, m1, m2, m3)
      n_p = j_max

      DO j=j_max-1,j_min+1,-1

         !Z.M. bugfix: denom wasn't defined here
         denom=Y_3j(j, j2, j3, m1, m2, m3)+X_3j(j, j2, j3, m1, m2, m3)*r(j+1)

         IF (denom .ne. 0) THEN
            r(j) = -Z_3j(j, j2, j3, m1, m2, m3)/denom
         ELSE
            n_p = j+1
            EXIT
         ENDIF

         IF (r(j) > 1.0_cfp .or. r(j) .le. 0.0_cfp) THEN !classical region boundary
            n_p = j+1
            EXIT
         ENDIF

      ENDDO

      !calculate the value on the boundary of the 'outer' nonclassical region
      psi_p(:)=1.0_cfp
      DO k=1,j_max-n_p
         DO p=1,k
            psi_p(n_p+k)=psi_p(n_p+k)*r(n_p+p)
         ENDDO
      ENDDO

      IF (n_p .eq. j_max) THEN
         phase = .true.
         psi_p(j_max) = 1
         psi_p(j_max-1) = -Y_3j(j_max, j2, j3, m1, m2, m3)/Z_3j(j_max, j2, j3, m1, m2, m3)
         n_p = j_max-1
      ENDIF

      IF (n_p .eq. j_min) THEN !all coefficients have been evaluated

         f(:)=psi_p(:)

         !calculate the normalization of the whole thing
         norm=0.0_cfp
         DO j=j_min,j_max
            norm=norm+(2.0_cfp*j+1)*f(j)**2.0_cfp
         ENDDO
         !normalize and correct the phase
         tmp=sign(1.0_cfp,f(j_max)) !this makes sure that I get rid of a possibly wrong sign of f(j_max)
         f(:)=tmp*f(:)/sqrt(norm)*(-1)**(j2-j3-m1)
         RETURN

      ENDIF

   ENDIF

   !forward propagation in the 'inner' nonclassical region:
   IF (m1 .eq. 0 .and. m2 .eq.0 .and. m3 .eq. 0) THEN

      n_m = j_min+1

      !if 0=m1=m2=m3 then w(j) = 0 for j+j2+j3 = odd!! j_min+j2+j3 = even (always), because j_min=abs(j2-j3); (m1=0)
      psi_m(j_min) = 1
      psi_m(j_min+1) = 0

   ELSEIF (Y_3j(j_min, j2, j3, m1, m2, m3) .eq. 0.0_cfp) THEN !the two-term recurrence relation cannot be used

      IF (X_3j(j_min, j2, j3, m1, m2, m3) .eq. 0.0_cfp) THEN !the three term forward recurrence relation cannot be used from this point
         n_m = j_min
         !set some logical variable saying that the three-term recurrence must be iterated from above:
         from_a = .true.
      ELSE
         n_m = j_min + 1
         psi_m(j_min) = 1
         psi_m(n_m) = 0
      ENDIF

   ELSE

      !print *,'here'

      s(j_min)=-X_3j(j_min, j2, j3, m1, m2, m3)/Y_3j(j_min, j2, j3, m1, m2, m3)
      n_m = j_min

      DO j=j_min+1,n_p

         denom = (Y_3j(j, j2, j3, m1, m2, m3)+Z_3j(j, j2, j3, m1, m2, m3)*s(j-1))

         IF (denom .ne. 0) THEN
            s(j) = -X_3j(j, j2, j3, m1, m2, m3)/denom
         ELSE
            n_m = j-1
            EXIT
         ENDIF

         IF (s(j) > 1.0_cfp .or. s(j) .le. 0.0_cfp) THEN !classical region boundary
            n_m = j-1
            EXIT
         ENDIF

      ENDDO

      !calculate the value on the 'inner' boundary of the classical region
      psi_m(:)=1.0_cfp
      DO k=1,n_m-j_min
         DO p=1,k
            psi_m(n_m-k)=psi_m(n_m-k)*s(n_m-p)
         ENDDO
      ENDDO

      IF (n_m .eq. j_min) THEN
         psi_m(j_min) = 1
         psi_m(j_min+1) = -Y_3j(j_min, j2, j3, m1, m2, m3)/X_3j(j_min, j2, j3, m1, m2, m3)
         n_m = j_min+1
      ENDIF

   ENDIF

   IF (from_b .and. .not.(from_a)) THEN
      !forward propagation from the 'inner' boundary of the nonclassical region using the three term recurrence relation
      !PRINT *,'from_b true, from_a false',j_min,j_max
      s_tmp=psi_m(n_m)
      sm1=psi_m(n_m-1)
      DO j=n_m,n_p
      
         psi_m(j) = s_tmp
         sp1 = (-Y_3j(j, j2, j3, m1, m2, m3)*s_tmp-Z_3j(j, j2, j3, m1, m2, m3)*sm1)/X_3j(j, j2, j3, m1, m2, m3)
         sm1=s_tmp
         s_tmp=sp1
      
      ENDDO
   
      f(j_min:j_max) = psi_m(j_min:j_max)
   
   ELSEIF (.not.(from_b) .and. from_a) THEN
      !backward propagation from the 'outer' boundary of the nonclassical region using the three term recurrence relation
      !PRINT *,'from_b false, from_a .true.',j_min,j_max
      s_tmp=psi_p(n_p)
      sp1=psi_p(n_p+1)
      DO j=n_p,n_m,-1

         psi_p(j) = s_tmp

         nom = -X_3j(j, j2, j3, m1, m2, m3)*sp1-Y_3j(j, j2, j3, m1, m2, m3)*s_tmp
         denom = Z_3j(j, j2, j3, m1, m2, m3)

         !this situation would have raised a floating point exception if -fpe0 was used for compilation, hence nom and denom checks
         if (nom .eq. 0.0_cfp .and. denom .eq. 0.0_cfp) then
            sm1 = 0.0_cfp
         else
            sm1 = (-X_3j(j, j2, j3, m1, m2, m3)*sp1-Y_3j(j, j2, j3, m1, m2, m3)*s_tmp)/Z_3j(j, j2, j3, m1, m2, m3)
         endif
         sp1=s_tmp
         s_tmp=sm1

      ENDDO
 
      f(j_min:j_max) = psi_p(j_min:j_max)

   ELSEIF (.not.(from_b) .and. .not.(from_a)) THEN !from n_m to n_p or n_p to n_m using the three-term recurrence relation.
      !forward propagation from the 'inner' boundary of the nonclassical region using the three term recurrence relation
      !PRINT *,'from_* false',n_m,n_p
      s_tmp=psi_m(n_m)
      sm1=psi_m(n_m-1)
      DO j=n_m,n_p

         psi_m(j) = s_tmp
         sp1 = (-Y_3j(j, j2, j3, m1, m2, m3)*s_tmp-Z_3j(j, j2, j3, m1, m2, m3)*sm1)/X_3j(j, j2, j3, m1, m2, m3)
         sm1=s_tmp
         s_tmp=sp1

      ENDDO

      IF (n_p+1 .le. j_max .and. n_p .ne. j_min) THEN !this condition is a bit clumsy, but it works. The whole IF can be rearranged.
         psi_m(n_p+1)=sp1
         !match the solutions based on the value on the boundary of the 'outer' nonclassical region
         tmp=psi_m(n_p+1)
         psi_m(j_min:n_p)=psi_m(:)/tmp*psi_p(n_p+1) !psi_p(n_p+1) can never be zero

         !print *,'ns',n_p,n_m,j_max,j_min

         !if (n_p+1 .eq. j_max .and. n_p-n_m>1) psi_p(n_p+1:j_max) = sign(1.0_cfp,tmp)*psi_p(n_p+1:j_max)

         f(j_min:n_p) = psi_m(j_min:n_p)
         f(n_p+1:j_max) = psi_p(n_p+1:j_max)
      ELSEIF (n_p == j_max) THEN !we have psi_m for all j
         f(j_min:j_max) = psi_m(j_min:j_max)
      ELSEIF (n_p == j_min) THEN !we have psi_p up to j_min+1
         psi_p(j_min) = (-X_3j(j_min+1, j2, j3, m1, m2, m3) * psi_p(j_min+2) - Y_3j(j_min+1, j2, j3, m1, m2, m3) * psi_p(j_min+1)) &
                        / Z_3j(j_min+1, j2, j3, m1, m2, m3)
         f(j_min:j_max) = psi_p(j_min:j_max)
      ENDIF

   ELSE !from_b .and. from_a .eq..true.

      STOP "Error in Wigner_3j. Can't calculate the given coefficient using the simplified algorithm. &
           &Try three term recurrence only at this point?"

   ENDIF
   
   !calculate the normalization of the whole thing
   norm=0.0_cfp
   DO j=j_min,j_max
      norm=norm+(2.0_cfp*j+1)*f(j)**2.0_cfp
   ENDDO

   !normalize and correct the phase
   tmp=sign(1.0_cfp,f(j_max)) !this makes sure that I get rid of a possibly wrong sign of f(j_max)
   f(:)=tmp*f(:)/sqrt(norm)*(-1)**(j2-j3-m1)
   
   DEALLOCATE(r,s,psi_m,psi_p)

END SUBROUTINE Wigner_3j

SUBROUTINE Wigner_6j(f,j2,j3,l1,l2,l3,j_min,j_max)
!calculates (recursively) the Wigner 6j symbol using the method of J. H. Luscombe and M. Luban: Phys. Rev. E, vol. 57, no. 6, 7274, 1998.
!the array f(:) is allocated like this f(j_min:j_max), where j_min and j_max are the limits of allowed j1 for which the 6j symbol is calculated.
use utils, only: xermsg
IMPLICIT NONE
INTEGER :: err
INTEGER, INTENT(IN) :: j2, j3, l1, l2, l3 !assumed to be positive integers
INTEGER, INTENT(OUT) :: j_min, j_max
REAL(kind=cfp), ALLOCATABLE, INTENT(OUT) :: f(:)

INTEGER :: j, p, k, n_m, n_p
REAL(kind=cfp) :: norm, tmp, sp1, sm1, s_tmp, denom
REAL(kind=cfp), ALLOCATABLE :: r(:), s(:), psi_m(:), psi_p(:)
LOGICAL :: from_b, from_a

   call xermsg ('wigner_cf', 'Wigner_6j', 'Use Wigner3j instead. This has some phase problems for very large L > 14.', 1, 1)

   j_min = max(abs(j2-j3),abs(l2-l3))
   j_max = min(j2+j3,l2+l3)

   IF (j_min > j_max) STOP "Error in Wigner_6j; Unphysical values of j_min and j_max."

   ALLOCATE(f(j_min:j_max),stat=err)
   if (err .ne. 0) call xermsg ('wigner_cf', 'Wigner_6j', 'Memory allocation failed; see Wigner_6j', err, 1)
   f(:) = 1.0_cfp
   
!selection rules below. Please note that we are not checking if sums of each row are integers as required,i.e. this is an assuption about sanity of the user!
   IF (j2 < abs(l1-l3) .or. j2 > l1+l3 .or. l2 < abs(l1-j3) .or. l2 > l1+j3) THEN !triangular inequalities
      f(:) = 0.0_cfp
      RETURN
   ENDIF
   
   IF (j_max == j_min) THEN !only one allowed value
      f(:) = 1.0_cfp/sqrt((2.0_cfp*j_max+1.0_cfp)*(2.0_cfp*l1+1))
      IF ((-1)**(j2+j3+l2+l3) < 0) f(:)=-f(:)
      RETURN
   ENDIF
   
   ALLOCATE(r(j_min:j_max),s(j_min:j_max),psi_m(j_min:j_max),psi_p(j_min:j_max),stat=err)
   if (err .ne. 0) call xermsg ('wigner_cf', 'Wigner_6j', 'Memory allocation 2 failed; see Wigner_6j', err, 1)
   r(:)=0; s(:)=0

   from_b = .false.; from_a = .false.

   !backward propagation in the 'outer' nonclassical region
   !if denom = 0 at the start then check if the backward three term recurrence can be started from here. If not then start it from below. 
   IF (j_max > j2+j3 .or. j_max < abs(j2-j3) .or. j_max > l2+l3 .or. j_max < abs(l2-l3)) THEN

      n_p = j_max-1

      !if m1=0=m2=m3 then w(j) = 0 for j+j2+j3 = odd. For j = j_max and j_max = j2+j3 this means that j+j2+j3 is even, i.e. psi_p(j_max) .ne. 0
      psi_p(j_max) = 1
      psi_p(j_max-1) = 0
 
   ELSEIF (Y_6j(j_max, j2, j3, l1, l2, l3) .eq. 0.0_cfp) THEN !the two-term recurrence relation cannot be used (Y_6j is the denominator of the recurrence).

      IF (Z_6j(j_max, j2, j3, l1, l2, l3) .eq. 0.0_cfp) THEN !the three term backward recurrence relation cannot be used from this point
         n_p = j_max
         !set some logical variable saying that the three-term recurrence must be iterated from below:
         from_b = .true.
      ELSE
         n_p = j_max -1
         psi_p(j_max) = 1
         psi_p(n_p) = 0 !(= -Y/Z, but Y=0)
      ENDIF

   ELSE !no problem - we can iterate the two term recurrence relation until we hit either the classical region or one of the r(j) is 0.

      r(j_max) = -Z_6j(j_max, j2, j3, l1, l2, l3)/Y_6j(j_max, j2, j3, l1, l2, l3)
      n_p = j_max

      DO j=j_max-1,j_min+1,-1

         !Z.M. bugfix: denom wasn't defined here
         denom=Y_6j(j, j2, j3, l1, l2, l3)+X_6j(j, j2, j3, l1, l2, l3)*r(j+1);

         IF (denom .ne. 0) THEN
            r(j) = -Z_6j(j, j2, j3, l1, l2, l3)/denom
         ELSE
            n_p = j+1
            EXIT
         ENDIF

         IF (r(j) > 1.0_cfp .or. r(j) .le. 0.0_cfp) THEN !classical region boundary
            n_p = j+1
            EXIT
         ENDIF

      ENDDO

      !calculate the value on the boundary of the 'outer' nonclassical region
      psi_p(:)=1.0_cfp
      DO k=1,j_max-n_p
         DO p=1,k
            psi_p(n_p+k)=psi_p(n_p+k)*r(n_p+p)
         ENDDO
      ENDDO

      IF (n_p .eq. j_max) THEN
         psi_p(j_max) = 1
         psi_p(j_max-1) = -Y_6j(j_max, j2, j3, l1, l2, l3)/Z_6j(j_max, j2, j3, l1, l2, l3)
         n_p = j_max-1
      ENDIF

      IF (n_p .eq. j_min) THEN !all coefficients have been evaluated

         f(:)=psi_p(:)

         !calculate the normalization of the whole thing
         norm=0.0_cfp
         DO j=j_min,j_max
            norm=norm+(2.0_cfp*j+1)*f(j)**2.0_cfp
         ENDDO
         norm=(2*l1+1)*norm

         !normalize and correct the phase
         tmp=sign(1.0_cfp,f(j_max)) !this makes sure that I get rid of a possibly wrong sign of f(j_max)
         IF (tmp .ne. (-1)**(j2+j3+l2+l3)) THEN
            f(:)=-f(:)/sqrt(norm)
         ELSE
            f(:)=f(:)/sqrt(norm)
         ENDIF

      ENDIF

   ENDIF

   !forward propagation in the 'inner' nonclassical region:
   IF (j_min > j2+j3 .or. j_min < abs(j2-j3) .or. j_min > l2+l3 .or. j_min < abs(l2-l3)) THEN

      n_m = j_min+1

      !if w(j) = 0
      psi_m(j_min) = 1
      psi_m(j_min+1) = 0

   ELSEIF (Y_6j(j_min, j2, j3, l1, l2, l3) .eq. 0.0_cfp) THEN !the two-term recurrence relation cannot be used

      IF (X_6j(j_min, j2, j3, l1, l2, l3) .eq. 0.0_cfp) THEN !the three term forward recurrence relation cannot be used from this point
         n_m = j_min
         !set some logical variable saying that the three-term recurrence must be iterated from above:
         from_a = .true.
      ELSE
         n_m = j_min + 1
         psi_m(j_min) = 1
         psi_m(n_m) = 0
      ENDIF

   ELSE

      s(j_min)=-X_6j(j_min, j2, j3, l1, l2, l3)/Y_6j(j_min, j2, j3, l1, l2, l3)
      n_m = j_min

      DO j=j_min+1,n_p

         denom = (Y_6j(j, j2, j3, l1, l2, l3)+Z_6j(j, j2, j3, l1, l2, l3)*s(j-1))

         IF (denom .ne. 0) THEN
            s(j) = -X_6j(j, j2, j3, l1, l2, l3)/denom
         ELSE
            n_m = j-1
            EXIT
         ENDIF

         IF (s(j) > 1.0_cfp .or. s(j) .le. 0.0_cfp) THEN !classical region boundary
            n_m = j-1
            EXIT
         ENDIF

      ENDDO

      !calculate the value on the 'inner' boundary of the classical region
      psi_m(:)=1.0_cfp
      DO k=1,n_m-j_min
         DO p=1,k
            psi_m(n_m-k)=psi_m(n_m-k)*s(n_m-p)
         ENDDO
      ENDDO

      IF (n_m .eq. j_min) THEN
         psi_m(j_min) = 1
         psi_m(j_min+1) = -Y_6j(j_min, j2, j3, l1, l2, l3)/X_6j(j_min, j2, j3, l1, l2, l3)
         n_m = j_min+1
      ENDIF

   ENDIF

   IF (from_b .and. .not.(from_a)) THEN
      !forward propagation from the 'inner' boundary of the nonclassical region using the three term recurrence relation
      !PRINT *,'from_b true, from_a false',j_min,j_max
      s_tmp=psi_m(n_m)
      sm1=psi_m(n_m-1)
      DO j=n_m,n_p
      
         psi_m(j) = s_tmp
         sp1 = (-Y_6j(j, j2, j3, l1, l2, l3)*s_tmp-Z_6j(j, j2, j3, l1, l2, l3)*sm1)/X_6j(j, j2, j3, l1, l2, l3)
         sm1=s_tmp
         s_tmp=sp1
      
      ENDDO
   
      f(j_min:j_max) = psi_m(j_min:j_max)
   
   ELSEIF (.not.(from_b) .and. from_a) THEN
      !backward propagation from the 'outer' boundary of the nonclassical region using the three term recurrence relation
      !PRINT *,'from_b false, from_a .true.',j_min,j_max
      s_tmp=psi_p(n_p)
      sp1=psi_p(n_p+1)
      DO j=n_p,n_m,-1

         psi_p(j) = s_tmp
         sm1 = (-X_6j(j, j2, j3, l1, l2, l3)*sp1-Y_6j(j, j2, j3, l1, l2, l3)*s_tmp)/Z_6j(j, j2, j3, l1, l2, l3)
         sp1=s_tmp
         s_tmp=sm1

      ENDDO
 
      f(j_min:j_max) = psi_p(j_min:j_max)

   ELSEIF (.not.(from_b) .and. .not.(from_a)) THEN !from n_m to n_p or n_p to n_m using the three-term recurrence relation.
      !forward propagation from the 'inner' boundary of the nonclassical region using the three term recurrence relation
      !PRINT *,'from_* false',n_m,n_p
      s_tmp=psi_m(n_m)
      sm1=psi_m(n_m-1)
      DO j=n_m,n_p

         psi_m(j) = s_tmp
         sp1 = (-Y_6j(j, j2, j3, l1, l2, l3)*s_tmp-Z_6j(j, j2, j3, l1, l2, l3)*sm1)/X_6j(j, j2, j3, l1, l2, l3)
         sm1=s_tmp
         s_tmp=sp1

      ENDDO

      IF (n_p+1 .le. j_max .and. n_p .ne. j_min) THEN !this condition is a bit clumsy, but it works. The whole IF can be rearranged.
         psi_m(n_p+1)=sp1
         !match the solutions based on the value on the boundary of the 'outer' nonclassical region
         tmp=psi_m(n_p+1)
         psi_m(j_min:n_p)=psi_m(:)/tmp*psi_p(n_p+1) !psi_p(n_p+1) can never be zero

         f(j_min:n_p) = psi_m(j_min:n_p)
         f(n_p+1:j_max) = psi_p(n_p+1:j_max)
      ELSEIF (n_p == j_max) THEN !we have psi_m for all j
         f(j_min:j_max) = psi_m(j_min:j_max)
      ELSEIF (n_p == j_min) THEN !we have psi_p up to j_min+1
         psi_p(j_min) = (-X_6j(j_min+1, j2, j3, l1, l2, l3) * psi_p(j_min+2) - Y_6j(j_min+1, j2, j3, l1, l2, l3) * psi_p(j_min+1)) &
                        / Z_6j(j_min+1, j2, j3, l1, l2, l3)
         f(j_min:j_max) = psi_p(j_min:j_max)
      ENDIF

   ELSE !from_b .and. from_a .eq..true.

      STOP "Error in Wigner_6j. Can't calculate the given coefficient using the simplified algorithm. &
           &Try three term recurrence only at this point?"

   ENDIF
   
   !calculate the normalization of the whole thing
   norm=0.0_cfp
   DO j=j_min,j_max
      norm=norm+(2.0_cfp*j+1)*f(j)**2.0_cfp
   ENDDO
   norm=(2*l1+1)*norm
   
   !normalize and correct the phase
   tmp=sign(1.0_cfp,f(j_max)) !this makes sure that I get rid of a possibly wrong sign of f(j_max)
!   PRINT *,'sign',tmp,j_min,j_max
   IF (tmp .ne. (-1)**(j2+j3+l2+l3)) THEN
      f(:)=-f(:)/sqrt(norm)
   ELSE
      f(:)=f(:)/sqrt(norm)
   ENDIF
!   PRINT *,f(j_min:j_max)
   
   DEALLOCATE(r,s,psi_m,psi_p)

END SUBROUTINE Wigner_6j

FUNCTION X_3j(j, j2, j3, m1, m2, m3)
!this is X for general f(j)
IMPLICIT NONE
INTEGER, INTENT(IN) :: j, j2, j3, m1, m2, m3
REAL(kind=cfp) :: X_3j

X_3j = j*sqrt(((j+1)**2.0_cfp-(j2-j3)**2)*((j2+j3+1)**2-(j+1)**2)*((j+1)**2-m1**2))

END FUNCTION X_3j

FUNCTION Y_3j(j, j2, j3, m1, m2, m3)
!this is Y for general f(j)
IMPLICIT NONE
INTEGER, INTENT(IN) :: j, j2, j3, m1, m2, m3
REAL(kind=cfp) :: Y_3j

Y_3j = -(2.0_cfp*j+1)*(m1*(j2*(j2+1)-j3*(j3+1))-(m3-m2)*j*(j+1))

END FUNCTION Y_3j

FUNCTION Z_3j(j, j2, j3, m1, m2, m3)
!this is Z for general f(j)
IMPLICIT NONE
INTEGER, INTENT(IN) :: j, j2, j3, m1, m2, m3
REAL(kind=cfp) :: Z_3j

Z_3j = (j+1)*sqrt((j**2.0_cfp-(j2-j3)**2)*((j2+j3+1)**2-j**2)*(j**2-m1**2))

END FUNCTION Z_3j

FUNCTION X_6j(j, j2, j3, l1, l2, l3)
!this is X for general f(j)
IMPLICIT NONE
INTEGER, INTENT(IN) :: j, j2, j3, l1, l2, l3
REAL(kind=cfp) :: X_6j

X_6j = j*sqrt(((j+1)**2.0_cfp-(j2-j3)**2)*((j2+j3+1)**2-(j+1)**2)*((j+1)**2-(l2-l3)**2)*((l2+l3+1)**2-(j+1)**2))

END FUNCTION X_6j

FUNCTION Y_6j(j, j2, j3, l1, l2, l3)
!this is Y for general f(j)
IMPLICIT NONE
INTEGER, INTENT(IN) :: j, j2, j3, l1, l2, l3
REAL(kind=cfp) :: Y_6j

Y_6j = (2.0_cfp*j+1)*( j * ( j + 1) * (-j * (j + 1) + j2 * (j2 + 1) + j3 * (j3 + 1) - 2 * l1 * (l1 + 1)) &
                    + l2 * (l2 + 1) * ( j * (j + 1) + j2 * (j2 + 1) - j3 * (j3 + 1)) &
                    + l3 * (l3 + 1) * ( j * (j + 1) - j2 * (j2 + 1) + j3 * (j3 + 1)))

END FUNCTION Y_6j

FUNCTION Z_6j(j, j2, j3, l1, l2, l3)
!this is Z for general f(j)
IMPLICIT NONE
INTEGER, INTENT(IN) :: j, j2, j3, l1, l2, l3
REAL(kind=cfp) :: Z_6j

Z_6j = (j+1)*sqrt((j**2.0_cfp-(j2-j3)**2)*((j2+j3+1)**2-j**2)*(j**2-(l2-l3)**2)*((l2+l3+1)**2-j**2))

END FUNCTION Z_6j

END MODULE wigner_cf
