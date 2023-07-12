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
module special_functions
use utils, only: xermsg
use precisn, only: cfp, wp, ep1, cfp_dummy, wp_dummy, ep_dummy, f1mach, i1mach

  implicit none

  interface cfp_gamma_fun
     module procedure wp_gamma_fun, ep_gamma_fun
  end interface

  !> \warning Quad precision version not implemented yet. 
  interface cfp_besi
     module procedure wp_besi, ep_besi
  end interface

  !> \warning Quad precision version not implemented yet.
  interface cfp_besj
     module procedure wp_besj, ep_besj
  end interface

  !> \warning Quad precision version not implemented yet.
  interface cfp_asyik
     module procedure wp_asyik, ep_asyik
  end interface

  !> \warning Quad precision version not implemented yet.
  interface cfp_asyjy
     module procedure wp_asyjy, ep_asyjy
  end interface

  !> \warning Quad precision version not implemented yet.
  interface cfp_jairy
     module procedure wp_jairy, ep_jairy
  end interface

  interface cfp_lngam
     module procedure wp_lngam, ep_lngam
  end interface

  interface cfp_9lgmc
     module procedure wp_9lgmc, ep_9lgmc
  end interface

  interface cfp_initds
     module procedure wp_initds, ep_initds
  end interface

  interface cfp_csevl
     module procedure wp_csevl, ep_csevl
  end interface

  interface cfp_binom
     module procedure wp_binom, ep_binom
  end interface

  interface cfp_gamlm
     module procedure wp_gamlm, ep_gamlm
  end interface

  interface cfp_gamma_slatec
     module procedure wp_gamma, ep_gamma
  end interface

  interface cfp_lnrel
     module procedure wp_lnrel, ep_lnrel
  end interface

  interface cfp_9gmic
     module procedure wp_9gmic, ep_9gmic
  end interface

  interface cfp_9gmit
     module procedure wp_9gmit, ep_9gmit
  end interface

  interface cfp_9lgic
     module procedure wp_9lgic, ep_9lgic
  end interface

  interface cfp_9lgit
     module procedure wp_9lgit, ep_9lgit
  end interface

  interface cfp_gamic
     module procedure wp_gamic, ep_gamic
  end interface

  interface cfp_lgams
     module procedure wp_lgams, ep_lgams
  end interface

  interface cfp_eval_poly_horner
     module procedure cfp_eval_poly_horner_single
  end interface

  private boys_function

  private

  public cfp_gamma_fun, cfp_besi, cfp_lngam, cfp_9lgmc, cfp_initds, cfp_csevl, cfp_binom, cfp_gamlm, cfp_lnrel
  public cfp_9gmic, cfp_9gmit, cfp_9lgic
  public cfp_9lgit, cfp_gamic, cfp_lgams, cfp_nlm, cfp_resh, cfp_solh, cfp_zhar, cfp_sph_to_cart_mapping
  public cfp_sph_shell_to_cart_shell, cfp_sph_shell_to_cart_lshells, cfp_eval_poly_horner
  public boys_function_quad, cfp_besj, ipair, unpack_pq
  public cfp_eval_poly_horner_many

  private cfp_asyjy, cfp_jairy, cfp_asyik

contains

   !> This function is used to index an ordered pair of values (i,j), where i .ge. j, i .ge. 1.
   !> The index of the ordered pair is: ipair(i) + j. This process of indexing can be used in a nested way to index
   !> quartets of integers, say (i,j,k,l). We compute ij = ipair(i)+j and kl = ipair(k)+l.
   !> The index of the quartet (i,j,k,l) is: ipair(ij)+kl. We assume that i.ge.j, k.ge.l, ij.ge.kl. This nesting
   !> (triangularization) is used heavilly to index the 2-particle symmetric integrals.
   elemental function ipair(i)
      implicit none
      integer, intent(in) :: i
      integer :: ipair

         ipair = i*(i-1)/2

   end function ipair

   !> Assuming \f$ pq = p + (q-1)*n \f$ this routine returns p,q given pq,n.
   subroutine unpack_pq(pq,n,p,q)
      implicit none
      integer, intent(in) :: pq, n
      integer, intent(out) :: p,q

         q = pq/n+1
         p = mod(pq,n)
         if (p .eq. 0) then
            q = q-1
            p = pq - n*(q-1)
         endif

   end subroutine unpack_pq

  !> The routine that calculates the gamma function; by default this is the intrinsic gamma(x) function, but it does not have
  !> to be intrinsic depending on the compiler used. Therefore if the intrisic gamma function is not available,
  !> this routine should be replaced by a call to a different implementation of this function, e.g. the cfp_gamma_slatec function
  !> from the slatec library (included). However, we also need for some programs the quad precision gamma, so this should be
  !> provided as well. The special functions from the netlib library (i.e. the ones in the dprec_routines directory) have not
  !> been converted to use the gamma_fun. The only routines that use it are my own routines.
  elemental function wp_gamma_fun(x)
     real(kind=wp), intent(in) :: x
     real(kind=wp) :: wp_gamma_fun

        wp_gamma_fun = gamma(x)

  end function wp_gamma_fun

  elemental function ep_gamma_fun(x)
     real(kind=ep1), intent(in) :: x
     real(kind=ep1) :: ep_gamma_fun

        ep_gamma_fun = gamma(x)

  end function ep_gamma_fun

!> \verbatim
!>***SUBSIDIARY
!>***PURPOSE  Compute the log Gamma correction factor so that
!>            LOG(cfp_gamma(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
!>            + wp_9lgmc(X).
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7E
!>***TYPE      REAL(kind=wp) (R9LGMC-S, wp_9lgmc-D, C9LGMC-C)
!>***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
!>             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> Compute the log gamma correction factor for X .GE. 10. so that
!> LOG (wp_gamma(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + wp_9lgmc(X)
!>
!> Series for ALGM       on the interval  0.          to  1.00000E-02
!>                                        with weighted error   1.28E-31
!>                                         log weighted error  30.89
!>                               significant figures required  29.81
!>                                    decimal places required  31.48
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, cfp_csevl, cfp_initds, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770601  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900720  Routine changed from user-callable to subsidiary.  (WRB)
!> \endverbatim
!
      REAL(kind=wp) FUNCTION wp_9lgmc (X)
      IMPLICIT NONE
!
      INTEGER :: NALGM
      REAL(kind=wp) X, ALGMCS(15), XBIG, XMAX
      LOGICAL FIRST
      SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
      DATA ALGMCS(  1) / +.1666389480451863247205729650822E+0_wp      /
      DATA ALGMCS(  2) / -.1384948176067563840732986059135E-4_wp      /
      DATA ALGMCS(  3) / +.9810825646924729426157171547487E-8_wp      /
      DATA ALGMCS(  4) / -.1809129475572494194263306266719E-10_wp     /
      DATA ALGMCS(  5) / +.6221098041892605227126015543416E-13_wp     /
      DATA ALGMCS(  6) / -.3399615005417721944303330599666E-15_wp     /
      DATA ALGMCS(  7) / +.2683181998482698748957538846666E-17_wp     /
      DATA ALGMCS(  8) / -.2868042435334643284144622399999E-19_wp     /
      DATA ALGMCS(  9) / +.3962837061046434803679306666666E-21_wp     /
      DATA ALGMCS( 10) / -.6831888753985766870111999999999E-23_wp     /
      DATA ALGMCS( 11) / +.1429227355942498147573333333333E-24_wp     /
      DATA ALGMCS( 12) / -.3547598158101070547199999999999E-26_wp     /
      DATA ALGMCS( 13) / +.1025680058010470912000000000000E-27_wp     /
      DATA ALGMCS( 14) / -.3401102254316748799999999999999E-29_wp     /
      DATA ALGMCS( 15) / +.1276642195630062933333333333333E-30_wp     /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_9lgmc
      IF (FIRST) THEN
         NALGM = cfp_initds (ALGMCS, 15, REAL(F1MACH(3,wp_dummy)) )
         XBIG = 1.0_wp/SQRT(F1MACH(3,wp_dummy))
         XMAX = EXP (MIN(LOG(F1MACH(2,wp_dummy)/12.0_wp), -LOG(12.0_wp*F1MACH(1,wp_dummy))))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LT. 10.0_wp) CALL XERMSG ('SLATEC', 'wp_9lgmc', 'X MUST BE GE 10', 1, 2)
      IF (X.GE.XMAX) GO TO 20
!
      wp_9lgmc = 1.0_wp/(12.0_wp*X)
      IF (X.LT.XBIG) wp_9lgmc = cfp_csevl (2.00_wp*(10.0_wp/X)**2-1.0_wp, ALGMCS, NALGM) / X
      RETURN
!
 20   wp_9lgmc = 0.0_wp
      CALL XERMSG ('SLATEC', 'wp_9lgmc', 'X SO BIG wp_9lgmc UNDERFLOWS', 2, 1)
      RETURN
!
      END FUNCTION
!
      !> Quad precision version of wp_9lgmc.
      REAL(kind=ep1) FUNCTION ep_9lgmc (X)
      IMPLICIT NONE
!
      INTEGER, parameter :: algmcs_len = 36
      REAL(kind=ep1), parameter :: ALGMCS(1:algmcs_len) = (/ &
        1.666389480451863247205729650822634E-01_ep1,-1.384948176067563840732986059135685E-05_ep1, &
        9.810825646924729426157171547494646E-09_ep1,-1.809129475572494194263306266228199E-11_ep1, &
        6.221098041892605227126015425026470E-14_ep1,-3.399615005417721944303239312834238E-16_ep1, &
        2.683181998482698748954316939042279E-18_ep1,-2.868042435334643284529047659564974E-20_ep1, &
        3.962837061046434902918518354024164E-22_ep1,-6.831888753985773031044337598163979E-24_ep1, &
        1.429227355942400147750478163327568E-25_ep1,-3.547598158100687420273243440310758E-27_ep1, &
        1.025680057924723603180821548541475E-28_ep1,-3.401102235191009965917634926326868E-30_ep1, &
        1.276642119142403037183509657588180E-31_ep1,-5.363992996332298948751802756889457E-33_ep1, &
        2.498419748615974137373733023216102E-34_ep1,-1.279180004435981239791780867519034E-35_ep1, &
        7.146084883638207209714256656556028E-37_ep1,-4.327468162848449421867897638599174E-38_ep1, &
        2.824278083799969375380195561854690E-39_ep1,-1.976242554861509518781751528816163E-40_ep1, &
        1.475767033968755183800242001951694E-41_ep1,-1.171191914267258725710752170105842E-42_ep1, &
        9.840900425245182857831150498791323E-44_ep1,-8.724801246834106040983336496233320E-45_ep1, &
        8.136589380331286133540745054384699E-46_ep1,-7.959145144741623801395941296680164E-47_ep1, &
        8.145237413333173981221657997512119E-48_ep1,-8.700092109564861725122326075938473E-49_ep1, &
        9.677861764016593924355835046119559E-50_ep1,-1.118916060710528009572203837277006E-50_ep1, &
        1.342059788020631580157540218713364E-51_ep1,-1.667068284935961661354922920118027E-52_ep1, &
        2.140417469185841158842409837058527E-53_ep1,-2.784806354790998375237622438158237E-54_ep1  &
      /)
      INTEGER :: NALGM
      REAL(kind=ep1) X, XBIG, XMAX
      LOGICAL FIRST
      SAVE NALGM, XBIG, XMAX, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_9lgmc
      IF (FIRST) THEN
         NALGM = cfp_initds (ALGMCS, algmcs_len, REAL(F1MACH(3,ep_dummy)) )
         XBIG = 1.0_ep1/SQRT(F1MACH(3,ep_dummy))
         XMAX = EXP (MIN(LOG(F1MACH(2,ep_dummy)/12.0_ep1), -LOG(12.0_ep1*F1MACH(1,ep_dummy))))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LT. 10.0_ep1) CALL XERMSG ('SLATEC', 'ep_9lgmc', 'X MUST BE GE 10', 1, 2)
      IF (X.GE.XMAX) GO TO 20
!
      ep_9lgmc = 1.0_ep1/(12.0_ep1*X)
      IF (X.LT.XBIG) ep_9lgmc = cfp_csevl (2.00_ep1*(10.0_ep1/X)**2-1.0_ep1, ALGMCS, NALGM) / X
      RETURN
!
 20   ep_9lgmc = 0.0_ep1
      CALL XERMSG ('SLATEC', 'ep_9lgmc', 'X SO BIG ep_9lgmc UNDERFLOWS', 2, 1)
      RETURN
!
      END FUNCTION
!> \verbatim
!>***SUBSIDIARY
!>***PURPOSE  Subsidiary to cfp_besi and DBESK
!>***LIBRARY   SLATEC
!>***TYPE      REAL(kind=wp) (ASYIK-S, wp_asyik-D)
!>***AUTHOR  Amos, D. E., (SNLA)
!>***DESCRIPTION
!>
!>                    wp_asyik computes Bessel functions I and K
!>                  for arguments X.GT.0.0 and orders FNU.GE.35
!>                  on FLGIK = 1 and FLGIK = -1 respectively.
!>
!>                                    INPUT
!>
!>      X    - Argument, X.GT.0.00_wp
!>      FNU  - Order of first Bessel function
!>      KODE - A parameter to indicate the scaling option
!>             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN
!>                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN
!>                    on FLGIK = 1.00_wp or FLGIK = -1.00_wp
!>             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN
!>                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN
!>                    on FLGIK = 1.00_wp or FLGIK = -1.00_wp
!>     FLGIK - Selection parameter for I or K FUNCTION
!>             FLGIK =  1.00_wp gives the I function
!>             FLGIK = -1.00_wp gives the K function
!>        RA - SQRT(1.+Z*Z), Z=X/FNU
!>       ARG - Argument of the leading exponential
!>        IN - Number of functions desired, IN=1 or 2
!>
!>                                    OUTPUT
!>
!>         Y - A vector whose first IN components contain the sequence
!>
!>     Abstract  **** A REAL(kind=wp) routine ****
!>         wp_asyik implements the uniform asymptotic expansion of
!>         the I and K Bessel functions for FNU.GE.35 and real
!>         X.GT.0.00_wp. The forms are identical except for a change
!>         in sign of some of the terms. This change in sign is
!>         accomplished by means of the FLAG FLGIK = 1 or -1.
!>
!>***SEE ALSO  cfp_besi, DBESK
!>***ROUTINES CALLED  F1MACH
!>***REVISION HISTORY  (YYMMDD)
!>   750101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>   910408  Updated the AUTHOR section.  (WRB)
!> \endverbatim
!
      SUBROUTINE wp_asyik (X, FNU, KODE, FLGIK, RA, ARG, IN, Y)
      IMPLICIT NONE
!
!
      INTEGER IN, J, JN, K, KK, KODE, L
      REAL(kind=wp) AK,AP,ARG,C,COEF,CON,ETX,FLGIK,FN,FNU,GLN,RA, S1, S2, T, TOL, T2, X, Y, Z
      DIMENSION Y(*), C(65), CON(2)
      SAVE CON, C
      DATA CON(1), CON(2)  /3.98942280401432678E-01_wp,    1.25331413731550025E+00_wp/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),&
     &     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),&
     &     C(19), C(20), C(21), C(22), C(23), C(24)/&
     &       -2.08333333333333E-01_wp,        1.25000000000000E-01_wp,&
     &        3.34201388888889E-01_wp,       -4.01041666666667E-01_wp,&
     &        7.03125000000000E-02_wp,       -1.02581259645062E+00_wp,&
     &        1.84646267361111E+00_wp,       -8.91210937500000E-01_wp,&
     &        7.32421875000000E-02_wp,        4.66958442342625E+00_wp,&
     &       -1.12070026162230E+01_wp,        8.78912353515625E+00_wp,&
     &       -2.36408691406250E+00_wp,        1.12152099609375E-01_wp,&
     &       -2.82120725582002E+01_wp,        8.46362176746007E+01_wp,&
     &       -9.18182415432400E+01_wp,        4.25349987453885E+01_wp,&
     &       -7.36879435947963E+00_wp,        2.27108001708984E-01_wp,&
     &        2.12570130039217E+02_wp,       -7.65252468141182E+02_wp,&
     &        1.05999045252800E+03_wp,       -6.99579627376133E+02_wp/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),&
     &     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),&
     &     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/&
     &        2.18190511744212E+02_wp,       -2.64914304869516E+01_wp,&
     &        5.72501420974731E-01_wp,       -1.91945766231841E+03_wp,&
     &        8.06172218173731E+03_wp,       -1.35865500064341E+04_wp,&
     &        1.16553933368645E+04_wp,       -5.30564697861340E+03_wp,&
     &        1.20090291321635E+03_wp,       -1.08090919788395E+02_wp,&
     &        1.72772750258446E+00_wp,        2.02042913309661E+04_wp,&
     &       -9.69805983886375E+04_wp,        1.92547001232532E+05_wp,&
     &       -2.03400177280416E+05_wp,        1.22200464983017E+05_wp,&
     &       -4.11926549688976E+04_wp,        7.10951430248936E+03_wp,&
     &       -4.93915304773088E+02_wp,        6.07404200127348E+00_wp,&
     &       -2.42919187900551E+05_wp,        1.31176361466298E+06_wp,&
     &       -2.99801591853811E+06_wp,        3.76327129765640E+06_wp/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),&
     &     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),&
     &     C(65)/&
     &       -2.81356322658653E+06_wp,        1.26836527332162E+06_wp,&
     &       -3.31645172484564E+05_wp,        4.52187689813627E+04_wp,&
     &       -2.49983048181121E+03_wp,        2.43805296995561E+01_wp,&
     &        3.28446985307204E+06_wp,       -1.97068191184322E+07_wp,&
     &        5.09526024926646E+07_wp,       -7.41051482115327E+07_wp,&
     &        6.63445122747290E+07_wp,       -3.75671766607634E+07_wp,&
     &        1.32887671664218E+07_wp,       -2.78561812808645E+06_wp,&
     &        3.08186404612662E+05_wp,       -1.38860897537170E+04_wp,&
     &        1.10017140269247E+02_wp/
!***FIRST EXECUTABLE STATEMENT  wp_asyik
      TOL = F1MACH(3,wp_dummy)
      TOL = MAX(TOL,1.0E-15_wp) !1.0E-15_wp should be replaced by 1.0E-33_wp in quad precision code
      FN = FNU
      Z  = (3.00_wp-FLGIK)/2.00_wp
      KK = INT(Z)
      DO 50 JN=1,IN
        IF (JN.EQ.1) GO TO 10
        FN = FN - FLGIK
        Z = X/FN
        RA = SQRT(1.00_wp+Z*Z)
        GLN = LOG((1.00_wp+RA)/Z)
        ETX = KODE - 1
        T = RA*(1.00_wp-ETX) + ETX/(Z+RA)
        ARG = FN*(T-GLN)*FLGIK
   10   COEF = EXP(ARG)
        T = 1.00_wp/RA
        T2 = T*T
        T = T/FN
        T = SIGN(T,FLGIK)
        S2 = 1.00_wp
        AP = 1.00_wp
        L = 0
        DO 30 K=2,11
          L = L + 1
          S1 = C(L)
          DO 20 J=2,K
            L = L + 1
            S1 = S1*T2 + C(L)
   20     CONTINUE
          AP = AP*T
          AK = AP*S1
          S2 = S2 + AK
          IF (MAX(ABS(AK),ABS(AP)) .LT.TOL) GO TO 40
   30   CONTINUE
   40   CONTINUE
      T = ABS(T)
      Y(JN) = S2*COEF*SQRT(T)*CON(KK)
   50 CONTINUE
      RETURN
      END SUBROUTINE
      !> Quad prec version of wp_asyik tranlated to F90 standard.
      SUBROUTINE ep_asyik (X, FNU, KODE, FLGIK, RA, ARG, IN, Y)
      IMPLICIT NONE
!
!
      !CON(1) = sqrt(1/(2*pi)); CON(2) = sqrt(pi/2)
      REAL(kind=ep1), parameter :: CON(1:2) = (/0.3989422804014326779399460599343819_ep1, 1.253314137315500251207882642405523_ep1/)
!     Coefficients for the uniform asymptotic expansion of I,K (see eq. 9.7.7. in Abramowitz, Stegun and the paper by D.E.Amos references in the wp version of this routine).
      INTEGER, parameter :: kmax = 22 !Maximum k for which the coefficients below are given. Note that for wp precision this was 11. The quad version value 22 was found using experiments in Mathematica.
      REAL(kind=ep1), parameter :: C(1:275) = (/ &
        -2.083333333333333333333333333333333E-01_ep1, 1.250000000000000000000000000000000E-01_ep1, &
         3.342013888888888888888888888888889E-01_ep1,-4.010416666666666666666666666666667E-01_ep1, &
         7.031250000000000000000000000000000E-02_ep1,-1.025812596450617283950617283950617E+00_ep1, &
         1.846462673611111111111111111111111E+00_ep1,-8.912109375000000000000000000000000E-01_ep1, &
         7.324218750000000000000000000000000E-02_ep1, 4.669584423426247427983539094650206E+00_ep1, &
        -1.120700261622299382716049382716049E+01_ep1, 8.789123535156250000000000000000000E+00_ep1, &
        -2.364086914062500000000000000000000E+00_ep1, 1.121520996093750000000000000000000E-01_ep1, &
        -2.821207255820024487740054869684499E+01_ep1, 8.463621767460073463220164609053498E+01_ep1, &
        -9.181824154324001736111111111111111E+01_ep1, 4.253499874538845486111111111111111E+01_ep1, &
        -7.368794359479631696428571428571429E+00_ep1, 2.271080017089843750000000000000000E-01_ep1, &
         2.125701300392171228609694120560890E+02_ep1,-7.652524681411816422994898834019204E+02_ep1, &
         1.059990452527999877929687500000000E+03_ep1,-6.995796273761325412326388888888889E+02_ep1, &
         2.181905117442115904792906746031746E+02_ep1,-2.649143048695155552455357142857143E+01_ep1, &
         5.725014209747314453125000000000000E-01_ep1,-1.919457662318406996310063083863613E+03_ep1, &
         8.061722181737309384502264952227176E+03_ep1,-1.358655000643413743855040750385802E+04_ep1, &
         1.165539333686453324777108651620370E+04_ep1,-5.305646978613403108384874131944444E+03_ep1, &
         1.200902913216352462768554687500000E+03_ep1,-1.080909197883946555001395089285714E+02_ep1, &
         1.727727502584457397460937500000000E+00_ep1, 2.020429133096614864345123694004355E+04_ep1, &
        -9.698059838863751348856593731220906E+04_ep1, 1.925470012325315323590578202193984E+05_ep1, &
        -2.034001772804155342781658198771326E+05_ep1, 1.222004649830174597877043264883536E+05_ep1, &
        -4.119265496889755129814147949218750E+04_ep1, 7.109514302489363721438816615513393E+03_ep1, &
        -4.939153047730880124228341238839286E+02_ep1, 6.074042001273483037948608398437500E+00_ep1, &
        -2.429191879005513334585317700615422E+05_ep1, 1.311763614662977200676071558332328E+06_ep1, &
        -2.998015918538106750091346203054420E+06_ep1, 3.763271297656403996402105622276303E+06_ep1, &
        -2.813563226586534110707868355618910E+06_ep1, 1.268365273321624781625966231028239E+06_ep1, &
        -3.316451724845635778315010524931408E+05_ep1, 4.521876898136272627328123365129743E+04_ep1, &
        -2.499830481811209624125198884443803E+03_ep1, 2.438052969955606386065483093261719E+01_ep1, &
         3.284469853072037821137231641040435E+06_ep1,-1.970681911843222692682338984624261E+07_ep1, &
         5.095260249266464220638182198049918E+07_ep1,-7.410514821153265774833562096441470E+07_ep1, &
         6.634451227472902666479879845432838E+07_ep1,-3.756717666076335130816319796406193E+07_ep1, &
         1.328876716642181832943741163169896E+07_ep1,-2.785618128086454688959444562594096E+06_ep1, &
         3.081864046126623984803907842771021E+05_ep1,-1.388608975371704053197225386446173E+04_ep1, &
         1.100171402692467381712049245834351E+02_ep1,-4.932925366450996197276183127547471E+07_ep1, &
         3.255730741857657490202280864181331E+08_ep1,-9.394623596815784025462443009203891E+08_ep1, &
         1.553596899570580056158121044387996E+09_ep1,-1.621080552108337075248175882636769E+09_ep1, &
         1.106842816823014468259666669096246E+09_ep1,-4.958897842750303092546362453742585E+08_ep1, &
         1.420629077975330951856532785179162E+08_ep1,-2.447406272573872846781300815601586E+07_ep1, &
         2.243768177922449429230737780239813E+06_ep1,-8.400543360302408528867828125668156E+04_ep1, &
         5.513358961220205856079701334238052E+02_ep1, 8.147890961183121149459306645049764E+08_ep1, &
        -5.866481492051847227610700784435830E+09_ep1, 1.868820750929582492236591930279163E+10_ep1, &
        -3.463204338815877792290241333559551E+10_ep1, 4.128018557975397395513147102709743E+10_ep1, &
        -3.302659974980072314009099267577854E+10_ep1, 1.795421373115560008015220585382804E+10_ep1, &
        -6.563293792619284332035016850974713E+09_ep1, 1.559279864879257513349646204741952E+09_ep1, &
        -2.251056618894152778040714269630530E+08_ep1, 1.739510755397816453810439631423674E+07_ep1, &
        -5.498423275722886871349019329372950E+05_ep1, 3.038090510922384268610585422720760E+03_ep1, &
        -1.467926124769561666061242392686690E+10_ep1, 1.144982377320258099527769066295618E+11_ep1, &
        -3.990961752244664979552346232418551E+11_ep1, 8.192186695485773286413033254921622E+11_ep1, &
        -1.098375156081223306827064535415196E+12_ep1, 1.008158106865382094769125164830028E+12_ep1, &
        -6.453648692453765032808836899474692E+11_ep1, 2.879006499061505887229132920572018E+11_ep1, &
        -8.786707217802326567663590419006899E+10_ep1, 1.763473060683496938315197395758436E+10_ep1, &
        -2.167164983223795093518416127711201E+09_ep1, 1.431578767188889812910572701178294E+08_ep1, &
        -3.871833442572612620626626666931147E+06_ep1, 1.825775547429317469116938355000457E+04_ep1, &
         2.864640357176790429870109038347210E+11_ep1,-2.406297900028503961090891592211656E+12_ep1, &
         9.109341185239898955907876541296591E+12_ep1,-2.051689941093443739076047796919586E+13_ep1, &
         3.056512551993532061172003688280095E+13_ep1,-3.166708858478515840255256788568800E+13_ep1, &
         2.334836404458184093765746780272027E+13_ep1,-1.232049130559828715978770065317654E+13_ep1, &
         4.612725780849131966803816033578737E+12_ep1,-1.196552880196181598974160686063231E+12_ep1, &
         2.059145032324100156890817255326430E+11_ep1,-2.182292775752922372939877564974734E+10_ep1, &
         1.247009293512710324825868373805287E+09_ep1,-2.918838812222081340342732031993895E+07_ep1, &
         1.188384262567832531237721482852976E+05_ep1,-6.019723417234005444990937465304623E+12_ep1, &
         5.417751075510604900491843718774161E+13_ep1,-2.213496387025251959655937979407546E+14_ep1, &
         5.427396649876597227020591239685814E+14_ep1,-8.894969398810264418128257191774092E+14_ep1, &
         1.026955196082762488813740580594022E+15_ep1,-8.574610329828950513961987089197640E+14_ep1, &
         5.230548825784446555790535196170039E+14_ep1,-2.326048311889399252321748606017756E+14_ep1, &
         7.437312290867914494114728953717076E+13_ep1,-1.663482472489248051865692506018648E+13_ep1, &
         2.485000928034085323647452365637966E+12_ep1,-2.296193729682464681659534773231422E+11_ep1, &
         1.146575489944823715692235895954002E+10_ep1,-2.345579635222515247762632483406070E+08_ep1, &
         8.328593040162892989757698058994606E+05_ep1, 1.355221587030936902915277458009335E+14_ep1, &
        -1.301012723549699426798666359688962E+15_ep1, 5.705782159023670809618694505693765E+15_ep1, &
        -1.512982632245768118084636116034478E+16_ep1, 2.705471130619708124101419805879396E+16_ep1, &
        -3.444722600648514469779708309880293E+16_ep1, 3.213827526858624120000619279550618E+16_ep1, &
        -2.226822513391114256219382687736836E+16_ep1, 1.148670697844975210969241162584644E+16_ep1, &
        -4.379325838364015437780098517928750E+15_ep1, 1.212675804250347416525907257338235E+15_ep1, &
        -2.366525304516492516817769490479771E+14_ep1, 3.100743647289646141719069923627260E+13_ep1, &
        -2.521558474912854621312538497418401E+12_ep1, 1.109974051391790127937406706965332E+11_ep1, &
        -2.001646928191776331529938815671142E+09_ep1, 6.252951493434797002466521745854544E+06_ep1, &
        -3.254192619642668832809062072577808E+15_ep1, 3.319276472035522209465243314029364E+16_ep1, &
        -1.555298350431390256212648930015349E+17_ep1, 4.434795461417190406002566704378017E+17_ep1, &
        -8.592577980317547990581328866810682E+17_ep1, 1.196199114275630785068459026679820E+18_ep1, &
        -1.233611693196069502238697805757509E+18_ep1, 9.575335098169138663533895519883678E+17_ep1, &
        -5.626317880746360283949116996609716E+17_ep1, 2.496036512616042570994262490257148E+17_ep1, &
        -8.270945651585064278725937951200775E+16_ep1, 2.006427147630953080010052087140459E+16_ep1, &
        -3.450385511846272492011831924964665E+15_ep1, 4.000444570430362415133450812654904E+14_ep1, &
        -2.886383763141476025414316301471395E+13_ep1, 1.128709145410874078578624908477651E+12_ep1, &
        -1.807822038465806371713485435332554E+10_ep1, 5.006958953198892599769148662673234E+07_ep1, &
         8.301957606731910464441822477287042E+16_ep1,-8.966114215270463301597168275470005E+17_ep1, &
         4.470200964012310169294212039718634E+18_ep1,-1.363942041057159065682587128977512E+19_ep1, &
         2.846521225167657097650533573844960E+19_ep1,-4.301555703831443742343849559927269E+19_ep1, &
         4.859942729324835775153498734270616E+19_ep1,-4.178861444656838881754858156644904E+19_ep1, &
         2.757282981650518864947605602180672E+19_ep1,-1.397080351644337385472472411321626E+19_ep1, &
         5.402894876715981887221861297046391E+18_ep1,-1.573643476518959871900805130320975E+18_ep1, &
         3.376676249790609622988679489577904E+17_ep1,-5.160509319348522743652109163378487E+16_ep1, &
         5.335106978708838675506690952489470E+15_ep1,-3.439653047430759474698419168725856E+14_ep1, &
         1.203011582641919172809950344616877E+13_ep1,-1.722832387173504987359310146915852E+11_ep1, &
         4.259392165047669051886949383176883E+08_ep1,-2.242438856186775026108112444139134E+18_ep1, &
         2.556380296052923529763248186318612E+19_ep1,-1.351217503435996111683396148700986E+20_ep1, &
         4.392792200888712002497385246154331E+20_ep1,-9.824438427689858246661460629133877E+20_ep1, &
         1.601689857369359736514880523609124E+21_ep1,-1.967724707705312458948384730248364E+21_ep1, &
         1.857108932146345179545529847779226E+21_ep1,-1.360203777284994087313165869102075E+21_ep1, &
         7.756704953461136792953564435694892E+20_ep1,-3.434621399768416893167721964613248E+20_ep1, &
         1.170749053579725885376371165883123E+20_ep1,-3.025566598990372035718148898948393E+19_ep1, &
         5.789887667664653131092223684204065E+18_ep1,-7.921651119323832137067359486449903E+17_ep1, &
         7.351663610930970405128460344287513E+16_ep1,-4.261935510426898338177749237549401E+15_ep1, &
         1.341241691518063854324417782641181E+14_ep1,-1.727704012352999522442090987607040E+12_ep1, &
         3.836255180230433507916601122084969E+09_ep1, 6.393286613940836715060316416259176E+19_ep1, &
        -7.671943936729004058072379699511011E+20_ep1, 4.285296082829493950777900479572493E+21_ep1, &
        -1.478774352843361445883955678341307E+22_ep1, 3.528435843903409379223597935378111E+22_ep1, &
        -6.173206302884414597368837221789481E+22_ep1, 8.194331005435129643139474666294128E+22_ep1, &
        -8.423222750084322624731938520667074E+22_ep1, 6.783661642951883229678547204155410E+22_ep1, &
        -4.302534303482378471023824962121826E+22_ep1, 2.148741481505588275526310883859541E+22_ep1, &
        -8.405915817108350448584740956999061E+21_ep1, 2.548961114664971585268545305016171E+21_ep1, &
        -5.891794135069496380504705156059513E+20_ep1, 1.012677416953659245416131770457788E+20_ep1, &
        -1.248370099504723315233152656418892E+19_ep1, 1.046172113113434395507698714401890E+18_ep1, &
        -5.484033603883289655520138802579367E+16_ep1, 1.561312393048467278412079969121827E+15_ep1, &
        -1.818726203851103723856933168276792E+13_ep1, 3.646840080706555853463218941682024E+10_ep1, &
        -1.918620238806649907049350908649207E+21_ep1, 2.417461500896378882882182144898001E+22_ep1, &
        -1.422839482332141380896819974811198E+23_ep1, 5.194289094766812226908110353617543E+23_ep1, &
        -1.317096961809238583721848995108535E+24_ep1, 2.461506085403875122901808631476913E+24_ep1, &
        -3.511096528332644078960709357969561E+24_ep1, 3.905264103536984928826687105359672E+24_ep1, &
        -3.430872898515745847695634879628826E+24_ep1, 2.396723774435168338766083562293517E+24_ep1, &
        -1.333717890779830224713540284497184E+24_ep1, 5.896543461978244771497081903622794E+23_ep1, &
        -2.056614913627154329822225073511004E+23_ep1, 5.591591380366263143498749296465660E+22_ep1, &
        -1.164024646146536927974080437089376E+22_ep1, 1.808159405713194358474659756990244E+21_ep1, &
        -2.019733541930087336814753176424288E+20_ep1, 1.536502521844337298059386429447268E+19_ep1, &
        -7.319501491566133145629484752049250E+17_ep1, 1.894406984252143386289224691111005E+16_ep1, &
        -2.005244012362711215413005056906471E+14_ep1, 3.649010818849833565280756572004454E+11_ep1, &
         6.045470627467089868102282399090314E+22_ep1,-7.980021228256558625895012766799214E+23_ep1, &
         4.936185283790662299213549866598326E+24_ep1,-1.900680753566443321252444955261614E+25_ep1, &
         5.103920268388801657606762578711132E+25_ep1,-1.014804898276639585373698271981776E+26_ep1, &
         1.548092083577385108248351501295410E+26_ep1,-1.852673104154991739253362836347395E+26_ep1, &
         1.763571327232664474625195159063481E+26_ep1,-1.345919399455641577189736905850979E+26_ep1, &
         8.262585357989550245211694677096963E+25_ep1,-4.077501349206541341009894674990116E+25_ep1, &
         1.610312854113731522960944257174675E+25_ep1,-5.046359865254400339044848638105650E+24_ep1, &
         1.238524103792451951435606681582903E+24_ep1,-2.336107524486965003556852475613820E+23_ep1, &
         3.297557757461477698552299619447756E+22_ep1,-3.354468912222678442775948942599043E+21_ep1, &
         2.327534625808941314633500392129141E+20_ep1,-1.012181837994208883276041022357394E+19_ep1, &
         2.392028012026999584409482284997005E+17_ep1,-2.310915976132356556501437801947091E+15_ep1, &
         3.833534661393944467161431194111497E+12_ep1 /)
!
      REAL(kind=ep1), intent(in) :: X, FNU, FLGIK
      INTEGER, intent(in) :: KODE, IN
      REAL(kind=ep1), intent(out) :: RA, ARG, Y(*)
!
      INTEGER J, JN, K, KK, L
      REAL(kind=ep1) AK,AP,COEF,ETX,FN,GLN, S1, S2, T, TOL, T2, Z
!***FIRST EXECUTABLE STATEMENT  ep_asyik
      TOL = F1MACH(3,ep_dummy)
      TOL = MAX(TOL,1.0E-33_ep1)
      FN = FNU
      Z  = (3.00_ep1-FLGIK)/2.00_ep1
      KK = INT(Z)
      DO JN=1,IN
        IF (JN.EQ.1) GO TO 10
        FN = FN - FLGIK
        Z = X/FN
        RA = SQRT(1.00_ep1+Z*Z)
        GLN = LOG((1.00_ep1+RA)/Z)
        ETX = KODE - 1
        T = RA*(1.00_ep1-ETX) + ETX/(Z+RA)
        ARG = FN*(T-GLN)*FLGIK
   10   COEF = EXP(ARG)
        T = 1.00_ep1/RA
        T2 = T*T
        T = T/FN
        T = SIGN(T,FLGIK)
        S2 = 1.00_ep1
        AP = 1.00_ep1
        L = 0
        DO K=2,kmax
          L = L + 1
          S1 = C(L)
          DO J=2,K
            L = L + 1
            S1 = S1*T2 + C(L)
          ENDDO
          AP = AP*T
          AK = AP*S1
          S2 = S2 + AK
          IF (MAX(ABS(AK),ABS(AP)) .LT.TOL) EXIT
        ENDDO
      T = ABS(T)
      Y(JN) = S2*COEF*SQRT(T)*CON(KK)
      ENDDO
      END SUBROUTINE
!> \verbatim
!>***SUBSIDIARY
!>***PURPOSE  Subsidiary to cfp_besj and DBESY
!>***LIBRARY   SLATEC
!>***TYPE      REAL(kind=wp) (ASYJY-S, wp_asyjy-D)
!>***AUTHOR  Amos, D. E., (SNLA)
!>***DESCRIPTION
!>
!>                 wp_asyjy computes Bessel functions J and Y
!>               for arguments X.GT.0.0 and orders FNU .GE. 35.0
!>               on FLGJY = 1 and FLGJY = -1 respectively
!>
!>                                  INPUT
!>
!>      FUNJY - External subroutine JAIRY or YAIRY
!>          X - Argument, X.GT.0.00_wp
!>        FNU - Order of the first Bessel function
!>      FLGJY - Selection flag
!>              FLGJY =  1.00_wp gives the J function
!>              FLGJY = -1.00_wp gives the Y function
!>         IN - Number of functions desired, IN = 1 or 2
!>
!>                                  OUTPUT
!>
!>         Y  - A vector whose first IN components contain the sequence
!>       IFLW - A flag indicating underflow or overflow
!>                    return variables for BESJ only
!>      WK(1) = 1 - (X/FNU)**2 = W**2
!>      WK(2) = SQRT(ABS(WK(1)))
!>      WK(3) = ABS(WK(2) - ATAN(WK(2)))  or
!>              ABS(LN((1 + WK(2))/(X/FNU)) - WK(2))
!>            = ABS((2/3)*ZETA**(3/2))
!>      WK(4) = FNU*WK(3)
!>      WK(5) = (1.5*WK(3)*FNU)**(1/3) = SQRT(ZETA)*FNU**(1/3)
!>      WK(6) = SIGN(1.,W**2)*WK(5)**2 = SIGN(1.,W**2)*ZETA*FNU**(2/3)
!>      WK(7) = FNU**(1/3)
!>
!>     Abstract   **** A REAL(kind=wp) Routine ****
!>         wp_asyjy implements the uniform asymptotic expansion of
!>         the J and Y Bessel functions for FNU.GE.35 and real
!>         X.GT.0.00_wp. The forms are identical except for a change
!>         in sign of some of the terms. This change in sign is
!>         accomplished by means of the flag FLGJY = 1 or -1. On
!>         FLGJY = 1 the Airy functions AI(X) and DAI(X) are
!>         supplied by the external function JAIRY, and on
!>         FLGJY = -1 the Airy functions BI(X) and DBI(X) are
!>         supplied by the external function YAIRY.
!>
!>***SEE ALSO  cfp_besj, DBESY
!>***ROUTINES CALLED  F1MACH, I1MACH
!>***REVISION HISTORY  (YYMMDD)
!>   750101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   891004  Correction computation of ELIM.  (WRB)
!>   891009  Removed unreferenced variable.  (WRB)
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>   910408  Updated the AUTHOR section.  (WRB)
!> \endverbatim
!
      SUBROUTINE wp_asyjy (FUNJY, X, FNU, FLGJY, IN, Y, WK, IFLW)
      IMPLICIT NONE
!
      EXTERNAL FUNJY
      INTEGER I, IFLW, IN, J, JN,JR,JU,K, KB,KLAST,KMAX,KP1, KS, KSP1, KSTEMP, L, LR, LRP1, ISETA, ISETB
      REAL(kind=wp) ABW2, AKM, ALFA, ALFA1, ALFA2, AP, AR, ASUM, AZ, &
     & BETA, BETA1, BETA2, BETA3, BR, BSUM, C, CON1, CON2,              &
     & CON548,CR,CRZ32, DFI,ELIM, DR,FI, FLGJY, FN, FNU,                &
     & FN2, GAMA, PHI,  RCZ, RDEN, RELB, RFN2,  RTZ, RZDEN,             &
     & SA, SB, SUMA, SUMB, S1, TA, TAU, TB, TFN, TOL, TOLS, T2, UPOL,   &
     &  WK, X, XX, Y, Z, Z32
      DIMENSION Y(*), WK(*), C(65)
      DIMENSION ALFA(26,4), BETA(26,5)
      DIMENSION ALFA1(26,2), ALFA2(26,2)
      DIMENSION BETA1(26,2), BETA2(26,2), BETA3(26,1)
      DIMENSION GAMA(26), KMAX(5), AR(8), BR(10), UPOL(10)
      DIMENSION CR(10), DR(10)
      EQUIVALENCE (ALFA(1,1),ALFA1(1,1))
      EQUIVALENCE (ALFA(1,3),ALFA2(1,1))
      EQUIVALENCE (BETA(1,1),BETA1(1,1))
      EQUIVALENCE (BETA(1,3),BETA2(1,1))
      EQUIVALENCE (BETA(1,5),BETA3(1,1))
      SAVE TOLS, CON1, CON2, CON548, AR, BR, C, ALFA1, ALFA2, BETA1, BETA2, BETA3, GAMA
      DATA TOLS            /-6.90775527898214E+00_wp/
      DATA CON1,CON2,CON548/6.66666666666667E-01_wp, 3.33333333333333E-01_wp, 1.04166666666667E-01_wp/
      DATA  AR(1),  AR(2),  AR(3),  AR(4),  AR(5),  AR(6),  AR(7), AR(8)&
     & / 8.35503472222222E-02_wp, 1.28226574556327E-01_wp,&
     & 2.91849026464140E-01_wp, 8.81627267443758E-01_wp, 3.32140828186277E+00_wp,&
     & 1.49957629868626E+01_wp, 7.89230130115865E+01_wp, 4.74451538868264E+02_wp/
      DATA  BR(1), BR(2), BR(3), BR(4), BR(5), BR(6), BR(7), BR(8),&
     &      BR(9), BR(10)  /-1.45833333333333E-01_wp,-9.87413194444444E-02_wp,&
     &-1.43312053915895E-01_wp,-3.17227202678414E-01_wp,-9.42429147957120E-01_wp,&
     &-3.51120304082635E+00_wp,-1.57272636203680E+01_wp,-8.22814390971859E+01_wp,&
     &-4.92355370523671E+02_wp,-3.31621856854797E+03_wp/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10),&
     &     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18),&
     &     C(19), C(20), C(21), C(22), C(23), C(24)/&
     &       -2.08333333333333E-01_wp,        1.25000000000000E-01_wp,&
     &        3.34201388888889E-01_wp,       -4.01041666666667E-01_wp,&
     &        7.03125000000000E-02_wp,       -1.02581259645062E+00_wp,&
     &        1.84646267361111E+00_wp,       -8.91210937500000E-01_wp,&
     &        7.32421875000000E-02_wp,        4.66958442342625E+00_wp,&
     &       -1.12070026162230E+01_wp,        8.78912353515625E+00_wp,&
     &       -2.36408691406250E+00_wp,        1.12152099609375E-01_wp,&
     &       -2.82120725582002E+01_wp,        8.46362176746007E+01_wp,&
     &       -9.18182415432400E+01_wp,        4.25349987453885E+01_wp,&
     &       -7.36879435947963E+00_wp,        2.27108001708984E-01_wp,&
     &        2.12570130039217E+02_wp,       -7.65252468141182E+02_wp,&
     &        1.05999045252800E+03_wp,       -6.99579627376133E+02_wp/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32),&
     &     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40),&
     &     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/&
     &        2.18190511744212E+02_wp,       -2.64914304869516E+01_wp,&
     &        5.72501420974731E-01_wp,       -1.91945766231841E+03_wp,&
     &        8.06172218173731E+03_wp,       -1.35865500064341E+04_wp,&
     &        1.16553933368645E+04_wp,       -5.30564697861340E+03_wp,&
     &        1.20090291321635E+03_wp,       -1.08090919788395E+02_wp,&
     &        1.72772750258446E+00_wp,        2.02042913309661E+04_wp,&
     &       -9.69805983886375E+04_wp,        1.92547001232532E+05_wp,&
     &       -2.03400177280416E+05_wp,        1.22200464983017E+05_wp,&
     &       -4.11926549688976E+04_wp,        7.10951430248936E+03_wp,&
     &       -4.93915304773088E+02_wp,        6.07404200127348E+00_wp,&
     &       -2.42919187900551E+05_wp,        1.31176361466298E+06_wp,&
     &       -2.99801591853811E+06_wp,        3.76327129765640E+06_wp/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56),&
     &     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64),&
     &     C(65)/&
     &       -2.81356322658653E+06_wp,        1.26836527332162E+06_wp,&
     &       -3.31645172484564E+05_wp,        4.52187689813627E+04_wp,&
     &       -2.49983048181121E+03_wp,        2.43805296995561E+01_wp,&
     &        3.28446985307204E+06_wp,       -1.97068191184322E+07_wp,&
     &        5.09526024926646E+07_wp,       -7.41051482115327E+07_wp,&
     &        6.63445122747290E+07_wp,       -3.75671766607634E+07_wp,&
     &        1.32887671664218E+07_wp,       -2.78561812808645E+06_wp,&
     &        3.08186404612662E+05_wp,       -1.38860897537170E+04_wp,&
     &        1.10017140269247E+02_wp/
      DATA ALFA1(1,1), ALFA1(2,1), ALFA1(3,1), ALFA1(4,1), ALFA1(5,1),&
     &     ALFA1(6,1), ALFA1(7,1), ALFA1(8,1), ALFA1(9,1), ALFA1(10,1),&
     &     ALFA1(11,1),ALFA1(12,1),ALFA1(13,1),ALFA1(14,1),ALFA1(15,1),&
     &     ALFA1(16,1),ALFA1(17,1),ALFA1(18,1),ALFA1(19,1),ALFA1(20,1),&
     &     ALFA1(21,1),ALFA1(22,1),ALFA1(23,1),ALFA1(24,1),ALFA1(25,1),&
     &     ALFA1(26,1)     /-4.44444444444444E-03_wp,-9.22077922077922E-04_wp,&
     &-8.84892884892885E-05_wp, 1.65927687832450E-04_wp, 2.46691372741793E-04_wp,&
     & 2.65995589346255E-04_wp, 2.61824297061501E-04_wp, 2.48730437344656E-04_wp,&
     & 2.32721040083232E-04_wp, 2.16362485712365E-04_wp, 2.00738858762752E-04_wp,&
     & 1.86267636637545E-04_wp, 1.73060775917876E-04_wp, 1.61091705929016E-04_wp,&
     & 1.50274774160908E-04_wp, 1.40503497391270E-04_wp, 1.31668816545923E-04_wp,&
     & 1.23667445598253E-04_wp, 1.16405271474738E-04_wp, 1.09798298372713E-04_wp,&
     & 1.03772410422993E-04_wp, 9.82626078369363E-05_wp, 9.32120517249503E-05_wp,&
     & 8.85710852478712E-05_wp, 8.42963105715700E-05_wp, 8.03497548407791E-05_wp/
      DATA ALFA1(1,2), ALFA1(2,2), ALFA1(3,2), ALFA1(4,2), ALFA1(5,2),&
     &     ALFA1(6,2), ALFA1(7,2), ALFA1(8,2), ALFA1(9,2), ALFA1(10,2),&
     &     ALFA1(11,2),ALFA1(12,2),ALFA1(13,2),ALFA1(14,2),ALFA1(15,2),&
     &     ALFA1(16,2),ALFA1(17,2),ALFA1(18,2),ALFA1(19,2),ALFA1(20,2),&
     &     ALFA1(21,2),ALFA1(22,2),ALFA1(23,2),ALFA1(24,2),ALFA1(25,2),&
     &     ALFA1(26,2)     / 6.93735541354589E-04_wp, 2.32241745182922E-04_wp,&
     &-1.41986273556691E-05_wp,-1.16444931672049E-04_wp,-1.50803558053049E-04_wp,&
     &-1.55121924918096E-04_wp,-1.46809756646466E-04_wp,-1.33815503867491E-04_wp,&
     &-1.19744975684254E-04_wp,-1.06184319207974E-04_wp,-9.37699549891194E-05_wp,&
     &-8.26923045588193E-05_wp,-7.29374348155221E-05_wp,-6.44042357721016E-05_wp,&
     &-5.69611566009369E-05_wp,-5.04731044303562E-05_wp,-4.48134868008883E-05_wp,&
     &-3.98688727717599E-05_wp,-3.55400532972042E-05_wp,-3.17414256609022E-05_wp,&
     &-2.83996793904175E-05_wp,-2.54522720634871E-05_wp,-2.28459297164725E-05_wp,&
     &-2.05352753106481E-05_wp,-1.84816217627666E-05_wp,-1.66519330021394E-05_wp/
      DATA ALFA2(1,1), ALFA2(2,1), ALFA2(3,1), ALFA2(4,1), ALFA2(5,1),&
     &     ALFA2(6,1), ALFA2(7,1), ALFA2(8,1), ALFA2(9,1), ALFA2(10,1),&
     &     ALFA2(11,1),ALFA2(12,1),ALFA2(13,1),ALFA2(14,1),ALFA2(15,1),&
     &     ALFA2(16,1),ALFA2(17,1),ALFA2(18,1),ALFA2(19,1),ALFA2(20,1),&
     &     ALFA2(21,1),ALFA2(22,1),ALFA2(23,1),ALFA2(24,1),ALFA2(25,1),&
     &     ALFA2(26,1)     /-3.54211971457744E-04_wp,-1.56161263945159E-04_wp,&
     & 3.04465503594936E-05_wp, 1.30198655773243E-04_wp, 1.67471106699712E-04_wp,&
     & 1.70222587683593E-04_wp, 1.56501427608595E-04_wp, 1.36339170977445E-04_wp,&
     & 1.14886692029825E-04_wp, 9.45869093034688E-05_wp, 7.64498419250898E-05_wp,&
     & 6.07570334965197E-05_wp, 4.74394299290509E-05_wp, 3.62757512005344E-05_wp,&
     & 2.69939714979225E-05_wp, 1.93210938247939E-05_wp, 1.30056674793963E-05_wp,&
     & 7.82620866744497E-06_wp, 3.59257485819352E-06_wp, 1.44040049814252E-07_wp,&
     &-2.65396769697939E-06_wp,-4.91346867098486E-06_wp,-6.72739296091248E-06_wp,&
     &-8.17269379678658E-06_wp,-9.31304715093561E-06_wp,-1.02011418798016E-05_wp/
      DATA ALFA2(1,2), ALFA2(2,2), ALFA2(3,2), ALFA2(4,2), ALFA2(5,2),&
     &     ALFA2(6,2), ALFA2(7,2), ALFA2(8,2), ALFA2(9,2), ALFA2(10,2),&
     &     ALFA2(11,2),ALFA2(12,2),ALFA2(13,2),ALFA2(14,2),ALFA2(15,2),&
     &     ALFA2(16,2),ALFA2(17,2),ALFA2(18,2),ALFA2(19,2),ALFA2(20,2),&
     &     ALFA2(21,2),ALFA2(22,2),ALFA2(23,2),ALFA2(24,2),ALFA2(25,2),&
     &     ALFA2(26,2)     / 3.78194199201773E-04_wp, 2.02471952761816E-04_wp,&
     &-6.37938506318862E-05_wp,-2.38598230603006E-04_wp,-3.10916256027362E-04_wp,&
     &-3.13680115247576E-04_wp,-2.78950273791323E-04_wp,-2.28564082619141E-04_wp,&
     &-1.75245280340847E-04_wp,-1.25544063060690E-04_wp,-8.22982872820208E-05_wp,&
     &-4.62860730588116E-05_wp,-1.72334302366962E-05_wp, 5.60690482304602E-06_wp,&
     & 2.31395443148287E-05_wp, 3.62642745856794E-05_wp, 4.58006124490189E-05_wp,&
     & 5.24595294959114E-05_wp, 5.68396208545815E-05_wp, 5.94349820393104E-05_wp,&
     & 6.06478527578422E-05_wp, 6.08023907788436E-05_wp, 6.01577894539460E-05_wp,&
     & 5.89199657344698E-05_wp, 5.72515823777593E-05_wp, 5.52804375585853E-05_wp/
      DATA BETA1(1,1), BETA1(2,1), BETA1(3,1), BETA1(4,1), BETA1(5,1),&
     &     BETA1(6,1), BETA1(7,1), BETA1(8,1), BETA1(9,1), BETA1(10,1),&
     &     BETA1(11,1),BETA1(12,1),BETA1(13,1),BETA1(14,1),BETA1(15,1),&
     &     BETA1(16,1),BETA1(17,1),BETA1(18,1),BETA1(19,1),BETA1(20,1),&
     &     BETA1(21,1),BETA1(22,1),BETA1(23,1),BETA1(24,1),BETA1(25,1),&
     &     BETA1(26,1)     / 1.79988721413553E-02_wp, 5.59964911064388E-03_wp,&
     & 2.88501402231133E-03_wp, 1.80096606761054E-03_wp, 1.24753110589199E-03_wp,&
     & 9.22878876572938E-04_wp, 7.14430421727287E-04_wp, 5.71787281789705E-04_wp,&
     & 4.69431007606482E-04_wp, 3.93232835462917E-04_wp, 3.34818889318298E-04_wp,&
     & 2.88952148495752E-04_wp, 2.52211615549573E-04_wp, 2.22280580798883E-04_wp,&
     & 1.97541838033063E-04_wp, 1.76836855019718E-04_wp, 1.59316899661821E-04_wp,&
     & 1.44347930197334E-04_wp, 1.31448068119965E-04_wp, 1.20245444949303E-04_wp,&
     & 1.10449144504599E-04_wp, 1.01828770740567E-04_wp, 9.41998224204238E-05_wp,&
     & 8.74130545753834E-05_wp, 8.13466262162801E-05_wp, 7.59002269646219E-05_wp/
      DATA BETA1(1,2), BETA1(2,2), BETA1(3,2), BETA1(4,2), BETA1(5,2),&
     &     BETA1(6,2), BETA1(7,2), BETA1(8,2), BETA1(9,2), BETA1(10,2),&
     &     BETA1(11,2),BETA1(12,2),BETA1(13,2),BETA1(14,2),BETA1(15,2),&
     &     BETA1(16,2),BETA1(17,2),BETA1(18,2),BETA1(19,2),BETA1(20,2),&
     &     BETA1(21,2),BETA1(22,2),BETA1(23,2),BETA1(24,2),BETA1(25,2),&
     &     BETA1(26,2)     /-1.49282953213429E-03_wp,-8.78204709546389E-04_wp,&
     &-5.02916549572035E-04_wp,-2.94822138512746E-04_wp,-1.75463996970783E-04_wp,&
     &-1.04008550460816E-04_wp,-5.96141953046458E-05_wp,-3.12038929076098E-05_wp,&
     &-1.26089735980230E-05_wp,-2.42892608575730E-07_wp, 8.05996165414274E-06_wp,&
     & 1.36507009262147E-05_wp, 1.73964125472926E-05_wp, 1.98672978842134E-05_wp,&
     & 2.14463263790823E-05_wp, 2.23954659232457E-05_wp, 2.28967783814713E-05_wp,&
     & 2.30785389811178E-05_wp, 2.30321976080909E-05_wp, 2.28236073720349E-05_wp,&
     & 2.25005881105292E-05_wp, 2.20981015361991E-05_wp, 2.16418427448104E-05_wp,&
     & 2.11507649256221E-05_wp, 2.06388749782171E-05_wp, 2.01165241997082E-05_wp/
      DATA BETA2(1,1), BETA2(2,1), BETA2(3,1), BETA2(4,1), BETA2(5,1),&
     &     BETA2(6,1), BETA2(7,1), BETA2(8,1), BETA2(9,1), BETA2(10,1),&
     &     BETA2(11,1),BETA2(12,1),BETA2(13,1),BETA2(14,1),BETA2(15,1),&
     &     BETA2(16,1),BETA2(17,1),BETA2(18,1),BETA2(19,1),BETA2(20,1),&
     &     BETA2(21,1),BETA2(22,1),BETA2(23,1),BETA2(24,1),BETA2(25,1),&
     &     BETA2(26,1)     / 5.52213076721293E-04_wp, 4.47932581552385E-04_wp,&
     & 2.79520653992021E-04_wp, 1.52468156198447E-04_wp, 6.93271105657044E-05_wp,&
     & 1.76258683069991E-05_wp,-1.35744996343269E-05_wp,-3.17972413350427E-05_wp,&
     &-4.18861861696693E-05_wp,-4.69004889379141E-05_wp,-4.87665447413787E-05_wp,&
     &-4.87010031186735E-05_wp,-4.74755620890087E-05_wp,-4.55813058138628E-05_wp,&
     &-4.33309644511266E-05_wp,-4.09230193157750E-05_wp,-3.84822638603221E-05_wp,&
     &-3.60857167535411E-05_wp,-3.37793306123367E-05_wp,-3.15888560772110E-05_wp,&
     &-2.95269561750807E-05_wp,-2.75978914828336E-05_wp,-2.58006174666884E-05_wp,&
     &-2.41308356761280E-05_wp,-2.25823509518346E-05_wp,-2.11479656768913E-05_wp/
      DATA BETA2(1,2), BETA2(2,2), BETA2(3,2), BETA2(4,2), BETA2(5,2),&
     &     BETA2(6,2), BETA2(7,2), BETA2(8,2), BETA2(9,2), BETA2(10,2),&
     &     BETA2(11,2),BETA2(12,2),BETA2(13,2),BETA2(14,2),BETA2(15,2),&
     &     BETA2(16,2),BETA2(17,2),BETA2(18,2),BETA2(19,2),BETA2(20,2),&
     &     BETA2(21,2),BETA2(22,2),BETA2(23,2),BETA2(24,2),BETA2(25,2),&
     &     BETA2(26,2)     /-4.74617796559960E-04_wp,-4.77864567147321E-04_wp,&
     &-3.20390228067038E-04_wp,-1.61105016119962E-04_wp,-4.25778101285435E-05_wp,&
     & 3.44571294294968E-05_wp, 7.97092684075675E-05_wp, 1.03138236708272E-04_wp,&
     & 1.12466775262204E-04_wp, 1.13103642108481E-04_wp, 1.08651634848774E-04_wp,&
     & 1.01437951597662E-04_wp, 9.29298396593364E-05_wp, 8.40293133016090E-05_wp,&
     & 7.52727991349134E-05_wp, 6.69632521975731E-05_wp, 5.92564547323195E-05_wp,&
     & 5.22169308826976E-05_wp, 4.58539485165361E-05_wp, 4.01445513891487E-05_wp,&
     & 3.50481730031328E-05_wp, 3.05157995034347E-05_wp, 2.64956119950516E-05_wp,&
     & 2.29363633690998E-05_wp, 1.97893056664022E-05_wp, 1.70091984636413E-05_wp/
      DATA BETA3(1,1), BETA3(2,1), BETA3(3,1), BETA3(4,1), BETA3(5,1),&
     &     BETA3(6,1), BETA3(7,1), BETA3(8,1), BETA3(9,1), BETA3(10,1),&
     &     BETA3(11,1),BETA3(12,1),BETA3(13,1),BETA3(14,1),BETA3(15,1),&
     &     BETA3(16,1),BETA3(17,1),BETA3(18,1),BETA3(19,1),BETA3(20,1),&
     &     BETA3(21,1),BETA3(22,1),BETA3(23,1),BETA3(24,1),BETA3(25,1),&
     &     BETA3(26,1)     / 7.36465810572578E-04_wp, 8.72790805146194E-04_wp,&
     & 6.22614862573135E-04_wp, 2.85998154194304E-04_wp, 3.84737672879366E-06_wp,&
     &-1.87906003636972E-04_wp,-2.97603646594555E-04_wp,-3.45998126832656E-04_wp,&
     &-3.53382470916038E-04_wp,-3.35715635775049E-04_wp,-3.04321124789040E-04_wp,&
     &-2.66722723047613E-04_wp,-2.27654214122820E-04_wp,-1.89922611854562E-04_wp,&
     &-1.55058918599094E-04_wp,-1.23778240761874E-04_wp,-9.62926147717644E-05_wp,&
     &-7.25178327714425E-05_wp,-5.22070028895634E-05_wp,-3.50347750511901E-05_wp,&
     &-2.06489761035552E-05_wp,-8.70106096849767E-06_wp, 1.13698686675100E-06_wp,&
     & 9.16426474122779E-06_wp, 1.56477785428873E-05_wp, 2.08223629482467E-05_wp/
      DATA GAMA(1),   GAMA(2),   GAMA(3),   GAMA(4),   GAMA(5),&
     &     GAMA(6),   GAMA(7),   GAMA(8),   GAMA(9),   GAMA(10),&
     &     GAMA(11),  GAMA(12),  GAMA(13),  GAMA(14),  GAMA(15),&
     &     GAMA(16),  GAMA(17),  GAMA(18),  GAMA(19),  GAMA(20),&
     &     GAMA(21),  GAMA(22),  GAMA(23),  GAMA(24),  GAMA(25),&
     &     GAMA(26)        / 6.29960524947437E-01_wp, 2.51984209978975E-01_wp,&
     & 1.54790300415656E-01_wp, 1.10713062416159E-01_wp, 8.57309395527395E-02_wp,&
     & 6.97161316958684E-02_wp, 5.86085671893714E-02_wp, 5.04698873536311E-02_wp,&
     & 4.42600580689155E-02_wp, 3.93720661543510E-02_wp, 3.54283195924455E-02_wp,&
     & 3.21818857502098E-02_wp, 2.94646240791158E-02_wp, 2.71581677112934E-02_wp,&
     & 2.51768272973862E-02_wp, 2.34570755306079E-02_wp, 2.19508390134907E-02_wp,&
     & 2.06210828235646E-02_wp, 1.94388240897881E-02_wp, 1.83810633800683E-02_wp,&
     & 1.74293213231963E-02_wp, 1.65685837786612E-02_wp, 1.57865285987918E-02_wp,&
     & 1.50729501494096E-02_wp, 1.44193250839955E-02_wp, 1.38184805735342E-02_wp/
!***FIRST EXECUTABLE STATEMENT  wp_asyjy
      TA = F1MACH(3,wp_dummy)
      TOL = MAX(TA,1.0E-15_wp)
      TB = F1MACH(5,wp_dummy)
      JU = I1MACH(15)
      IF(FLGJY.EQ.1.00_wp) GO TO 6
      JR = I1MACH(14)
      ELIM = -2.3030_wp*TB*(JU+JR)
      GO TO 7
    6 CONTINUE
      ELIM = -2.3030_wp*(TB*JU+3.00_wp)
    7 CONTINUE
      FN = FNU
      IFLW = 0
      DO 170 JN=1,IN
        XX = X/FN
        WK(1) = 1.00_wp - XX*XX
        ABW2 = ABS(WK(1))
        WK(2) = SQRT(ABW2)
        WK(7) = FN**CON2
        IF (ABW2.GT.0.277500_wp) GO TO 80
!
!     ASYMPTOTIC EXPANSION
!     CASES NEAR X=FN, ABS(1.-(X/FN)**2).LE.0.2775
!     COEFFICIENTS OF ASYMPTOTIC EXPANSION BY SERIES
!
!     ZETA AND TRUNCATION FOR A(ZETA) AND B(ZETA) SERIES
!
!     KMAX IS TRUNCATION INDEX FOR A(ZETA) AND B(ZETA) SERIES=MAX(2,SA)
!
        SA = 0.00_wp
        IF (ABW2.EQ.0.00_wp) GO TO 10
        SA = TOLS/LOG(ABW2)
   10   SB = SA
        DO 20 I=1,5
          AKM = MAX(SA,2.00_wp)
          KMAX(I) = INT(AKM)
          SA = SA + SB
   20   CONTINUE
        KB = KMAX(5)
        KLAST = KB - 1
        SA = GAMA(KB)
        DO 30 K=1,KLAST
          KB = KB - 1
          SA = SA*WK(1) + GAMA(KB)
   30   CONTINUE
        Z = WK(1)*SA
        AZ = ABS(Z)
        RTZ = SQRT(AZ)
        WK(3) = CON1*AZ*RTZ
        WK(4) = WK(3)*FN
        WK(5) = RTZ*WK(7)
        WK(6) = -WK(5)*WK(5)
        IF(Z.LE.0.00_wp) GO TO 35
        IF(WK(4).GT.ELIM) GO TO 75
        WK(6) = -WK(6)
   35   CONTINUE
        PHI = SQRT(SQRT(SA+SA+SA+SA))
!
!     B(ZETA) FOR S=0
!
        KB = KMAX(5)
        KLAST = KB - 1
        SB = BETA(KB,1)
        DO 40 K=1,KLAST
          KB = KB - 1
          SB = SB*WK(1) + BETA(KB,1)
   40   CONTINUE
        KSP1 = 1
        FN2 = FN*FN
        RFN2 = 1.00_wp/FN2
        RDEN = 1.00_wp
        ASUM = 1.00_wp
        RELB = TOL*ABS(SB)
        BSUM = SB
        DO 60 KS=1,4
          KSP1 = KSP1 + 1
          RDEN = RDEN*RFN2
!
!     A(ZETA) AND B(ZETA) FOR S=1,2,3,4
!
          KSTEMP = 5 - KS
          KB = KMAX(KSTEMP)
          KLAST = KB - 1
          SA = ALFA(KB,KS)
          SB = BETA(KB,KSP1)
          DO 50 K=1,KLAST
            KB = KB - 1
            SA = SA*WK(1) + ALFA(KB,KS)
            SB = SB*WK(1) + BETA(KB,KSP1)
   50     CONTINUE
          TA = SA*RDEN
          TB = SB*RDEN
          ASUM = ASUM + TA
          BSUM = BSUM + TB
          IF (ABS(TA).LE.TOL .AND. ABS(TB).LE.RELB) GO TO 70
   60   CONTINUE
   70   CONTINUE
        BSUM = BSUM/(FN*WK(7))
        GO TO 160
!
   75   CONTINUE
        IFLW = 1
        RETURN
!
   80   CONTINUE
        UPOL(1) = 1.00_wp
        TAU = 1.00_wp/WK(2)
        T2 = 1.00_wp/WK(1)
        IF (WK(1).GE.0.00_wp) GO TO 90
!
!     CASES FOR (X/FN).GT.SQRT(1.2775)
!
        WK(3) = ABS(WK(2)-ATAN(WK(2)))
        WK(4) = WK(3)*FN
        RCZ = -CON1/WK(4)
        Z32 = 1.50_wp*WK(3)
        RTZ = Z32**CON2
        WK(5) = RTZ*WK(7)
        WK(6) = -WK(5)*WK(5)
        GO TO 100
   90   CONTINUE
!
!     CASES FOR (X/FN).LT.SQRT(0.7225)
!
        WK(3) = ABS(LOG((1.00_wp+WK(2))/XX)-WK(2))
        WK(4) = WK(3)*FN
        RCZ = CON1/WK(4)
        IF(WK(4).GT.ELIM) GO TO 75
        Z32 = 1.50_wp*WK(3)
        RTZ = Z32**CON2
        WK(7) = FN**CON2
        WK(5) = RTZ*WK(7)
        WK(6) = WK(5)*WK(5)
  100   CONTINUE
        PHI = SQRT((RTZ+RTZ)*TAU)
        TB = 1.00_wp
        ASUM = 1.00_wp
        TFN = TAU/FN
        RDEN=1.00_wp/FN
        RFN2=RDEN*RDEN
        RDEN=1.00_wp
        UPOL(2) = (C(1)*T2+C(2))*TFN
        CRZ32 = CON548*RCZ
        BSUM = UPOL(2) + CRZ32
        RELB = TOL*ABS(BSUM)
        AP = TFN
        KS = 0
        KP1 = 2
        RZDEN = RCZ
        L = 2
        ISETA=0
        ISETB=0
        DO 140 LR=2,8,2
!
!     COMPUTE TWO U POLYNOMIALS FOR NEXT A(ZETA) AND B(ZETA)
!
          LRP1 = LR + 1
          DO 120 K=LR,LRP1
            KS = KS + 1
            KP1 = KP1 + 1
            L = L + 1
            S1 = C(L)
            DO 110 J=2,KP1
              L = L + 1
              S1 = S1*T2 + C(L)
  110       CONTINUE
            AP = AP*TFN
            UPOL(KP1) = AP*S1
            CR(KS) = BR(KS)*RZDEN
            RZDEN = RZDEN*RCZ
            DR(KS) = AR(KS)*RZDEN
  120     CONTINUE
          SUMA = UPOL(LRP1)
          SUMB = UPOL(LR+2) + UPOL(LRP1)*CRZ32
          JU = LRP1
          DO 130 JR=1,LR
            JU = JU - 1
            SUMA = SUMA + CR(JR)*UPOL(JU)
            SUMB = SUMB + DR(JR)*UPOL(JU)
  130     CONTINUE
          RDEN=RDEN*RFN2
          TB = -TB
          IF (WK(1).GT.0.00_wp) TB = ABS(TB)
          IF(RDEN.LT.TOL) GO TO 131
          ASUM = ASUM + SUMA*TB
          BSUM = BSUM + SUMB*TB
          GO TO 140
  131     IF(ISETA.EQ.1) GO TO 132
          IF(ABS(SUMA).LT.TOL) ISETA=1
          ASUM=ASUM+SUMA*TB
  132     IF(ISETB.EQ.1) GO TO 133
          IF(ABS(SUMB).LT.RELB) ISETB=1
          BSUM=BSUM+SUMB*TB
  133     IF(ISETA.EQ.1 .AND. ISETB.EQ.1) GO TO 150
  140   CONTINUE
  150   TB = WK(5)
        IF (WK(1).GT.0.00_wp) TB = -TB
        BSUM = BSUM/TB
!
  160   CONTINUE
        CALL FUNJY(WK(6), WK(5), WK(4), FI, DFI)
        TA=1.00_wp/TOL
        TB=F1MACH(1,wp_dummy)*TA*1.0E+3_wp
        IF(ABS(FI).GT.TB) GO TO 165
        FI=FI*TA
        DFI=DFI*TA
        PHI=PHI*TOL
  165   CONTINUE
        Y(JN) = FLGJY*PHI*(FI*ASUM+DFI*BSUM)/WK(7)
        FN = FN - FLGJY
  170 CONTINUE
      RETURN
      END SUBROUTINE
!
      !> Quad precision version of wp_asyjy.
      SUBROUTINE ep_asyjy (FUNJY, X, FNU, FLGJY, IN, Y, WK, IFLW)
      IMPLICIT NONE
!
      EXTERNAL FUNJY
      INTEGER IFLW, IN
      REAL(kind=ep1) X, FNU, FLGJY, Y, WK
      DIMENSION Y(*), WK(*)

         CALL XERMSG('SLATEC','ep_asyjy','Quad precision version not implemented yet.',1,1)

      END SUBROUTINE
!> \verbatim
!>***PURPOSE  Compute an N member sequence of I Bessel functions
!>            I/SUB(ALPHA+K-1)/(X), K=1,...,N or scaled Bessel functions
!>            EXP(-X)*I/SUB(ALPHA+K-1)/(X), K=1,...,N for nonnegative
!>            ALPHA and X.
!>***LIBRARY   SLATEC
!>***CATEGORY  C10B3
!>***TYPE      REAL(kind=wp) (BESI-S, wp_besi-D)
!>***KEYWORDS  I BESSEL FUNCTION, SPECIAL FUNCTIONS
!>***AUTHOR  Amos, D. E., (SNLA)
!>           Daniel, S. L., (SNLA)
!>***DESCRIPTION
!>
!>     Abstract  **** a REAL(kind=wp) routine ****
!>         wp_besi computes an N member sequence of I Bessel functions
!>         I/sub(ALPHA+K-1)/(X), K=1,...,N or scaled Bessel functions
!>         EXP(-X)*I/sub(ALPHA+K-1)/(X), K=1,...,N for nonnegative ALPHA
!>         and X.  A combination of the power series, the asymptotic
!>         expansion for X to infinity, and the uniform asymptotic
!>         expansion for NU to infinity are applied over subdivisions of
!>         the (NU,X) plane.  For values not covered by one of these
!>         formulae, the order is incremented by an integer so that one
!>         of these formulae apply.  Backward recursion is used to reduce
!>         orders by integer values.  The asymptotic expansion for X to
!>         infinity is used only when the entire sequence (specifically
!>         the last member) lies within the region covered by the
!>         expansion.  Leading terms of these expansions are used to test
!>         for over or underflow where appropriate.  If a sequence is
!>         requested and the last member would underflow, the result is
!>         set to zero and the next lower order tried, etc., until a
!>         member comes on scale or all are set to zero.  An overflow
!>         cannot occur with scaling.
!>
!>         The maximum number of significant digits obtainable
!>         is the smaller of 14 and the number of digits carried in
!>         REAL(kind=wp) arithmetic.
!>
!>     Description of Arguments
!>
!>         Input      X,ALPHA are REAL(kind=wp)
!>           X      - X .GE. 0.00_wp
!>           ALPHA  - order of first member of the sequence,
!>                    ALPHA .GE. 0.00_wp
!>           KODE   - a parameter to indicate the scaling option
!>                    KODE=1 returns
!>                           Y(K)=        I/sub(ALPHA+K-1)/(X),
!>                                K=1,...,N
!>                    KODE=2 returns
!>                           Y(K)=EXP(-X)*I/sub(ALPHA+K-1)/(X),
!>                                K=1,...,N
!>           N      - number of members in the sequence, N .GE. 1
!>
!>         Output     Y is REAL(kind=wp)
!>           Y      - a vector whose first N components contain
!>                    values for I/sub(ALPHA+K-1)/(X) or scaled
!>                    values for EXP(-X)*I/sub(ALPHA+K-1)/(X),
!>                    K=1,...,N depending on KODE
!>           NZ     - number of components of Y set to zero due to
!>                    underflow,
!>                    NZ=0   , normal return, computation completed
!>                    NZ .NE. 0, last NZ components of Y set to zero,
!>                             Y(K)=0.00_wp, K=N-NZ+1,...,N.
!>
!>     Error Conditions
!>         Improper input arguments - a fatal error
!>         Overflow with KODE=1 - a fatal error
!>         Underflow - a non-fatal error(NZ .NE. 0)
!>
!>***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
!>                 subroutines IBESS and JBESS for Bessel functions
!>                 I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0, ACM
!>                 Transactions on Mathematical Software 3, (1977),
!>                 pp. 76-92.
!>               F. W. J. Olver, Tables of Bessel Functions of Moderate
!>                 or Large Orders, NPL Mathematical Tables 6, Her
!>                 Majesty's Stationery Office, London, 1962.
!>***ROUTINES CALLED  F1MACH, cfp_asyik, cfp_lngam, I1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   750101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   890911  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900326  Removed duplicate information from DESCRIPTION section.
!>           (WRB)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!> \endverbatim
!
      SUBROUTINE wp_besi (X, ALPHA, KODE, N, Y, NZ)
      IMPLICIT NONE
!
!
      INTEGER I, IALP, IN, INLIM, IS, I1, K, KK, KM, KODE, KT, N, NN, NS, NZ
      REAL(kind=wp) AIN,AK,AKM,ALPHA,ANS,AP,ARG,ATOL,TOLLN,DFN,&
     & DTM, DX, EARG, ELIM, ETX, FLGIK,FN, FNF, FNI,FNP1,FNU,GLN,RA,&
     & RTTPI, S, SX, SXO2, S1, S2, T, TA, TB, TEMP, TFN, TM, TOL,&
     & TRX, T2, X, XO2, XO2L, Y, Z
      DIMENSION Y(*), TEMP(3)
      SAVE RTTPI, INLIM
      DATA RTTPI           / 3.98942280401433E-01_wp/ != sqrt(1/(2*pi))
      DATA INLIM           /          80         / !this may need to be modified to get quad precision efficient routine
!***FIRST EXECUTABLE STATEMENT  wp_besi
      NZ = 0
      KT = 1
!     I1MACH(15) REPLACES I1MACH(12) IN A REAL(kind=wp) CODE
!     I1MACH(14) REPLACES I1MACH(11) IN A REAL(kind=wp) CODE
      RA = F1MACH(3,wp_dummy)
      TOL = MAX(RA,1.0E-15_wp) !1.0E-15_wp should be replaced by 1.0E-33_wp in quad prec code
      I1 = -I1MACH(15)
      GLN = F1MACH(5,wp_dummy)
      ELIM = 2.3030_wp*(I1*GLN-3.00_wp) !2.3030 is probably universal - a factor in under/overflow tests
!     TOLLN = -LN(TOL)
      I1 = I1MACH(14)+1
      TOLLN = 2.3030_wp*GLN*I1
      TOLLN = MIN(TOLLN,34.53880_wp) !34.53880 is -Log(epsilon); this is -75.9853080688 in quad precision
      IF (N < 1) THEN
         GO TO 590
      ELSE IF (N == 1) THEN
         GO TO 10
      ELSE
         GO TO 20
      END IF
   10 KT = 2
   20 NN = N
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 570
      IF (X < 0) THEN
         GO TO 600
      ELSE IF (X == 0) THEN
         GO TO 30
      ELSE
         GO TO 80
      END IF
   30 IF (ALPHA < 0) THEN
         GO TO 580
      ELSE IF (ALPHA == 0) THEN
         GO TO 40
      ELSE
         GO TO 50
      END IF
   40 Y(1) = 1.00_wp
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.00_wp
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.00_wp) GO TO 580
!
      IALP = INT(ALPHA)
      FNI = IALP + N - 1
      FNF = ALPHA - IALP
      DFN = FNI + FNF
      FNU = DFN
      IN = 0
      XO2 = X*0.50_wp
      SXO2 = XO2*XO2
      ETX = KODE - 1
      SX = ETX*X
!
!     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
!     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
!     APPLIED.
!
      IF (SXO2.LE.(FNU+1.00_wp)) GO TO 90
      IF (X.LE.12.00_wp) GO TO 110
      FN = 0.550_wp*FNU*FNU !0.55 is a universal factor coming from the requirement 8*x > 1.1* 4*FN^2 => x > 0.55*FN^2
      FN = MAX(17.00_wp,FN)
      IF (X.GE.FN) GO TO 430 !asymptotic expansion: this rule is the same regardless of precision used.
      ANS = MAX(36.00_wp-FNU,0.00_wp)
      NS = INT(ANS)
      FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      KM = N - 1 + NS
      IF (KM.GT.0) IS = 3
      GO TO 120 !uniform asymptotic expansion
   90 FN = FNU
      FNP1 = FN + 1.00_wp
      XO2L = LOG(XO2)
      IS = KT
      IF (X.LE.0.50_wp) GO TO 230 !power series
      NS = 0
  100 FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.00_wp
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 230
  110 XO2L = LOG(XO2)
      NS = INT(SXO2-FNU)
      GO TO 100
  120 CONTINUE
!
!     OVERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
!
      IF (KODE.EQ.2) GO TO 130
      IF (ALPHA.LT.1.00_wp) GO TO 150
      Z = X/ALPHA
      RA = SQRT(1.00_wp+Z*Z)
      GLN = LOG((1.00_wp+RA)/Z)
      T = RA*(1.00_wp-ETX) + ETX/(Z+RA)
      ARG = ALPHA*(T-GLN)
      IF (ARG.GT.ELIM) GO TO 610
      IF (KM.EQ.0) GO TO 140
  130 CONTINUE
!
!     UNDERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
!
      Z = X/FN
      RA = SQRT(1.00_wp+Z*Z)
      GLN = LOG((1.00_wp+RA)/Z)
      T = RA*(1.00_wp-ETX) + ETX/(Z+RA)
      ARG = FN*(T-GLN)
  140 IF (ARG.LT.(-ELIM)) GO TO 280
      GO TO 190
  150 IF (X.GT.ELIM) GO TO 610
      GO TO 130
!
!     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
!
  160 IF (KM.NE.0) GO TO 170
      Y(1) = TEMP(3)
      RETURN
  170 TEMP(1) = TEMP(3)
      IN = NS
      KT = 1
      I1 = 0
  180 CONTINUE
      IS = 2
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 350
      Z = X/FN
      RA = SQRT(1.00_wp+Z*Z)
      GLN = LOG((1.00_wp+RA)/Z)
      T = RA*(1.00_wp-ETX) + ETX/(Z+RA)
      ARG = FN*(T-GLN)
  190 CONTINUE
      I1 = ABS(3-IS)
      I1 = MAX(I1,1)
      FLGIK = 1.00_wp
      CALL wp_asyik(X,FN,KODE,FLGIK,RA,ARG,I1,TEMP(IS))
      GO TO (180, 350, 510), IS
!
!     SERIES FOR (X/2)**2.LE.NU+1
!
  230 CONTINUE
      GLN = cfp_lngam(FNP1)
      ARG = FN*XO2L - GLN - SX
      IF (ARG.LT.(-ELIM)) GO TO 300
      EARG = EXP(ARG)
  240 CONTINUE
      S = 1.00_wp
      IF (X.LT.TOL) GO TO 260
      AK = 3.00_wp
      T2 = 1.00_wp
      T = 1.00_wp
      S1 = FN
      DO 250 K=1,17 !17 should be replaced by a higher number ~29-30 for quad precision code
        S2 = T2 + S1
        T = T*SXO2/S2
        S = S + T
        IF (ABS(T).LT.TOL) GO TO 260
        T2 = T2 + AK
        AK = AK + 2.00_wp
        S1 = S1 + FN
  250 CONTINUE
  260 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (270, 350, 500), IS
  270 EARG = EARG*FN/XO2
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 240
!
!     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
!
  280 Y(NN) = 0.00_wp
      NN = NN - 1
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      IF (NN < 1) THEN
         GO TO 340
      ELSE IF (NN == 1) THEN
         GO TO 290
      ELSE
         GO TO 130
      END IF
  290 KT = 2
      IS = 2
      GO TO 130
  300 Y(NN) = 0.00_wp
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      IF (NN < 1) THEN
         GO TO 340
      ELSE IF (NN == 1) THEN
         GO TO 310
      ELSE
         GO TO 320
      END IF
  310 KT = 2
      IS = 2
  320 IF (SXO2.LE.FNP1) GO TO 330
      GO TO 130
  330 ARG = ARG - XO2L + LOG(FNP1)
      IF (ARG.LT.(-ELIM)) GO TO 300
      GO TO 230
  340 NZ = N - NN
      RETURN
!
!     BACKWARD RECURSION SECTION
!
  350 CONTINUE
      NZ = N - NN
  360 CONTINUE
      IF(KT.EQ.2) GO TO 420
      S1 = TEMP(1)
      S2 = TEMP(2)
      TRX = 2.00_wp/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      IF (IN.EQ.0) GO TO 390
!     BACKWARD RECUR TO INDEX ALPHA+NN-1
      DO 380 I=1,IN
        S = S2
        S2 = TM*S2 + S1
        S1 = S
        DTM = DTM - 1.00_wp
        TM = (DTM+FNF)*TRX
  380 CONTINUE
      Y(NN) = S1
      IF (NN.EQ.1) RETURN
      Y(NN-1) = S2
      IF (NN.EQ.2) RETURN
      GO TO 400
  390 CONTINUE
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = S1
      Y(NN-1) = S2
      IF (NN.EQ.2) RETURN
  400 K = NN + 1
      DO 410 I=3,NN
        K = K - 1
        Y(K-2) = TM*Y(K-1) + Y(K)
        DTM = DTM - 1.00_wp
        TM = (DTM+FNF)*TRX
  410 CONTINUE
      RETURN
  420 Y(1) = TEMP(2)
      RETURN
!
!     ASYMPTOTIC EXPANSION FOR X TO INFINITY
!
  430 CONTINUE
      EARG = RTTPI/SQRT(X)
      IF (KODE.EQ.2) GO TO 440
      IF (X.GT.ELIM) GO TO 610
      EARG = EARG*EXP(X)
  440 ETX = 8.00_wp*X
      IS = KT
      IN = 0
      FN = FNU
  450 DX = FNI + FNI
      TM = 0.00_wp
      IF (FNI.EQ.0.00_wp .AND. ABS(FNF).LT.TOL) GO TO 460
      TM = 4.00_wp*FNF*(FNI+FNI+FNF)
  460 CONTINUE
      DTM = DX*DX
      S1 = ETX
      TRX = DTM - 1.00_wp
      DX = -(TRX+TM)/ETX
      T = DX
      S = 1.00_wp + DX
      ATOL = TOL*ABS(S)
      S2 = 1.00_wp
      AK = 8.00_wp
      DO 470 K=1,25
        S1 = S1 + ETX
        S2 = S2 + AK
        DX = DTM - S2
        AP = DX + TM
        T = -T*AP/S1
        S = S + T
        IF (ABS(T).LE.ATOL) GO TO 480
        AK = AK + 8.00_wp
  470 CONTINUE
  480 TEMP(IS) = S*EARG
      IF(IS.EQ.2) GO TO 360
      IS = 2
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      GO TO 450
!
!     BACKWARD RECURSION WITH NORMALIZATION BY
!     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
!
  500 CONTINUE
!     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = MAX(3.00_wp-FN,0.00_wp)
      KM = INT(AKM)
      TFN = FN + KM
      TA = (GLN+TFN-0.91893853320_wp-0.08333333330_wp/TFN)/(TFN+0.50_wp) !universal factors: 0.918938 = Log(sqrt(2*pi)); 0.0833333=1/12
      TA = XO2L - TA
      TB = -(1.00_wp-1.00_wp/TFN)/TFN
      AIN = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.50_wp
      IN = INT(AIN)
      IN = IN + KM
      GO TO 520
  510 CONTINUE
!     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      T = 1.00_wp/(FN*RA)
      AIN = TOLLN/(GLN+SQRT(GLN*GLN+T*TOLLN)) + 1.50_wp
      IN = INT(AIN)
      IF (IN.GT.INLIM) GO TO 160
  520 CONTINUE
      TRX = 2.00_wp/X
      DTM = FNI + IN
      TM = (DTM+FNF)*TRX
      TA = 0.00_wp
      TB = TOL
      KK = 1
  530 CONTINUE
!
!     BACKWARD RECUR UNINDEXED
!
      DO 540 I=1,IN
        S = TB
        TB = TM*TB + TA
        TA = S
        DTM = DTM - 1.00_wp
        TM = (DTM+FNF)*TRX
  540 CONTINUE
!     NORMALIZATION
      IF (KK.NE.1) GO TO 550
      TA = (TA/TB)*TEMP(3)
      TB = TEMP(3)
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 530
  550 Y(NN) = TB
      NZ = N - NN
      IF (NN.EQ.1) RETURN
      TB = TM*TB + TA
      K = NN - 1
      Y(K) = TB
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.00_wp
      TM = (DTM+FNF)*TRX
      KM = K - 1
!
!     BACKWARD RECUR INDEXED
!
      DO 560 I=1,KM
        Y(K-1) = TM*Y(K) + Y(K+1)
        DTM = DTM - 1.00_wp
        TM = (DTM+FNF)*TRX
        K = K - 1
  560 CONTINUE
      RETURN
!
!
!
  570 CONTINUE
      CALL XERMSG ('SLATEC', 'wp_besi', 'SCALING OPTION, KODE, NOT 1 OR 2.', 2, 1)
      RETURN
  580 CONTINUE
      CALL XERMSG ('SLATEC', 'wp_besi', 'ORDER, ALPHA, LESS THAN ZERO.', 2, 1)
      RETURN
  590 CONTINUE
      CALL XERMSG ('SLATEC', 'wp_besi', 'N LESS THAN ONE.', 2, 1)
      RETURN
  600 CONTINUE
      CALL XERMSG ('SLATEC', 'wp_besi', 'X LESS THAN ZERO.', 2, 1)
      RETURN
  610 CONTINUE
      CALL XERMSG ('SLATEC', 'wp_besi', 'OVERFLOW, X TOO LARGE FOR KODE = 1.', 6, 1)
      RETURN
      END SUBROUTINE
      !> Quad precision equivalent of wp_besi.
      !> Precision of the results for X=0,1,...,999,1000 and NU=0,1/2,...,99,99+1/2 was compared with Mathematica. The relative error is at most 2*EPS, where EPS = 10^-33.
      !> For very large X: X >~ 2000 the results have relative precision around EPS.
      SUBROUTINE ep_besi (X, ALPHA, KODE, N, Y, NZ)
      IMPLICIT NONE
!
!
      INTEGER, parameter :: numin = 68 !Order from which the uniform asymptotic expansion gives results in full quad precision accuracy
      REAL(kind=ep1), parameter :: numinp1 = numin + 1.0_ep1
      INTEGER, parameter :: kmaxseries = 30 !Maximum order required in the power series expansion to obtain full quad precision accuracy in the parameter space (x/2)^2 = nu + 1, nu >= 0.
      INTEGER, parameter :: kmaxae = 30 !Maximum order required in the asymptotic expansion to obtain full quad precision accuracy in the space x > max(76,0.55*nu^2).
      REAL(kind=ep1), parameter :: xmax = 2.0_ep1*sqrt(real(numin,ep1)) !~ 16.6
      REAL(kind=ep1), parameter :: xmin = 76.0_ep1 !minimum x value for which the asymptotic expansion can be used.
      INTEGER I, IALP, IN, INLIM, IS, I1, K, KK, KM, KODE, KT, N, NN, NS, NZ, tmpi
      REAL(kind=ep1) AIN,AK,AKM,ALPHA,ANS,AP,ARG,ATOL,TOLLN,DFN,&
     & DTM, DX, EARG, ELIM, ETX, FLGIK,FN, FNF, FNI,FNP1,FNU,GLN,RA,&
     & RTTPI, S, SX, SXO2, S1, S2, T, TA, TB, TEMP, TFN, TM, TOL,&
     & TRX, T2, X, XO2, XO2L, Y, Z
      DIMENSION Y(*), TEMP(3)
      SAVE RTTPI, INLIM
      DATA RTTPI           / 0.3989422804014326779399460599343819_ep1/ != sqrt(1/(2*pi))
      DATA INLIM           /         160         / !this has been doubled: see the line INT(2.3*AIN) which causes the recursion to be ~2 times longer compared to the double precision case.
!***FIRST EXECUTABLE STATEMENT  ep_besi
      NZ = 0
      KT = 1
      RA = F1MACH(3,ep_dummy)
      TOL = MAX(RA,1.0E-33_ep1)
      I1 = -I1MACH(18)
      GLN = F1MACH(5,ep_dummy)
      ELIM = 2.3030_ep1*(I1*GLN-3.00_ep1) !2.3030 is probably universal - a scaling factor in under/overflow tests
!     TOLLN = -LN(TOL)
      I1 = I1MACH(17)+1
      TOLLN = 2.3030_ep1*GLN*I1
      TOLLN = MIN(TOLLN,75.9853080688_ep1) != -Log(epsilon); for epsilon = 10E-33.
      IF (N < 1) THEN
         GO TO 590
      ELSE IF (N == 1) THEN
         GO TO 10
      ELSE
         GO TO 20
      END IF
   10 KT = 2
   20 NN = N
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 570
      IF (X < 0) THEN
         GO TO 600
      ELSE IF (X == 0) THEN
         GO TO 30
      ELSE
         GO TO 80
      END IF
   30 IF (ALPHA < 0) THEN
         GO TO 580
      ELSE IF (ALPHA == 0) THEN
         GO TO 40
      ELSE
         GO TO 50
      END IF
   40 Y(1) = 1.00_ep1
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.00_ep1
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.00_ep1) GO TO 580
!
      IALP = INT(ALPHA)
      FNI = IALP + N - 1
      FNF = ALPHA - IALP
      DFN = FNI + FNF
      FNU = DFN
      IN = 0
      XO2 = X*0.50_ep1
      SXO2 = XO2*XO2
      ETX = KODE - 1
      SX = ETX*X
!
!     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
!     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
!     APPLIED.
!
      IF (SXO2.LE.(FNU+1.00_ep1)) GO TO 90
      IF (X.LE.XMAX) GO TO 110
      FN = 0.550_ep1*FNU*FNU !0.55 is a universal factor coming from the requirement 8*x > 1.1* 4*FN^2 => x > 0.55*FN^2
      FN = MAX(XMIN,FN)
      IF (X.GE.FN) GO TO 430 !asymptotic expansion: this rule is the same regardless of precision used.
      ANS = MAX(NUMINP1-FNU,0.00_ep1)
      NS = INT(ANS)
      FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      KM = N - 1 + NS
      IF (KM.GT.0) IS = 3
      GO TO 120 !uniform asymptotic expansion
   90 FN = FNU
      FNP1 = FN + 1.00_ep1
      XO2L = LOG(XO2)
      IS = KT
      IF (X.LE.0.50_ep1) GO TO 230 !power series
      NS = 0
  100 FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.00_ep1
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 230
  110 XO2L = LOG(XO2)
      NS = INT(SXO2-FNU)
      GO TO 100
  120 CONTINUE
!
!     OVERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
!
      IF (KODE.EQ.2) GO TO 130
      IF (ALPHA.LT.1.00_ep1) GO TO 150
      Z = X/ALPHA
      RA = SQRT(1.00_ep1+Z*Z)
      GLN = LOG((1.00_ep1+RA)/Z)
      T = RA*(1.00_ep1-ETX) + ETX/(Z+RA)
      ARG = ALPHA*(T-GLN)
      IF (ARG.GT.ELIM) GO TO 610
      IF (KM.EQ.0) GO TO 140
  130 CONTINUE
!
!     UNDERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
!
      Z = X/FN
      RA = SQRT(1.00_ep1+Z*Z)
      GLN = LOG((1.00_ep1+RA)/Z)
      T = RA*(1.00_ep1-ETX) + ETX/(Z+RA)
      ARG = FN*(T-GLN)
  140 IF (ARG.LT.(-ELIM)) GO TO 280
      GO TO 190
  150 IF (X.GT.ELIM) GO TO 610
      GO TO 130
!
!     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
!
  160 IF (KM.NE.0) GO TO 170
      Y(1) = TEMP(3)
      RETURN
  170 TEMP(1) = TEMP(3)
      IN = NS
      KT = 1
      I1 = 0
  180 CONTINUE
      IS = 2
      FNI = FNI - 1.00_ep1
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 350
      Z = X/FN
      RA = SQRT(1.00_ep1+Z*Z)
      GLN = LOG((1.00_ep1+RA)/Z)
      T = RA*(1.00_ep1-ETX) + ETX/(Z+RA)
      ARG = FN*(T-GLN)
  190 CONTINUE
      I1 = ABS(3-IS)
      I1 = MAX(I1,1)
      FLGIK = 1.00_ep1
      CALL ep_asyik(X,FN,KODE,FLGIK,RA,ARG,I1,TEMP(IS))
      GO TO (180, 350, 510), IS
!
!     SERIES FOR (X/2)**2.LE.NU+1
!
  230 CONTINUE
      GLN = cfp_lngam(FNP1)
      ARG = FN*XO2L - GLN - SX
      IF (ARG.LT.(-ELIM)) GO TO 300
      EARG = EXP(ARG)
  240 CONTINUE
      S = 1.00_ep1
      IF (X.LT.TOL) GO TO 260
      AK = 3.00_ep1
      T2 = 1.00_ep1
      T = 1.00_ep1
      S1 = FN
      DO K=1,KMAXSERIES !17 replaced by 30 for quad precision code
        S2 = T2 + S1
        T = T*SXO2/S2
        S = S + T
        IF (ABS(T).LT.TOL) GO TO 260
        T2 = T2 + AK
        AK = AK + 2.00_ep1
        S1 = S1 + FN
      ENDDO
  260 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (270, 350, 500), IS
  270 EARG = EARG*FN/XO2
      FNI = FNI - 1.00_ep1
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 240
!
!     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
!
  280 Y(NN) = 0.00_ep1
      NN = NN - 1
      FNI = FNI - 1.00_ep1
      DFN = FNI + FNF
      FN = DFN
      IF (NN < 1) THEN
         GO TO 340
      ELSE IF (NN == 1) THEN
         GO TO 290
      ELSE
         GO TO 130
      END IF
  290 KT = 2
      IS = 2
      GO TO 130
  300 Y(NN) = 0.00_ep1
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.00_ep1
      DFN = FNI + FNF
      FN = DFN
      IF (NN < 1) THEN
         GO TO 340
      ELSE IF (NN == 1) THEN
         GO TO 310
      ELSE
         GO TO 320
      END IF
  310 KT = 2
      IS = 2
  320 IF (SXO2.LE.FNP1) GO TO 330
      GO TO 130
  330 ARG = ARG - XO2L + LOG(FNP1)
      IF (ARG.LT.(-ELIM)) GO TO 300
      GO TO 230
  340 NZ = N - NN
      RETURN
!
!     BACKWARD RECURSION SECTION
!
  350 CONTINUE
      NZ = N - NN
  360 CONTINUE
      IF(KT.EQ.2) GO TO 420
      S1 = TEMP(1)
      S2 = TEMP(2)
      TRX = 2.00_ep1/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      IF (IN.EQ.0) GO TO 390
!     BACKWARD RECUR TO INDEX ALPHA+NN-1
      DO I=1,IN
        S = S2
        S2 = TM*S2 + S1
        S1 = S
        DTM = DTM - 1.00_ep1
        TM = (DTM+FNF)*TRX
      ENDDO
      Y(NN) = S1
      IF (NN.EQ.1) RETURN
      Y(NN-1) = S2
      IF (NN.EQ.2) RETURN
      GO TO 400
  390 CONTINUE
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = S1
      Y(NN-1) = S2
      IF (NN.EQ.2) RETURN
  400 K = NN + 1
      DO 410 I=3,NN
        K = K - 1
        Y(K-2) = TM*Y(K-1) + Y(K)
        DTM = DTM - 1.00_ep1
        TM = (DTM+FNF)*TRX
  410 CONTINUE
      RETURN
  420 Y(1) = TEMP(2)
      RETURN
!
!     ASYMPTOTIC EXPANSION FOR X TO INFINITY
!
  430 CONTINUE
      EARG = RTTPI/SQRT(X)
      IF (KODE.EQ.2) GO TO 440
      IF (X.GT.ELIM) GO TO 610
      EARG = EARG*EXP(X)
  440 ETX = 8.00_ep1*X
      IS = KT
      IN = 0
      FN = FNU
  450 DX = FNI + FNI
      TM = 0.00_ep1
      IF (FNI.EQ.0.00_ep1 .AND. ABS(FNF).LT.TOL) GO TO 460
      TM = 4.00_ep1*FNF*(FNI+FNI+FNF)
  460 CONTINUE
      DTM = DX*DX
      S1 = ETX
      TRX = DTM - 1.00_ep1
      DX = -(TRX+TM)/ETX
      T = DX
      S = 1.00_ep1 + DX
      ATOL = TOL*ABS(S)
      S2 = 1.00_ep1
      AK = 8.00_ep1
      DO 470 K=1,KMAXAE !changed to 30 for quad precision code
        S1 = S1 + ETX
        S2 = S2 + AK
        DX = DTM - S2
        AP = DX + TM
        T = -T*AP/S1
        S = S + T
        IF (ABS(T).LE.ATOL) GO TO 480
        AK = AK + 8.00_ep1
  470 CONTINUE
  480 TEMP(IS) = S*EARG
      IF(IS.EQ.2) GO TO 360
      IS = 2
      FNI = FNI - 1.00_ep1
      DFN = FNI + FNF
      FN = DFN
      GO TO 450
!
!     BACKWARD RECURSION WITH NORMALIZATION BY
!     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
!     ZM: FOR QUAD PRECISION THE BACKWARD RECURRENCE MUST START AT A POINT ROUGHLY 2 TIMES FURTHER AWAY IN ORDER TO ACHIEVE FULL QUAD PRECISION AT NU=KM.
!
  500 CONTINUE
!     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = MAX(3.00_ep1-FN,0.00_ep1)
      KM = INT(AKM)
      TFN = FN + KM
      TA = (GLN+TFN-0.9189385332046727417803297364056176_ep1-0.08333333333333333333333333333333333_ep1/TFN)/(TFN+0.50_ep1) !universal factors: Log(sqrt(2*pi)); 1/12
      TA = XO2L - TA
      TB = -(1.00_ep1-1.00_ep1/TFN)/TFN
      AIN = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.50_ep1
      IN = INT(2.3*AIN)
      IN = IN + KM
      GO TO 520
  510 CONTINUE
!     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      T = 1.00_ep1/(FN*RA)
      AIN = TOLLN/(GLN+SQRT(GLN*GLN+T*TOLLN)) + 1.50_ep1
      IN = INT(2.3*AIN)
      IF (IN.GT.INLIM) GO TO 160
  520 CONTINUE
      TRX = 2.00_ep1/X
      DTM = FNI + IN
      TM = (DTM+FNF)*TRX
      TA = 0.00_ep1
      TB = TOL
      KK = 1
  530 CONTINUE
!
!     BACKWARD RECUR UNINDEXED
!
      DO 540 I=1,IN
        S = TB
        TB = TM*TB + TA
        TA = S
        DTM = DTM - 1.00_ep1
        TM = (DTM+FNF)*TRX
  540 CONTINUE
!     NORMALIZATION
      IF (KK.NE.1) GO TO 550
      TA = (TA/TB)*TEMP(3)
      TB = TEMP(3)
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 530
  550 Y(NN) = TB
      NZ = N - NN
      IF (NN.EQ.1) RETURN
      TB = TM*TB + TA
      K = NN - 1
      Y(K) = TB
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.00_ep1
      TM = (DTM+FNF)*TRX
      KM = K - 1
!
!     BACKWARD RECUR INDEXED
!
      DO 560 I=1,KM
        Y(K-1) = TM*Y(K) + Y(K+1)
        DTM = DTM - 1.00_ep1
        TM = (DTM+FNF)*TRX
        K = K - 1
  560 CONTINUE
      RETURN
!
!
!
  570 CONTINUE
      CALL XERMSG ('SLATEC', 'ep_besi', 'SCALING OPTION, KODE, NOT 1 OR 2.', 2, 1)
      RETURN
  580 CONTINUE
      CALL XERMSG ('SLATEC', 'ep_besi', 'ORDER, ALPHA, LESS THAN ZERO.', 2, 1)
      RETURN
  590 CONTINUE
      CALL XERMSG ('SLATEC', 'ep_besi', 'N LESS THAN ONE.', 2, 1)
      RETURN
  600 CONTINUE
      CALL XERMSG ('SLATEC', 'ep_besi', 'X LESS THAN ZERO.', 2, 1)
      RETURN
  610 CONTINUE
      CALL XERMSG ('SLATEC', 'ep_besi', 'OVERFLOW, X TOO LARGE FOR KODE = 1.', 6, 1)
      RETURN
      END SUBROUTINE ep_besi
!> \verbatim
!>***PURPOSE  Compute an N member sequence of J Bessel functions
!>            J/SUB(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA
!>            and X.
!>***LIBRARY   SLATEC
!>***CATEGORY  C10A3
!>***TYPE      REAL(kind=wp) (BESJ-S, wp_besj-D)
!>***KEYWORDS  J BESSEL FUNCTION, SPECIAL FUNCTIONS
!>***AUTHOR  Amos, D. E., (SNLA)
!>           Daniel, S. L., (SNLA)
!>           Weston, M. K., (SNLA)
!>***DESCRIPTION
!>
!>     Abstract  **** a REAL(kind=wp) routine ****
!>         wp_besj computes an N member sequence of J Bessel functions
!>         J/sub(ALPHA+K-1)/(X), K=1,...,N for non-negative ALPHA and X.
!>         A combination of the power series, the asymptotic expansion
!>         for X to infinity and the uniform asymptotic expansion for
!>         NU to infinity are applied over subdivisions of the (NU,X)
!>         plane.  For values of (NU,X) not covered by one of these
!>         formulae, the order is incremented or decremented by integer
!>         values into a region where one of the formulae apply. Backward
!>         recursion is applied to reduce orders by integer values except
!>         where the entire sequence lies in the oscillatory region.  In
!>         this case forward recursion is stable and values from the
!>         asymptotic expansion for X to infinity start the recursion
!>         when it is efficient to do so. Leading terms of the series and
!>         uniform expansion are tested for underflow.  If a sequence is
!>         requested and the last member would underflow, the result is
!>         set to zero and the next lower order tried, etc., until a
!>         member comes on scale or all members are set to zero.
!>         Overflow cannot occur.
!>
!>         The maximum number of significant digits obtainable
!>         is the smaller of 14 and the number of digits carried in
!>         REAL(kind=wp) arithmetic.
!>
!>     Description of Arguments
!>
!>         Input      X,ALPHA are REAL(kind=wp)
!>           X      - X .GE. 0.00_wp
!>           ALPHA  - order of first member of the sequence,
!>                    ALPHA .GE. 0.00_wp
!>           N      - number of members in the sequence, N .GE. 1
!>
!>         Output     Y is REAL(kind=wp)
!>           Y      - a vector whose first N components contain
!>                    values for J/sub(ALPHA+K-1)/(X), K=1,...,N
!>           NZ     - number of components of Y set to zero due to
!>                    underflow,
!>                    NZ=0   , normal return, computation completed
!>                    NZ .NE. 0, last NZ components of Y set to zero,
!>                             Y(K)=0.00_wp, K=N-NZ+1,...,N.
!>
!>     Error Conditions
!>         Improper input arguments - a fatal error
!>         Underflow  - a non-fatal error (NZ .NE. 0)
!>
!>***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
!>                 subroutines IBESS and JBESS for Bessel functions
!>                 I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0, ACM
!>                 Transactions on Mathematical Software 3, (1977),
!>                 pp. 76-92.
!>               F. W. J. Olver, Tables of Bessel Functions of Moderate
!>                 or Large Orders, NPL Mathematical Tables 6, Her
!>                 Majesty's Stationery Office, London, 1962.
!>***ROUTINES CALLED  F1MACH, wp_asyjy, cfp_jairy, cfp_lngam, I1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   750101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   890911  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900326  Removed duplicate information from DESCRIPTION section.
!>           (WRB)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!> \endverbatim
!
      SUBROUTINE wp_besj (X, ALPHA, N, Y, NZ)
      IMPLICIT NONE
!
      !EXTERNAL wp_jairy
      INTEGER I,IALP,IDALP,IFLW,IN,INLIM,IS,I1,I2,K,KK,KM,KT,N,NN, NS,NZ
      REAL(kind=wp) AK,AKM,ALPHA,ANS,AP,ARG,COEF,DALPHA,DFN,DTM,    &
     &           EARG,ELIM1,ETX,FIDAL,FLGJY,FN,FNF,FNI,FNP1,FNU,       &
     &           FNULIM,GLN,PDF,PIDT,PP,RDEN,RELB,RTTP,RTWO,RTX,RZDEN, &
     &           S,SA,SB,SXO2,S1,S2,T,TA,TAU,TB,TEMP,TFN,TM,TOL,       &
     &           TOLLN,TRX,TX,T1,T2,WK,X,XO2,XO2L,Y,SLIM,RTOL          
      SAVE RTWO, PDF, RTTP, PIDT, PP, INLIM, FNULIM
      DIMENSION Y(*), TEMP(3), FNULIM(2), PP(4), WK(7)
      DATA RTWO, PDF,  RTTP, PIDT / 1.34839972492648E+00_wp,7.85398163397448E-01_wp,7.97884560802865E-01_wp,1.57079632679490E+00_wp/
      DATA PP(1),PP(2),PP(3),PP(4)/ 8.72909153935547E+00_wp,2.65693932265030E-01_wp,1.24578576865586E-01_wp,7.70133747430388E-04_wp/
      DATA INLIM           /      150            /
      DATA FNULIM(1), FNULIM(2) /      100.00_wp,     60.00_wp     /
!***FIRST EXECUTABLE STATEMENT  wp_besj
      NZ = 0
      KT = 1
      NS=0
!     I1MACH(14) REPLACES I1MACH(11) IN A REAL(kind=wp) CODE
!     I1MACH(15) REPLACES I1MACH(12) IN A REAL(kind=wp) CODE
      TA = F1MACH(3,wp_dummy)
      TOL = MAX(TA,1.0E-15_wp)
      I1 = I1MACH(14) + 1
      I2 = I1MACH(15)
      TB = F1MACH(5,wp_dummy)
      ELIM1 = -2.3030_wp*(I2*TB+3.00_wp)
      RTOL=1.00_wp/TOL
      SLIM=F1MACH(1,wp_dummy)*RTOL*1.0E+3
!     TOLLN = -LN(TOL)
      TOLLN = 2.3030_wp*TB*I1
      TOLLN = MIN(TOLLN,34.53880_wp)
      IF (N < 1) THEN
         GO TO 720
      ELSE IF (N == 1) THEN
         GO TO 10
      ELSE
         GO TO 20
      END IF
   10 KT = 2
   20 NN = N
      IF (X < 0) THEN
         GO TO 730
      ELSE IF (X == 0) THEN
         GO TO 30
      ELSE
         GO TO 80
      END IF
   30 IF (ALPHA < 0) THEN
         GO TO 710
      ELSE IF (ALPHA == 0) THEN
         GO TO 40
      ELSE
         GO TO 50
      END IF
   40 Y(1) = 1.00_wp
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.00_wp
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.00_wp) GO TO 710
!
      IALP = INT(ALPHA)
      FNI = IALP + N - 1
      FNF = ALPHA - IALP
      DFN = FNI + FNF
      FNU = DFN
      XO2 = X*0.50_wp
      SXO2 = XO2*XO2
!
!     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
!     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
!     APPLIED.
!
      IF (SXO2.LE.(FNU+1.00_wp)) GO TO 90
      TA = MAX(20.00_wp,FNU)
      IF (X.GT.TA) GO TO 120
      IF (X.GT.12.00_wp) GO TO 110
      XO2L = LOG(XO2)
      NS = INT(SXO2-FNU) + 1
      GO TO 100
   90 FN = FNU
      FNP1 = FN + 1.00_wp
      XO2L = LOG(XO2)
      IS = KT
      IF (X.LE.0.500_wp) GO TO 330
      NS = 0
  100 FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.00_wp
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 330
  110 ANS = MAX(36.00_wp-FNU,0.00_wp)
      NS = INT(ANS)
      FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 130
  120 CONTINUE
      RTX = SQRT(X)
      TAU = RTWO*RTX
      TA = TAU + FNULIM(KT)
      IF (FNU.LE.TA) GO TO 480
      FN = FNU
      IS = KT
!
!     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
!
  130 CONTINUE
      I1 = ABS(3-IS)
      I1 = MAX(I1,1)
      FLGJY = 1.00_wp
      CALL wp_asyjy(wp_jairy,X,FN,FLGJY,I1,TEMP(IS),WK,IFLW)
      IF(IFLW.NE.0) GO TO 380
      GO TO (320, 450, 620), IS
  310 TEMP(1) = TEMP(3)
      KT = 1
  320 IS = 2
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 450
      GO TO 130
!
!     SERIES FOR (X/2)**2.LE.NU+1
!
  330 CONTINUE
      GLN = cfp_lngam(FNP1)
      ARG = FN*XO2L - GLN
      IF (ARG.LT.(-ELIM1)) GO TO 400
      EARG = EXP(ARG)
  340 CONTINUE
      S = 1.00_wp
      IF (X.LT.TOL) GO TO 360
      AK = 3.00_wp
      T2 = 1.00_wp
      T = 1.00_wp
      S1 = FN
      DO 350 K=1,17
        S2 = T2 + S1
        T = -T*SXO2/S2
        S = S + T
        IF (ABS(T).LT.TOL) GO TO 360
        T2 = T2 + AK
        AK = AK + 2.00_wp
        S1 = S1 + FN
  350 CONTINUE
  360 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (370, 450, 610), IS
  370 EARG = EARG*FN/XO2
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 340
!
!     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
!     UNDERFLOW CAN ONLY OCCUR FOR NS=0 SINCE THE ORDER MUST BE LARGER
!     THAN 36. THEREFORE, NS NEE NOT BE TESTED.
!
  380 Y(NN) = 0.00_wp
      NN = NN - 1
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      IF (NN < 1) THEN
         GO TO 440
      ELSE IF (NN == 1) THEN
         GO TO 390
      ELSE
         GO TO 130
      END IF
  390 KT = 2
      IS = 2
      GO TO 130
  400 Y(NN) = 0.00_wp
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.00_wp
      DFN = FNI + FNF
      FN = DFN
      IF (NN < 1) THEN
         GO TO 440
      ELSE IF (NN == 1) THEN
         GO TO 410
      ELSE
         GO TO 420
      END IF
  410 KT = 2
      IS = 2
  420 IF (SXO2.LE.FNP1) GO TO 430
      GO TO 130
  430 ARG = ARG - XO2L + LOG(FNP1)
      IF (ARG.LT.(-ELIM1)) GO TO 400
      GO TO 330
  440 NZ = N - NN
      RETURN
!
!     BACKWARD RECURSION SECTION
!
  450 CONTINUE
      IF(NS.NE.0) GO TO 451
      NZ = N - NN
      IF (KT.EQ.2) GO TO 470
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = TEMP(1)
      Y(NN-1) = TEMP(2)
      IF (NN.EQ.2) RETURN
  451 CONTINUE
      TRX = 2.00_wp/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      AK=1.00_wp
      TA=TEMP(1)
      TB=TEMP(2)
      IF(ABS(TA).GT.SLIM) GO TO 455
      TA=TA*RTOL
      TB=TB*RTOL
      AK=TOL
  455 CONTINUE
      KK=2
      IN=NS-1
      IF(IN.EQ.0) GO TO 690
      IF(NS.NE.0) GO TO 670
      K=NN-2
      DO 460 I=3,NN
        S=TB
        TB = TM*TB - TA
        TA=S
        Y(K)=TB*AK
        DTM = DTM - 1.00_wp
        TM = (DTM+FNF)*TRX
        K = K - 1
  460 CONTINUE
      RETURN
  470 Y(1) = TEMP(2)
      RETURN
!
!     ASYMPTOTIC EXPANSION FOR X TO INFINITY WITH FORWARD RECURSION IN
!     OSCILLATORY REGION X.GT.MAX(20, NU), PROVIDED THE LAST MEMBER
!     OF THE SEQUENCE IS ALSO IN THE REGION.
!
  480 CONTINUE
      IN = INT(ALPHA-TAU+2.00_wp)
      IF (IN.LE.0) GO TO 490
      IDALP = IALP - IN - 1
      KT = 1
      GO TO 500
  490 CONTINUE
      IDALP = IALP
      IN = 0
  500 IS = KT
      FIDAL = IDALP
      DALPHA = FIDAL + FNF
      ARG = X - PIDT*DALPHA - PDF
      SA = SIN(ARG)
      SB = COS(ARG)
      COEF = RTTP/RTX
      ETX = 8.00_wp*X
  510 CONTINUE
      DTM = FIDAL + FIDAL
      DTM = DTM*DTM
      TM = 0.00_wp
      IF (FIDAL.EQ.0.00_wp .AND. ABS(FNF).LT.TOL) GO TO 520
      TM = 4.00_wp*FNF*(FIDAL+FIDAL+FNF)
  520 CONTINUE
      TRX = DTM - 1.00_wp
      T2 = (TRX+TM)/ETX
      S2 = T2
      RELB = TOL*ABS(T2)
      T1 = ETX
      S1 = 1.00_wp
      FN = 1.00_wp
      AK = 8.00_wp
      DO 530 K=1,13
        T1 = T1 + ETX
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = -T2*AP/T1
        S1 = S1 + T2
        T1 = T1 + ETX
        AK = AK + 8.00_wp
        FN = FN + AK
        TRX = DTM - FN
        AP = TRX + TM
        T2 = T2*AP/T1
        S2 = S2 + T2
        IF (ABS(T2).LE.RELB) GO TO 540
        AK = AK + 8.00_wp
  530 CONTINUE
  540 TEMP(IS) = COEF*(S1*SB-S2*SA)
      IF(IS.EQ.2) GO TO 560
      FIDAL = FIDAL + 1.00_wp
      DALPHA = FIDAL + FNF
      IS = 2
      TB = SA
      SA = -SB
      SB = TB
      GO TO 510
!
!     FORWARD RECURSION SECTION
!
  560 IF (KT.EQ.2) GO TO 470
      S1 = TEMP(1)
      S2 = TEMP(2)
      TX = 2.00_wp/X
      TM = DALPHA*TX
      IF (IN.EQ.0) GO TO 580
!
!     FORWARD RECUR TO INDEX ALPHA
!
      DO 570 I=1,IN
        S = S2
        S2 = TM*S2 - S1
        TM = TM + TX
        S1 = S
  570 CONTINUE
      IF (NN.EQ.1) GO TO 600
      S = S2
      S2 = TM*S2 - S1
      TM = TM + TX
      S1 = S
  580 CONTINUE
!
!     FORWARD RECUR FROM INDEX ALPHA TO ALPHA+N-1
!
      Y(1) = S1
      Y(2) = S2
      IF (NN.EQ.2) RETURN
      DO 590 I=3,NN
        Y(I) = TM*Y(I-1) - Y(I-2)
        TM = TM + TX
  590 CONTINUE
      RETURN
  600 Y(1) = S2
      RETURN
!
!     BACKWARD RECURSION WITH NORMALIZATION BY
!     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
!
  610 CONTINUE
!     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = MAX(3.00_wp-FN,0.00_wp)
      KM = INT(AKM)
      TFN = FN + KM
      TA = (GLN+TFN-0.91893853320_wp-0.08333333330_wp/TFN)/(TFN+0.50_wp)
      TA = XO2L - TA
      TB = -(1.00_wp-1.50_wp/TFN)/TFN
      AKM = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.50_wp
      IN = KM + INT(AKM)
      GO TO 660
  620 CONTINUE
!     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      GLN = WK(3) + WK(2)
      IF (WK(6).GT.30.00_wp) GO TO 640
      RDEN = (PP(4)*WK(6)+PP(3))*WK(6) + 1.00_wp
      RZDEN = PP(1) + PP(2)*WK(6)
      TA = RZDEN/RDEN
      IF (WK(1).LT.0.100_wp) GO TO 630
      TB = GLN/WK(5)
      GO TO 650
  630 TB=(1.2599210490_wp+(0.16798947300_wp+0.08879443580_wp*WK(1))*WK(1))/WK(7)
      GO TO 650
  640 CONTINUE
      TA = 0.50_wp*TOLLN/WK(4)
      TA=((0.04938271600_wp*TA-0.11111111110_wp)*TA+0.66666666670_wp)*TA*WK(6)
      IF (WK(1).LT.0.100_wp) GO TO 630
      TB = GLN/WK(5)
  650 IN = INT(TA/TB+1.50_wp)
      IF (IN.GT.INLIM) GO TO 310
  660 CONTINUE
      DTM = FNI + IN
      TRX = 2.00_wp/X
      TM = (DTM+FNF)*TRX
      TA = 0.00_wp
      TB = TOL
      KK = 1
      AK=1.00_wp
  670 CONTINUE
!
!     BACKWARD RECUR UNINDEXED
!
      DO 680 I=1,IN
        S = TB
        TB = TM*TB - TA
        TA = S
        DTM = DTM - 1.00_wp
        TM = (DTM+FNF)*TRX
  680 CONTINUE
!     NORMALIZATION
      IF (KK.NE.1) GO TO 690
      S=TEMP(3)
      SA=TA/TB
      TA=S
      TB=S
      IF(ABS(S).GT.SLIM) GO TO 685
      TA=TA*RTOL
      TB=TB*RTOL
      AK=TOL
  685 CONTINUE
      TA=TA*SA
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 670
  690 Y(NN) = TB*AK
      NZ = N - NN
      IF (NN.EQ.1) RETURN
      K = NN - 1
      S=TB
      TB = TM*TB - TA
      TA=S
      Y(K)=TB*AK
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.00_wp
      TM = (DTM+FNF)*TRX
      K=NN-2
!
!     BACKWARD RECUR INDEXED
!
      DO 700 I=3,NN
        S=TB
        TB = TM*TB - TA
        TA=S
        Y(K)=TB*AK
        DTM = DTM - 1.00_wp
        TM = (DTM+FNF)*TRX
        K = K - 1
  700 CONTINUE
      RETURN
!
!
!
  710 CONTINUE
      CALL XERMSG ('SLATEC', 'wp_besj', 'ORDER, ALPHA, LESS THAN ZERO.', 2, 1)
      RETURN
  720 CONTINUE
      CALL XERMSG ('SLATEC', 'wp_besj', 'N LESS THAN ONE.', 2, 1)
      RETURN
  730 CONTINUE
      CALL XERMSG ('SLATEC', 'wp_besj', 'X LESS THAN ZERO.', 2, 1)
      RETURN
      END SUBROUTINE
      !> Quad precision version of wp_besj.
      !> \warning Not implemented yet.
      SUBROUTINE ep_besj (X, ALPHA, N, Y, NZ)
      IMPLICIT NONE
!
      !EXTERNAL ep1_jairy
      INTEGER N,NZ
      REAL(kind=ep1) X,ALPHA,Y
      DIMENSION Y(*)

         CALL XERMSG('SLATEC','ep_besj','Quad precision version not implemented yet.',1,1)

      END SUBROUTINE
!> \verbatim
!>***PURPOSE  Compute the binomial coefficients.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C1
!>***TYPE      REAL(kind=wp) (BINOM-S, wp_binom-D)
!>***KEYWORDS  BINOMIAL COEFFICIENTS, FNLIB, SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> wp_binom(N,M) calculates the REAL(kind=wp) binomial coefficient
!> for integer arguments N and M.  The result is (N!)/((M!)(N-M)!).
!> The argument 'p' is a dummy variable of the same type as the result, i.e.
!> if the double precision result is required then p must be a dummy variable 
!> of type kind=wp. For quad precision p must be kind=ep1.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, wp_9lgmc, cfp_lnrel, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770601  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!> \endverbatim
!
      REAL(kind=wp) FUNCTION wp_binom (N, M, p)
      IMPLICIT NONE
!
      INTEGER N, M
      INTEGER K, I
      REAL(kind=wp) CORR, FINTMX, SQ2PIL, XK, XN, XNK, BILNMX, p
      LOGICAL FIRST
      SAVE SQ2PIL, BILNMX, FINTMX, FIRST
      DATA SQ2PIL / 0.918938533204672741780329736405620_wp /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_binom
      IF (FIRST) THEN
         BILNMX = LOG(F1MACH(2,wp_dummy)) - 0.00010_wp
         FINTMX = 0.90_wp/F1MACH(3,wp_dummy)
      ENDIF
      FIRST = .FALSE.
!
      IF (N .LT. 0 .OR. M .LT. 0) CALL XERMSG ('SLATEC', 'wp_binom', 'N OR M LT ZERO', 1, 2)
      IF (N .LT. M) CALL XERMSG ('SLATEC', 'wp_binom', 'N LT M', 2, 2)
!
      K = MIN (M, N-M)
      IF (K.GT.20) GO TO 30
      IF (K*LOG(AMAX0(N,1)).GT.BILNMX) GO TO 30
!
      wp_binom = 1.00_wp
      IF (K.EQ.0) RETURN
      DO 20 I=1,K
        XN = N - I + 1
        XK = I
        wp_binom = wp_binom * (XN/XK)
 20   CONTINUE
!
      IF (wp_binom.LT.FINTMX) wp_binom = AINT (wp_binom+0.50_wp)
      RETURN
!
! IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
 30   IF (K .LT. 9) CALL XERMSG ('SLATEC', 'wp_binom', 'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
!
      XN = N + 1
      XK = K + 1
      XNK = N - K + 1
!
      CORR = cfp_9lgmc(XN) - cfp_9lgmc(XK) - cfp_9lgmc(XNK)
      wp_binom = XK*LOG(XNK/XK) - XN*cfp_lnrel(-(XK-1.00_wp)/XN)-0.50_wp*LOG(XN*XNK/XK) + 1.00_wp - SQ2PIL + CORR
!
      IF (wp_binom .GT. BILNMX) CALL XERMSG ('SLATEC', 'wp_binom', 'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
!
      wp_binom = EXP (wp_binom)
      IF (wp_binom.LT.FINTMX) wp_binom = AINT (wp_binom+0.50_wp)
!
      RETURN
      END FUNCTION
!
      !> Quad precision version of wp_binom.
      REAL(kind=ep1) FUNCTION ep_binom (N, M, p)
      IMPLICIT NONE
!
      INTEGER N, M
      INTEGER K, I
      REAL(kind=ep1) CORR, FINTMX, SQ2PIL, XK, XN, XNK, BILNMX, p
      LOGICAL FIRST
      SAVE SQ2PIL, BILNMX, FINTMX, FIRST
      DATA SQ2PIL / 0.918938533204672741780329736405620_ep1 /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_binom
      IF (FIRST) THEN
         BILNMX = LOG(F1MACH(2,ep_dummy)) - 0.00010_ep1
         FINTMX = 0.90_ep1/F1MACH(3,ep_dummy)
      ENDIF
      FIRST = .FALSE.
!
      IF (N .LT. 0 .OR. M .LT. 0) CALL XERMSG ('SLATEC', 'ep_binom', 'N OR M LT ZERO', 1, 2)
      IF (N .LT. M) CALL XERMSG ('SLATEC', 'ep_binom', 'N LT M', 2, 2)
!
      K = MIN (M, N-M)
      IF (K.GT.20) GO TO 30
      IF (K*LOG(AMAX0(N,1)).GT.BILNMX) GO TO 30
!
      ep_binom = 1.00_ep1
      IF (K.EQ.0) RETURN
      DO 20 I=1,K
        XN = N - I + 1
        XK = I
        ep_binom = ep_binom * (XN/XK)
 20   CONTINUE
!
      IF (ep_binom.LT.FINTMX) ep_binom = AINT (ep_binom+0.50_ep1)
      RETURN
!
! IF K.LT.9, APPROX IS NOT VALID AND ANSWER IS CLOSE TO THE OVERFLOW LIM
 30   IF (K .LT. 9) CALL XERMSG ('SLATEC', 'ep_binom', 'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
!
      XN = N + 1
      XK = K + 1
      XNK = N - K + 1
!
      CORR = cfp_9lgmc(XN) - cfp_9lgmc(XK) - cfp_9lgmc(XNK)
      ep_binom = XK*LOG(XNK/XK) - XN*cfp_lnrel(-(XK-1.00_ep1)/XN)-0.50_ep1*LOG(XN*XNK/XK) + 1.00_ep1 - SQ2PIL + CORR
!
      IF (ep_binom .GT. BILNMX) CALL XERMSG ('SLATEC', 'ep_binom', 'RESULT OVERFLOWS BECAUSE N AND/OR M TOO BIG', 3, 2)
!
      ep_binom = EXP (ep_binom)
      IF (ep_binom.LT.FINTMX) ep_binom = AINT (ep_binom+0.50_ep1)
!
      RETURN
      END FUNCTION
!> \verbatim
!>***PURPOSE  Evaluate a Chebyshev series.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C3A2
!>***TYPE      REAL(kind=wp) (CSEVL-S, wp_csevl-D)
!>***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!>  Evaluate the N-term Chebyshev series CS at X.  Adapted from
!>  a method presented in the paper by Broucke referenced below.
!>
!>       Input Arguments --
!>  X    value at which the series is to be evaluated.
!>  CS   array of N terms of a Chebyshev series.  In evaluating
!>       CS, only half the first coefficient is summed.
!>  N    number of terms in array CS.
!>
!>***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
!>                 Chebyshev series, Algorithm 446, Communications of
!>                 the A.C.M. 16, (1973) pp. 254-256.
!>               L. Fox and I. B. Parker, Chebyshev Polynomials in
!>                 Numerical Analysis, Oxford University Press, 1968,
!>                 page 56.
!>***ROUTINES CALLED  F1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770401  DATE WRITTEN
!>   890831  Modified array declarations.  (WRB)
!>   890831  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900329  Prologued revised extensively and code rewritten to allow
!>           X to be slightly outside interval (-1,+1).  (WRB)
!>   920501  Reformatted the REFERENCES section.  (WRB)
!> \endverbatim
!
      REAL(kind=wp) FUNCTION wp_csevl (X, CS, N)
      IMPLICIT NONE
!
      INTEGER :: N
      REAL(kind=wp) B0, B1, B2, CS(*), ONEPL, TWOX, X 
      INTEGER I, NI
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_csevl
      IF (FIRST) ONEPL = 1.00_wp + F1MACH(4,wp_dummy)
      FIRST = .FALSE.
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'wp_csevl', 'NUMBER OF TERMS .LE. 0', 2, 2)
      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'wp_csevl', 'NUMBER OF TERMS .GT. 1000', 3, 2)
      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'wp_csevl', 'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
!
      B1 = 0.00_wp
      B0 = 0.00_wp
      TWOX = 2.00_wp*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
!
      wp_csevl = 0.50_wp*(B0-B2)
!
      RETURN
      END FUNCTION
!
      !> Quad precision version of wp_csevl.
      REAL(kind=ep1) FUNCTION ep_csevl (X, CS, N)
      IMPLICIT NONE
!
      INTEGER :: N
      REAL(kind=ep1) B0, B1, B2, CS(*), ONEPL, TWOX, X 
      INTEGER I, NI
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_csevl
      IF (FIRST) ONEPL = 1.00_ep1 + F1MACH(4,ep_dummy)
      FIRST = .FALSE.
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'ep_csevl', 'NUMBER OF TERMS .LE. 0', 2, 2)
      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'ep_csevl', 'NUMBER OF TERMS .GT. 1000', 3, 2)
      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'ep_csevl', 'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
!
      B1 = 0.00_ep1
      B0 = 0.00_ep1
      TWOX = 2.00_ep1*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
!
      ep_csevl = 0.50_ep1*(B0-B2)
!
      RETURN
      END FUNCTION
!> \verbatim
!>***PURPOSE  Compute the minimum and maximum bounds for the argument in
!>            the Gamma function.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7A, R2
!>***TYPE      REAL(kind=wp) (GAMLIM-S, wp_gamlm-D)
!>***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> Calculate the minimum and maximum legal bounds for X in gamma(X).
!> XMIN and XMAX are not the only bounds, but they are the only non-
!> trivial ones to calculate.
!>
!>             Output Arguments --
!> XMIN   REAL(kind=wp) minimum legal value of X in gamma(X).  Any
!>        smaller value of X might result in underflow.
!> XMAX   REAL(kind=wp) maximum legal value of X in gamma(X).  Any
!>        larger value of X might cause overflow.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770601  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!> \endverbatim
!
      SUBROUTINE wp_gamlm (XMIN, XMAX)
      IMPLICIT NONE
!
      REAL(kind=wp) XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD 
      INTEGER I
!***FIRST EXECUTABLE STATEMENT  wp_gamlm
      ALNSML = LOG(F1MACH(1,wp_dummy))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.50_wp)*XLN - XMIN - 0.22580_wp + ALNSML) / (XMIN*XLN+0.50_wp)
        IF (ABS(XMIN-XOLD).LT.0.0050_wp) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'wp_gamlm', 'UNABLE TO FIND XMIN', 1, 2)
!
 20   XMIN = -XMIN + 0.010_wp
!
      ALNBIG = LOG (F1MACH(2,wp_dummy))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.50_wp)*XLN - XMAX + 0.91890_wp - ALNBIG) / (XMAX*XLN-0.50_wp)
        IF (ABS(XMAX-XOLD).LT.0.0050_wp) GO TO 40
 30   CONTINUE
      CALL XERMSG ('SLATEC', 'wp_gamlm', 'UNABLE TO FIND XMAX', 2, 2)
!
 40   XMAX = XMAX - 0.010_wp
      XMIN = MAX (XMIN, -XMAX+1.0_wp)
!
      RETURN
      END SUBROUTINE
!
      !> Quad precision version of wp_gamlm.
      SUBROUTINE ep_gamlm (XMIN, XMAX)
      IMPLICIT NONE
!
      REAL(kind=ep1) XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD 
      INTEGER I
!***FIRST EXECUTABLE STATEMENT  ep_gamlm
      ALNSML = LOG(F1MACH(1,ep_dummy))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.50_ep1)*XLN - XMIN - 0.22580_ep1 + ALNSML) / (XMIN*XLN+0.50_ep1)
        IF (ABS(XMIN-XOLD).LT.0.0050_ep1) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'ep_gamlm', 'UNABLE TO FIND XMIN', 1, 2)
!
 20   XMIN = -XMIN + 0.010_ep1
!
      ALNBIG = LOG (F1MACH(2,ep_dummy))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.50_ep1)*XLN - XMAX + 0.91890_ep1 - ALNBIG) / (XMAX*XLN-0.50_ep1)
        IF (ABS(XMAX-XOLD).LT.0.0050_ep1) GO TO 40
 30   CONTINUE
      CALL XERMSG ('SLATEC', 'ep_gamlm', 'UNABLE TO FIND XMAX', 2, 2)
!
 40   XMAX = XMAX - 0.010_ep1
      XMIN = MAX (XMIN, -XMAX+1.0_ep1)
!
      RETURN
      END SUBROUTINE
!> \verbatim
!>***PURPOSE  Compute the complete Gamma function.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7A
!>***TYPE      REAL(kind=wp) (GAMMA-S, wp_gamma-D, CGAMMA-C)
!>***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> wp_gamma(X) calculates the REAL(kind=wp) complete Gamma function
!> for REAL(kind=wp) argument X.
!>
!> Series for GAM        on the interval  0.          to  1.00000E+00_wp
!>                                        with weighted error   5.79E-32
!>                                         log weighted error  31.24
!>                               significant figures required  30.00
!>                                    decimal places required  32.05
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, cfp_9lgmc, cfp_csevl, cfp_gamlm, cfp_initds, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770601  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   890911  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   920618  Removed space from variable name.  (RWC, WRB)
!> \endverbatim
!
      REAL(kind=wp) FUNCTION wp_gamma (X)
      IMPLICIT NONE
!
      REAL(kind=wp) X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX, XMIN, Y 
      LOGICAL FIRST
!
      INTEGER I, N
      SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
      DATA GAMCS(  1) / +.8571195590989331421920062399942E-2_wp      /
      DATA GAMCS(  2) / +.4415381324841006757191315771652E-2_wp      /
      DATA GAMCS(  3) / +.5685043681599363378632664588789E-1_wp      /
      DATA GAMCS(  4) / -.4219835396418560501012500186624E-2_wp      /
      DATA GAMCS(  5) / +.1326808181212460220584006796352E-2_wp      /
      DATA GAMCS(  6) / -.1893024529798880432523947023886E-3_wp      /
      DATA GAMCS(  7) / +.3606925327441245256578082217225E-4_wp      /
      DATA GAMCS(  8) / -.6056761904460864218485548290365E-5_wp      /
      DATA GAMCS(  9) / +.1055829546302283344731823509093E-5_wp      /
      DATA GAMCS( 10) / -.1811967365542384048291855891166E-6_wp      /
      DATA GAMCS( 11) / +.3117724964715322277790254593169E-7_wp      /
      DATA GAMCS( 12) / -.5354219639019687140874081024347E-8_wp      /
      DATA GAMCS( 13) / +.9193275519859588946887786825940E-9_wp      /
      DATA GAMCS( 14) / -.1577941280288339761767423273953E-9_wp      /
      DATA GAMCS( 15) / +.2707980622934954543266540433089E-10_wp     /
      DATA GAMCS( 16) / -.4646818653825730144081661058933E-11_wp     /
      DATA GAMCS( 17) / +.7973350192007419656460767175359E-12_wp     /
      DATA GAMCS( 18) / -.1368078209830916025799499172309E-12_wp     /
      DATA GAMCS( 19) / +.2347319486563800657233471771688E-13_wp     /
      DATA GAMCS( 20) / -.4027432614949066932766570534699E-14_wp     /
      DATA GAMCS( 21) / +.6910051747372100912138336975257E-15_wp     /
      DATA GAMCS( 22) / -.1185584500221992907052387126192E-15_wp     /
      DATA GAMCS( 23) / +.2034148542496373955201026051932E-16_wp     /
      DATA GAMCS( 24) / -.3490054341717405849274012949108E-17_wp     /
      DATA GAMCS( 25) / +.5987993856485305567135051066026E-18_wp     /
      DATA GAMCS( 26) / -.1027378057872228074490069778431E-18_wp     /
      DATA GAMCS( 27) / +.1762702816060529824942759660748E-19_wp     /
      DATA GAMCS( 28) / -.3024320653735306260958772112042E-20_wp     /
      DATA GAMCS( 29) / +.5188914660218397839717833550506E-21_wp     /
      DATA GAMCS( 30) / -.8902770842456576692449251601066E-22_wp     /
      DATA GAMCS( 31) / +.1527474068493342602274596891306E-22_wp     /
      DATA GAMCS( 32) / -.2620731256187362900257328332799E-23_wp     /
      DATA GAMCS( 33) / +.4496464047830538670331046570666E-24_wp     /
      DATA GAMCS( 34) / -.7714712731336877911703901525333E-25_wp     /
      DATA GAMCS( 35) / +.1323635453126044036486572714666E-25_wp     /
      DATA GAMCS( 36) / -.2270999412942928816702313813333E-26_wp     /
      DATA GAMCS( 37) / +.3896418998003991449320816639999E-27_wp     /
      DATA GAMCS( 38) / -.6685198115125953327792127999999E-28_wp     /
      DATA GAMCS( 39) / +.1146998663140024384347613866666E-28_wp     /
      DATA GAMCS( 40) / -.1967938586345134677295103999999E-29_wp     /
      DATA GAMCS( 41) / +.3376448816585338090334890666666E-30_wp     /
      DATA GAMCS( 42) / -.5793070335782135784625493333333E-31_wp     /
      DATA PI / 3.141592653589793238462643383279500_wp /
      DATA SQ2PIL / 0.918938533204672741780329736405620_wp /
      DATA FIRST /.TRUE./
      INTEGER NGAM
!***FIRST EXECUTABLE STATEMENT  wp_gamma
      IF (FIRST) THEN
         NGAM = cfp_initds (GAMCS, 42, 0.1*REAL(F1MACH(3,wp_dummy)) )
!
         CALL cfp_gamlm (XMIN, XMAX)
         DXREL = SQRT(F1MACH(4,wp_dummy))
      ENDIF
      FIRST = .FALSE.
!
      Y = ABS(X)
      IF (Y.GT.10.0_wp) GO TO 50
!
! COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
!
      N = X
      IF (X.LT.0.0_wp) N = N - 1
      Y = X - N
      N = N - 1
      wp_gamma = 0.93750_wp + cfp_csevl (2.0_wp*Y-1.0_wp, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
!
      IF (N.GT.0) GO TO 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.0
!
      N = -N
      IF (X .EQ. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_gamma', 'X IS 0', 4, 2)
      IF (X .LT. 0.0 .AND. X+N-2 .EQ. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_gamma', 'X IS A NEGATIVE INTEGER', 4, 2)
      IF (X .LT. (-0.50_wp) .AND. ABS((X-AINT(X-0.50_wp))/X) .LT. DXREL) &
        CALL XERMSG ('SLATEC', 'wp_gamma', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      DO 20 I=1,N
        wp_gamma = wp_gamma/(X+I-1 )
 20   CONTINUE
      RETURN
!
! GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
!
 30   DO 40 I=1,N
        wp_gamma = (Y+I) * wp_gamma
 40   CONTINUE
      RETURN
!
! GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
 50   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'wp_gamma', 'X SO BIG GAMMA OVERFLOWS', 3, 2)
!
      wp_gamma = 0.0_wp
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'wp_gamma', 'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
      IF (X.LT.XMIN) RETURN
!
      wp_gamma = EXP ((Y-0.50_wp)*LOG(Y) - Y + SQ2PIL + cfp_9lgmc(Y) )
      IF (X.GT.0.0_wp) RETURN
!
      IF (ABS((X-AINT(X-0.50_wp))/X) .LT. DXREL) &
        CALL XERMSG ('SLATEC', 'wp_gamma', 'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      SINPIY = SIN (PI*Y)
      IF (SINPIY .EQ. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_gamma', 'X IS A NEGATIVE INTEGER', 4, 2)
!
      wp_gamma = -PI/(Y*SINPIY*wp_gamma)
!
      RETURN
      END FUNCTION
!
      !> Quad precision version of wp_gamma.
      REAL(kind=ep1) FUNCTION ep_gamma (X)
      IMPLICIT NONE
!
      REAL(kind=ep1) X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX, XMIN, Y 
      LOGICAL FIRST
!
      INTEGER I, N
      SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
      DATA GAMCS(  1) / +.8571195590989331421920062399942E-2_ep1      /
      DATA GAMCS(  2) / +.4415381324841006757191315771652E-2_ep1      /
      DATA GAMCS(  3) / +.5685043681599363378632664588789E-1_ep1      /
      DATA GAMCS(  4) / -.4219835396418560501012500186624E-2_ep1      /
      DATA GAMCS(  5) / +.1326808181212460220584006796352E-2_ep1      /
      DATA GAMCS(  6) / -.1893024529798880432523947023886E-3_ep1      /
      DATA GAMCS(  7) / +.3606925327441245256578082217225E-4_ep1      /
      DATA GAMCS(  8) / -.6056761904460864218485548290365E-5_ep1      /
      DATA GAMCS(  9) / +.1055829546302283344731823509093E-5_ep1      /
      DATA GAMCS( 10) / -.1811967365542384048291855891166E-6_ep1      /
      DATA GAMCS( 11) / +.3117724964715322277790254593169E-7_ep1      /
      DATA GAMCS( 12) / -.5354219639019687140874081024347E-8_ep1      /
      DATA GAMCS( 13) / +.9193275519859588946887786825940E-9_ep1      /
      DATA GAMCS( 14) / -.1577941280288339761767423273953E-9_ep1      /
      DATA GAMCS( 15) / +.2707980622934954543266540433089E-10_ep1     /
      DATA GAMCS( 16) / -.4646818653825730144081661058933E-11_ep1     /
      DATA GAMCS( 17) / +.7973350192007419656460767175359E-12_ep1     /
      DATA GAMCS( 18) / -.1368078209830916025799499172309E-12_ep1     /
      DATA GAMCS( 19) / +.2347319486563800657233471771688E-13_ep1     /
      DATA GAMCS( 20) / -.4027432614949066932766570534699E-14_ep1     /
      DATA GAMCS( 21) / +.6910051747372100912138336975257E-15_ep1     /
      DATA GAMCS( 22) / -.1185584500221992907052387126192E-15_ep1     /
      DATA GAMCS( 23) / +.2034148542496373955201026051932E-16_ep1     /
      DATA GAMCS( 24) / -.3490054341717405849274012949108E-17_ep1     /
      DATA GAMCS( 25) / +.5987993856485305567135051066026E-18_ep1     /
      DATA GAMCS( 26) / -.1027378057872228074490069778431E-18_ep1     /
      DATA GAMCS( 27) / +.1762702816060529824942759660748E-19_ep1     /
      DATA GAMCS( 28) / -.3024320653735306260958772112042E-20_ep1     /
      DATA GAMCS( 29) / +.5188914660218397839717833550506E-21_ep1     /
      DATA GAMCS( 30) / -.8902770842456576692449251601066E-22_ep1     /
      DATA GAMCS( 31) / +.1527474068493342602274596891306E-22_ep1     /
      DATA GAMCS( 32) / -.2620731256187362900257328332799E-23_ep1     /
      DATA GAMCS( 33) / +.4496464047830538670331046570666E-24_ep1     /
      DATA GAMCS( 34) / -.7714712731336877911703901525333E-25_ep1     /
      DATA GAMCS( 35) / +.1323635453126044036486572714666E-25_ep1     /
      DATA GAMCS( 36) / -.2270999412942928816702313813333E-26_ep1     /
      DATA GAMCS( 37) / +.3896418998003991449320816639999E-27_ep1     /
      DATA GAMCS( 38) / -.6685198115125953327792127999999E-28_ep1     /
      DATA GAMCS( 39) / +.1146998663140024384347613866666E-28_ep1     /
      DATA GAMCS( 40) / -.1967938586345134677295103999999E-29_ep1     /
      DATA GAMCS( 41) / +.3376448816585338090334890666666E-30_ep1     /
      DATA GAMCS( 42) / -.5793070335782135784625493333333E-31_ep1     /
      DATA PI / 3.141592653589793238462643383279500_ep1 /
      DATA SQ2PIL / 0.918938533204672741780329736405620_ep1 /
      DATA FIRST /.TRUE./
      INTEGER NGAM
!***FIRST EXECUTABLE STATEMENT  ep_gamma
      IF (FIRST) THEN
         NGAM = cfp_initds (GAMCS, 42, 0.1*REAL(F1MACH(3,ep_dummy)) )
!
         CALL cfp_gamlm (XMIN, XMAX)
         DXREL = SQRT(F1MACH(4,ep_dummy))
      ENDIF
      FIRST = .FALSE.
!
      Y = ABS(X)
      IF (Y.GT.10.0_ep1) GO TO 50
!
! COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
!
      N = X
      IF (X.LT.0.0_ep1) N = N - 1
      Y = X - N
      N = N - 1
      ep_gamma = 0.93750_ep1 + cfp_csevl (2.0_ep1*Y-1.0_ep1, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
!
      IF (N.GT.0) GO TO 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.0
!
      N = -N
      IF (X .EQ. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_gamma', 'X IS 0', 4, 2)
      IF (X .LT. 0.0 .AND. X+N-2 .EQ. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_gamma', 'X IS A NEGATIVE INTEGER', 4, 2)
      IF (X .LT. (-0.50_ep1) .AND. ABS((X-AINT(X-0.50_ep1))/X) .LT. DXREL) &
        CALL XERMSG ('SLATEC', 'ep_gamma', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      DO 20 I=1,N
        ep_gamma = ep_gamma/(X+I-1 )
 20   CONTINUE
      RETURN
!
! GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
!
 30   DO 40 I=1,N
        ep_gamma = (Y+I) * ep_gamma
 40   CONTINUE
      RETURN
!
! GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
 50   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'ep_gamma', 'X SO BIG GAMMA OVERFLOWS', 3, 2)
!
      ep_gamma = 0.0_ep1
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'ep_gamma', 'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
      IF (X.LT.XMIN) RETURN
!
      ep_gamma = EXP ((Y-0.50_ep1)*LOG(Y) - Y + SQ2PIL + cfp_9lgmc(Y) )
      IF (X.GT.0.0_ep1) RETURN
!
      IF (ABS((X-AINT(X-0.50_ep1))/X) .LT. DXREL) &
        CALL XERMSG ('SLATEC', 'ep_gamma', 'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      SINPIY = SIN (PI*Y)
      IF (SINPIY .EQ. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_gamma', 'X IS A NEGATIVE INTEGER', 4, 2)
!
      ep_gamma = -PI/(Y*SINPIY*ep_gamma)
!
      RETURN
      END FUNCTION
!> \verbatim
!>***SUBSIDIARY
!>***PURPOSE  Subsidiary to cfp_besj and DBESY
!>***LIBRARY   SLATEC
!>***TYPE      REAL(kind=wp) (JAIRY-S, wp_jairy-D)
!>***AUTHOR  Amos, D. E., (SNLA)
!>           Daniel, S. L., (SNLA)
!>           Weston, M. K., (SNLA)
!>***DESCRIPTION
!>
!>                  wp_jairy computes the Airy function AI(X)
!>                   and its derivative DAI(X) for wp_asyjy
!>
!>                                   INPUT
!>
!>         X - Argument, computed by wp_asyjy, X unrestricted
!>        RX - RX=SQRT(ABS(X)), computed by wp_asyjy
!>         C - C=2.*(ABS(X)**1.5)/3., computed by wp_asyjy
!>
!>                                  OUTPUT
!>
!>        AI - Value of function AI(X)
!>       DAI - Value of the derivative DAI(X)
!>
!>***SEE ALSO  cfp_besj, DBESY
!>***ROUTINES CALLED  (NONE)
!>***REVISION HISTORY  (YYMMDD)
!>   750101  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   891009  Removed unreferenced variable.  (WRB)
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900328  Added TYPE section.  (WRB)
!>   910408  Updated the AUTHOR section.  (WRB)
!> \endverbatim
!
      SUBROUTINE wp_jairy (X, RX, C, AI, DAI)
      IMPLICIT NONE
!
!
      INTEGER I, J, M1, M1D, M2, M2D, M3, M3D, M4, M4D, N1, N1D, N2, N2D, N3, N3D, N4, N4D
      REAL(kind=wp) A,AI,AJN,AJP,AK1,AK2,AK3,B,C,CCV,CON2, CON3, CON4, CON5, CV, DA, DAI, DAJN, DAJP, DAK1
      REAL(kind=wp) DAK2, DAK3, DB, EC, E1, E2, FPI12, F1, F2, RTRX, RX, SCV, T, TEMP1, TEMP2, TT, X
      DIMENSION AJP(19), AJN(19), A(15), B(15)
      DIMENSION AK1(14), AK2(23), AK3(14)
      DIMENSION DAJP(19), DAJN(19), DA(15), DB(15)
      DIMENSION DAK1(14), DAK2(24), DAK3(14)
      SAVE N1, N2, N3, N4, M1, M2, M3, M4, FPI12, CON2, CON3, CON4, CON5, AK1, AK2, AK3, AJP, AJN, A, B, N1D, N2D, N3D, &
           N4D, M1D, M2D, M3D, M4D, DAK1, DAK2, DAK3, DAJP, DAJN, DA, DB
      DATA N1,N2,N3,N4/14,23,19,15/
      DATA M1,M2,M3,M4/12,21,17,13/
      DATA FPI12,CON2,CON3,CON4,CON5/&
     & 1.30899693899575E+00_wp, 5.03154716196777E+00_wp, 3.80004589867293E-01_wp,&
     & 8.33333333333333E-01_wp, 8.66025403784439E-01_wp/
      DATA AK1(1), AK1(2), AK1(3), AK1(4), AK1(5), AK1(6), AK1(7),&
     &     AK1(8), AK1(9), AK1(10),AK1(11),AK1(12),AK1(13),&
     &     AK1(14)         / 2.20423090987793E-01_wp,-1.25290242787700E-01_wp,&
     & 1.03881163359194E-02_wp, 8.22844152006343E-04_wp,-2.34614345891226E-04_wp,&
     & 1.63824280172116E-05_wp, 3.06902589573189E-07_wp,-1.29621999359332E-07_wp,&
     & 8.22908158823668E-09_wp, 1.53963968623298E-11_wp,-3.39165465615682E-11_wp,&
     & 2.03253257423626E-12_wp,-1.10679546097884E-14_wp,-5.16169497785080E-15_wp/
      DATA AK2(1), AK2(2), AK2(3), AK2(4), AK2(5), AK2(6), AK2(7),&
     &     AK2(8), AK2(9), AK2(10),AK2(11),AK2(12),AK2(13),AK2(14),&
     &     AK2(15),AK2(16),AK2(17),AK2(18),AK2(19),AK2(20),AK2(21),&
     &     AK2(22),AK2(23) / 2.74366150869598E-01_wp, 5.39790969736903E-03_wp,&
     &-1.57339220621190E-03_wp, 4.27427528248750E-04_wp,-1.12124917399925E-04_wp,&
     & 2.88763171318904E-05_wp,-7.36804225370554E-06_wp, 1.87290209741024E-06_wp,&
     &-4.75892793962291E-07_wp, 1.21130416955909E-07_wp,-3.09245374270614E-08_wp,&
     & 7.92454705282654E-09_wp,-2.03902447167914E-09_wp, 5.26863056595742E-10_wp,&
     &-1.36704767639569E-10_wp, 3.56141039013708E-11_wp,-9.31388296548430E-12_wp,&
     & 2.44464450473635E-12_wp,-6.43840261990955E-13_wp, 1.70106030559349E-13_wp,&
     &-4.50760104503281E-14_wp, 1.19774799164811E-14_wp,-3.19077040865066E-15_wp/
      DATA AK3(1), AK3(2), AK3(3), AK3(4), AK3(5), AK3(6), AK3(7),&
     &     AK3(8), AK3(9), AK3(10),AK3(11),AK3(12),AK3(13),&
     &     AK3(14)         / 2.80271447340791E-01_wp,-1.78127042844379E-03_wp,&
     & 4.03422579628999E-05_wp,-1.63249965269003E-06_wp, 9.21181482476768E-08_wp,&
     &-6.52294330229155E-09_wp, 5.47138404576546E-10_wp,-5.24408251800260E-11_wp,&
     & 5.60477904117209E-12_wp,-6.56375244639313E-13_wp, 8.31285761966247E-14_wp,&
     &-1.12705134691063E-14_wp, 1.62267976598129E-15_wp,-2.46480324312426E-16_wp/
      DATA AJP(1), AJP(2), AJP(3), AJP(4), AJP(5), AJP(6), AJP(7),&
     &     AJP(8), AJP(9), AJP(10),AJP(11),AJP(12),AJP(13),AJP(14),&
     &     AJP(15),AJP(16),AJP(17),AJP(18),&
     &     AJP(19)         / 7.78952966437581E-02_wp,-1.84356363456801E-01_wp,&
     & 3.01412605216174E-02_wp, 3.05342724277608E-02_wp,-4.95424702513079E-03_wp,&
     &-1.72749552563952E-03_wp, 2.43137637839190E-04_wp, 5.04564777517082E-05_wp,&
     &-6.16316582695208E-06_wp,-9.03986745510768E-07_wp, 9.70243778355884E-08_wp,&
     & 1.09639453305205E-08_wp,-1.04716330588766E-09_wp,-9.60359441344646E-11_wp,&
     & 8.25358789454134E-12_wp, 6.36123439018768E-13_wp,-4.96629614116015E-14_wp,&
     &-3.29810288929615E-15_wp, 2.35798252031104E-16_wp/
      DATA AJN(1), AJN(2), AJN(3), AJN(4), AJN(5), AJN(6), AJN(7),&
     &     AJN(8), AJN(9), AJN(10),AJN(11),AJN(12),AJN(13),AJN(14),&
     &     AJN(15),AJN(16),AJN(17),AJN(18),&
     &     AJN(19)         / 3.80497887617242E-02_wp,-2.45319541845546E-01_wp,&
     & 1.65820623702696E-01_wp, 7.49330045818789E-02_wp,-2.63476288106641E-02_wp,&
     &-5.92535597304981E-03_wp, 1.44744409589804E-03_wp, 2.18311831322215E-04_wp,&
     &-4.10662077680304E-05_wp,-4.66874994171766E-06_wp, 7.15218807277160E-07_wp,&
     & 6.52964770854633E-08_wp,-8.44284027565946E-09_wp,-6.44186158976978E-10_wp,&
     & 7.20802286505285E-11_wp, 4.72465431717846E-12_wp,-4.66022632547045E-13_wp,&
     &-2.67762710389189E-14_wp, 2.36161316570019E-15_wp/
      DATA A(1),   A(2),   A(3),   A(4),   A(5),   A(6),   A(7),&
     &     A(8),   A(9),   A(10),  A(11),  A(12),  A(13),  A(14),&
     &     A(15)           / 4.90275424742791E-01_wp, 1.57647277946204E-03_wp,&
     &-9.66195963140306E-05_wp, 1.35916080268815E-07_wp, 2.98157342654859E-07_wp,&
     &-1.86824767559979E-08_wp,-1.03685737667141E-09_wp, 3.28660818434328E-10_wp,&
     &-2.57091410632780E-11_wp,-2.32357655300677E-12_wp, 9.57523279048255E-13_wp,&
     &-1.20340828049719E-13_wp,-2.90907716770715E-15_wp, 4.55656454580149E-15_wp,&
     &-9.99003874810259E-16_wp/
      DATA B(1),   B(2),   B(3),   B(4),   B(5),   B(6),   B(7),&
     &     B(8),   B(9),   B(10),  B(11),  B(12),  B(13),  B(14),&
     &     B(15)           / 2.78593552803079E-01_wp,-3.52915691882584E-03_wp,&
     &-2.31149677384994E-05_wp, 4.71317842263560E-06_wp,-1.12415907931333E-07_wp,&
     &-2.00100301184339E-08_wp, 2.60948075302193E-09_wp,-3.55098136101216E-11_wp,&
     &-3.50849978423875E-11_wp, 5.83007187954202E-12_wp,-2.04644828753326E-13_wp,&
     &-1.10529179476742E-13_wp, 2.87724778038775E-14_wp,-2.88205111009939E-15_wp,&
     &-3.32656311696166E-16_wp/
      DATA N1D,N2D,N3D,N4D/14,24,19,15/
      DATA M1D,M2D,M3D,M4D/12,22,17,13/
      DATA DAK1(1), DAK1(2), DAK1(3), DAK1(4), DAK1(5), DAK1(6),&
     &     DAK1(7), DAK1(8), DAK1(9), DAK1(10),DAK1(11),DAK1(12),&
     &    DAK1(13),DAK1(14)/ 2.04567842307887E-01_wp,-6.61322739905664E-02_wp,&
     &-8.49845800989287E-03_wp, 3.12183491556289E-03_wp,-2.70016489829432E-04_wp,&
     &-6.35636298679387E-06_wp, 3.02397712409509E-06_wp,-2.18311195330088E-07_wp,&
     &-5.36194289332826E-10_wp, 1.13098035622310E-09_wp,-7.43023834629073E-11_wp,&
     & 4.28804170826891E-13_wp, 2.23810925754539E-13_wp,-1.39140135641182E-14_wp/
      DATA DAK2(1), DAK2(2), DAK2(3), DAK2(4), DAK2(5), DAK2(6),&
     &     DAK2(7), DAK2(8), DAK2(9), DAK2(10),DAK2(11),DAK2(12),&
     &     DAK2(13),DAK2(14),DAK2(15),DAK2(16),DAK2(17),DAK2(18),&
     &     DAK2(19),DAK2(20),DAK2(21),DAK2(22),DAK2(23),&
     &     DAK2(24)        / 2.93332343883230E-01_wp,-8.06196784743112E-03_wp,&
     & 2.42540172333140E-03_wp,-6.82297548850235E-04_wp, 1.85786427751181E-04_wp,&
     &-4.97457447684059E-05_wp, 1.32090681239497E-05_wp,-3.49528240444943E-06_wp,&
     & 9.24362451078835E-07_wp,-2.44732671521867E-07_wp, 6.49307837648910E-08_wp,&
     &-1.72717621501538E-08_wp, 4.60725763604656E-09_wp,-1.23249055291550E-09_wp,&
     & 3.30620409488102E-10_wp,-8.89252099772401E-11_wp, 2.39773319878298E-11_wp,&
     &-6.48013921153450E-12_wp, 1.75510132023731E-12_wp,-4.76303829833637E-13_wp,&
     & 1.29498241100810E-13_wp,-3.52679622210430E-14_wp, 9.62005151585923E-15_wp,&
     &-2.62786914342292E-15_wp/
      DATA DAK3(1), DAK3(2), DAK3(3), DAK3(4), DAK3(5), DAK3(6),&
     &     DAK3(7), DAK3(8), DAK3(9), DAK3(10),DAK3(11),DAK3(12),&
     &    DAK3(13),DAK3(14)/ 2.84675828811349E-01_wp, 2.53073072619080E-03_wp,&
     &-4.83481130337976E-05_wp, 1.84907283946343E-06_wp,-1.01418491178576E-07_wp,&
     & 7.05925634457153E-09_wp,-5.85325291400382E-10_wp, 5.56357688831339E-11_wp,&
     &-5.90889094779500E-12_wp, 6.88574353784436E-13_wp,-8.68588256452194E-14_wp,&
     & 1.17374762617213E-14_wp,-1.68523146510923E-15_wp, 2.55374773097056E-16_wp/
      DATA DAJP(1), DAJP(2), DAJP(3), DAJP(4), DAJP(5), DAJP(6),&
     &     DAJP(7), DAJP(8), DAJP(9), DAJP(10),DAJP(11),DAJP(12),&
     &     DAJP(13),DAJP(14),DAJP(15),DAJP(16),DAJP(17),DAJP(18),&
     &     DAJP(19)        / 6.53219131311457E-02_wp,-1.20262933688823E-01_wp,&
     & 9.78010236263823E-03_wp, 1.67948429230505E-02_wp,-1.97146140182132E-03_wp,&
     &-8.45560295098867E-04_wp, 9.42889620701976E-05_wp, 2.25827860945475E-05_wp,&
     &-2.29067870915987E-06_wp,-3.76343991136919E-07_wp, 3.45663933559565E-08_wp,&
     & 4.29611332003007E-09_wp,-3.58673691214989E-10_wp,-3.57245881361895E-11_wp,&
     & 2.72696091066336E-12_wp, 2.26120653095771E-13_wp,-1.58763205238303E-14_wp,&
     &-1.12604374485125E-15_wp, 7.31327529515367E-17_wp/
      DATA DAJN(1), DAJN(2), DAJN(3), DAJN(4), DAJN(5), DAJN(6),&
     &     DAJN(7), DAJN(8), DAJN(9), DAJN(10),DAJN(11),DAJN(12),&
     &     DAJN(13),DAJN(14),DAJN(15),DAJN(16),DAJN(17),DAJN(18),&
     &     DAJN(19)        / 1.08594539632967E-02_wp, 8.53313194857091E-02_wp,&
     &-3.15277068113058E-01_wp,-8.78420725294257E-02_wp, 5.53251906976048E-02_wp,&
     & 9.41674060503241E-03_wp,-3.32187026018996E-03_wp,-4.11157343156826E-04_wp,&
     & 1.01297326891346E-04_wp, 9.87633682208396E-06_wp,-1.87312969812393E-06_wp,&
     &-1.50798500131468E-07_wp, 2.32687669525394E-08_wp, 1.59599917419225E-09_wp,&
     &-2.07665922668385E-10_wp,-1.24103350500302E-11_wp, 1.39631765331043E-12_wp,&
     & 7.39400971155740E-14_wp,-7.32887475627500E-15_wp/
      DATA DA(1),  DA(2),  DA(3),  DA(4),  DA(5),  DA(6),  DA(7),&
     &     DA(8),  DA(9),  DA(10), DA(11), DA(12), DA(13), DA(14),&
     &     DA(15)          / 4.91627321104601E-01_wp, 3.11164930427489E-03_wp,&
     & 8.23140762854081E-05_wp,-4.61769776172142E-06_wp,-6.13158880534626E-08_wp,&
     & 2.87295804656520E-08_wp,-1.81959715372117E-09_wp,-1.44752826642035E-10_wp,&
     & 4.53724043420422E-11_wp,-3.99655065847223E-12_wp,-3.24089119830323E-13_wp,&
     & 1.62098952568741E-13_wp,-2.40765247974057E-14_wp, 1.69384811284491E-16_wp,&
     & 8.17900786477396E-16_wp/
      DATA DB(1),  DB(2),  DB(3),  DB(4),  DB(5),  DB(6),  DB(7),&
     &     DB(8),  DB(9),  DB(10), DB(11), DB(12), DB(13), DB(14),&
     &     DB(15)          /-2.77571356944231E-01_wp, 4.44212833419920E-03_wp,&
     &-8.42328522190089E-05_wp,-2.58040318418710E-06_wp, 3.42389720217621E-07_wp,&
     &-6.24286894709776E-09_wp,-2.36377836844577E-09_wp, 3.16991042656673E-10_wp,&
     &-4.40995691658191E-12_wp,-5.18674221093575E-12_wp, 9.64874015137022E-13_wp,&
     &-4.90190576608710E-14_wp,-1.77253430678112E-14_wp, 5.55950610442662E-15_wp,&
     &-7.11793337579530E-16_wp/
!***FIRST EXECUTABLE STATEMENT  wp_jairy
      IF (X.LT.0.00_wp) GO TO 90
      IF (C.GT.5.00_wp) GO TO 60
      IF (X.GT.1.200_wp) GO TO 30
      T = (X+X-1.20_wp)*CON4
      TT = T + T
      J = N1
      F1 = AK1(J)
      F2 = 0.00_wp
      DO 10 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK1(J)
        F2 = TEMP1
   10 CONTINUE
      AI = T*F1 - F2 + AK1(1)
!
      J = N1D
      F1 = DAK1(J)
      F2 = 0.00_wp
      DO 20 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK1(J)
        F2 = TEMP1
   20 CONTINUE
      DAI = -(T*F1-F2+DAK1(1))
      RETURN
!
   30 CONTINUE
      T = (X+X-CON2)*CON3
      TT = T + T
      J = N2
      F1 = AK2(J)
      F2 = 0.00_wp
      DO 40 I=1,M2
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK2(J)
        F2 = TEMP1
   40 CONTINUE
      RTRX = SQRT(RX)
      EC = EXP(-C)
      AI = EC*(T*F1-F2+AK2(1))/RTRX
      J = N2D
      F1 = DAK2(J)
      F2 = 0.00_wp
      DO 50 I=1,M2D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK2(J)
        F2 = TEMP1
   50 CONTINUE
      DAI = -EC*(T*F1-F2+DAK2(1))*RTRX
      RETURN
!
   60 CONTINUE
      T = 10.00_wp/C - 1.00_wp
      TT = T + T
      J = N1
      F1 = AK3(J)
      F2 = 0.00_wp
      DO 70 I=1,M1
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + AK3(J)
        F2 = TEMP1
   70 CONTINUE
      RTRX = SQRT(RX)
      EC = EXP(-C)
      AI = EC*(T*F1-F2+AK3(1))/RTRX
      J = N1D
      F1 = DAK3(J)
      F2 = 0.00_wp
      DO 80 I=1,M1D
        J = J - 1
        TEMP1 = F1
        F1 = TT*F1 - F2 + DAK3(J)
        F2 = TEMP1
   80 CONTINUE
      DAI = -RTRX*EC*(T*F1-F2+DAK3(1))
      RETURN
!
   90 CONTINUE
      IF (C.GT.5.00_wp) GO TO 120
      T = 0.40_wp*C - 1.00_wp
      TT = T + T
      J = N3
      F1 = AJP(J)
      E1 = AJN(J)
      F2 = 0.00_wp
      E2 = 0.00_wp
      DO 100 I=1,M3
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + AJP(J)
        E1 = TT*E1 - E2 + AJN(J)
        F2 = TEMP1
        E2 = TEMP2
  100 CONTINUE
      AI = (T*E1-E2+AJN(1)) - X*(T*F1-F2+AJP(1))
      J = N3D
      F1 = DAJP(J)
      E1 = DAJN(J)
      F2 = 0.00_wp
      E2 = 0.00_wp
      DO 110 I=1,M3D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DAJP(J)
        E1 = TT*E1 - E2 + DAJN(J)
        F2 = TEMP1
        E2 = TEMP2
  110 CONTINUE
      DAI = X*X*(T*F1-F2+DAJP(1)) + (T*E1-E2+DAJN(1))
      RETURN
!
  120 CONTINUE
      T = 10.00_wp/C - 1.00_wp
      TT = T + T
      J = N4
      F1 = A(J)
      E1 = B(J)
      F2 = 0.00_wp
      E2 = 0.00_wp
      DO 130 I=1,M4
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + A(J)
        E1 = TT*E1 - E2 + B(J)
        F2 = TEMP1
        E2 = TEMP2
  130 CONTINUE
      TEMP1 = T*F1 - F2 + A(1)
      TEMP2 = T*E1 - E2 + B(1)
      RTRX = SQRT(RX)
      CV = C - FPI12
      CCV = COS(CV)
      SCV = SIN(CV)
      AI = (TEMP1*CCV-TEMP2*SCV)/RTRX
      J = N4D
      F1 = DA(J)
      E1 = DB(J)
      F2 = 0.00_wp
      E2 = 0.00_wp
      DO 140 I=1,M4D
        J = J - 1
        TEMP1 = F1
        TEMP2 = E1
        F1 = TT*F1 - F2 + DA(J)
        E1 = TT*E1 - E2 + DB(J)
        F2 = TEMP1
        E2 = TEMP2
  140 CONTINUE
      TEMP1 = T*F1 - F2 + DA(1)
      TEMP2 = T*E1 - E2 + DB(1)
      E1 = CCV*CON5 + 0.50_wp*SCV
      E2 = SCV*CON5 - 0.50_wp*CCV
      DAI = (TEMP1*E1-TEMP2*E2)*RTRX
      RETURN
      END SUBROUTINE
!
      !> Quad precision version of wp_jairy.
      !> \warning Not implemented yet.
      SUBROUTINE ep_jairy (X, RX, C, AI, DAI)
      IMPLICIT NONE
!
!
      REAL(kind=ep1) RX, X, C, AI, DAI

         CALL XERMSG('SLATEC','ep_jairy','Quad precision version not implemented yet.',1,1)

      END SUBROUTINE
!> \verbatim
!>***PURPOSE  Compute the logarithm of the absolute value of the Gamma
!>            function.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7A
!>***TYPE      REAL(kind=cfp) (ALNGAM-S, cfp_lngam-D, CLNGAM-C)
!>***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
!>             SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> cfp_lngam(X) calculates the REAL(kind=cfp) logarithm of the
!> absolute value of the Gamma function for REAL(kind=cfp)
!> argument X.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, cfp_9lgmc, cfp_gamma, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770601  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900727  Added EXTERNAL statement.  (WRB)
!> \endverbatim
!
      REAL(kind=wp) FUNCTION wp_lngam (X)
      IMPLICIT NONE
!
      REAL(kind=wp) X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX, Y, TEMP
      LOGICAL FIRST
      !EXTERNAL cfp_gamma_fun
      SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
      DATA SQ2PIL / 0.918938533204672741780329736405620_wp /
      DATA SQPI2L / 0.225791352644727432363097614947441_wp    /
      DATA PI / 3.141592653589793238462643383279500_wp /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_lngam
      IF (FIRST) THEN
         TEMP = 1.0_wp/LOG(F1MACH(2,wp_dummy))
         XMAX = TEMP*F1MACH(2,wp_dummy)
         DXREL = SQRT(F1MACH(4,wp_dummy))
      ENDIF
      FIRST = .FALSE.
!
      Y = ABS (X)
      IF (Y.GT.10.0_wp) GO TO 20
!
! LOG (ABS (wp_gamma(X)) ) FOR ABS(X) .LE. 10.0
!
      wp_lngam = LOG (ABS (cfp_gamma_fun(X)) )
      RETURN
!
! LOG ( ABS (wp_gamma(X)) ) FOR ABS(X) .GT. 10.0
!
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'wp_lngam', 'ABS(X) SO BIG wp_lngam OVERFLOWS', 2, 2)
!
      IF (X.GT.0.0_wp) wp_lngam = SQ2PIL + (X-0.50_wp)*LOG(X) - X + cfp_9lgmc(Y)
      IF (X.GT.0.0_wp) RETURN
!
      SINPIY = ABS (SIN(PI*Y))
      IF (SINPIY .EQ. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_lngam', 'X IS A NEGATIVE INTEGER', 3, 2)
!
      IF (ABS((X-AINT(X-0.50_wp))/X) .LT. DXREL) &
        CALL XERMSG ('SLATEC', 'wp_lngam', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      wp_lngam = SQPI2L + (X-0.50_wp)*LOG(Y) - X - LOG(SINPIY) - cfp_9lgmc(Y)
      RETURN
!
      END FUNCTION
!
      !> Quad precision version of wp_lngam.
      REAL(kind=ep1) FUNCTION ep_lngam (X)
      IMPLICIT NONE
!
      REAL(kind=ep1) X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX, Y, TEMP
      LOGICAL FIRST
      !EXTERNAL cfp_gamma_fun
      SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
      DATA SQ2PIL / 0.918938533204672741780329736405620_ep1 /
      DATA SQPI2L / 0.225791352644727432363097614947441_ep1    /
      DATA PI / 3.141592653589793238462643383279500_ep1 /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_lngam
      IF (FIRST) THEN
         TEMP = 1.0_ep1/LOG(F1MACH(2,ep_dummy))
         XMAX = TEMP*F1MACH(2,ep_dummy)
         DXREL = SQRT(F1MACH(4,ep_dummy))
      ENDIF
      FIRST = .FALSE.
!
      Y = ABS (X)
      IF (Y.GT.10.0_ep1) GO TO 20
!
! LOG (ABS (cfp_gamma_fun(X)) ) FOR ABS(X) .LE. 10.0
!
      ep_lngam = LOG (ABS (cfp_gamma_fun(X)) )
      RETURN
!
! LOG ( ABS (cfp_gamma_fun(X)) ) FOR ABS(X) .GT. 10.0
!
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'ep_lngam', 'ABS(X) SO BIG ep_lngam OVERFLOWS', 2, 2)
!
      IF (X.GT.0.0_ep1) ep_lngam = SQ2PIL + (X-0.50_ep1)*LOG(X) - X + cfp_9lgmc(Y)
      IF (X.GT.0.0_ep1) RETURN
!
      SINPIY = ABS (SIN(PI*Y))
      IF (SINPIY .EQ. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_lngam', 'X IS A NEGATIVE INTEGER', 3, 2)
!
      IF (ABS((X-AINT(X-0.50_ep1))/X) .LT. DXREL) &
        CALL XERMSG ('SLATEC', 'ep_lngam', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      ep_lngam = SQPI2L + (X-0.50_ep1)*LOG(Y) - X - LOG(SINPIY) - cfp_9lgmc(Y)
      RETURN
!
      END FUNCTION
!> \verbatim
!>***PURPOSE  Evaluate ln(1+X) accurate in the sense of relative error.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C4B
!>***TYPE      REAL(kind=wp) (ALNREL-S, wp_lnrel-D, CLNREL-C)
!>***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> wp_lnrel(X) calculates the REAL(kind=wp) natural logarithm of
!> (1.0+X) for REAL(kind=wp) argument X.  This routine should
!> be used when X is small and accurate to calculate the logarithm
!> accurately (in the relative error sense) in the neighborhood
!> of 1.0.
!>
!> Series for ALNR       on the interval -3.75000E-01_wp to  3.75000E-01_wp
!>                                        with weighted error   6.35E-32
!>                                         log weighted error  31.20
!>                               significant figures required  30.93
!>                                    decimal places required  32.01
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, cfp_csevl, cfp_initds, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770601  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!> \endverbatim
!
      REAL(kind=wp) FUNCTION wp_lnrel (X)
      IMPLICIT NONE
!
      REAL(kind=wp) ALNRCS(43), X, XMIN
      INTEGER NLNREL
      LOGICAL FIRST
      SAVE ALNRCS, NLNREL, XMIN, FIRST
      DATA ALNRCS(  1) / +.10378693562743769800686267719098E+1_wp     /
      DATA ALNRCS(  2) / -.13364301504908918098766041553133E+0_wp     /
      DATA ALNRCS(  3) / +.19408249135520563357926199374750E-1_wp     /
      DATA ALNRCS(  4) / -.30107551127535777690376537776592E-2_wp     /
      DATA ALNRCS(  5) / +.48694614797154850090456366509137E-3_wp     /
      DATA ALNRCS(  6) / -.81054881893175356066809943008622E-4_wp     /
      DATA ALNRCS(  7) / +.13778847799559524782938251496059E-4_wp     /
      DATA ALNRCS(  8) / -.23802210894358970251369992914935E-5_wp     /
      DATA ALNRCS(  9) / +.41640416213865183476391859901989E-6_wp     /
      DATA ALNRCS( 10) / -.73595828378075994984266837031998E-7_wp     /
      DATA ALNRCS( 11) / +.13117611876241674949152294345011E-7_wp     /
      DATA ALNRCS( 12) / -.23546709317742425136696092330175E-8_wp     /
      DATA ALNRCS( 13) / +.42522773276034997775638052962567E-9_wp     /
      DATA ALNRCS( 14) / -.77190894134840796826108107493300E-10_wp    /
      DATA ALNRCS( 15) / +.14075746481359069909215356472191E-10_wp    /
      DATA ALNRCS( 16) / -.25769072058024680627537078627584E-11_wp    /
      DATA ALNRCS( 17) / +.47342406666294421849154395005938E-12_wp    /
      DATA ALNRCS( 18) / -.87249012674742641745301263292675E-13_wp    /
      DATA ALNRCS( 19) / +.16124614902740551465739833119115E-13_wp    /
      DATA ALNRCS( 20) / -.29875652015665773006710792416815E-14_wp    /
      DATA ALNRCS( 21) / +.55480701209082887983041321697279E-15_wp    /
      DATA ALNRCS( 22) / -.10324619158271569595141333961932E-15_wp    /
      DATA ALNRCS( 23) / +.19250239203049851177878503244868E-16_wp    /
      DATA ALNRCS( 24) / -.35955073465265150011189707844266E-17_wp    /
      DATA ALNRCS( 25) / +.67264542537876857892194574226773E-18_wp    /
      DATA ALNRCS( 26) / -.12602624168735219252082425637546E-18_wp    /
      DATA ALNRCS( 27) / +.23644884408606210044916158955519E-19_wp    /
      DATA ALNRCS( 28) / -.44419377050807936898878389179733E-20_wp    /
      DATA ALNRCS( 29) / +.83546594464034259016241293994666E-21_wp    /
      DATA ALNRCS( 30) / -.15731559416479562574899253521066E-21_wp    /
      DATA ALNRCS( 31) / +.29653128740247422686154369706666E-22_wp    /
      DATA ALNRCS( 32) / -.55949583481815947292156013226666E-23_wp    /
      DATA ALNRCS( 33) / +.10566354268835681048187284138666E-23_wp    /
      DATA ALNRCS( 34) / -.19972483680670204548314999466666E-24_wp    /
      DATA ALNRCS( 35) / +.37782977818839361421049855999999E-25_wp    /
      DATA ALNRCS( 36) / -.71531586889081740345038165333333E-26_wp    /
      DATA ALNRCS( 37) / +.13552488463674213646502024533333E-26_wp    /
      DATA ALNRCS( 38) / -.25694673048487567430079829333333E-27_wp    /
      DATA ALNRCS( 39) / +.48747756066216949076459519999999E-28_wp    /
      DATA ALNRCS( 40) / -.92542112530849715321132373333333E-29_wp    /
      DATA ALNRCS( 41) / +.17578597841760239233269760000000E-29_wp    /
      DATA ALNRCS( 42) / -.33410026677731010351377066666666E-30_wp    /
      DATA ALNRCS( 43) / +.63533936180236187354180266666666E-31_wp    /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_lnrel
      IF (FIRST) THEN
         NLNREL = cfp_initds (ALNRCS, 43, 0.1*REAL(F1MACH(3,wp_dummy)))
         XMIN = -1.00_wp + SQRT(F1MACH(4,wp_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LE. (-1.0_wp)) CALL XERMSG ('SLATEC', 'wp_lnrel', 'X IS LE -1', 2, 2)
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'wp_lnrel', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 1, 1)
!
      IF (ABS(X).LE.0.3750_wp) wp_lnrel = X*(1.0_wp -X*cfp_csevl (X/.3750_wp, ALNRCS, NLNREL))
!
      IF (ABS(X).GT.0.3750_wp) wp_lnrel = LOG (1.00_wp+X)
!
      RETURN
      END FUNCTION
!
      !> Quad precision version of wp_lnrel.
      REAL(kind=ep1) FUNCTION ep_lnrel (X)
      IMPLICIT NONE
!
      REAL(kind=ep1) ALNRCS(43), X, XMIN
      INTEGER NLNREL
      LOGICAL FIRST
      SAVE ALNRCS, NLNREL, XMIN, FIRST
      DATA ALNRCS(  1) / +.10378693562743769800686267719098E+1_ep1     /
      DATA ALNRCS(  2) / -.13364301504908918098766041553133E+0_ep1     /
      DATA ALNRCS(  3) / +.19408249135520563357926199374750E-1_ep1     /
      DATA ALNRCS(  4) / -.30107551127535777690376537776592E-2_ep1     /
      DATA ALNRCS(  5) / +.48694614797154850090456366509137E-3_ep1     /
      DATA ALNRCS(  6) / -.81054881893175356066809943008622E-4_ep1     /
      DATA ALNRCS(  7) / +.13778847799559524782938251496059E-4_ep1     /
      DATA ALNRCS(  8) / -.23802210894358970251369992914935E-5_ep1     /
      DATA ALNRCS(  9) / +.41640416213865183476391859901989E-6_ep1     /
      DATA ALNRCS( 10) / -.73595828378075994984266837031998E-7_ep1     /
      DATA ALNRCS( 11) / +.13117611876241674949152294345011E-7_ep1     /
      DATA ALNRCS( 12) / -.23546709317742425136696092330175E-8_ep1     /
      DATA ALNRCS( 13) / +.42522773276034997775638052962567E-9_ep1     /
      DATA ALNRCS( 14) / -.77190894134840796826108107493300E-10_ep1    /
      DATA ALNRCS( 15) / +.14075746481359069909215356472191E-10_ep1    /
      DATA ALNRCS( 16) / -.25769072058024680627537078627584E-11_ep1    /
      DATA ALNRCS( 17) / +.47342406666294421849154395005938E-12_ep1    /
      DATA ALNRCS( 18) / -.87249012674742641745301263292675E-13_ep1    /
      DATA ALNRCS( 19) / +.16124614902740551465739833119115E-13_ep1    /
      DATA ALNRCS( 20) / -.29875652015665773006710792416815E-14_ep1    /
      DATA ALNRCS( 21) / +.55480701209082887983041321697279E-15_ep1    /
      DATA ALNRCS( 22) / -.10324619158271569595141333961932E-15_ep1    /
      DATA ALNRCS( 23) / +.19250239203049851177878503244868E-16_ep1    /
      DATA ALNRCS( 24) / -.35955073465265150011189707844266E-17_ep1    /
      DATA ALNRCS( 25) / +.67264542537876857892194574226773E-18_ep1    /
      DATA ALNRCS( 26) / -.12602624168735219252082425637546E-18_ep1    /
      DATA ALNRCS( 27) / +.23644884408606210044916158955519E-19_ep1    /
      DATA ALNRCS( 28) / -.44419377050807936898878389179733E-20_ep1    /
      DATA ALNRCS( 29) / +.83546594464034259016241293994666E-21_ep1    /
      DATA ALNRCS( 30) / -.15731559416479562574899253521066E-21_ep1    /
      DATA ALNRCS( 31) / +.29653128740247422686154369706666E-22_ep1    /
      DATA ALNRCS( 32) / -.55949583481815947292156013226666E-23_ep1    /
      DATA ALNRCS( 33) / +.10566354268835681048187284138666E-23_ep1    /
      DATA ALNRCS( 34) / -.19972483680670204548314999466666E-24_ep1    /
      DATA ALNRCS( 35) / +.37782977818839361421049855999999E-25_ep1    /
      DATA ALNRCS( 36) / -.71531586889081740345038165333333E-26_ep1    /
      DATA ALNRCS( 37) / +.13552488463674213646502024533333E-26_ep1    /
      DATA ALNRCS( 38) / -.25694673048487567430079829333333E-27_ep1    /
      DATA ALNRCS( 39) / +.48747756066216949076459519999999E-28_ep1    /
      DATA ALNRCS( 40) / -.92542112530849715321132373333333E-29_ep1    /
      DATA ALNRCS( 41) / +.17578597841760239233269760000000E-29_ep1    /
      DATA ALNRCS( 42) / -.33410026677731010351377066666666E-30_ep1    /
      DATA ALNRCS( 43) / +.63533936180236187354180266666666E-31_ep1    /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_lnrel
      IF (FIRST) THEN
         NLNREL = cfp_initds (ALNRCS, 43, 0.1*REAL(F1MACH(3,ep_dummy)))
         XMIN = -1.00_ep1 + SQRT(F1MACH(4,ep_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LE. (-1.0_ep1)) CALL XERMSG ('SLATEC', 'ep_lnrel', 'X IS LE -1', 2, 2)
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'ep_lnrel', 'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 1, 1)
!
      IF (ABS(X).LE.0.3750_ep1) ep_lnrel = X*(1.0_ep1 -X*cfp_csevl (X/.3750_ep1, ALNRCS, NLNREL))
!
      IF (ABS(X).GT.0.3750_ep1) ep_lnrel = LOG (1.00_ep1+X)
!
      RETURN
      END FUNCTION
!> \verbatim
!>***PURPOSE  Determine the number of terms needed in an orthogonal
!>            polynomial series so that it meets a specified accuracy.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C3A2
!>***TYPE      REAL(kind=wp) (INITS-S, wp_initds-D)
!>***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!>             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!>  Initialize the orthogonal series, represented by the array OS, so
!>  that wp_initds is the number of terms needed to insure the error is no
!>  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
!>  machine precision.
!>
!>             Input Arguments --
!>   OS     REAL(kind=wp) array of NOS coefficients in an orthogonal
!>          series.
!>   NOS    number of coefficients in OS.
!>   ETA    single precision scalar containing requested accuracy of
!>          series.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770601  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890831  Modified array declarations.  (WRB)
!>   891115  Modified error message.  (WRB)
!>   891115  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!> \endverbatim
!
      FUNCTION wp_initds (OS, NOS, ETA)
      IMPLICIT NONE
!
      INTEGER :: wp_initds, NOS
      INTEGER I, II
      REAL :: ETA, ERR
      REAL(kind=wp) OS(*)
!***FIRST EXECUTABLE STATEMENT  wp_initds
      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'wp_initds', 'Number of coefficients is less than 1', 2, 1)
!
      ERR = 0.0
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
!
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'wp_initds', 'Chebyshev series too short for specified accuracy', 1, 1)
      wp_initds = I
!
      RETURN
      END FUNCTION
!     
      !> Quad prec version of wp_initds.
      FUNCTION ep_initds (OS, NOS, ETA)
      IMPLICIT NONE
!
      INTEGER :: ep_initds, NOS
      INTEGER I, II
      REAL :: ETA
      REAL(kind=wp) :: ERR
      REAL(kind=ep1) OS(*)
!***FIRST EXECUTABLE STATEMENT  ep_initds
      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'ep_initds', 'Number of coefficients is less than 1', 2, 1)
!
      ERR = 0.0_wp
      DO II = 1,NOS
         I = NOS + 1 - II
         ERR = ERR + ABS(REAL(OS(I),kind=wp))
         IF (ERR.GT.ETA) GO TO 20
      ENDDO
!
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'ep_initds', 'Chebyshev series too short for specified accuracy', 1, 1)
      ep_initds = I
!
      RETURN
      END FUNCTION
!>DECK wp_9gmic
!>***BEGIN PROLOGUE  wp_9gmic
!>***SUBSIDIARY
!>***PURPOSE  Compute the complementary incomplete Gamma function for A
!>            near a negative integer and X small.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7E
!>***TYPE      REAL(kind=wp) (R9GMIC-S, wp_9gmic-D)
!>***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
!>             SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> Compute the complementary incomplete gamma function for A near
!> a negative integer and for small X.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, cfp_lngam, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770701  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   890911  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900720  Routine changed from user-callable to subsidiary.  (WRB)
!>***END PROLOGUE  wp_9gmic
      REAL(kind=wp) FUNCTION wp_9gmic (A, X, ALX)
      IMPLICIT NONE
      REAL(kind=wp) A, X, ALX, ALNG, BOT, EPS, EULER, FK, FKP1, FM, S, SGNG, T, TE
      INTEGER MM1, K, M
      LOGICAL FIRST
      SAVE EULER, EPS, BOT, FIRST
      DATA EULER / 0.577215664901532860606512090082400_wp /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_9gmic
      IF (FIRST) THEN
         EPS = 0.50_wp*F1MACH(3,wp_dummy)
         BOT = LOG (F1MACH(1,wp_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (A .GT. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_9gmic', 'A MUST BE NEAR A NEGATIVE INTEGER', 2, 2)
      IF (X .LE. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_9gmic', 'X MUST BE GT ZERO', 3, 2)
!
      M = -(A - 0.50_wp)
      FM = M
!
      TE = 1.00_wp
      T = 1.00_wp
      S = T
      DO 20 K=1,200
        FKP1 = K + 1
        TE = -X*TE/(FM+FKP1)
        T = TE/FKP1
        S = S + T
        IF (ABS(T).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'wp_9gmic', 'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 4, 2)
!
 30   wp_9gmic = -ALX - EULER + X*S/(FM+1.00_wp)
      IF (M.EQ.0) RETURN
!
      IF (M.EQ.1) wp_9gmic = -wp_9gmic - 1.0_wp + 1.0_wp/X
      IF (M.EQ.1) RETURN
!
      TE = FM
      T = 1.0_wp
      S = T
      MM1 = M - 1
      DO 40 K=1,MM1
        FK = K
        TE = -X*TE/FK
        T = TE/(FM-FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
!
 50   DO 60 K=1,M
        wp_9gmic = wp_9gmic + 1.00_wp/K
 60   CONTINUE
!
      SGNG = 1.00_wp
      IF (MOD(M,2).EQ.1) SGNG = -1.00_wp
      ALNG = LOG(wp_9gmic) - cfp_lngam(FM+1.0_wp)
!
      wp_9gmic = 0.0_wp
      IF (ALNG.GT.BOT) wp_9gmic = SGNG * EXP(ALNG)
      IF (S.NE.0.0_wp) wp_9gmic = wp_9gmic + SIGN (EXP(-FM*ALX+LOG(ABS(S)/FM)), S)
!
      IF (wp_9gmic .EQ. 0.0_wp .AND. S .EQ. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_9gmic', 'RESULT UNDERFLOWS', 1, 1)
      RETURN
!
      END FUNCTION wp_9gmic
!
      !> Quad precision version of wp_9gmic.
      REAL(kind=ep1) FUNCTION ep_9gmic (A, X, ALX)
      IMPLICIT NONE
      REAL(kind=ep1) A, X, ALX, ALNG, BOT, EPS, EULER, FK, FKP1, FM, S, SGNG, T, TE
      INTEGER MM1, K, M
      LOGICAL FIRST
      SAVE EULER, EPS, BOT, FIRST
      DATA EULER / 0.577215664901532860606512090082400_ep1 /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_9gmic
      IF (FIRST) THEN
         EPS = 0.50_ep1*F1MACH(3,ep_dummy)
         BOT = LOG (F1MACH(1,ep_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (A .GT. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_9gmic', 'A MUST BE NEAR A NEGATIVE INTEGER', 2, 2)
      IF (X .LE. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_9gmic', 'X MUST BE GT ZERO', 3, 2)
!
      M = -(A - 0.50_ep1)
      FM = M
!
      TE = 1.00_ep1
      T = 1.00_ep1
      S = T
      DO 20 K=1,200
        FKP1 = K + 1
        TE = -X*TE/(FM+FKP1)
        T = TE/FKP1
        S = S + T
        IF (ABS(T).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'ep_9gmic', 'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 4, 2)
!
 30   ep_9gmic = -ALX - EULER + X*S/(FM+1.00_ep1)
      IF (M.EQ.0) RETURN
!
      IF (M.EQ.1) ep_9gmic = -ep_9gmic - 1.0_ep1 + 1.0_ep1/X
      IF (M.EQ.1) RETURN
!
      TE = FM
      T = 1.0_ep1
      S = T
      MM1 = M - 1
      DO 40 K=1,MM1
        FK = K
        TE = -X*TE/FK
        T = TE/(FM-FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
!
 50   DO 60 K=1,M
        ep_9gmic = ep_9gmic + 1.00_ep1/K
 60   CONTINUE
!
      SGNG = 1.00_ep1
      IF (MOD(M,2).EQ.1) SGNG = -1.00_ep1
      ALNG = LOG(ep_9gmic) - cfp_lngam(FM+1.0_ep1)
!
      ep_9gmic = 0.0_ep1
      IF (ALNG.GT.BOT) ep_9gmic = SGNG * EXP(ALNG)
      IF (S.NE.0.0_ep1) ep_9gmic = ep_9gmic + SIGN (EXP(-FM*ALX+LOG(ABS(S)/FM)), S)
!
      IF (ep_9gmic .EQ. 0.0_ep1 .AND. S .EQ. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_9gmic', 'RESULT UNDERFLOWS', 1, 1)
      RETURN
!
      END FUNCTION ep_9gmic
!>DECK wp_9gmit
!>***BEGIN PROLOGUE  wp_9gmit
!>***SUBSIDIARY
!>***PURPOSE  Compute Tricomi's incomplete Gamma function for small
!>            arguments.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7E
!>***TYPE      REAL(kind=wp) (R9GMIT-S, wp_9gmit-D)
!>***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, SMALL X,
!>             SPECIAL FUNCTIONS, TRICOMI
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> Compute Tricomi's incomplete gamma function for small X.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, wp_lngam, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770701  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890911  Removed unnecessary intrinsics.  (WRB)
!>   890911  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900720  Routine changed from user-callable to subsidiary.  (WRB)
!>***END PROLOGUE  wp_9gmit
      REAL(kind=wp) FUNCTION wp_9gmit (A, X, ALGAP1, SGNGAM, ALX)
      IMPLICIT NONE
      REAL(kind=wp) A, X, ALGAP1, SGNGAM, ALX, AE, AEPS, ALGS, ALG2, BOT, EPS, FK, S, SGNG2, T, TE
      INTEGER M, K, MA
      LOGICAL FIRST
      SAVE EPS, BOT, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_9gmit
      IF (FIRST) THEN
         EPS = 0.50_wp*F1MACH(3,wp_dummy)
         BOT = LOG (F1MACH(1,wp_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LE. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_9gmit', 'X SHOULD BE GT 0', 1, 2)
!
      MA = A + 0.50_wp
      IF (A.LT.0.0_wp) MA = A - 0.50_wp
      AEPS = A - MA
!
      AE = A
      IF (A.LT.(-0.50_wp)) AE = AEPS
!
      T = 1.0_wp
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'wp_9gmit', 'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES', 2, 2)
!
 30   IF (A.GE.(-0.50_wp)) ALGS = -ALGAP1 + LOG(S)
      IF (A.GE.(-0.50_wp)) GO TO 60
!
      ALGS = -cfp_lngam(1.0_wp+AEPS) + LOG(S)
      S = 1.00_wp
      M = -MA - 1
      IF (M.EQ.0) GO TO 50
      T = 1.00_wp
      DO 40 K=1,M
        T = X*T/(AEPS-(M+1-K))
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
!
 50   wp_9gmit = 0.00_wp
      ALGS = -MA*LOG(X) + ALGS
      IF (S.EQ.0.0_wp .OR. AEPS.EQ.0.0_wp) GO TO 60
!
      SGNG2 = SGNGAM * SIGN (1.00_wp, S)
      ALG2 = -X - ALGAP1 + LOG(ABS(S))
!
      IF (ALG2.GT.BOT) wp_9gmit = SGNG2 * EXP(ALG2)
      IF (ALGS.GT.BOT) wp_9gmit = wp_9gmit + EXP(ALGS)
      RETURN
!
 60   wp_9gmit = EXP (ALGS)
      RETURN
!
      END FUNCTION
!
      !> Quad precision version of wp_9gmit.
      REAL(kind=ep1) FUNCTION ep_9gmit (A, X, ALGAP1, SGNGAM, ALX)
      IMPLICIT NONE
      REAL(kind=ep1) A, X, ALGAP1, SGNGAM, ALX, AE, AEPS, ALGS, ALG2, BOT, EPS, FK, S, SGNG2, T, TE
      INTEGER M, K, MA
      LOGICAL FIRST
      SAVE EPS, BOT, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_9gmit
      IF (FIRST) THEN
         EPS = 0.50_ep1*F1MACH(3,ep_dummy)
         BOT = LOG (F1MACH(1,ep_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LE. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_9gmit', 'X SHOULD BE GT 0', 1, 2)
!
      MA = A + 0.50_ep1
      IF (A.LT.0.0_ep1) MA = A - 0.50_ep1
      AEPS = A - MA
!
      AE = A
      IF (A.LT.(-0.50_ep1)) AE = AEPS
!
      T = 1.0_ep1
      TE = AE
      S = T
      DO 20 K=1,200
        FK = K
        TE = -X*TE/FK
        T = TE/(AE+FK)
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'ep_9gmit', 'NO CONVERGENCE IN 200 TERMS OF TAYLOR-S SERIES', 2, 2)
!
 30   IF (A.GE.(-0.50_ep1)) ALGS = -ALGAP1 + LOG(S)
      IF (A.GE.(-0.50_ep1)) GO TO 60
!
      ALGS = -cfp_lngam(1.0_ep1+AEPS) + LOG(S)
      S = 1.00_ep1
      M = -MA - 1
      IF (M.EQ.0) GO TO 50
      T = 1.00_ep1
      DO 40 K=1,M
        T = X*T/(AEPS-(M+1-K))
        S = S + T
        IF (ABS(T).LT.EPS*ABS(S)) GO TO 50
 40   CONTINUE
!
 50   ep_9gmit = 0.00_ep1
      ALGS = -MA*LOG(X) + ALGS
      IF (S.EQ.0.0_ep1 .OR. AEPS.EQ.0.0_ep1) GO TO 60
!
      SGNG2 = SGNGAM * SIGN (1.00_ep1, S)
      ALG2 = -X - ALGAP1 + LOG(ABS(S))
!
      IF (ALG2.GT.BOT) ep_9gmit = SGNG2 * EXP(ALG2)
      IF (ALGS.GT.BOT) ep_9gmit = ep_9gmit + EXP(ALGS)
      RETURN
!
 60   ep_9gmit = EXP (ALGS)
      RETURN
!
      END FUNCTION
!>DECK wp_9lgic
!>***BEGIN PROLOGUE  wp_9lgic
!>***SUBSIDIARY
!>***PURPOSE  Compute the log complementary incomplete Gamma function
!>            for large X and for A .LE. X.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7E
!>***TYPE      REAL(kind=wp) (R9LGIC-S, wp_9lgic-D)
!>***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB, LARGE X,
!>             LOGARITHM, SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> Compute the log complementary incomplete gamma function for large X
!> and for A .LE. X.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770701  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900720  Routine changed from user-callable to subsidiary.  (WRB)
!>***END PROLOGUE  wp_9lgic
      REAL(kind=wp) FUNCTION wp_9lgic (A, X, ALX)
      IMPLICIT NONE
      REAL(kind=wp) A, X, ALX, EPS, FK, P, R, S, T, XMA, XPA
      INTEGER K
      SAVE EPS
      DATA EPS / 0.0_wp /
!***FIRST EXECUTABLE STATEMENT  wp_9lgic
      IF (EPS.EQ.0.0_wp) EPS = 0.50_wp*F1MACH(3,wp_dummy)
!
      XPA = X + 1.00_wp - A
      XMA = X - 1.0_wp - A
!
      R = 0.0_wp
      P = 1.0_wp
      S = P
      DO 10 K=1,300
        FK = K
        T = FK*(A-FK)*(1.0_wp+R)
        R = -T/((XMA+2.0_wp*FK)*(XPA+2.0_wp*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'wp_9lgic', 'NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION', 1, 2)
!
 20   wp_9lgic = A*ALX - X + LOG(S/XPA)
!
      RETURN
      END FUNCTION
!
      !> Quad precision version of wp_9lgic.
      REAL(kind=ep1) FUNCTION ep_9lgic (A, X, ALX)
      IMPLICIT NONE
      REAL(kind=ep1) A, X, ALX, EPS, FK, P, R, S, T, XMA, XPA
      INTEGER K
      SAVE EPS
      DATA EPS / 0.0_ep1 /
!***FIRST EXECUTABLE STATEMENT  ep_9lgic
      IF (EPS.EQ.0.0_ep1) EPS = 0.50_ep1*F1MACH(3,ep_dummy)
!
      XPA = X + 1.00_ep1 - A
      XMA = X - 1.0_ep1 - A
!
      R = 0.0_ep1
      P = 1.0_ep1
      S = P
      DO 10 K=1,300
        FK = K
        T = FK*(A-FK)*(1.0_ep1+R)
        R = -T/((XMA+2.0_ep1*FK)*(XPA+2.0_ep1*FK)+T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'ep_9lgic', 'NO CONVERGENCE IN 300 TERMS OF CONTINUED FRACTION', 1, 2)
!
 20   ep_9lgic = A*ALX - X + LOG(S/XPA)
!
      RETURN
      END FUNCTION

!>DECK wp_9lgit
!>***BEGIN PROLOGUE  wp_9lgit
!>***SUBSIDIARY
!>***PURPOSE  Compute the logarithm of Tricomi's incomplete Gamma
!>            function with Perron's continued fraction for large X and
!>            A .GE. X.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7E
!>***TYPE      REAL(kind=wp) (R9LGIT-S, wp_9lgit-D)
!>***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
!>             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> Compute the log of Tricomi's incomplete gamma function with Perron's
!> continued fraction for large X and for A .GE. X.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  F1MACH, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770701  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   900720  Routine changed from user-callable to subsidiary.  (WRB)
!>***END PROLOGUE  wp_9lgit
      REAL(kind=wp) FUNCTION wp_9lgit (A, X, ALGAP1)
      IMPLICIT NONE
      REAL(kind=wp) A, X, ALGAP1, AX, A1X, EPS, FK, HSTAR, P, R, S, SQEPS, T
      INTEGER K
      LOGICAL FIRST
      SAVE EPS, SQEPS, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_9lgit
      IF (FIRST) THEN
         EPS = 0.50_wp*F1MACH(3,wp_dummy)
         SQEPS = SQRT(F1MACH(4,wp_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LE. 0.0_wp .OR. A .LT. X) CALL XERMSG ('SLATEC', 'wp_9lgit', 'X SHOULD BE GT 0.0 AND LE A', 2, 2)
!
      AX = A + X
      A1X = AX + 1.00_wp
      R = 0.0_wp
      P = 1.0_wp
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.0_wp+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'wp_9lgit', 'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 3, 2)
!
 30   HSTAR = 1.00_wp - X*S/A1X
      IF (HSTAR .LT. SQEPS) CALL XERMSG ('SLATEC', 'wp_9lgit', 'RESULT LESS THAN HALF PRECISION', 1, 1)
!
      wp_9lgit = -X - ALGAP1 - LOG(HSTAR)
      RETURN
!
      END FUNCTION
!
      !> Quad precision version of wp_9lgit.
      REAL(kind=ep1) FUNCTION ep_9lgit (A, X, ALGAP1)
      IMPLICIT NONE
      REAL(kind=ep1) A, X, ALGAP1, AX, A1X, EPS, FK, HSTAR, P, R, S, SQEPS, T
      INTEGER K
      LOGICAL FIRST
      SAVE EPS, SQEPS, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_9lgit
      IF (FIRST) THEN
         EPS = 0.50_ep1*F1MACH(3,ep_dummy)
         SQEPS = SQRT(F1MACH(4,ep_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LE. 0.0_ep1 .OR. A .LT. X) CALL XERMSG ('SLATEC', 'ep_9lgit', 'X SHOULD BE GT 0.0 AND LE A', 2, 2)
!
      AX = A + X
      A1X = AX + 1.00_ep1
      R = 0.0_ep1
      P = 1.0_ep1
      S = P
      DO 20 K=1,200
        FK = K
        T = (A+FK)*X*(1.0_ep1+R)
        R = T/((AX+FK)*(A1X+FK)-T)
        P = R*P
        S = S + P
        IF (ABS(P).LT.EPS*S) GO TO 30
 20   CONTINUE
      CALL XERMSG ('SLATEC', 'ep_9lgit', 'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 3, 2)
!
 30   HSTAR = 1.00_ep1 - X*S/A1X
      IF (HSTAR .LT. SQEPS) CALL XERMSG ('SLATEC', 'ep_9lgit', 'RESULT LESS THAN HALF PRECISION', 1, 1)
!
      ep_9lgit = -X - ALGAP1 - LOG(HSTAR)
      RETURN
!
      END FUNCTION
!>DECK wp_gamic
!>***BEGIN PROLOGUE  wp_gamic
!>***PURPOSE  Calculate the complementary incomplete Gamma function.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7E
!>***TYPE      REAL(kind=wp) (GAMIC-S, wp_gamic-D)
!>***KEYWORDS  COMPLEMENTARY INCOMPLETE GAMMA FUNCTION, FNLIB,
!>             SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!>   Evaluate the complementary incomplete Gamma function
!>
!>   wp_gamic = integral from X to infinity of EXP(-T) * T**(A-1.)  .
!>
!>   wp_gamic is evaluated for arbitrary real values of A and for non-
!>   negative values of X (even though wp_gamic is defined for X .LT.
!>   0.0), except that for X = 0 and A .LE. 0.0, wp_gamic is undefined.
!>
!>   wp_gamic, A, and X are REAL(kind=wp).
!>
!>   A slight deterioration of 2 or 3 digits accuracy will occur when
!>   wp_gamic is very large or very small in absolute value, because log-
!>   arithmic variables are used.  Also, if the parameter A is very close
!>   to a negative INTEGER (but not a negative integer), there is a loss
!>   of accuracy, which is reported if the result is less than half
!>   machine precision.
!>
!>***REFERENCES  W. Gautschi, A computational procedure for incomplete
!>                 gamma functions, ACM Transactions on Mathematical
!>                 Software 5, 4 (December 1979), pp. 466-481.
!>               W. Gautschi, Incomplete gamma functions, Algorithm 542,
!>                 ACM Transactions on Mathematical Software 5, 4
!>                 (December 1979), pp. 482-489.
!>***ROUTINES CALLED  F1MACH, cfp_9gmic, cfp_9gmit, cfp_9lgic, cfp_9lgit, cfp_lgams,
!>                    cfp_lngam, XERCLR, XERMSG
!>***REVISION HISTORY  (YYMMDD)
!>   770701  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!>   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
!>***END PROLOGUE  wp_gamic
      REAL(kind=wp) FUNCTION wp_gamic (A, X)
      IMPLICIT NONE
      REAL(kind=wp) A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNGS, ALX, BOT, E, EPS, GSTAR, H, SGA, SGNG, SGNGAM, SGNGS, SQEPS, T
      INTEGER IZERO
      LOGICAL FIRST
      SAVE EPS, SQEPS, ALNEPS, BOT, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  wp_gamic
      IF (FIRST) THEN
         EPS = 0.50_wp*F1MACH(3,wp_dummy)
         SQEPS = SQRT(F1MACH(4,wp_dummy))
         ALNEPS = -LOG (F1MACH(3,wp_dummy))
         BOT = LOG (F1MACH(1,wp_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LT. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_gamic', 'X IS NEGATIVE' , 2, 2)
!
      IF (X.GT.0.0_wp) GO TO 20
      IF (A .LE. 0.0_wp) CALL XERMSG ('SLATEC', 'wp_gamic', 'X = 0 AND A LE 0 SO wp_gamic IS UNDEFINED', 3, 2)
!
      wp_gamic = EXP (cfp_lngam(A+1.0_wp) - LOG(A))
      RETURN
!
 20   ALX = LOG (X)
      SGA = 1.00_wp
      IF (A.NE.0.0_wp) SGA = SIGN (1.00_wp, A)
      AINTA = AINT (A + 0.50_wp*SGA)
      AEPS = A - AINTA
!
      IZERO = 0
      IF (X.GE.1.00_wp) GO TO 40
!
      IF (A.GT.0.50_wp .OR. ABS(AEPS).GT.0.0010_wp) GO TO 30
      E = 2.00_wp
      IF (-AINTA.GT.1.0_wp) E = 2.0_wp*(-AINTA+2.0_wp)/(AINTA*AINTA-1.00_wp)
      E = E - ALX * X**(-0.0010_wp)
      IF (E*ABS(AEPS).GT.EPS) GO TO 30
!
      wp_gamic = cfp_9gmic (A, X, ALX)
      RETURN
!
 30   CALL cfp_lgams (A+1.00_wp, ALGAP1, SGNGAM)
      GSTAR = cfp_9gmit (A, X, ALGAP1, SGNGAM, ALX)
      IF (GSTAR.EQ.0.0_wp) IZERO = 1
      IF (GSTAR.NE.0.0_wp) ALNGS = LOG (ABS(GSTAR))
      IF (GSTAR.NE.0.0_wp) SGNGS = SIGN (1.00_wp, GSTAR)
      GO TO 50
!
 40   IF (A.LT.X) wp_gamic = EXP (cfp_9lgic(A, X, ALX))
      IF (A.LT.X) RETURN
!
      SGNGAM = 1.00_wp
      ALGAP1 = cfp_lngam (A+1.00_wp)
      SGNGS = 1.00_wp
      ALNGS = cfp_9lgit (A, X, ALGAP1)
!
! EVALUATION OF wp_gamic(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
!
 50   H = 1.0_wp
      IF (IZERO.EQ.1) GO TO 60
!
      T = A*ALX + ALNGS
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.00_wp - SGNGS*EXP(T)
!
!      IF (ABS(H).LT.SQEPS) CALL XERCLR
      IF (ABS(H) .LT. SQEPS) CALL XERMSG ('SLATEC', 'wp_gamic', 'RESULT LT HALF PRECISION', 1, 1)
!
 60   SGNG = SIGN (1.00_wp, H) * SGA * SGNGAM
      T = LOG(ABS(H)) + ALGAP1 - LOG(ABS(A))
!      IF (T.LT.BOT) CALL XERCLR
      wp_gamic = SGNG * EXP(T)
      RETURN
!
 70   SGNG = -SGNGS * SGA * SGNGAM
      T = T + ALGAP1 - LOG(ABS(A))
!      IF (T.LT.BOT) CALL XERCLR
      wp_gamic = SGNG * EXP(T)
      RETURN
!
      END FUNCTION
!
      !> Quad precision version of wp_gamic.
      REAL(kind=ep1) FUNCTION ep_gamic (A, X)
      IMPLICIT NONE
      REAL(kind=ep1) A, X, AEPS, AINTA, ALGAP1, ALNEPS, ALNGS, ALX, BOT, E, EPS, GSTAR, H, SGA, SGNG, SGNGAM, SGNGS, SQEPS, T
      INTEGER IZERO
      LOGICAL FIRST
      SAVE EPS, SQEPS, ALNEPS, BOT, FIRST
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ep_gamic
      IF (FIRST) THEN
         EPS = 0.50_ep1*F1MACH(3,ep_dummy)
         SQEPS = SQRT(F1MACH(4,ep_dummy))
         ALNEPS = -LOG (F1MACH(3,ep_dummy))
         BOT = LOG (F1MACH(1,ep_dummy))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LT. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_gamic', 'X IS NEGATIVE' , 2, 2)
!
      IF (X.GT.0.0_ep1) GO TO 20
      IF (A .LE. 0.0_ep1) CALL XERMSG ('SLATEC', 'ep_gamic', 'X = 0 AND A LE 0 SO ep_gamic IS UNDEFINED', 3, 2)
!
      ep_gamic = EXP (cfp_lngam(A+1.0_ep1) - LOG(A))
      RETURN
!
 20   ALX = LOG (X)
      SGA = 1.00_ep1
      IF (A.NE.0.0_ep1) SGA = SIGN (1.00_ep1, A)
      AINTA = AINT (A + 0.50_ep1*SGA)
      AEPS = A - AINTA
!
      IZERO = 0
      IF (X.GE.1.00_ep1) GO TO 40
!
      IF (A.GT.0.50_ep1 .OR. ABS(AEPS).GT.0.0010_ep1) GO TO 30
      E = 2.00_ep1
      IF (-AINTA.GT.1.0_ep1) E = 2.0_ep1*(-AINTA+2.0_ep1)/(AINTA*AINTA-1.00_ep1)
      E = E - ALX * X**(-0.0010_ep1)
      IF (E*ABS(AEPS).GT.EPS) GO TO 30
!
      ep_gamic = cfp_9gmic (A, X, ALX)
      RETURN
!
 30   CALL cfp_lgams (A+1.00_ep1, ALGAP1, SGNGAM)
      GSTAR = cfp_9gmit (A, X, ALGAP1, SGNGAM, ALX)
      IF (GSTAR.EQ.0.0_ep1) IZERO = 1
      IF (GSTAR.NE.0.0_ep1) ALNGS = LOG (ABS(GSTAR))
      IF (GSTAR.NE.0.0_ep1) SGNGS = SIGN (1.00_ep1, GSTAR)
      GO TO 50
!
 40   IF (A.LT.X) ep_gamic = EXP (cfp_9lgic(A, X, ALX))
      IF (A.LT.X) RETURN
!
      SGNGAM = 1.00_ep1
      ALGAP1 = cfp_lngam (A+1.00_ep1)
      SGNGS = 1.00_ep1
      ALNGS = cfp_9lgit (A, X, ALGAP1)
!
! EVALUATION OF ep_gamic(A,X) IN TERMS OF TRICOMI-S INCOMPLETE GAMMA FN.
!
 50   H = 1.0_ep1
      IF (IZERO.EQ.1) GO TO 60
!
      T = A*ALX + ALNGS
      IF (T.GT.ALNEPS) GO TO 70
      IF (T.GT.(-ALNEPS)) H = 1.00_ep1 - SGNGS*EXP(T)
!
!      IF (ABS(H).LT.SQEPS) CALL XERCLR
      IF (ABS(H) .LT. SQEPS) CALL XERMSG ('SLATEC', 'ep_gamic', 'RESULT LT HALF PRECISION', 1, 1)
!
 60   SGNG = SIGN (1.00_ep1, H) * SGA * SGNGAM
      T = LOG(ABS(H)) + ALGAP1 - LOG(ABS(A))
!      IF (T.LT.BOT) CALL XERCLR
      ep_gamic = SGNG * EXP(T)
      RETURN
!
 70   SGNG = -SGNGS * SGA * SGNGAM
      T = T + ALGAP1 - LOG(ABS(A))
!      IF (T.LT.BOT) CALL XERCLR
      ep_gamic = SGNG * EXP(T)
      RETURN
!
      END FUNCTION
!>DECK wp_lgams
!>***BEGIN PROLOGUE  wp_lgams
!>***PURPOSE  Compute the logarithm of the absolute value of the Gamma
!>            function.
!>***LIBRARY   SLATEC (FNLIB)
!>***CATEGORY  C7A
!>***TYPE      REAL(kind=wp) (ALGAMS-S, wp_lgams-D)
!>***KEYWORDS  ABSOLUTE VALUE OF THE LOGARITHM OF THE GAMMA FUNCTION,
!>             FNLIB, SPECIAL FUNCTIONS
!>***AUTHOR  Fullerton, W., (LANL)
!>***DESCRIPTION
!>
!> wp_lgams(X,DLGAM,SGNGAM) calculates the REAL(kind=wp) natural
!> logarithm of the absolute value of the Gamma function for
!> REAL(kind=wp) argument X and stores the result in double
!> precision argument DLGAM.
!>
!>***REFERENCES  (NONE)
!>***ROUTINES CALLED  cfp_lngam
!>***REVISION HISTORY  (YYMMDD)
!>   770701  DATE WRITTEN
!>   890531  Changed all specific intrinsics to generic.  (WRB)
!>   890531  REVISION DATE from Version 3.2
!>   891214  Prologue converted to Version 4.0 format.  (BAB)
!>***END PROLOGUE  wp_lgams
      SUBROUTINE wp_lgams (X, DLGAM, SGNGAM)
      IMPLICIT NONE
      REAL(kind=wp) X, DLGAM, SGNGAM
      INTEGER INT
!***FIRST EXECUTABLE STATEMENT  wp_lgams
      DLGAM = cfp_lngam(X)
      SGNGAM = 1.00_wp
      IF (X.GT.0.0_wp) RETURN
!
      INT = MOD (-AINT(X), 2.00_wp) + 0.10_wp
      IF (INT.EQ.0) SGNGAM = -1.00_wp
!
      RETURN
      END SUBROUTINE
!
      !> Quad precision version of wp_lgams.
      SUBROUTINE ep_lgams (X, DLGAM, SGNGAM)
      IMPLICIT NONE
      REAL(kind=ep1) X, DLGAM, SGNGAM
      INTEGER INT
!***FIRST EXECUTABLE STATEMENT  ep_lgams
      DLGAM = cfp_lngam(X)
      SGNGAM = 1.00_ep1
      IF (X.GT.0.0_ep1) RETURN
!
      INT = MOD (-AINT(X), 2.00_ep1) + 0.10_ep1
      IF (INT.EQ.0) SGNGAM = -1.00_ep1
!
      RETURN
      END SUBROUTINE
  !
  !> \par Purpose:
  !> \verbatim
  !> Calculate the full set of inverse normalization factors Ilm for the real unnormalized solid spherical harmonics with l,m .le. L using recursion.
  !> The unnormalized solid spherical harmonics are defined in Rico, Lopez, et al. Int. J. Q. Chem.,2012,DOI:10.1002/qua.24356.
  !> \endverbatim
  !> \f[
  !>   I_{lm} = \frac{1}{N_{lm}} = \sqrt{\frac{2\pi(1+\delta_{m,0})}{2l+1}\frac{(l+|m|)!}{(l-|m|)!}}.
  !> \f]
  !> \todo {add error checking for case l > size(INorm(,:))}
  !> \verbatim
  !>  Nlm is the normalization factor for the real unnormalized solid spherical harmonics zlm.
  !> \endverbatim
  !> \param[out] INorm
  !> \verbatim
  !>  Real array (-L:L,0:L) containing the values of the inverse normalization factors Ilm for all l .le. L and the corresponding values of m.
  !> \endverbatim
  !
  !> \param[in] L
  !> \verbatim
  !>  Integer angular momentum of the last Ilm to be generated.
  !> \endverbatim
  subroutine cfp_nlm(INorm,L)
    use precisn
    use phys_const, only: twopi, fourpi, rtwo
    implicit none
    integer, intent(in) :: L
    real(kind=cfp), intent(out) :: INorm(-L:L,0:L)
  
    integer :: l_it, m_it
    real(kind=cfp) :: f
  
    INorm = 0.0_cfp
  
    INorm(0,0) = sqrt(fourpi)
   
    !use recursion to calculate the INorm
    do l_it = 1,L
       INorm(0,l_it) = 1.0_cfp
       do m_it = 1,l_it
          INorm(m_it,l_it) = INorm(m_it-1,l_it)*sqrt(1.0_cfp*(l_it+m_it)*(l_it-m_it+1))
          INorm(-m_it,l_it) = INorm(m_it,l_it)
       enddo
       f = sqrt(twopi/(2*l_it+1.0_cfp))
       INorm(0,l_it) = rtwo
       INorm(-l_it:l_it,l_it) = f*INorm(-l_it:l_it,l_it)
    enddo
  
  end subroutine cfp_nlm
  !
  !> \par Purpose:
  !> \verbatim
  !> Calculate the full set of real normalized spherical harmonics with l,m .le. L using recursion. See Helgaker, p.210, p.215-218 (Section 6.4) for details.
  !> The real spherical harmonics are related to the solid ones using the formula:
  !> \endverbatim
  !> \f[
  !>   X_{lm}(\Omega) = S_{lm}(x/r,y/r,z/r)\sqrt{\frac{2l+1}{4\pi}}.
  !> \f]
  !> We define the normalization factor:
  !> \f[
  !>   n_{lm} = \sqrt{\frac{2l+1}{4\pi}}
  !> \f]
  !> and then use Helgaker's method to calculate the 'scaled' solid harmonics:
  !> \f[
  !>   X_{lm}(\Omega) = {\overline{S}_{lm}}(x/r,y/r,z/r)=n_{lm}S_{lm}(x/r,y/r,z/r).
  !> \f]
  !> \verbatim
  !> The real spherical harmonics defined in this way are equivalent to the real spherical harmonics defined using the standard formula:
  !> \endverbatim
  !> \f[
  !>     X_{lm}(\Omega) = 
  !>     \begin{cases}
  !>       {\sqrt{2}}(-1)^{m}{\cal{R}}\left(Y_{l\vert m\vert}(\Omega)\right), m>0 \\
  !>       {\sqrt{2}}(-1)^{m}{\cal{I}}\left(Y_{l\vert m\vert}(\Omega)\right), m<0 \\
  !>       Y_{l0}(\Omega), m=0,
  !>     \end{cases}
  !> \f]
  !> where the complex spherical harmonics \$Y_{lm}(\Omega)\$ satisfy the Condon-Shortley phase convention.
  !> Note that the convention of Homeier and Steinborn is not to include the factor (-1)**m in the definition of the real spherical harmonics. However, here we follow the std. definition and include the factor
  !> (-1)**m in the definition of the real spherical harmonics.
  !> \todo {add error checking for case L > size(SH(,:))}
  !
  !> \param[out] SH
  !> \verbatim
  !>  Real array (-L:L,0:L) containing the values of the spherical harmonics for all l .le. L and the corresponding values of m.
  !>  According to my tests the results are in vast majority of cases accurate to full REAL(kind=cfp). In only a few cases the precision in the last 1 or 2 digits is lost.
  !> \endverbatim
  !
  !> \param[in] X, Y, Z
  !> \verbatim
  !>  Real input coordinates for which the real harmonics will be evaluated. 
  !>  The radius vector r={X,Y,Z} need not be normalized to 1 since this is done automatically during the calculation.
  !> \endverbatim
  !
  !> \param[in] L
  !> \verbatim
  !>  Integer angular momentum of the last spherical harmonic to be generated. For L = 0; SH = 1/sqrt(4*pi)
  !> \endverbatim
  !rsdef
  subroutine cfp_resh(SH,X,Y,Z,L)
    use precisn
    use phys_const, only: fourpi
    implicit none
    integer, intent(in) :: L
    real(kind=cfp), intent(out) :: SH(-L:L,0:L)
    real(kind=cfp), intent(in) :: x, y, z

    integer :: l_it, m_it, lp1
    real(kind=cfp) :: r, xn, yn, zn, fac, l2p1, l2m1, rl2p1
    real(kind=cfp), parameter :: norm = 1.0_cfp/sqrt(fourpi)

    SH = 0.0_cfp !vectorized
    SH(0,0) = norm !ensure the S(l=0,m=0) is always returned even if r .eq. 0 below

    if (L .eq. 0) return

    !normalize the radius vector
    r = sqrt(x**2.0_cfp+y**2.0_cfp+z**2.0_cfp)
    if (r == 0.0_cfp) return
    xn = x/r
    yn = y/r
    zn = z/r

    !initialize the starting values
    SH(0,0) = norm
    SH(-1,0) = 0.0_cfp

    !recursion
    !l_it = 0 case: L=1
    fac = sqrt(3.0_cfp/fourpi)
    SH(1,1) = fac*xn
    SH(0,1) = fac*zn
    SH(-1,1)= fac*yn

    !diagonal recursions
    do l_it = 1, L-1
       fac = sqrt((2.0_cfp*l_it+3.0_cfp)/(2.0_cfp*l_it+2.0_cfp))
       SH(l_it+1, l_it+1) = fac*(xn*SH(l_it,l_it)-yn*SH(-l_it,l_it))
       SH(-l_it-1,l_it+1) = fac*(yn*SH(l_it,l_it)+xn*SH(-l_it,l_it))

       l2p1 = 2*l_it + 1.0_cfp
       rl2p1 = sqrt(l2p1)
       l2m1 = 2*l_it - 1.0_cfp
       lp1 = l_it + 1
       !vertical recursions (vectorized)
       do m_it = -l_it,l_it
          SH(m_it,l_it+1) = sqrt((l2p1+2.0_cfp)/((lp1+m_it)*(lp1-m_it))) &
                * (rl2p1*zn*SH(m_it,l_it)-sqrt((l_it+m_it)*(l_it-m_it)/l2m1)*SH(m_it,l_it-1))
       end do
    end do
  
  end subroutine cfp_resh
  !
  !> \par Purpose:
  !> \verbatim
  !> Calculate the full set of real normalized solid spherical harmonics with l,m .le. L using recursion. See Helgaker, p.218 (Section 6.4) for details.
  !> For l=m=0: Slm=1.
  !> \endverbatim
  !> \todo {add error checking for case L > size(SH(,:))}
  !
  !> \param[out] SH
  !> \verbatim
  !>  Real array (-L:L,0:L) containing the values of the solid spherical harmonics for all l .le. L and the corresponding values of m.
  !> \endverbatim
  !
  !> \param[in] X, Y, Z
  !> \verbatim
  !>  real input coordinates at which the solid harmonics will be evaluated.
  !> \endverbatim
  !
  !> \param[in] L
  !> \verbatim
  !>  Integer angular momentum of the last spherical harmonic to be generated.
  !> \endverbatim
   
  subroutine cfp_solh(SH,x,y,z,L)
    use precisn
    implicit none
    integer, intent(in) :: L
    real(kind=cfp), intent(out) :: SH(-L:L,0:L)
    real(kind=cfp), intent(in) :: x, y, z
  
    integer :: l_it, m_it, lp1, l2p1
    real(kind=cfp) :: rsq, fac
  
    SH = 0.0_cfp !vectorized
   
    !initialize the starting values
    SH(0,0) = 1.0_cfp
    SH(-1,0) = 0.0_cfp
    rsq = x*x + y*y + z*z
  
    !recursion
    !l_it = 0 case: L=1
    SH(1,1) = x
    SH(0,1) = z
    SH(-1,1)= y
  
    !diagonal recursions
    do l_it = 1, L-1
       fac = sqrt((2.0_cfp*l_it+1.0_cfp)/(2.0_cfp*l_it+2.0_cfp))
       SH(l_it+1, l_it+1) = fac*(x*SH(l_it,l_it)-y*SH(-l_it,l_it))
       SH(-l_it-1,l_it+1) = fac*(y*SH(l_it,l_it)+x*SH(-l_it,l_it))
       l2p1 = 2*l_it + 1
       lp1 = l_it + 1
       !vertical recursions (vectorized)
       do m_it = -l_it,l_it
          SH(m_it,l_it+1) = (l2p1*z*SH(m_it,l_it)-sqrt(1.0_cfp*(l_it+m_it)*(l_it-m_it))*rsq*SH(m_it,l_it-1)) &
                            / sqrt(1.0_cfp*(lp1+m_it)*(lp1-m_it))
       end do
    end do
  
  end subroutine cfp_solh
  !
  !> \par Purpose:
  !> \verbatim
  !> Calculate the full set of real unnormalized solid spherical harmonics with l,m .le. L using recursion. See Helgaker, p.210, p.218 (Section 6.4) for details.
  !> The unnormalized solid spherical harmonics are defined in Rico, Lopez, et al. Int. J. Q. Chem.,2012,DOI:10.1002/qua.24356.
  !> The mehod of calculation is to calculate the real solid spherical harmonics and then renormalize them to obtain the surface 
  !> spherical harmonics using the formula:
  !> \endverbatim
  !> \f[
  !>   z_{lm}(x,y,z) = S_{lm}(x,y,z)\sqrt{\frac{1+\delta_{m,0}}{2}\frac{(l+|m|)!}{(l-|m|)!}}.
  !> \f]
  !> \todo {rewrite the algorithm to calculate zlm directly using Helgaker's algorithm for the scaled solid harmonics rather than by calculating Slm first and then renormalizing}
  !> \todo {add error checking for case L > size(SH(,:))}
  !
  !> \param[out] SH
  !> \verbatim
  !>  Real array (-L:L,0:L) containing the values of the real unnormalized solidspherical harmonics for all l .le. L and the corresponding values of m.
  !> \endverbatim
  !
  !> \param[in] X, Y, Z
  !> \verbatim
  !>  Real input coordinates for which the harmonics will be evaluated.
  !> \endverbatim
  !
  !> \param[in] L
  !> \verbatim
  !>  Integer angular momentum of the last harmonic to be generated.
  !> \endverbatim
   
  subroutine cfp_zhar(SH,X,Y,Z,L)
    use precisn
    use phys_const, only: roneh
    implicit none
    integer, intent(in) :: L
    real(kind=cfp), intent(out) :: SH(-L:L,0:L)
    real(kind=cfp), intent(in) :: x, y, z
  
    integer :: l_it, m_it
    real(kind=cfp) :: Nlm
  
    !calculate the solid spherical harmonics
    call cfp_solh(SH,x,y,z,L)
  
    !unnormalize
    do l_it = 2,L
       Nlm = 1.0_cfp
       do m_it = 1,l_it
          Nlm = Nlm*sqrt(1.0_cfp*(l_it+m_it)*(l_it-m_it+1))
          SH(m_it,l_it) = Nlm*roneh*SH(m_it,l_it)
          SH(-m_it,l_it)= Nlm*roneh*SH(-m_it,l_it)
       enddo
    enddo
  
  end subroutine cfp_zhar

  !> Returns the string of Boys functions \f$F_{m}(T)\f$ for \f$m=0,\dots,m_{max}\f$.
  !> \f[
  !>    F_{m}(T) = \int_{0}^{1}x^{2m}\exp[-Tx^2]dx.
  !> \f]
  !> The result is a double precision array boys_function(1:mmax+1) containing the value of the Boys function for \f$m=0,\dots,m_{max}\f$.
  !> The Boys function is calculated first for \f$m = m_{max}\f$ using the series expansion:
  !> \f[
  !>    F_{m}(T) = \exp[-T]\sum_{i=0}^{\infty}\frac{(2m-1)!!(2T)^{i}}{(2m+2i+1)!!}
  !> \f]
  !> and then for \f$m = m_{max}-1,\dots,0\f$ using the recurrence:
  !> \f[
  !>    F_{m}(T) = \frac{2TF_{m+1}(T)+\exp[-T]}{2m+1}.
  !> \f]
  !> For T > 85 we use the asymptotic formula for the Boys function:
  !> \f[
  !>    F_{m}(T) \approx \frac{1}{2T^{m+1/2}}\Gamma(m+1/2).
  !> \f]
  !> \param[in] T Real value corresponding to T in \f$F_{m}(T)\f$.
  !> \param[in] mmax Integer value corresponding to \f$m_{max}\f$ in \f$F_{m}(T)\f$, \f$m=0,\dots,m_{max}\f$.
  !> boys_function Double precision array boys_function(1:mmax+1) containing the values: \f$F_{m}(T)\f$, \f$m=0,\dots,m_{max}\f$.
  !> \warning It is vital for calculation of the two electron integrals that the Boys function is calculated as accurately as possible since it is used to start the recurrent evaluation of the integrals.
  !> It is important to understand that the limits on the maximum GTO L used are connected with the precision of calculation of the corresponding Boys functions. Given the largest GTO L in the basis the 
  !> largest required mmax is 4L. The parameters involved in the calculation of the Boys function (see the module const) have been tuned for the case L_max=6, i.e. mmax=24.
  function boys_function(T,mmax)
     use precisn
     use utils, only: xermsg
     use const, only: boys_f_dprec_asym_thr, boys_tol
     implicit none
     real(kind=wp), intent(in) :: T
     integer, intent(in) :: mmax
     real(kind=wp) :: boys_function(1:mmax+1)

     integer :: m, prev
     real(kind=wp) :: T_pow, two_T, exp_T, s, s_prev, term, fac

        if (T < 0.0_wp .or. mmax < 0) then
           call xermsg('special_functions','boys_function','Invalid input parameters.',1,1)
        endif

        two_T = 2.0_wp*T
        exp_T = exp(-T)

        if (T .ge. boys_f_dprec_asym_thr) then !use the asymptotic formula; 85 is the smallest exponent for which the asymptotic formula gives fully accurate results up to mmax = 24, but we can afford 60.

           boys_function(mmax+1) = cfp_gamma_fun(mmax+0.5_wp)*0.5_wp/(T**(mmax+0.5_wp))
           
        else !use the power series for F_{m}(T)

           !i=0: starting term for the sum
           fac = 2*mmax+1.0_wp !the ratio of the factorials: (2*m-1)!!/(2*m+2*i+1)!! reduces to 1/(polynomial_in_m); fac = polynomial_in_m
           prev = 2*mmax+1
           T_pow = 1.0_wp
           s = T_pow/fac
           s_prev = 0.0_wp
   
           !sum over i=1,...,until convergence
           !note that the terms responsible for numerical problems for T > 65 are T_pow and fac; if these are evaluated in quad precision then the limitations on T and mmax can be dropped.
           do
              prev = prev + 2
              fac = fac*prev        !1/fac = (2*m-1)!!/(2*m+2*i+1)!!
              T_pow = T_pow*two_T   !(2*T)**i
              !todo precalculate 1/fac terms - that should result in a significant speed up
              term = T_pow/fac      !the next term in the sum: (2*m-1)!!/(2*m+2*i+1)!!*(2*T)**i
              s = s + term
              if (term .le. s*boys_tol) exit !convergence criterion
              s_prev = s
           enddo
   
           boys_function(mmax+1) = exp_T*s !Boys function for m=mmax

        endif

        !Boys function for m=mmax-1,...,0
        do m=mmax-1,0,-1
           boys_function(m+1) = (boys_function(m+2)*two_T+exp_T)/real(2*m+1,kind=wp)
        enddo

  end function boys_function

  !> Quadruple precision version of boys_function.
  !> See boys_function for details on the method of evaluation.
  function boys_function_quad(T,mmax)
     use precisn
     use utils, only: xermsg
     use const, only: boys_tol_qprec
     implicit none
     real(kind=ep1), intent(in) :: T
     integer, intent(in) :: mmax
     real(kind=ep1) :: boys_function_quad(1:mmax+1)

     integer :: m, prev
     real(kind=ep1) :: qT, T_pow, two_T, exp_T, s, s_prev, term, fac

        if (T < 0.0_ep1 .or. mmax < 0) then
           call xermsg('special_functions','boys_function_quad','Invalid input parameters.',1,1)
        endif

        qT = real(T,ep1) !make T in the quad precision
        two_T = 2.0_ep1*qT
        exp_T = exp(-qT)

        !use the power series for F_{m}(T)

        !i=0: starting term for the sum
        fac = 2*mmax+1.0_ep1 !the ratio of the factorials: (2*m-1)!!/(2*m+2*i+1)!! reduces to 1/(polynomial_in_m); fac = polynomial_in_m
        prev = 2*mmax+1
        T_pow = 1.0_ep1
        s = T_pow/fac
        s_prev = 0.0_ep1
   
        !sum over i=1,...,until convergence
        do
           prev = prev + 2
           !calculating fac and T_pow (and term) in quad precision is the step which removes the numerical instabilities present in boys_function_quad
           fac = fac*prev        !1/fac = (2*m-1)!!/(2*m+2*i+1)!!
           T_pow = T_pow*two_T   !(2*T)**i
           term = T_pow/fac      !the next term in the sum: (2*m-1)!!/(2*m+2*i+1)!!*(2*T)**i
           s = s + term
           if (term .le. s*boys_tol_qprec) exit !convergence criterion
           s_prev = s
        enddo
   
        boys_function_quad(mmax+1) = exp_T*s !Boys function for m=mmax

        !Boys function for m=mmax-1,...,0
        do m=mmax-1,0,-1
           boys_function_quad(m+1) = (boys_function_quad(m+2)*two_T+exp_T)/(2*m+1.0_ep1)
        enddo

  end function boys_function_quad

  !> Computes the double factorial \f$(n)!!\f$ of the integer number n.
  !> n must be .ge. -1.
  elemental function dfact(n)
     use precisn
     implicit none
     integer, intent(in) :: n
     integer :: dfact

     integer :: i

        if (n < -1) then
           dfact = 0
           return
        endif

        if (n .eq. 0 .or. n .eq. -1) then
           dfact = 1
           return
        endif

        dfact = n
        do i=n-2,1,-2
           dfact = dfact*i
        enddo

  end function dfact

  !> This function evaluates a polynomial of degree \f$n\f$ at point \f$x\f$ using the Horner form for polynomials.
  !> \param[in] n Order of the polynomial
  !> \param[in] x The point at which to evaluate.
  !> \param[in] a The array of the coefficients a(1:n+1): \f$a_{0},a_{1},\dots,a_{n}\f$.
  function cfp_eval_poly_horner_single(n,x,a)
     implicit none
     integer, intent(in) :: n
     real(kind=cfp), intent(in) :: x
     real(kind=cfp), intent(in) :: a(1:n+1)
     real(kind=cfp) :: cfp_eval_poly_horner_single

     integer :: i
     real(kind=cfp) :: inv_x

        !Implement a stable version of the scheme for both x small and x large.
        if (abs(x) .le. 1.0_cfp) then
           cfp_eval_poly_horner_single = a(n+1)
           do i=n,1,-1
              cfp_eval_poly_horner_single = x*cfp_eval_poly_horner_single + a(i)
           enddo
        else !abs(x) > 1
           inv_x = 1.0_cfp/x
           cfp_eval_poly_horner_single = a(1)
           do i=2,n+1
              cfp_eval_poly_horner_single = cfp_eval_poly_horner_single*inv_x + a(i)
           enddo
           cfp_eval_poly_horner_single = cfp_eval_poly_horner_single*x**n
        endif

  end function cfp_eval_poly_horner_single

  !> This subroutine evaluates a set of polynomials of degree \f$n\f$ at a number of points \f$x\f$ using the Horner form for polynomials.
  !> \param[in] n Order of the polynomials.
  !> \param[in] x Array of points at which to evaluate.
  !> \param[in] n_x The number of points (dimension of x).
  !> \param[in] a The array of the coefficients a(1:n_x,1:n+1): different for each point x.
  !> \param[out] res The polynomials evaluated at points x.
  subroutine cfp_eval_poly_horner_many(n,x,n_x,a,res)
     implicit none
     integer, intent(in) :: n, n_x
     real(kind=cfp), intent(in) :: x(n_x)
     real(kind=cfp), intent(in) :: a(n_x,n+1)
     real(kind=cfp), intent(out) :: res(n_x)

     integer :: i, i_1, err, j
     real(kind=cfp), allocatable :: inv_x(:)

        !We assume that the values in x are in non-decreasing order.
        i_1 = -1
        do i=1,n_x
           if (x(i) .le. 1.0_cfp) i_1 = i
        enddo !i

        !Implement a stable version of the scheme for both x small and x large.
        if (i_1 > 0) then
           res(1:i_1) = a(1:i_1,n+1)
           do i=n,1,-1
              res(1:i_1) = x(1:i_1)*res(1:i_1) + a(1:i_1,i)
           enddo
        endif

        if (i_1 < n_x) then
           j = max(i_1+1,1)

           allocate(inv_x(j:n_x),stat=err)
           if (err .ne. 0) call xermsg('cgto_pw_expansions_mod','cfp_eval_poly_horner_many','Memory allocation failed.',err,1)

           inv_x(j:n_x) = 1.0_cfp/x(j:n_x)
           res(j:n_x) = a(j:n_x,1)
           do i=2,n+1
              res(j:n_x) = res(j:n_x)*inv_x(j:n_x) + a(j:n_x,i)
           enddo
           res(j:n_x) = res(j:n_x)*x(j:n_x)**n
        endif

  end subroutine cfp_eval_poly_horner_many

  !> For given L,M values of the real solid harmonic this function returns the number of terms in the spherical -> cartesian mapping.
  function no_terms_sph_to_cart_mapping(l,m)
     implicit none
     integer, intent(in) :: l,m
     integer :: no_terms_sph_to_cart_mapping

     integer :: tmp_i
     real(kind=cfp) :: vm

        if (m .ge. 0) then
           vm = 0.0_cfp
        else
           vm = 0.5_cfp
        endif

        !calculate the number of terms in the expansion and allocate storage for the list of cartesian exponents and coefficients
        tmp_i = floor((l-abs(m))/2.0)
        no_terms_sph_to_cart_mapping = (floor(abs(m)/2.0-vm)+1)*(tmp_i+1+(tmp_i*(tmp_i+1))/2) !number of terms in the cartesian -> spherical GTO formula

  end function no_terms_sph_to_cart_mapping

  !> This routine constructs for a given real solid harmonic L,M the list of cartesian functions which build up this real solid harmonic.
  !> On output the array c(:) contains the coefficients with which the individual cartesians contribute to the solid harmonic. Exponents of the cartesians are then given in the arrays i_exp,j_exp,k_exp.
  !> \todo the construction of the coefficient list (the summations) should be included along the lines of sph_shell_to_cart_shell
  subroutine cfp_sph_to_cart_mapping(l,m,c,i_exp,j_exp,k_exp)
     implicit none
     integer, intent(in) :: l,m
     integer, allocatable, intent(out) :: i_exp(:),j_exp(:),k_exp(:)
     real(kind=cfp), allocatable, intent(out) :: c(:)

     integer :: tmp_i, t, u, v_it, no_terms, cnt, err
     real(kind=cfp) :: v, vm, nlm_s 

        if (l < 0 .or. abs(m) > l) call xermsg ('gto_function', 'sph_to_cart_mapping', 'The spherical GTO L,M are invalid.', 1, 1)

        if (allocated(c)) deallocate(c)
        if (allocated(i_exp)) deallocate(i_exp)
        if (allocated(j_exp)) deallocate(j_exp)
        if (allocated(k_exp)) deallocate(k_exp)

        if (m .ge. 0) then
           vm = 0.0_cfp
        else
           vm = 0.5_cfp
        endif

        !calculate the number of terms in the expansion and allocate storage for the list of cartesian exponents and coefficients
        tmp_i = floor((l-abs(m))/2.0)
        no_terms = no_terms_sph_to_cart_mapping(l,m) !(floor(abs(m)/2.0-vm)+1)*(tmp_i+1+(tmp_i*(tmp_i+1))/2) !number of terms in the cartesian -> spherical GTO formula
        allocate(i_exp(1:no_terms),j_exp(1:no_terms),k_exp(1:no_terms),c(1:no_terms),stat=err)
        if (err .ne. 0) call xermsg ('gto_function', 'sph_to_cart_mapping', 'Memory allocation failed.',2,1)

        !generate the list of cartesians forming this spherical GTO, see Helgaker Section 9.1.2
        nlm_s = 1.0_cfp / (2**abs(m)*cfp_gamma_fun(l+1.0_cfp)) &
                * sqrt(2.0_cfp*cfp_gamma_fun(l+abs(m)+1.0_cfp)*cfp_gamma_fun(l-abs(m)+1.0_cfp))
        if (m .eq. 0) nlm_s = nlm_s/sqrt(2.0_cfp)
        cnt = 0
        do t=0,tmp_i
           do u=0,t
              do v_it=0,floor(abs(m)/2.0-vm)
                 v = v_it + vm
                 cnt = cnt + 1
                 i_exp(cnt) = 2*t+abs(m)-nint(2*(u+v)) !x exponent of the cartesian GTO
                 j_exp(cnt) = nint(2*(u+v))            !y exponent of the cartesian GTO
                 k_exp(cnt) = l-2*t-abs(m)             !z exponent of the cartesian GTO
                 c(cnt) = nlm_s*(-1)**(t+v_it)*(0.25_cfp**t)*cfp_binom(l,t,cfp_dummy)*cfp_binom(l-t,abs(m)+t,cfp_dummy) &
                        * cfp_binom(t,u,cfp_dummy)*cfp_binom(abs(m),nint(2*v),cfp_dummy) !coefficient of this cartesian GTO
              enddo !v_it
           enddo !u
        enddo !t

  end subroutine cfp_sph_to_cart_mapping

  !> This routine constructs for a given real solid harmonic L the matrix of coefficients for all M of the transformation from the cartesian harmonics (~x^i*y^j*z^k, i+j+k=L).
  !> On output the linear array c(:) emulates 2D matrix with 2*L+1 rows and (L+1)*(L+2)/2 columns. The rows of c correspond to the M values: -L,-L+1,...,0,1,...,L.  
  !> The coefficients in each row are ordered so that the n-th coefficient in the given row of c corresponds to the coefficient for the cartesian harmonic with the shell canonical index n.
  subroutine cfp_sph_shell_to_cart_shell(l,c)
     implicit none
     integer, intent(in) :: l
     real(kind=cfp), intent(out) :: c(:)

     integer :: i_exp,j_exp,k_exp
     integer :: m, tmp_i, t, u, v_it, can, ncart, space, stride
     real(kind=cfp) :: v, vm, nlm_s

        if (l < 0) call xermsg ('gto_function', 'sph_shell_to_cart_shell', 'The spherical GTO L<0.', 1, 1)

        ncart = (l+1)*(l+2)/2 !the number of cartesian harmonics in this shell
        space = ncart*(2*l+1)

        stride = 2*l+1 !the memory stride in the array c

        if (size(c) < space) call xermsg ('gto_function', 'sph_shell_to_cart_shell', 'The output array c is too small.',2,1)

        c(1:space) = 0.0_cfp

        if (l .eq. 0) then
           c(1) = 1.0_cfp
           return
        endif

!        print *,'OLD'
        do m =-l,l

           if (m .ge. 0) then
              vm = 0.0_cfp
           else
              vm = 0.5_cfp
           endif

           !calculate the number of terms in the expansion and allocate storage for the list of cartesian exponents and coefficients
           tmp_i = floor((l-abs(m))/2.0)

           !generate the list of cartesians forming this spherical GTO, see Helgaker Section 9.1.2
           nlm_s = 1.0_cfp / (2**abs(m)*cfp_gamma_fun(l+1.0_cfp)) &
                    * sqrt(2.0_cfp*cfp_gamma_fun(l+abs(m)+1.0_cfp)*cfp_gamma_fun(l-abs(m)+1.0_cfp))
           if (m .eq. 0) nlm_s = nlm_s/sqrt(2.0_cfp)
           v = 0.0_cfp
           do t=0,tmp_i
              do u=0,t
                 do v_it=0,floor(abs(m)/2.0-vm)
                    v = v_it + vm
                    i_exp = 2*t+abs(m)-nint(2*(u+v)) !x exponent of the cartesian GTO
                    j_exp = nint(2*(u+v))            !y exponent of the cartesian GTO
                    k_exp = l-2*t-abs(m)             !z exponent of the cartesian GTO
                    can = (l-i_exp)*(l-i_exp+1)/2 + k_exp + 1 !the in-shell canonical index of this cartesian harmonic
                    !coefficient of this cartesian GTO: note that we add it; in some cases this algorithm generates only partial cartesian contributions for the current (i,j,k) triplet and it is these
                    !contributions that we are gradually summing here.
                    c(m+l+1+stride*(can-1)) = c(m+l+1+stride*(can-1)) + nlm_s*(-1)**(t+v_it)*(0.25_cfp**t) &
                                            * cfp_binom(l,t,cfp_dummy)*cfp_binom(l-t,abs(m)+t,cfp_dummy) &
                                            * cfp_binom(t,u,cfp_dummy)*cfp_binom(abs(m),nint(2*v),cfp_dummy)
                 enddo !v_it
              enddo !u
           enddo !t

!           print *,'l,m',l,m
!           t = 0
!           do can=1,ncart
!              if (c(m+l+1+stride*(can-1)) .ne. 0.0_cfp) then
!                 t = t + 1
!                 print *,can,c(m+l+1+stride*(can-1))
!              endif
!           enddo
!           print *,'nonzero/full',t,ncart

        enddo !m

  end subroutine cfp_sph_shell_to_cart_shell

  !> This routine constructs for all real solid harmonic with L .le. l the matrix of coefficients for all M of the transformation from the cartesian harmonics (~x^i*y^j*z^k, i+j+k=L).
  !> On output the linear array c(:) emulates 2D matrix with 2*L+1 rows and (L+1)*(L+2)/2 columns. The rows of c correspond to the M values: -L,-L+1,...,0,1,...,L.  
  !> The coefficients in each row are ordered so that the n-th coefficient in the given row of c corresponds to the coefficient for the cartesian harmonic with the shell canonical index n.
  subroutine cfp_sph_shell_to_cart_lshells(l,nz,c_nz,nz_can)
     implicit none
     integer, intent(in) :: l
     real(kind=cfp), intent(out) :: c_nz(:)
     integer, intent(out) :: nz(:), nz_can(:)

     integer :: i_exp,j_exp,k_exp, l_it
     integer :: m, tmp_i, t, u, v_it, can, ncart, space, stride, ind, ind_lm
     real(kind=cfp) :: v, vm, nlm_s
     real(kind=cfp), allocatable :: c(:)

        if (l < 0) call xermsg ('gto_function', 'sph_shell_to_cart_shell', 'The spherical GTO L<0.', 1, 1)

        ncart = (l+1)*(l+2)/2 !the number of cartesian harmonics in this shell
        space = ncart*(2*l+1)

        allocate(c(space),stat=m)
        if (m .ne. 0) call xermsg ('gto_function', 'sph_shell_to_cart_shell', 'Memory allocation failed.', m, 1)

        c_nz = 0.0_cfp
        nz = 0
        nz_can = 0

        ind = 0
        ind_lm = 0
          
        do l_it = 0,l

           stride = 2*l_it+1 !the memory stride in the array c (row length)

           do m =-l_it,l_it

              ind_lm = ind_lm + 1
              if (l_it .eq. 0) then
                 ind = ind + 1
                 c_nz(ind) = 1.0_cfp
                 nz(ind) = 1
                 nz_can(ind) = 1
                 cycle !go to next l
              endif

              c = 0.0_cfp
   
              if (m .ge. 0) then
                 vm = 0.0_cfp
              else
                 vm = 0.5_cfp
              endif
   
              !calculate the number of terms in the expansion and allocate storage for the list of cartesian exponents and coefficients
              tmp_i = floor((l_it-abs(m))/2.0)
   
              !generate the list of cartesians forming this spherical GTO, see Helgaker Section 9.1.2
              nlm_s = 1.0_cfp/(2**abs(m)*cfp_gamma_fun(l_it+1.0_cfp)) &
                        * sqrt(2.0_cfp*cfp_gamma_fun(l_it+abs(m)+1.0_cfp)*cfp_gamma_fun(l_it-abs(m)+1.0_cfp))
              if (m .eq. 0) nlm_s = nlm_s/sqrt(2.0_cfp)
              v = 0.0_cfp
              do t=0,tmp_i
                 do u=0,t
                    do v_it=0,floor(abs(m)/2.0-vm)
                       v = v_it + vm
                       i_exp = 2*t+abs(m)-nint(2*(u+v)) !x exponent of the cartesian GTO
                       j_exp = nint(2*(u+v))            !y exponent of the cartesian GTO
                       k_exp = l_it-2*t-abs(m)             !z exponent of the cartesian GTO
                       can = (l_it-i_exp)*(l_it-i_exp+1)/2 + k_exp + 1 !the in-shell canonical index of this cartesian harmonic
                       !coefficient of this cartesian GTO: note that we add it; this algorithm generates partial cartesian contributions for the current (i,j,k) triplet and it is these
                       !contributions that we are gradually summing here.
                       c(m+l_it+1+stride*(can-1)) = c(m+l_it+1+stride*(can-1)) + nlm_s*(-1)**(t+v_it)*(0.25_cfp**t) &
                                * cfp_binom(l_it,t,cfp_dummy)*cfp_binom(l_it-t,abs(m)+t,cfp_dummy) &
                                * cfp_binom(t,u,cfp_dummy)*cfp_binom(abs(m),nint(2*v),cfp_dummy)
                    enddo !v_it
                 enddo !u
              enddo !t
   
              !Save the nonzero coefficients and the canonical indices of the corresponding cartesians. All non-zero indices are saved in the c_nz array one-by-one.
              nz(ind_lm) = 0
              do can=1,ncart
                 if (c(m+l_it+1+stride*(can-1)) .ne. 0.0_cfp) then
                    nz(ind_lm) = nz(ind_lm) + 1
                    ind = ind + 1
                    nz_can(ind) = can
                    c_nz(ind) = c(m+l_it+1+stride*(can-1))
                 endif
              enddo
   
           enddo !m

        enddo !l_it

        deallocate(c)

  end subroutine cfp_sph_shell_to_cart_lshells

end module special_functions
