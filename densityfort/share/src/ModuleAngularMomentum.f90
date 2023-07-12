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
!> 
!> Defines a number of special functions and symbols of the
!> quantum theory of angular momentum. Unless otherwise
!> stated, this module follow the conventions of 
!>
!> \emph{Quantum Theory of Angular Momentum} by 
!> D. A. Varshalovich, A. N. Moskalev, and V. K. Kersonskii
!> World Scientific, Singapore 
!> \cite{Varshalovich}
!>
module ModuleAngularMomentum

  !implicit none

  private

  logical, private :: NOT_INITIALIZED_YET=.TRUE.
  
  real(kind(1d0)), parameter :: CGCONTRACTION_THRESHOLD = 1.d-5

  integer        , private, parameter :: DIMFACT=170
  real(kind(1d0)), private :: fact(0:DIMFACT)
  integer        , private, parameter :: DIMLOGFACT=500
  real(kind(1d0)), private :: logfact(0:DIMLOGFACT)

  private :: initsymba, delta

  public :: ClebschGordanCoefficient
  public :: ThreeJSymbol
  public :: SixJSymbol
  public :: NineJSymbol
  public :: ThreeJCoefficient

  public :: bispharm
  public :: RAssLegFun, RLegPol, LegPolF, LegPolS
  public :: Plm,Ylm,BipolarSphericalHarmonics,SymmetryAdaptedBipolarSphericalHarmonics,WRME
  public :: mtssabsh,smmtssabsh,spsabsh,angcoup2,angcouppm, cartesianToSpherical


  public :: TestClebschGordanCoefInteger
  public :: TestClebschGordanCoefHalfInteger
  public :: TestSixJSymbolsHalfInteger

  public :: CheckClebschGordanContractions_1st1B
  public :: CheckClebschGordanContractions_3rd1B
  public :: CheckClebschGordanContractions_1st2B
  public :: CheckClebschGordanContractions_2nd2B
  public :: CheckClebschGordanContractions_3rd2B
  public :: CheckClebschGordanContractions_4th2B
  public :: CheckClebschGordanContractions_6th2B

  public :: Check6jSpecialFormulas

!!$  public :: CheckClebschGordanFormulas
  public :: CG_SaS_10_SbS
  public :: CG_Tt_Jmt_K0
  public :: CG_T0_J0_K0
  public :: CG_12a_12b_Ttau
  
  public :: Clebsch
  public :: SixJSymbol_HalfHalf
  
contains

  subroutine CartesianToSpherical(Rvec,R,theta,phi)
    real(kind(1d0)), intent(in)  :: Rvec(3)
    real(kind(1d0)), intent(out) :: R, theta, phi
    R=sqrt(sum(Rvec*Rvec))
    if(abs(R)<1d-40)then
       theta=0.d0
       phi=0.d0
    endif
    theta=acos(Rvec(3)/R)
    phi=atan2(Rvec(2),Rvec(1))
  end subroutine CartesianToSpherical

  subroutine InitThisModule()
    integer :: i
    fact(0)=1.d0
    do i = 1,DIMFACT
       fact(i)=fact(i-1)*dble(i)
    enddo
    logfact(0:DIMFACT)=log(fact(0:DIMFACT))
    do i=DIMFACT+1,DIMLOGFACT
       logfact(i)=logfact(i-1)+log(1.d0*i)
    enddo
    NOT_INITIALIZED_YET=.FALSE.
    return
  end subroutine InitThisModule


  real(kind(1d0)) function ClebschGordanCoefficient(J1,J2,J3,M1,M2) result( cgc )

    IMPLICIT REAL*8(A-H,O-Z)

    if(Not_Initialized_Yet)call InitThisModule

    cgc = 0.d0

    IF(J1<ABS(M1).or.J2<ABS(M2).or.J3<ABS(M1+M2)) return
    IF(J3>(J1+J2).or.J3<ABS(J1-J2))return
    IF((ABS(M1)+ABS(M2))==0.and.BTEST(J1+J2+J3,0)) return 

    IA1=J3+J2-J1
    M3=M1+M2
    IA2=J3+M3
    IA3=J2+M3-J1
    IF(IA3.GT.0) THEN
       NI=IA3
    ELSE
       NI=0
    END IF
    IF(IA2.LE.IA1) THEN
       NM=IA2
    ELSE
       NM=IA1
    END IF
    L1=J3+J1-J2
    L2=J1+J2+J3+1
    L3=J1+J2-J3
    L4=J1-M1
    L5=J2-M2
    L6=J2+M2
    L7=J3-M3
    L8=J1+M1
    L9=2*J3+1
    CC=L9*fact(L1)/fact(L2)*fact(IA1)*fact(L3)/fact(L4)/fact(L5)*&
         fact(IA2)/fact(L6)*fact(L7)/fact(L8)
    CC=DSQRT(CC)
    IP1=J2+J3+M1-NI
    B1=fact(IP1)
    IP2=J1-M1+NI
    B2=fact(IP2)
    IP2=IP2+1
    D1=fact(NI)
    IR1=NI+1
    IR2=J3-J1+J2-NI
    D2=fact(IR2)
    IR3=J3+M3-NI
    D3=fact(IR3)
    IR4=J1-J2-M3+NI
    D4=fact(IR4)
    IR4=IR4+1
    FAC=1.d0
    IF(BTEST(NI+J2+M2,0)) FAC=-FAC
    S1=B1/D2*B2/(D1*D3*D4)*FAC
    N=NM-NI

    IF(N/=0)THEN
       FA=S1
       DO I=1,N
          FA=-FA*IP2*IR2/IP1*IR3/(IR1*IR4)
          S1=S1+FA
          IP1=IP1-1
          IP2=IP2+1
          IR1=IR1+1
          IR2=IR2-1
          IR3=IR3-1
          IR4=IR4+1
       ENDDO
    ENDIF

    CGC=CC*S1
    if(abs(CGC)<1.d-14)CGC=0.d0

    RETURN

  END FUNCTION ClebschGordanCoefficient


  real(kind(1d0)) function ThreeJCoefficient(j1,j2,j3,m1,m2)
    implicit none
    integer, intent(in) :: j1,j2,j3,m1,m2
    if(Not_Initialized_Yet)call InitThisModule
    ThreeJCoefficient = &
         ClebschGordanCoefficient(j1,j2,j3,m1,m2) / &
         sqrt(dble(2*j3+1))*(1-2*mod(abs(j1-j2+m1+m2),2))
    return
  end function ThreeJCoefficient


  real(kind(1d0)) function Racah(ja,jb,je,jd,jc,jf)
    integer, intent(in) :: ja,jb,je,jd,jc,jf
    Racah=SixJSymbol(ja,jb,je,jd,jc,jf)
    if(mod(ja+jb+jc+jd,2)==1)Racah=-Racah
    return
  end function Racah


  !> \internal
  !> If the triangular conditions for the 6j symbol
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    j_a & j_b & j_c \\
  !>    j_d & j_e & j_f
  !>    \end{array}\right\}
  !> \f], i.e.
  !> \f{eqnarray} 
  !>   |j_a-jb|\leq j_c\leq j_a+j_b, \\
  !>   |j_a-je|\leq j_f\leq j_a+j_e, \\
  !>   |j_d-jb|\leq j_f\leq j_d+j_b, \\
  !>   |j_d-je|\leq j_c\leq j_d+j_e,
  !> \f}
  !> are satisfied, then 
  !>
  !>     lSixJSymbol=.TRUE.
  !>
  !> otherwise
  !>
  !>     lSixJSymbol=.FALSE.
  !>
  logical function lSixJSymbol(ja,jb,jc,jd,je,jf)
    implicit none
    integer, intent(in) :: ja,jb,jc,jd,je,jf
    lSixJSymbol=.false.
    if(ja>jb+jc.or.ja<abs(jb-jc))return
    if(ja>je+jf.or.ja<abs(je-jf))return
    if(jb>jd+jf.or.jb<abs(jd-jf))return
    if(jc>jd+je.or.jc<abs(jd-je))return
    lSixJSymbol=.true.
    return
  end function lSixJSymbol


  !> Returns the 6j symbol
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    j_a & j_b & j_c \\
  !>    j_d & j_e & j_f
  !>    \end{array}\right\}
  !> \f].
  real(kind(1d0)) function SixJSymbol(ja,jb,jc,jd,je,jf) 
    !
    implicit real*8(a-h,o-z)
    !
    integer, intent(in) :: ja,jb,jc,jd,je,jf
    !
    if(Not_Initialized_Yet)call InitThisModule
    !
    SixJSymbol = 0.d0
    !
    if((ja+jb-jc)*(ja-jb+jc)*(jb+jc-ja).lt.0) return
    if((jd+jb-jf)*(jd-jb+jf)*(jb+jf-jd).lt.0) return
    if((je+jd-jc)*(je-jd+jc)*(jd+jc-je).lt.0) return
    if((ja+je-jf)*(ja-je+jf)*(je+jf-ja).lt.0) return
    !
    w = 0.d0
    !
    iabe=ja+jb+jc
    icde=je+jd+jc
    iacf=ja+je+jf
    ibdf=jb+jd+jf
    !
    jacdb=ja+je+jb+jd
    jadef=ja+jd+jc+jf
    jbcef=jb+je+jc+jf
    !
    ia=max(iabe,icde,iacf,ibdf)
    ib=min(jacdb,jadef,jbcef)
    !
    iin=ia
    ifin=ib
    !
    sig=1.d0
    sigf=sig
    !
    if(btest(iin,0)) sig=-sig
    !
    aa=0.d0
    do i=iin,ifin
       bb=fact(i-iabe)*fact(i-icde)*fact(i-iacf)*&
            fact(i-ibdf)*fact(jacdb-i)*fact(jadef-i)*fact(jbcef-i)
       aa=aa+sig*fact(i+1)/bb
       sig=-sig
    enddo
    !
    SixJSymbol=aa*delta(ja,jb,jc)*delta(jb,jd,jf)*&
         delta(je,jd,jc)*delta(ja,je,jf)*sigf
    !
    if(abs(SixJSymbol)<1.d-14)SixJSymbol=0.d0
    !
    return
    !
  end function SixJSymbol


  !> Computes the function
  !>\f[
  !>   \Delta = \sqrt{
  !>        \frac{ (j_1+j_2-j_3)! \, (j_1-j_2+j_3)! \, (j_2+j_3-j_1)! }{ 
  !>                                (j_1+j_2+j_3+1)! } }
  !>\f]
  real(kind(1d0)) function Delta(j1,j2,j3)
    implicit real*8 (a-h,o-z)
    integer, intent(in) :: j1, j2, j3
    if(Not_Initialized_Yet)call InitThisModule
    Delta=sqrt(fact(j1+j2-j3)*fact(j1-j2+j3)*fact(j2+j3-j1)/fact(j1+j2+j3+1))
    return
  end function Delta


  logical function TriangularCondition(j1,j2,j3)
    integer :: j1,j2,j3
    TriangularCondition = &
         (  j1 >= abs( j2 - j3 )  ) .and. &
         (  j1 <=      j2 + j3    )
  end function TriangularCondition


  !.. Check triangular conditions for the 9J Symbol
  !    
  !     | j11 j12 j13 |
  !    <  j21 j22 j23  >
  !     | j31 j32 j33 |
  !..
  logical function lNineJSymbol(j11,j12,j13,j21,j22,j23,j31,j32,j33)
    implicit none
    integer, intent(in) :: j11,j12,j13,j21,j22,j23,j31,j32,j33
    lNineJSymbol=.false.
    if(.not.TriangularCondition(j11,j12,j13))return
    if(.not.TriangularCondition(j21,j22,j23))return
    if(.not.TriangularCondition(j31,j32,j33))return
    if(.not.TriangularCondition(j11,j21,j31))return
    if(.not.TriangularCondition(j12,j22,j32))return
    if(.not.TriangularCondition(j13,j23,j33))return
    lNineJSymbol=.true.
  end function lNineJSymbol


  !> Returns the 9j symbol
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    j_{11} & j_{12} & j_{13} \\
  !>    j_{21} & j_{22} & j_{23} \\
  !>    j_{31} & j_{32} & j_{33}
  !>    \end{array}\right\}
  !> \f].
  !> 
  !> computed using Eq.5 on pag. 305 of Varshalovich 
  !>
  real(kind(1d0)) function NineJSymbol(j11,j12,j13,j21,j22,j23,j31,j32,j33)
    implicit none
    integer, intent(in) :: j11,j12,j13,j21,j22,j23,j31,j32,j33
    integer :: ell, ellmi, ellma

    if(Not_Initialized_Yet)call InitThisModule

    NineJSymbol=0.d0

    if(.not.TriangularCondition(j11,j12,j13))return
    if(.not.TriangularCondition(j21,j22,j23))return
    if(.not.TriangularCondition(j31,j32,j33))return
    if(.not.TriangularCondition(j11,j21,j31))return
    if(.not.TriangularCondition(j12,j22,j32))return
    if(.not.TriangularCondition(j13,j23,j33))return

    ellmi = max(abs(j11-j33),abs(j32-j21),abs(j23-j12))
    ellma = min(    j11+j33 ,    j32+j21 ,    j23+j12 )
    do ell = ellmi, ellma
       NineJSymbol = NineJSymbol + dble(2*ell+1) * &
            SixJSymbol(j11,j33,ell,j32,j21,j31) * &
            SixJSymbol(j32,j21,ell,j23,j12,j22) * &
            SixJSymbol(j23,j12,ell,j11,j33,j13)
    enddo
  end function NineJSymbol



  logical function delparf(j1,j2,j3,j4,ll)
    !VERIFICA PARITA' DI j1+j2+j3+j4 E LIMITI li e lm
    !DEI POSSIBILI VALORI DI l COMPATIBILI CON LE
    !CONDIZIONI DI TRIANGOLARITA' D(j1,j2,l) D(j3,j4,l)
    !log1=TRUE SE AMBEDUE I TESTS SONO POSITIVI
    integer, intent(in) :: j1,j2,j3,j4,ll
    delparf=.FALSE.
    if(btest(j1+j2+j3+j4,0))return
    if(btest(j1+j2+ll,0))return
    if(ll<max(abs(j1-j2),abs(j3-j4)))return
    if(ll>min(j1+j2,j3+j4))return
    delparf=.TRUE.
    return
  end function delparf

  logical function TRITPARF(j1,j2,j3,j4,j5,j6,ll)
    !CONTROLLA CONDIZIONI TRIANGOLARI D(j1,j2,j) D(j3,j4,j)
    !E D(j5,j6,j). SE LE CONDIZIONI SONO SODDISFATTE TUTTE
    !PER UN INTERVALLO DI VALORI DI j log1=true ed in ji e jm
    !SARANNO POSTI I VALORI MINIMO E MASSIMO DI j
    !ALTRIMENTI log1=false. Si assicura che la parita' di j5+j6 sia la
    !stessa che quella di j1+j2, aggiustando il range relativo.
    implicit none
    integer, intent(in) :: j1,j2,j3,j4,j5,j6,ll
    integer :: j12,l12,j34,l34,j56,l56,ji,jm
    tritparf=.FALSE.
    j12 =j1+j2
    l12=iabs(j1-j2)
    j34=j3+j4
    l34=iabs(j3-j4)
    j56=j5+j6
    l56=iabs(j5-j6)
    j56=j56-mod(j56+j12,2)
    l56=l56+mod(j56+j12,2)
    ji=max(l12,l34,l56)
    jm=min(j12,j34,j56)
    if(ll<ji.or.ll>jm)return
    tritparf=.TRUE.
    return
  end function TRITPARF

  logical function trif(j1,j2,j3)
    !CONTROLLA CONDIZIONE TRANGOLARE D(j1,j2,j3)
    !.TRUE. SE LA CONDIZIONE E' VERIFICATA
    implicit none
    integer, intent(in) :: j1,j2,j3
    trif=(j3>=abs(j1-j2)).and.(j3<=j1+j2)
    return
  end function trif

  logical function tritf(j1,j2,j3,j4,j5,j6,ji,jm)
    !Controlla l'intervallo [ji,jm] di j in cui sono soddisfatte 
    !le condizioni triangolari D(j1,j2,j) D(j3,j4,j) E D(j5,j6,j). 
    !Se l'intervallo esiste (ji<=jm) tritf=.TRUE.
    implicit none
    integer, intent(in) :: j1,j2,j3,j4,j5,j6
    integer, intent(out):: ji,jm
    jm=min(j1+j2,j3+j4,j5+j6)
    ji=max(abs(j1-j2),abs(j3-j4),abs(j5-j6))
    tritf=(jm>=ji)
    return
  end function tritf

  real(kind(1d0)) function sbcsn(j1,j2,j3)
    !CALCOLA IL COEFFICIENTE CS(J1,J2,J3)=C(J1,J2,J3;0,0)*[(2J1+1)*(2J2+1)/(2J3+1)]**1/2
    implicit none
    integer,intent(in) :: j1,j2,j3
    sbcsn=ClebschGordanCoefficient(j1,j2,j3,0,0)*sqrt(dble((2*j1+1)*(2*j2+1))/dble(2*j3+1))
    return
  end function sbcsn

  real(kind(1d0)) function mypar(n)
    implicit none
    integer, intent(in) :: n
    mypar=1.d0-2.d0*mod(n,2)
    return
  end function mypar

  complex(kind(1d0)) function BISPHARM(l1,l2,L,M,th1,ph1,th2,ph2)
    implicit none
    integer        , intent(in) :: l1,l2,L,M
    real(kind(1d0)), intent(in) :: th1,ph1,th2,ph2
    !Calcola l'armonica sferica bipolare assegnati i due angoli
    !e tutti i numeri quantici necessari: i momenti angolari
    !delle due armoniche sferiche accoppiate, quello dell'accoppiamento
    !complessivo e la sua proiezione sull'asse di quantizzazione.
    ![Y_l1(theta_1,phi_1)(x)Y_l2(theta_2,phi_2)]_{LM}=
    !sum_m C_{l_1m,l_2 M-m}^{LM}Y_{l_1m}(theta_1,phi_1)Y_{l_2M-m}(theta_2,phi_2)
    !(Varshalovic et al. "Quantum Theory of Angular Momentum" par. 15.16.1
    ! Ed. World Scientific (Singapore)

    complex(kind(1d0)) :: cf
    real(kind(1d0)) :: x1,x2,rf,w1,w2,w3,phi
    integer :: mu

    BISPHARM=(0.d0,0.d0)
    if(l1<0.or.l2<0.or.L<abs(l1-l2).or.L>l1+l2.or.abs(M)>L)return
    x1=cos(th1)
    x2=cos(th2)
    do mu=max(-l1,M-l2),min(l1,M+l2)
       w1=ClebschGordanCoefficient(l1,l2,L,mu,M-mu)
       w2=NORLEG(l1,mu  ,x1)
       w3=NORLEG(l2,M-mu,x2)
       rf=w1*w2*w3
       phi=dble(mu)*ph1+dble(M-mu)*ph2
       cf=rf*(cos(phi)*(1.d0,0.d0)+sin(phi)*(0.d0,1.d0))
       BISPHARM=BISPHARM+cf
    end do
    return 
  end function BISPHARM


  real(kind(1d0)) function NORLEG(l,m,c)
    implicit none
    integer        , intent(in) :: l, m
    real(kind(1d0)), intent(in) :: c
    !Funzione naive per il calcolo dei polinomi di 
    !legendre normalizzati, scritta per disperazione:
    !quelle di libreria non funzionano. In seguito 
    !porremo rimedio. In ogni caso questa subroutine
    !puo` essere utilizzata da una piu` generale nei
    !casi piu` semplici dal momento che e` prevedibilmente
    !piu` celere delle altre.
    !History (yymmdd):
    !051219: 1 error found.
    real(kind(1d0)), parameter :: M_PI=3.14159265358979323844d0
    real(kind(1d0)) :: s
    NORLEG=0.d0
    if(abs(c)>1.d0.or.l<0.or.l>5.or.abs(m)>l)return
    s=sqrt(1.d0-c*c)
    select case(l)
    case(0)
       NORLEG=1.d0/(2.d0*sqrt(M_PI))
    case(1)
       select case(m)
       case( 1) 
          NORLEG=-0.5d0*sqrt(1.5d0/M_PI)*s
       case( 0) 
          NORLEG= 0.5d0*sqrt(3.d0/M_PI)*c
       case(-1)
          NORLEG= 0.5d0*sqrt(1.5d0/M_PI)*s
       end select
    case(2)
       select case(m)
       case( 2) 
          NORLEG= 0.25d0*sqrt(7.5d0/M_PI)*s*s
       case( 1) 
          NORLEG=-0.50d0*sqrt(7.5d0/M_PI)*s*c
       case( 0)
          NORLEG= 0.25d0*sqrt(5.d0/M_PI)*(3.d0*c*c-1.d0)
       case(-1) 
          NORLEG= 0.50d0*sqrt(7.5d0/M_PI)*s*c
       case(-2) 
          NORLEG= 0.25d0*sqrt(7.5d0/M_PI)*s*s
       end select
    case( 3)
       select case(m)
       case( 3) 
          NORLEG=-0.125d0*sqrt(35.d0/M_PI)*s*s*s
       case( 2) 
          NORLEG= 0.250d0*sqrt(52.5d0/M_PI)*c*s*s
       case( 1) 
          NORLEG=-0.125d0*sqrt(21.d0/M_PI)*(5.d0*c*c-1.d0)*s
       case( 0) 
          NORLEG= 0.250d0*sqrt( 7.d0/M_PI)*(5.d0*c*c-3.d0)*c
       case(-1) 
          NORLEG= 0.125d0*sqrt(21.d0/M_PI)*(5.d0*c*c-1.d0)*s
       case(-2) 
          NORLEG= 0.250d0*sqrt(52.5d0/M_PI)*c*s*s
       case(-3) 
          NORLEG= 0.125d0*sqrt(35.d0/M_PI)*s*s*s
       end select
    case(4)
       select case(m)
       case( 4) 
          NORLEG= 0.1875d0*sqrt(17.5d0/M_PI)*s*s*s*s
       case( 3) 
          NORLEG=-0.3750d0*sqrt(35.0d0/M_PI)*c*s*s*s
       case( 2) 
          NORLEG= 0.3750d0*sqrt( 2.5d0/M_PI)*(7.d0*c*c-1.d0)*s*s
       case( 1) 
          NORLEG=-0.3750d0*sqrt( 5.0d0/M_PI)*(7.d0*c*c-3.d0)*c*s
       case( 0) 
          NORLEG= 0.1875d0*sqrt( 1.0d0/M_PI)*(35.d0*c*c*c*c-30.d0*c*c+3.d0)
       case(-1) 
          NORLEG= 0.3750d0*sqrt( 5.0d0/M_PI)*(7.d0*c*c-3.d0)*c*s
       case(-2) 
          NORLEG= 0.3750d0*sqrt( 2.5d0/M_PI)*(7.d0*c*c-1.d0)*s*s
       case(-3) 
          NORLEG= 0.3750d0*sqrt(35.0d0/M_PI)*c*s*s*s
       case(-4) 
          NORLEG= 0.1875d0*sqrt(17.5d0/M_PI)*s*s*s*s
       end select
    case(5)
       select case(m)
       case( 5) 
          NORLEG=-0.09375d0*sqrt( 77.d0/M_PI)*s*s*s*s*s
       case( 4) 
          NORLEG= 0.18750d0*sqrt(192.5d0/M_PI)*c*s*s*s*s
       case( 3) 
          NORLEG=-0.03125d0*sqrt(385.d0/M_PI)*(9.d0*c*c-1.d0)*s*s*s
       case( 2) 
          NORLEG= 0.12500d0*sqrt(577.5d0/M_PI)*(3.d0*c*c-1.d0)*c*s*s
       case( 1) 
          NORLEG=-0.06250d0*sqrt( 82.5d0/M_PI)*(21.d0*c*c*c*c-14.d0*c*c+ 1.d0)*s
       case( 0) 
          NORLEG= 0.06250d0*sqrt( 11.d0/M_PI)*(63.d0*c*c*c*c-70.d0*c*c+15.d0)*c
       case(-1) 
          NORLEG= 0.06250d0*sqrt( 82.5d0/M_PI)*(21.d0*c*c*c*c-14.d0*c*c+ 1.d0)*s
       case(-2) 
          NORLEG= 0.12500d0*sqrt(577.5d0/M_PI)*(3.d0*c*c-1.d0)*c*s*s
       case(-3) 
          NORLEG= 0.03125d0*sqrt(385.d0/M_PI)*(9.d0*c*c-1.d0)*s*s*s
       case(-4) 
          NORLEG= 0.18750d0*sqrt(192.5d0/M_PI)*c*s*s*s*s
       case(-5) 
          NORLEG= 0.09375d0*sqrt( 77.d0/M_PI)*s*s*s*s*s
       end select
    end select
    return
  end function NORLEG

!!$real(kind(1d0)) function legpolyfun(l,m_,x)
!!$  !This function follows the phase convention
!!$  !of Jackson  (Classical Electrodynamics)
!!$  !Abramovitz  (Handbook of Mathematical Functions)
!!$  !Varshalovic (Quantum Theory of Angular Momentum)
!!$  !while it differs by a factor (-1)^m from the
!!$  !conventions adopted by Messiah (Quantum Mechanics)
!!$  !and Friedrich (Theoretical Atomic Physics).
!!$  implicit none
!!$  integer        , intent(in) :: l,m_
!!$  real(kind(1d0)), intent(in) :: x
!!$  integer         :: m,j
!!$  real(kind(1d0)) :: c1,c2,c3,y
!!$  legpolyfun=0.d0
!!$  m=abs(m_)
!!$  !Check on input parameters
!!$  if(l<0.or.m>l.or.abs(x)>1.d0)return
!!$  !Case l==0
!!$  if(l==0)then
!!$     legpolyfun=1.d0
!!$     return
!!$  endif
!!$  !Case x = +/- 1
!!$  if(1.d0-abs(x)<=epsilon(1.d0))then
!!$     if(m/=0)return
!!$     legpolyfun=1.d0
!!$     if(x<0.and.mod(l,2)==1)legpolyfun=-1.d0
!!$     return
!!$  endif
!!$  !Regular cases
!!$  y=exp(0.5d0*dble(m)*log(1.d0-x*x))
!!$  do j=2*m-1,0,-2
!!$     y=y*dble(j)
!!$  enddo
!!$  if(mod(m,2)==1)y=-y
!!$  c1=y
!!$  c2=0.d0
!!$  do j=m+1,l
!!$     c3=c2
!!$     c2=c1
!!$     c1=dble(2*j-1)/dble(j-m)*x*c2-dble(j+m-1)/dble(j-m)*c3
!!$  enddo
!!$  !Anche se e` asimmetrico rispetto al cambio
!!$  !di segno di m, ce lo teniamo cosi`.
!!$  if(m_<0)then
!!$     y=1.d0
!!$     do j=l-m_+1,l+m_
!!$        y=y*dble(j)
!!$     enddo
!!$     c1=c1/y
!!$  endif
!!$  legpolyfun=c1
!!$  return
!!$end function legpolyfun

  recursive real(kind(1d0)) function RAssLegFun(l,m,x) result (res)
    !Associated Legendre Function
    implicit none
    integer        , intent(in) :: l, m
    real(kind(1d0)), intent(in) :: x
    real(kind(1d0))             :: y
    integer :: j
    res=0.d0
    if(l<0.or.abs(m)>l.or.abs(x)>1.d0)return
    if(abs(m)<l)then
       res =dble(2*l-1)/dble(l-m)*x*RAssLegFun(l-1,m,x)-&
            dble(l+m-1)/dble(l-m)*  RAssLegFun(l-2,m,x)
       return
    endif
    if(l==0)then
       res=1.d0
       return
    endif
    if(1.d0-abs(x)<=epsilon(1.d0))return
    res=exp(0.5d0*dble(l)*log(1.d0-x*x))
    if(mod(l,2)==1)res=-res
    if(m==l)then
       do j=2*l-1,1,-2
          res=res*dble(j)
       enddo
    else
       do j=2*l,1,-2
          res=res/dble(j)
       enddo
    endif
    return
  end function RAssLegFun

  real(kind(1d0)) function LegPolF(l,x) result(res)
    !Legendre Polynomials
    implicit none
    integer        , intent(in) :: l
    real(kind(1d0)), intent(in) :: x
    real(kind(1d0)) :: c_1,c_2
    integer :: ll
    res=0.d0;if(l<0.or.abs(x)>1.d0)return
    res=1.d0;if(l==0)return
    c_1=0
    do ll=1,l
       c_2=c_1
       c_1=res
       res=(dble(2*ll-1)*x*c_1-dble(ll-1)*c_2)/dble(ll)
    enddo
    return
  end function LegPolF

  subroutine LegPolS(l,x,P,dP,P_1,dP_1)
    !Legendre Polynomials
    !calcola P_l(x),dP_l(x)/dx,P_{l-1}(x),dP_{l-1}(x)/dx
    !L'ho testato e mi pare funzionare correttamente
    implicit none
    integer        , intent(in) :: l
    real(kind(1d0)), intent(in) :: x
    real(kind(1d0)), intent(out):: P,dP,P_1,dP_1
    integer :: ll
    real(kind(1d0)) :: c_1,c_2
    real(kind(1d0)) :: d_1,d_2
    P =0.d0
    dP=0.d0
    P_1=0.d0
    dP_1=0.d0
    if(l<0.or.abs(x)>1.d0)return
    P=1.d0;if(l==0)return
    c_1=0
    d_1=0
    do ll=1,l
       c_2=c_1
       d_2=d_1
       c_1=P
       d_1=dP
       P   =(dble(2*ll-1)*x*c_1-dble(ll-1)*c_2)/dble(ll)
       dP  =(dble(2*ll-1)*c_1+dble(2*ll-1)*x*d_1-dble(ll-1)*d_2)/dble(ll)
       P_1 =c_1
       dP_1=d_1
    enddo
    return
  end subroutine LegPolS

  real(kind(1d0)) function Plm(x,l,m) result(res)
    !Compute the Associated Legendre functions defined as 
    !Plm(x)=(-1)^m (1-x^2)^{m/2} d^m/dx^m Pl(x)
    !where Pl(x) is the Legendre polynomials
    implicit none
    real(kind(1d0)), intent(in) :: x
    integer        , intent(in) :: l,m
    real(kind(1d0)) :: st,Pmm,c0,c1,c2
    integer         :: j
    res=0.d0
    if(abs(x)>1.d0.or.abs(m)>l.or.l<0)return
    !Compute P|m|m according to the formulas
    !P|m|  |m|(x)  =(-1)^|m| (2|m|)!/ [2^|m| |m|!] sin^|m|(x)
    !P|m|,-|m|(x)  =              1 / [2^|m| |m|!] sin^|m|(x)
    st=sqrt(1-x*x)
    if(m>=0)then
       Pmm=1.d0-2.d0*mod(m,2)
       do j=1,m
          Pmm=Pmm*0.5d0*dble(m+j)*st
       enddo
    else
       Pmm=1.d0
       do j=1,abs(m)
          Pmm=Pmm*0.5d0/dble(j)*st
       enddo
    endif
    if(l==abs(m))then
       res=Pmm
       return
    endif
    !Compute Plm with the recursive formula
    !(l-m)P_{l,m}=(2l-1)xP_{l-1,m}-(l-1+m)P_{l-2,m}
    c0=Pmm
    c1=Pmm*x*dble(2*m+1)
    do j=abs(m)+2,l
       c2=(dble(2*j-1)*x*c1-dble(j-1+m)*c0)/dble(j-m)
       c0=c1
       c1=c2
    enddo
    res=c1
    return
  end function Plm

  complex(kind(1d0)) function Ylm(theta,phi,l,m) result(res)
    !Computes the spherical harmonics with the conventions of 
    !Varshalovich et al:
    ! Y_{lm}*(theta,phi)=Y_{lm}(theta,-phi)=(-1)^m Y_{l,-m}(theta,phi)
    !Parrebbe funzionare
    implicit none
    real   (kind(1d0)), intent(in) :: theta,phi
    integer           , intent(in) :: l,m
    complex(kind(1d0)), parameter  :: IU=(0.d0,1.d0)
    real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
    integer         :: j
    real(kind(1d0)) :: st,norm,x,Plm_
    res=(0.d0,0.d0)
    if(theta<0.d0.or.theta>PI.or.l<0.or.abs(m)>l)return
    x =cos(theta)
    Plm_=Plm(x,l,m)
    !Normalization factor Sqrt[(2l+1)/(4 Pi) * (l-m)!/(l+m)!]
    norm=1.d0
    do j=-abs(m)+1,abs(m)
       norm=norm*dble(l+j)
    enddo
    if(m>0)norm=1.d0/norm
    norm=sqrt(dble(2*l+1)/(4.d0*PI)*norm)
    res = norm * Plm_ * exp( IU * dble(m) * phi )
    return
  end function Ylm

  !> Compute the bipolar spherical harmonics
  !> \f[
  !>   \mathcal{Y}_{l_1l_2}^{LM}(\Omega_1,\Omega_2)\,=\,
  !>   \sum_{m_1m_2}\,C_{l_1m_1,l_2m_2}^{LM}\,Y_{l_1m_1}(\Omega_1)\,Y_{l_2m_2}(\Omega_2)
  !> \f]
  complex(kind(1d0)) function BipolarSphericalHarmonics(theta1,phi1,theta2,phi2,l1,l2,L,M) result(res)
    implicit none
    real(kind(1d0)), intent(in) :: theta1,phi1,theta2,phi2
    integer        , intent(in) :: l1,l2,L,M
    real(kind(1d0)), parameter  :: PI=3.14159265358979323844d0
    integer            :: mu
    complex(kind(1d0)) :: z1,z2
    res=(0.d0,0.d0)
    if(abs(M)>L.or.max(theta1,theta2)>PI.or.min(theta1,theta2)<0.d0&
         .or.l1<0.or.l2<0.or.L<abs(l1-l2).or.L>(l1+l2))return
    do mu=max(-l1,M-l2),min(l1,M+l2)
       z1=Ylm(theta1,phi1,l1, mu )
       z2=Ylm(theta2,phi2,l2,M-mu)
       res=res+ClebschGordanCoefficient(l1,l2,L,mu,M-mu)*z1*z2
    enddo
    return
  end function BipolarSphericalHarmonics

  !> Compute the symmetry adapted bipolar spherical harmonics
  !> Yl1l2^LM+/-(Omega1,Omega2)=
  !>  =[Y_{l1,l2}^{LM}(Omega1,Omega2)+/-Y_{l1,l2}^{LM}(Omega2,Omega1)]/sqrt(2)
  !> If s is even => Y+, if s is odd => Y-
  complex(kind(1d0)) function SymmetryAdaptedBipolarSphericalHarmonics(theta1,phi1,theta2,phi2,l1,l2,L,M,s) result(res)
    implicit none
    real(kind(1d0)), intent(in) :: theta1,phi1,theta2,phi2
    integer        , intent(in) :: l1,l2,L,M,s
    complex(kind(1d0)) :: z1,z2
    z1=BipolarSphericalHarmonics(theta1,phi1,theta2,phi2,l1,l2,L,M)
    z2=BipolarSphericalHarmonics(theta2,phi2,theta1,phi1,l1,l2,L,M)*(1.d0-2.d0*mod(s,2))
    res=(z1+z2)/sqrt(2.d0)
    return
  end function SymmetryAdaptedBipolarSphericalHarmonics

  complex(kind(1d0)) function WRME(alpha,beta,gamma,m1,m2,j) result(res)
    !Wigner Rotation Matrix Elements
    !not checked yet: be carful!
    implicit none
    real(kind(1d0)), intent(in) :: alpha,beta,gamma
    integer        , intent(in) :: j,m1,m2
    complex(kind(1d0)), parameter  :: IU=(0.d0,1.d0)
    real(kind(1d0)) :: cbm,sbm
    integer         :: k
    res=(0.d0,0.d0)
    if(max(abs(m1),abs(m2))>j)return
    !Attenzione: estende gli elementi di matrice anche oltre 
    !il dominio alpha [0:2 Pi], beta [0,Pi], gamma[0:2 Pi]
    res=(1.d0,0.d0)
    cbm=cos(beta/2.d0)
    sbm=sin(beta/2.d0)
    !Calcolo di djmm
    do k=max(0,-m1-m2),min(j-m1,j-m2)
       res=res+(1.d0,0.d0)*(1.d0-2.d0*mod(k,2))*&
            mypow(cbm,m1+m2+2*k)*mypow(sbm,2*j-m1-m2-2*k)/&
            (fact(k)*fact(j-m1-k)*fact(j-m2-k)*fact(m1+m2+k))
    enddo
    res=res*sqrt(fact(j+m1)*fact(j-m1)*fact(j+m2)*fact(j-m2))*(1.d0-2.d0*mod(j-m2,2))
    res=exp(-IU*(m1*alpha+m2*gamma))*res
    return
  end function WRME

  real(kind(1d0)) function mypow(x,n) result(res)
    real(kind(1d0)), intent(in) :: x
    integer        , intent(in) :: n
    res=1.d0;if(n==0)return
    res=0.d0;if(x<=0)return
    res=exp(dble(n)*log(x))
    return
  end function mypow

  real(kind(1d0)) function spsabsh(x,J,S,l1,l2,l3,l4) result(res)
    !Compute the
    !Scalar Product between Symmetry Adapted Bipolar Spherical Harmonic
    implicit none
    real(kind(1d0)), intent(in) :: x
    integer        , intent(in) :: J,S,l1,l2,l3,l4
    real(kind(1d0)), parameter  :: PI = 3.14159265358979323844d0
    integer                     :: nat,sym,par
    integer                     :: l,lmi,lma

    res=0.d0
    if(  min(l1,l2,l3,l4)<0          .or.&
         J<max(abs(l1-l2),abs(l3-l4)).or.&
         J>min(l1+l2,l3+l4)          .or.&
         abs(x)>1.d0                .or.&
         mod(l1+l2+l3+l4,2)/=0 ) return

    par=mod(l1+l2,2)
    nat=mod(J+par,2)
    sym=mod(S,2)

    lmi=min(max(abs(l1-l3),abs(l2-l4)),max(abs(l1-l4),abs(l2-l3)))
    lma=max(min(l1+l3,l2+l4),min(l1+l4,l2+l3))
    do l=lmi,lma
       if(mod(l1+l3+l,2)==1.and.mod(l1+l4+l,2)==1)cycle
       res=res+Plm(x,l,0)*(1.d0-2.d0*mod(l,2))*(&
            ClebschGordanCoefficient(l1,l3,l,0,0)*&
            ClebschGordanCoefficient(l2,l4,l,0,0)*&
            SixJSymbol(l1,l2,J,l4,l3,l)+&
            ClebschGordanCoefficient(l1,l4,l,0,0)*&
            ClebschGordanCoefficient(l2,l3,l,0,0)*&
            SixJSymbol(l1,l2,J,l3,l4,l)*&
            (1.d0-2.d0*mod(nat+s,2)))
    enddo

    res=res*sqrt(dble((2*l1+1)*(2*l2+1)*(2*l3+1)*(2*l4+1))) *&
         dble(2*J+1)*(1.d0-2.d0*mod(J,2))/(4.d0*PI)**2

    return
  end function spsabsh

  subroutine mtssabsh(J,N,S,lv,DIM)
    !Compute the dimension DIM and angular momenta lv(1:DIM) 
    !of the pairs ( l, J + N - l ) with 2*l >= J + N 
    !of the minimal tensorial set of symmetry adapted
    !bipolar spherical harmonic with the assigned
    !quantum numbers 
    ! J : total angular momentum
    ! N : naturality (N odd => unnatural, N even => natural)
    ! S : parity upon interchange of solid angles
    !     (S even => even, S odd => odd)
    ! Checked
    integer, intent(in) :: J,N,S
    integer, intent(out):: lv(*)
    integer, intent(out):: DIM
    integer :: l
    DIM=0
    if(J<0)return
    l=int(dble(J+mod(N,2))/2.d0+0.6d0)
    if(mod(J+N,2)==0.and.(mod(N+S,2)==1))l=l+1
    do while(l<=J)
       DIM=DIM+1
       lv(DIM)=l
       l=l+1
    enddo
    return
  end subroutine mtssabsh

  subroutine smmtssabsh(x,J,N,S,A,LDA,lv,DIM)
    !Compute the
    !Superposition Matrix of a Minimal Tensorial Set of 
    !Symmetry Adapted Bipolar Spherical Harmonic
    implicit none
    real(kind(1d0)), intent(in) :: x
    integer        , intent(in) :: J,N,S
    real(kind(1d0)), intent(out):: A(LDA,*)
    integer        , intent(in) :: LDA
    integer        , intent(in) :: DIM,lv(*)
    !x     : cos(theta_{12}), -1 <= x <= 1
    ! J    : total angular momentum
    ! N    : naturality
    ! S    : parity upon exchange of the angular 
    !        variables:
    !        SABSA(Omega2,Omega1)=(-1)^s SABSA(Omega1,Omega2)
    !        Where SABSA stands for Symmetry Adapted Bipolar
    !        Spherical Harmonic, and is defined as
    !        SABSA_{l_1l_2}^{JMs}(Omega1,Omega2) =
    !          = ( BSA_{l_1l_2}^{JMs}(Omega1,Omega2) +
    !              BSA_{l_1l_2}^{JMs}(Omega2,Omega1) * (-1)^s )/sqrt(2)
    ! A    : superposition matrix
    ! lda  : leading dimension of A
    ! lv   : vector of highest angular momenta of the minimal basis
    ! dim  : dimension of the minimal tensorial set
    !        with angular momentum J, parity l1+l2 and
    !        simmetry upon exchange s
    integer :: i1,i2,j1,j2,j3,j4,nat
    if(J<0.or.abs(x)>1.d0.or.DIM<=0.or.DIM>LDA)return

    !Build the superposition matrix 
    nat=mod(N,2)
    A(1:LDA,1:dim)=0.d0
    do i1=1,dim
       j1=lv(i1)
       j2=J+nat-j1
       do i2=i1,dim
          j3=lv(i2)
          j4=J+nat-j3
          A(i1,i2)=spsabsh(x,J,S,j1,j2,j3,j4)
          A(i2,i1)=A(i1,i2)
       enddo
    enddo

    return
  end subroutine smmtssabsh

  subroutine angcoup2(LANG,N,S,LMAX,NAC,l1v,l2v)
    !Compute the allowed couples of angular momenta
    !compatible with the assigned quantum numbers
    implicit none
    integer, intent(in) :: LANG,N,S,LMAX
    integer, intent(out):: NAC,l1v(*),l2v(*)
    integer :: l1,l2
    NAC=0
    do l1=0,LMAX
       do l2=abs(LANG-l1)+mod(N,2),l1-mod(N+S,2),2
          NAC=NAC+1
          l1v(NAC)=l1
          l2v(NAC)=l2
       enddo
    enddo
    return
  end subroutine angcoup2


  subroutine angcouppm(L12,N12,SGN,LMAX,NCPL,l1v,l2v)
    !Compute the allowed couples of angular momenta
    !compatible with the assigned quantum numbers
    implicit none
    integer, intent(in) :: L12,N12,SGN,LMAX
    integer, intent(out):: NCPL,l1v(*),l2v(*)
    integer :: l1,l2,lambda,P12
    NCPL=0
    P12=mod(L12+N12,2)
    lambda=0
    if(P12==0.and.mod(N12+SGN,2)==1)lambda=1
    !SGN entra solo nel determinare se
    !l'accoppiamento ll e` accettabile o meno
    do l1=0,LMAX
       do l2=abs(L12-l1)+mod(N12,2),l1-lambda,2
          NCPL=NCPL+1
          l1v(NCPL)=l1
          l2v(NCPL)=l2
       enddo
    enddo
    return
  end subroutine angcouppm


  !===================================================================
  !   SUBROUTINES ACCEPTING HALF-INTEGER ANGULAR MOMENTA
  !===================================================================

  !> Here we use the formula number (5) of chapter 8.2.1 of
  !> \emph{Quantum Theory of Angular Momentum} by 
  !> D. A. Varshalovich, A. N. Moskalev, and V. K. Kersonskii
  !> World Scientific, Singapore 
  !> \cite{Varshalovich}
  subroutine WignerClebschGordan(AJ,BJ,CJ,AM,BM,CM,CGC)
    implicit none
    real(kind(1d0)) :: CGC
    real(kind(1d0)) :: AJ, BJ, CJ, AM, BM, CM
    real(kind(1d0)) :: delta, bracket, sum
    real(kind(1d0)) :: int1, int2, int3, int4, int5, int6, int7
    integer         :: z, zmin, zmax

    if(Not_Initialized_Yet)call InitThisModule

    CGC = 0.d0
    if(abs(CM - (AM + BM)) .gt. 1.d-5 )return
    if(AJ<ABS(AM)-1.d-5.or.BJ<ABS(BM)-1.d-5.or.CJ<ABS(AM+BM)-1.d-5) return
    !.. This condition below seems to apply only for integer values of J.
    !IF((ABS(BM)+ABS(BM))==0.and.BTEST(nint(AJ+BJ+CJ),0)) return  
    if((CJ.lt.abs(AJ-BJ)-1.d-5).or.(CJ.gt.AJ+BJ+1.d-5))return
    !.. Check consistency for half integers
    if( mod( nint( (AJ+BJ+CJ) * 2.d0 ), 2 ) == 1 ) return
    !.. Let's consider the integer and half integer cases only
    if(abs(dble(nint(AJ*2.d0))-(AJ*2.d0)).gt.0.00001d0)then !Repeat for rest of the input.
       write(*,*) "consitency of the routine is made only for integer and half integer angular momentum"
       stop
    end if

    !.. 

    !.. The quantities between parenthesis () are always integers,
    !   for any set of integer or half integer angular momentum quantum numbers.
    
    int1 = logfact( nint(   AJ + BJ - CJ     ) )
    int2 = logfact( nint(   AJ - BJ + CJ     ) )
    int3 = logfact( nint( - AJ + BJ + CJ     ) )
    int4 = logfact( nint(   AJ + BJ + CJ + 1 ) )
    
    delta = ( exp( int1 + int2 + int3 - int4 ) )**0.5d0
    !sqrt(fact(j1+j2-j3)*fact(j1-j2+j3)*fact(j2+j3-j1)/fact(j1+j2+j3+1))
    
    int1 = logfact( nint( CJ + CM ) )
    int2 = logfact( nint( CJ - CM ) )
    int3 = logfact( nint( CJ * 2.d0 )+1 )
    int4 = logfact( nint( AJ + AM ) )
    int5 = logfact( nint( AJ - AM ) )
    int6 = logfact( nint( BJ + BM ) )
    int7 = logfact( nint( BJ - BM ) )
    bracket = ( exp( int1 + int2 - int4 - int5 - int6 - int7)*(CJ * 2.d0 +1.d0) )**0.5d0

    zmin = 0
    zmin = max ( 0, nint( CM + BJ - AJ ) )  !(a-b-gamma+z)>=0
    zmax = nint( CJ + BJ + AM )
    zmax = min ( zmax, nint( CJ - AJ + BJ ) )
    zmax = min ( zmax, nint( CJ + CM ) )
    sum = 0.d0
    do z = zmin, zmax
       !
       int1 = logfact( nint( CJ + BJ + AM ) - z )
       int2 = logfact( nint( AJ - AM      ) + z )
       int3 = logfact(                        z )
       int4 = logfact( nint( CJ - AJ + BJ ) - z )
       int5 = logfact( nint( CJ + CM )      - z )
       int6 = logfact( nint( AJ - BJ - CM ) + z )
       sum  = sum + ((-1.d0)**( z*1.d0 )) * exp( int1 + int2 - int3 - int4 - int5 - int6 )
    end do
    sum = sum * ((-1.d0)**( BJ + BM ))

    CGC = delta * bracket * sum
    if(abs(CGC)<1.d-14)CGC=0.d0

    return

  end subroutine WignerClebschGordan


  !> Returns the 6j symbol accepting half integers as arguments (type real):
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    j_a & j_b & j_c \\
  !>    j_d & j_e & j_f
  !>    \end{array}\right\}
  !> \f].
  !> We use the formula (1) of chapter 9.2 of
    !> \emph{Quantum Theory of Angular Momentum} by 
  !> D. A. Varshalovich, A. N. Moskalev, and V. K. Kersonskii
  !> World Scientific, Singapore 
  !> \cite{Varshalovich}
  real(kind(1d0)) function SixJSymbol_dinp(dja,djb,djc,djd,dje,djf)
    !
    implicit none
    !
    real(kind(1d0)), intent(in) :: dja, djb, djc, djd, dje, djf
    integer                     :: n, nmin, nmax
    real(kind(1d0))             :: factor, sum
    !
    if(Not_Initialized_Yet)call InitThisModule
    !
    SixJSymbol_dinp = 0.d0

    !.. Triangle selection rules
    if((dja+djb-djc)*(dja-djb+djc)*(djb+djc-dja).lt.0) return
    if((djd+djb-djf)*(djd-djb+djf)*(djb+djf-djd).lt.0) return
    if((dje+djd-djc)*(dje-djd+djc)*(djd+djc-dje).lt.0) return
    if((dja+dje-djf)*(dja-dje+djf)*(dje+djf-dja).lt.0) return
    !.. Half-integer compatibility
    if(BTEST(nint((dja+djb+djc)*2.d0),0))return
    if(BTEST(nint((djc+djd+dje)*2.d0),0))return
    if(BTEST(nint((dja+dje+djf)*2.d0),0))return
    if(BTEST(nint((djb+djd+djf)*2.d0),0))return
    !
    nmin = max( nint( dja + djb + djc), nint( djc + djd + dje) )
    nmin = max( nint( dja + dje + djf), nmin )
    nmin = max( nint( djb + djd + djf), nmin )
    nmax = min( nint( dja + djb + djd + dje ), nint( dja + djc + djd + djf) )
    nmax = min( nint( djb + djc + dje + djf ), nmax )

    sum = 0.d0
    do n  = nmin, nmax
       !
       factor = logfact( n + 1 ) &
              - logfact(  n - nint( dja + djb + djc ) )&
              - logfact(  n - nint( djc + djd + dje ) )&
              - logfact(  n - nint( dja + dje + djf ) )&
              - logfact(  n - nint( djb + djd + djf ) )&
              - logfact( -n + nint( dja + djb + djd + dje ) )&
              - logfact( -n + nint( dja + djc + djd + djf ) )&
              - logfact( -n + nint( djb + djc + dje + djf ) )
       sum = sum + ( (-1.d0) ** dble(n) ) * exp( factor )
    enddo

    SixJSymbol_dinp = DeltaR(dja,djb,djc)*DeltaR(djb,djd,djf)*&
         DeltaR(dje,djd,djc)*DeltaR(dja,dje,djf)*sum
    !
    if(abs(SixJSymbol_dinp)<1.d-14)SixJSymbol_dinp=0.d0
    !
    return
    !
  contains
    real(kind(1d0)) function DeltaR(dj1,dj2,dj3)
      implicit none
      real(kind(1d0)), intent(in) :: dj1, dj2, dj3
      real(kind(1d0))             :: int1, int2, int3, int4
      if(Not_Initialized_Yet)call InitThisModule
      int1 = logfact( nint(   dj1 + dj2 - dj3     ) )
      int2 = logfact( nint(   dj1 - dj2 + dj3     ) )
      int3 = logfact( nint( - dj1 + dj2 + dj3     ) )
      int4 = logfact( nint(   dj1 + dj2 + dj3 + 1 ) )
      DeltaR= ( exp( int1 + int2 + int3 - int4 ) )**0.5d0
      return
    end function DeltaR
    
  end function SixJSymbol_dinp


  !> Computes the six-j symbol
  !!
  !!  | 1/2  1/2  c |
  !! <               >
  !!  |  d   e    f |
  !!
  real(kind(1d0)) function SixJSymbol_HalfHalf(c,d,e,f) result(res)
    real(kind(1d0)), intent(in) :: c,d,e,f
    res=0.d0
    if(abs(c)<0.1d0)then
       if(abs(d-e)>0.1d0)return
       if(abs(abs(d-f)-0.5d0)>0.1d0)return
       if(abs(f-d-0.5d0)<0.1d0)then
          res=SixJSymbol_a1(d)
       else
          res=SixJSymbol_a2(d)
       endif
    elseif(abs(c-1.d0)<0.1d0)then
       if(abs(d-e)<0.1d0)then
          if(abs(f-d-0.5d0)<0.1d0)then
             res=SixJSymbol_b1(d)
          else
             res=SixJSymbol_b2(d)
          endif
       elseif(abs(abs(d-e)-1.d0)<0.1d0)then
          res = SixJSymbol_c(min(d,e))
       endif
    endif
  end function SixJSymbol_HalfHalf
  

  !> Returns the 6j symbol for the correponding particular case:
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    1/2 & 1/2 &  0  \\
  !>     a  &  a  & a + 1/2 
  !>    \end{array}\right\}    = - (-1)^{2a}/(2\sqrt{a+1/2})
  !> \f].
  !> We use the formula (1) of chapter 9.5.1 of
  !> \emph{Quantum Theory of Angular Momentum} by 
  !> D. A. Varshalovich, A. N. Moskalev, and V. K. Kersonskii
  !> World Scientific, Singapore 
  !> \cite{Varshalovich}
  real(kind(1d0)) function SixJSymbol_a1(dja)
    !
    implicit none
    !
    real(kind(1d0)), intent(in) :: dja
    !
    SixJSymbol_a1 = - 0.5d0 * ( ( -1.d0 ) ** ( 2.d0 * dja ) )/sqrt( dja + 0.5d0 )
    !
    return
    
  end function SixJSymbol_a1
  
  !> Returns the 6j symbol for the correponding particular case:
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    1/2 & 1/2 &  0  \\
  !>     a  &  a  & a - 1/2 
  !>    \end{array}\right\}    =  (-1)^{2a}/(2\sqrt{a+1/2})
  !> \f].
  !> We use the formula (1) of chapter 9.5.1 of
  !> \emph{Quantum Theory of Angular Momentum} by 
  !> D. A. Varshalovich, A. N. Moskalev, and V. K. Kersonskii
  !> World Scientific, Singapore 
  !> \cite{Varshalovich}
  real(kind(1d0)) function SixJSymbol_a2(dja)
    !
    implicit none
    !
    real(kind(1d0)), intent(in) :: dja
    !
    SixJSymbol_a2 = 0.5d0 * ( ( -1.d0 ) ** ( 2.d0 * dja ) )/sqrt( dja + 0.5d0 )
    !
    return
    
  end function SixJSymbol_a2


  !> Returns the 6j symbol for the correponding particular case:
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    1/2 & 1/2 &  1  \\
  !>     a  &  a  & a + 1/2 
  !>    \end{array}\right\}    = - (-1)^{2a}(1/2)\times
  !>    [\frac{a + 1/2 - 1/2}{3(a+1/2)(a+1/2+ 1/2)}]^{1/2}
  !> \f].
  !> We use the formula (1) of chapter 9.5.1 of
  !> \emph{Quantum Theory of Angular Momentum} by 
  !> D. A. Varshalovich, A. N. Moskalev, and V. K. Kersonskii
  !> World Scientific, Singapore 
  !> \cite{Varshalovich}
  real(kind(1d0)) function SixJSymbol_b1(dja)
    !
    implicit none
    !
    real(kind(1d0)), intent(in) :: dja
    real(kind(1d0))             :: fun1, fun2
    !
    fun2 = 3.d0 * ( dja + 0.5d0 ) * ( dja + 1.d0 )
    fun1 = ( dja / fun2 )**0.5d0
    SixJSymbol_b1 = -0.5d0 * ( (-1.d0) ** ( 2.d0 * dja ) ) * fun1
    !
    return
    
  end function SixJSymbol_b1

  !> Returns the 6j symbol for the correponding particular case:
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    1/2 & 1/2 &  1  \\
  !>     a  &  a  & a - 1/2 
  !>    \end{array}\right\}    = - (-1)^{2a}(1/2)\times
  !>    [\frac{a + 1/2 + 1/2}{3(a+1/2)(a+1/2 - 1/2)}]^{1/2}
  !> \f].
  !> We use the formula (1) of chapter 9.5.1 of
  !> \emph{Quantum Theory of Angular Momentum} by 
  !> D. A. Varshalovich, A. N. Moskalev, and V. K. Kersonskii
  !> World Scientific, Singapore 
  !> \cite{Varshalovich}
  real(kind(1d0)) function SixJSymbol_b2(dja)
    !
    implicit none
    !
    real(kind(1d0)), intent(in) :: dja
    real(kind(1d0))             :: fun1, fun2
    !
    fun2 = 3.d0 * ( dja + 0.5d0 ) * ( dja )
    fun1 = ( ( dja + 1.d0 ) / fun2 )**0.5d0
    SixJSymbol_b2 = -0.5d0 * ( (-1.d0) ** ( 2.d0 * dja ) ) * fun1
    !
    return
    
  end function SixJSymbol_b2


  !> Returns the 6j symbol for the correponding particular case:
  !> \f[
  !>    \left\{\begin{array}{ccc}
  !>    1/2 & 1/2 &  1  \\
  !>     a+1  &  a  & a + 1/2 
  !>    \end{array}\right\}    = (-1)^{2a}/\sqrt{3(2 a + 2)}
  !> \f].
  !> We use the formula (1) of chapter 9.5.1 of
  !> \emph{Quantum Theory of Angular Momentum} by 
  !> D. A. Varshalovich, A. N. Moskalev, and V. K. Kersonskii
  !> World Scientific, Singapore 
  !> \cite{Varshalovich}
  real(kind(1d0)) function SixJSymbol_c(dja)
    !
    implicit none
    !
    real(kind(1d0)), intent(in) :: dja
    real(kind(1d0))             :: fun1
    !
    fun1 = sqrt( 3.d0 * ( 2.d0 * dja + 2.d0 ) )
    SixJSymbol_c = ( (-1.d0) ** ( 2.d0 * dja ) ) / fun1
    !
    return
    
  end function SixJSymbol_c


  !===================================================================
  !   EXPLICIT FORMULAS FOR THE CLEBSCH GORDAN COEFFICIENTS
  !===================================================================

  !> Returns C_{S_{B}\Sigma,10}^{S_{A}\Sigma}
  real(kind(1d0)) function CG_SaS_10_SbS(Sb2,Sa2,S2)
    !
    implicit none
    !
    integer, intent(in) :: Sa2, Sb2, S2
    !
    CG_SaS_10_SbS = 0.d0
    !.. Triangle rule
    if((Sa2+Sb2-2)*(Sa2-Sb2+2)*(Sb2-Sa2+2).lt.0) return
    if(Sa2.eq.(Sb2 + 2))then
       CG_SaS_10_SbS =  (dble((Sa2+S2  )*(Sa2-S2)  )/ dble(Sa2    *(2*Sa2-2)))**0.5d0
    endif
    if(Sb2.eq.Sa2)then
       CG_SaS_10_SbS =                        S2    /(dble(Sa2*(Sa2+2))**0.5d0)
    endif
    if(Sa2.eq.(Sb2 - 2))then
       CG_SaS_10_SbS = -(dble((Sa2+S2+2)*(Sa2-S2+2))/ dble((Sa2+3)*(Sa2+2)*2))**0.5d0
    endif

  end function CG_SaS_10_SbS

  !> Returns C_{T\tau,J-\tau}^{K0}
  !.. Tirangle rule is assumed
  real(kind(1d0)) function CG_Tt_Jmt_K0(T2,K2,tau2)
    !
    implicit none
    !
    integer, intent(in) :: T2, K2, tau2
    !
    CG_Tt_Jmt_K0 = 0.d0
    if(abs(tau2).ne.2) then
       write(*,*) "ModuleAngularMomentum: in CG_T1_Jm1_K0, tau2 must be either 2 or -2"
       stop
    endif
    if(K2.eq.(T2 + 2))then
       CG_Tt_Jmt_K0  =  (dble(K2-2)/(2.d0*dble(2*K2-2)))**0.5d0
    endif
    if(K2.eq.T2)then
       CG_Tt_Jmt_K0 =  tau2*0.5d0/(2.d0**0.5d0)
    endif
    if(K2.eq.(T2-2))then
       CG_Tt_Jmt_K0 =  (dble((K2+4)*(K2+2))/dble(4*(K2+3)*(K2+2)))**0.5d0
    endif
  end function CG_Tt_Jmt_K0


  !> Returns C_{T0,J0}^{K0}
  real(kind(1d0)) function CG_T0_J0_K0(T2,K2)
    !
    implicit none
    !
    integer, intent(in) :: T2, K2
    !
    CG_T0_J0_K0 = 0.d0
    if(K2.eq.(T2 + 2))then
      CG_T0_J0_K0  =  (dble(K2)/dble(2*K2-2))**0.5d0
   endif
    if(K2.eq.T2)then
       CG_T0_J0_K0 =  0.d0
    endif
    if(K2.eq.(T2-2))then
       CG_T0_J0_K0 = - (dble(K2+2)/dble(2*(K2+3)))**0.5d0
    endif
  end function CG_T0_J0_K0


  !> Returns C_{\frac{1}{2}\alpha,\frac{1}{2}\beta}^{T(\alpha+\beta)}
  real(kind(1d0)) function CG_12a_12b_Ttau(T2,beta2,tau2)
    !
    implicit none
    !
    integer, intent(in) :: T2, beta2, tau2
    integer             :: alpha2
    !
    if(abs(beta2).ne.1) then
       write(*,*) "ModuleAngularMomentum: in CG_12a_12b_Ttau, alpha2 must be either 1 or -1"
       stop
    endif
    if(abs(tau2).gt.T2) then
       write(*,*) "ModuleAngularMomentum: in CG_12a_12b_Ttau, inconsistent value of tau2 > T2"
       stop
    endif
    alpha2 = tau2 - beta2
    
    CG_12a_12b_Ttau = 0.d0
    select case (T2)
    case (2)
       CG_12a_12b_Ttau  =  (dble(T2+beta2*tau2)/dble(2*T2))**0.5d0
    case (0)
       CG_12a_12b_Ttau  = -beta2*(dble(T2-beta2*tau2+2)/dble(2*T2+4))**0.5d0     
    end select
    return
  end function CG_12A_12B_TTAU

!!$  !.. This routine checks the above formulas for the CG coefficients and compare them with the
!!$  !.. clebsch(AJ,BJ,CJ,AM,BM,CM,CG) routine
!!$  subroutine CheckClebschGordanFormulas( )
!!$
!!$    implicit none
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: fun0,fun1,fun2,fun3,fun4,fun6,fun5,sum1,sum2,sum3
!!$
!!$!    real(kind(1d0)), external  :: CG_SaS_10_SbS,CG_Tt_Jmt_K0
!!$!    real(kind(1d0)), external  ::CG_T0_J0_K0CG_12a_12b_Ttau
!!$    integer :: Sa2, Sb2, S2, tau2, T2, J2, K2
!!$    real(kind(1d0)) :: Sa,Sb,Sigma
!!$    real(kind(1d0)) :: T,K,J,tau
!!$
!!$
!!$
!!$    do Sa = 0.d0, 10.1d0, 0.5d0
!!$       do Sb = 0.d0, Sa + 1.1d0, 1.d0
!!$          do Sigma = -Sb, Sb*1.01d0, 1.d0
!!$             Sa2 = nint(2.d0*Sa)
!!$             Sb2 = nint(2.d0*Sb)
!!$             S2  = nint(2.d0*Sigma)
!!$             call clebsch(Sb,1.d0,Sa,Sigma,0.d0,Sigma,fun1)
!!$             !> Returns C_{S_{B}\Sigma,10}^{S_{A}\Sigma}
!!$             !CG_SaS_10_SbS(Sa2,Sb2,S2)
!!$             fun2 = CG_SaS_10_SbS(Sb2,Sa2,S2)
!!$             write(*,*) "1111111111111111111111111"
!!$             write(*,*) Sa, Sb, Sigma
!!$             write(*,*) fun1,fun2
!!$
!!$             if(abs(fun1-fun2).gt.0.000001d0)then
!!$                pause
!!$             endif
!!$             write(*,*) "1111111111111111111111111"
!!$          enddo
!!$       enddo
!!$    enddo
!!$write(*,*) "+++++++++++++++++++++++++++++++"
!!$    
!!$    do T = 0.d0, 1.1d0, 1.d0
!!$       do J = 0.d0, 1.1d0, 1.d0
!!$          do K = abs(T-J), (T+J)*1.1d0, 1.d0
!!$             do tau= -T, T*1.1d0, 1.d0
!!$                T2 = nint(2.d0*T)
!!$                J2 = nint(2.d0*J)
!!$                K2 = nint(2.d0*K)
!!$                tau2 = nint(2.d0*tau)
!!$                call clebsch(T,J,K,tau,-tau,0.d0,fun1)
!!$                !> Returns C_{T\tau,J-\tau}^{K0}
!!$                fun2 = CG_Tt_Jmt_K0(T2,K2,tau2)
!!$                write(*,*) "22222222222222222222222"
!!$                write(*,*) fun1,fun2
!!$                             if(abs(fun1-fun2).gt.0.000001d0)then
!!$                pause
!!$             endif
!!$                write(*,*) "22222222222222222222222"
!!$             enddo
!!$
!!$                call clebsch(T,J,K,0.d0,0.d0,0.d0,fun1)
!!$                !> Returns C_{T0,J0}^{K0}
!!$                fun2 = CG_T0_J0_K0(T2,K2)
!!$                write(*,*) "33333333333333333333333"
!!$                write(*,*) fun1,fun2
!!$                             if(abs(fun1-fun2).gt.0.000001d0)then
!!$                pause
!!$             endif
!!$                write(*,*) "33333333333333333333333"
!!$          enddo
!!$       enddo
!!$
!!$       do tau= -T, T*1.1d0, 1.d0
!!$          tau2 = nint(2.d0*tau)
!!$          call clebsch(0.5d0,0.5d0,T,0.5d0,-0.5d0,tau,fun1)
!!$          !> Returns C_{\frac{1}{2}\alpha,\frac{1}{2}\beta}^{T(\alpha+\beta)}
!!$          fun2=CG_12a_12b_Ttau(T2,-1,tau2)
!!$          write(*,*) "444444444444444444444444"
!!$          write(*,*) fun1,fun2
!!$                       if(abs(fun1-fun2).gt.0.000001d0)then
!!$                pause
!!$             endif
!!$          write(*,*) "444444444444444444444444"
!!$          call clebsch(0.5d0,0.5d0,T,-0.5d0,0.5d0,tau,fun1)
!!$          !> Returns C_{\frac{1}{2}\alpha,\frac{1}{2}\beta}^{T(\alpha+\beta)}
!!$          fun2=CG_12a_12b_Ttau(T2,1,tau2)
!!$          write(*,*) "55555555555555555555555"
!!$          write(*,*) fun1,fun2
!!$                       if(abs(fun1-fun2).gt.0.000001d0)then
!!$                pause
!!$             endif
!!$          write(*,*) "55555555555555555555555"
!!$       enddo
!!$
!!$    enddo
!!$
!!$    
!!$  end subroutine CheckClebschGordanFormulas
  
  
  !	ARTURO QUIRANTES SIERRA
  !	Department of Applied Physics, Faculty of Sciences
  !	University of Granada, 18071 Granada (SPAIN)
  !	http://www.ugr.es/local/aquiran/codigos.htm
  !	aquiran@ugr.es
  !	Last update: 20 May 2.003
  !	Subroutine NED
  !	to calculate Clebsch-Gordan coefficients
  !	You need to add a "NED(AJ,BJ,CJ,AM,BM,CM,CG)" in your main routine
  !	Input:
  !		AJ,BJ,CJ,AM,BM,CM (the usual Clebsch-Gordan indices)
  !	Output:
  !		CG=C-G(AJ,BJ,CJ,AM,BM,CM)
  !                          (J,J1,J2,M,M1,M2,CG)
  !
  SUBROUTINE clebsch(AJ,BJ,CJ,AM,BM,CM,CG)
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
    DIMENSION Q(500,500)
    DOUBLE PRECISION CG
    DOUBLE PRECISION AJ,BJ,CJ,AM,BM,CM
    INTEGER ZZ
    ZZ=MAX(2*AJ+1,2*BJ+1,2*CJ+1,AJ+BJ+CJ,AJ+AM,BJ+BM,CJ+CM)+2
    DO I=1,ZZ
       Q(I,1)=1.D0
       Q(I,I)=1.D0
    enddo
    DO I=2,ZZ-1
       DO K=2,I
          Q(I+1,K)=Q(I,K-1)+Q(I,K)
       enddo
    enddo
    CG=0.D0
    JA=AJ+AM+1.01D0
    MA=AJ-AM+1.01D0
    JB=BJ+BM+1.01D0
    MB=BJ-BM+1.01D0
    JC=CJ+CM+1.01D0
    MC=CJ-CM+1.01D0
    LA=BJ+CJ-AJ+1.01D0
    LB=CJ+AJ-BJ+1.01D0
    LC=AJ+BJ-CJ+1.01D0
    LT=AJ+BJ+CJ+1.01D0
    D=DABS(AM+BM-CM)-0.01D0
    IF ( D  >  0.d0 ) return
    LD=MIN0(JA,JB,JC,MA,MB,MC,LA,LB,LC)
    IF ( LD <= 0    ) return
    JA2=AJ+AJ+AM+AM
    JB2=BJ+BJ+BM+BM
    JC2=CJ+CJ-CM-CM
    I2=JA2+JB2+JC2-JA2/2*2-JB2/2*2-JC2/2*2
    IF ( I2 /= 0 ) return
    FN=Q(JA+MA-1,LC)/Q(LT,JC+MC-1)
    FN=FN*Q(JB+MB-1,LC)/Q(LT+1,2)
    FN=FN/Q(JA+MA-1,JA)
    FN=FN/Q(JB+MB-1,JB)
    FN=FN/Q(JC+MC-1,JC)
    K0=MAX(0,LC-JA,LC-MB)+1
    K1=MIN(LC,MA,JB)
    X=0.D0
    DO K = K0, K1
       X=-X-Q(LC,K)*Q(LB,MA-K+1)*Q(LA,JB-K+1)
    enddo
    IP=K1+LB+JC
    P=1-2*(IP-IP/2*2)
    CG=P*X*DSQRT(FN)
    !	What weŽve calculated is a Wigner 3-j coefficient
    !	Next, weŽll turn it into a Clebsch-Gordan coefficient
    CG=CG*DSQRT(2*CJ+1)*(-1)**IDNINT(AJ-BJ-CM)
  END SUBROUTINE clebsch


!!$  !===================================================================
!!$  !   TEST SUBROUTINES
!!$  !===================================================================
!!$
!!$  subroutine TestClebschGordanCoefInteger( lmax )
!!$
!!$    implicit none
!!$    integer, intent(in) :: lmax 
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: dVal1, dVal2, dVal3
!!$
!!$    !.. Integer CG
!!$    do l1 = 0, lmax
!!$       do m1 = -l1, l1
!!$          do l2 = 0, lmax
!!$             do m2 = -l2, l2
!!$                do l3 = abs(l1-l2), l1+l2
!!$                   m3 = m1 + m2
!!$                   dVal1 = ClebschGordanCoefficient(l1,l2,l3,m1,m2)
!!$                   call WignerClebschGordan(1.d0*l1,1.d0*l2,1.d0*l3,1.d0*m1,1.d0*m2,1.d0*m3,dVal2)
!!$                   call clebsch(1.d0*l1,1.d0*l2,1.d0*l3,1.d0*m1,1.d0*m2,1.d0*m3,dVal3)
!!$                   if(  abs(dVal2-dVal1) .gt. CGCONTRACTION_THRESHOLD .or. &
!!$                        abs(dVal3-dVal1) .gt. CGCONTRACTION_THRESHOLD ) then
!!$                      write(*,"(*(x,i4))") l1,l2,l3,m1,m2
!!$                      write(*,"(*(x,d24.16))") dVal1
!!$                      write(*,"(*(x,d24.16))") dVal2, abs(dVal2-dVal1)
!!$                      write(*,"(*(x,d24.16))") dVal3, abs(dVal3-dVal1), abs(dVal3-dVal2)
!!$                   end if
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine TestClebschGordanCoefInteger
!!$
!!$
!!$  subroutine TestClebschGordanCoefHalfInteger( lmax )
!!$    implicit none
!!$    integer, intent(in) :: lmax 
!!$    real(kind(1d0))     :: dl1,dl2,dl3,dm1,dm2,dm3
!!$    real(kind(1d0))     :: dVal1, dVal2, dVal3
!!$    !.. Real argument (integer and half integer) CG
!!$    do dl1 = 0.0d0, 1.d0*lmax+0.1d0, 0.5d0
!!$       do dm1 = -dl1, dl1+0.1d0, 1.d0
!!$          do dl2 = 0.d0, 1.d0*lmax+0.1d0, 0.5d0
!!$             do dm2 = -dl2, dl2+0.1d0, 1.d0
!!$                do dl3 = abs(dl1-dl2), dl1+dl2+0.1d0, 0.5d0
!!$                   dm3 = dm1 + dm2
!!$                   call WignerClebschGordan(dl1,dl2,dl3,dm1,dm2,dm3,dVal2)
!!$                   call clebsch(dl1,dl2,dl3,dm1,dm2,dm3,dVal1)
!!$                   if( abs(dVal2-dVal1) .gt. CGCONTRACTION_THRESHOLD )then
!!$                      write(*,"(*(x,f10.1 ))") dl1,dl2,dl3,dm1,dm2,dm3
!!$                      write(*,"(*(x,d24.16))") dVal1
!!$                      write(*,"(*(x,d24.16))") dVal2, abs(dVal1-dVal2)
!!$                   end if
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$  end subroutine TestClebschGordanCoefHalfInteger
!!$
!!$
!!$  !*********** WORKING NOW ******************
!!$  subroutine TestSixJSymbolsHalfInteger( lmax )
!!$    implicit none
!!$    integer, intent(in) :: lmax 
!!$
!!$    real(kind(1d0))     :: dl1,dl2,dl3,dl12,dl23,dl,dlmi,dlma, dm1,dm2,dm3,dm, dm12, dm23
!!$    real(kind(1d0))     :: dVal0, dVal1, dVal2, dVal3, Factor
!!$    real(kind(1d0))     :: dlmax, sumA, sumB, fa1,fa2,fa3,fa4,fb1,fb2,fb3,fb4
!!$    logical             :: AllIntegers
!!$
!!$    dlmax = dble(lmax)+0.1d0
!!$
!!$
!!$     !.. Comparison between the new 6j routine, the old one, and the definition in terms of CG
!!$     !.. Choose arbitrary values of l1,l2 and l3 and combine them to build l throug l12 and l23 in the different schemes
!!$    do dl1 = 0.d0, dlmax, 1.d0
!!$       do dl2 = 0.d0, dlmax, 1.d0
!!$          do dl12 = abs(dl1-dl2),dl1+dl2+0.1d0,1.d0
!!$             do dl3 = 0.d0, dlmax,1.d0
!!$                do dl23 = abs(dl3-dl2),dl3+dl2+0.1d0,1.d0
!!$                   dlmi = max( abs(dl1-dl23), abs(dl3-dl12) )
!!$                   dlma = min(     dl1+dl23 ,     dl3+dl12  ) + 0.1d0
!!$                   Factor = 1.d0 / sqrt((2.d0*dl12+1.d0)*(2.d0*dl23+1.d0))
!!$                   do dl = dlmi, dlma, 1.d0
!!$                      do dm = -dl, dl + 0.1d0, 1.d0
!!$
!!$                         AllIntegers = ( dl1 + dl2 + dl12 + dl3 + dl23 + dl &
!!$                              - int(dl1+0.1d0) - int(dl2+0.1d0)- int(dl12+0.1d0)- int(dl3+0.1d0)- int(dl23+0.1d0)- int(dl+0.1d0) ) < 0.2d0
!!$
!!$                         dVal0 = SixJSymbol_dinp(dl1,dl2,dl12,dl3,dl,dl23)
!!$                         if(AllIntegers) dVal1 = SixJSymbol(nint(dl1),nint(dl2),nint(dl12),nint(dl3),nint(dl),nint(dl23))
!!$
!!$                         sumA = 0.d0
!!$                         sumB = 0.d0
!!$                         do dm1 = -dl1, dl1+0.1d0, 1.d0
!!$                            do dm2 = -dl2, dl2+0.1d0, 1.d0
!!$                               dm12=dm1+dm2
!!$                               do dm3 = -dl3, dl3+0.1d0, 1.d0
!!$                                  dm23=dm2+dm3
!!$                                  !
!!$                                  call WignerClebschGordan(dl12, dl3, dl  , dm12, dm3 , dm   , fa1)
!!$                                  call WignerClebschGordan(dl1 , dl2, dl12, dm1 , dm2 , dm12 , fa2)
!!$                                  call WignerClebschGordan(dl1 , dl23, dl , dm1 , dm23, dm   , fa3)
!!$                                  call WignerClebschGordan(dl2 , dl3, dl23, dm2 , dm3 , dm23 , fa4)
!!$                                  !
!!$                                  call Clebsch(dl12, dl3, dl  , dm12, dm3 , dm   , fb1)
!!$                                  call Clebsch(dl1 , dl2, dl12, dm1 , dm2 , dm12 , fb2)
!!$                                  call Clebsch(dl1 , dl23, dl , dm1 , dm23, dm   , fb3)
!!$                                  call Clebsch(dl2 , dl3, dl23, dm2 , dm3 , dm23 , fb4)
!!$                                  !
!!$                                  suma = suma + fa1 * fa2 * fa3 * fa4
!!$                                  sumb = sumb + fb1 * fb2 * fb3 * fb4
!!$                               enddo
!!$                            enddo
!!$                         enddo
!!$                         if(mod(nint(dl1+dl2+dl3+dl),2)==1)then
!!$                            suma = -suma
!!$                            sumb = -sumb
!!$                         endif
!!$                         dVal2 = suma * Factor
!!$                         dVal3 = sumb * Factor
!!$
!!$
!!$                         if(  abs(dVal2-dVal0) .gt. CGCONTRACTION_THRESHOLD .or. &
!!$                              abs(dVal3-dVal0) .gt. CGCONTRACTION_THRESHOLD .or. &
!!$                              ( (abs(dVal1-dVal0) .gt. CGCONTRACTION_THRESHOLD ) .and. AllIntegers ) ) then
!!$                            write(*,*) "-----"
!!$                            write(*,"(*(x,f8.1))") dl1,dl2,dl3,dl12,dl23,dl,dm
!!$                            write(*,"(*(x,d24.16))") dVal0
!!$                            if(AllIntegers)then
!!$                               write(*,"(*(x,d24.16))") dVal1, abs(dVal1 - dVal0)
!!$                            else
!!$                               write(*,"(a)") "Integer formula not applicable"
!!$                            endif
!!$                            write(*,"(*(x,d24.16))") dVal2, abs(dVal2-dVal0)
!!$                            write(*,"(*(x,d24.16))") dVal3, abs(dVal3-dVal0)
!!$                         end if
!!$                      end do
!!$                   enddo
!!$                enddo
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$  end subroutine TestSixJSymbolsHalfInteger
!!$
!!$
!!$
!!$  
!!$  !.. Check the summation 
!!$  !   \[
!!$  !    - \Pi_{S_A}^{-1}\sum_{\sigma}\sum_{\Sigma_A\Sigma_B\pi\theta}\sum_{\mu\tau\kappa}(-1)^{-\sigma-\pi}
!!$  !    C_{S_A   \Sigma_A, 1/2  \pi    }^{S  \Sigma  }C_{S_B \Sigma_B, 1/2 \theta}^{S \Sigma } 
!!$  !    C_{S_B   \Sigma_B, K    \kappa }^{S_A\Sigma_A}C_{T       \tau, J  -\mu   }^{K \kappa }
!!$  !    C_{1/2  -\pi     , 1/2 -\sigma }^{J-  \mu    }C_{1/2   \theta, 1/2 \sigma}^{T \tau   }
!!$  !    = -(-1)^{S+S_B+1/2} (-1)^{J+K}\times
!!$  !    \Pi_{KTJ}\sjs{S_A}{K}{S_B}{1/2}{S}{1/2} \sjs{J}{T}{K}{1/2}{1/2}{1/2} 
!!$  !   \]
!!$  !
!!$  
!!$  subroutine CheckClebschGordanContractions_1st1B( lmax )
!!$
!!$    implicit none
!!$    integer, intent(in) :: lmax
!!$    real(kind(1d0)) :: dlmax
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: fun0,fun1,fun2,fun3,fun4,fun6,fun5,sum1,sum2,sum3
!!$    real(kind(1d0)) :: dl1,dl2,dl3,dm1,dm2,dm3
!!$    real(kind(1d0)) :: dl,dl12,dl23,dm,dm12,dm23
!!$    real(kind(1d0)) :: Sa,Sia,Sb,Sib,S,Sig,factor1,factor2
!!$    real(kind(1d0)) :: T,ta,K,ka,J,mu,te,pi,si
!!$
!!$    dlmax = dble( lmax ) + 0.1d0
!!$
!!$    do S = 0.d0, dlmax, 0.5d0
!!$       do Sig = -S, S + 0.1d0, 1.d0
!!$          do Sa = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$             do Sb = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$                !
!!$                do K = 0.d0, 2.1d0, 1.d0
!!$                   do T = 0.d0, 1.1d0, 1.d0
!!$                      do J = 0.d0, 1.1d0, 1.d0
!!$                         !
!!$                         sum1 = 0.d0
!!$                         do Sia = -Sa, Sa + 0.1d0, 1.d0
!!$                            do Sib = -Sb, Sb + 0.1d0, 1.d0
!!$                               do pi = -0.5d0, 0.51d0, 1.d0
!!$                                  do te = -0.5d0, 0.51d0, 1.d0
!!$                                     do si = -0.5d0, 0.51d0, 1.d0
!!$                                        do ka = -K, K + 0.1d0, 1.d0
!!$                                           do ta = -T, T + 0.1d0, 1.d0
!!$                                              do mu = -J, J + 0.1d0, 1.d0
!!$                                                 !
!!$                                                 call WignerClebschGordan(Sa   ,0.5d0, S , Sia , pi , Sig , fun1)
!!$                                                 call WignerClebschGordan(Sb   ,0.5d0, S , Sib , te , Sig , fun2)
!!$                                                 call WignerClebschGordan(Sb   ,  K  , Sa, Sib , ka , Sia , fun3)
!!$                                                 call WignerClebschGordan(  T  ,  J  , K , ta  ,-mu , ka  , fun4)
!!$                                                 call WignerClebschGordan(0.5d0,0.5d0, J ,-pi  ,-si ,-mu  , fun5)
!!$                                                 call WignerClebschGordan(0.5d0,0.5d0, T , te  , si , ta  , fun6)
!!$
!!$                                                 sum1 = sum1 + ((-1.d0)**(-si-pi))*fun1*fun2*fun3*fun4*fun5*fun6
!!$                                                 
!!$                                              enddo
!!$                                           enddo
!!$                                        enddo
!!$                                     end do
!!$                                  end do
!!$                               enddo
!!$                            enddo
!!$                         enddo
!!$
!!$                         factor1 = -1.d0/sqrt(2.d0*Sa+1.d0)
!!$                         sum1 = factor1 * sum1
!!$
!!$                         factor2 = -((-1.d0)**(S + Sb + 0.5d0 + J + K))*sqrt((2.d0*K + 1.d0)*(2.d0*T + 1.d0)*(2.d0*J + 1.d0))
!!$                         !factor2 = -factor2
!!$                         sum2 = factor2 * SixJSymbol_dinp( Sa , K , Sb , 0.5d0 ,   S   , 0.5d0 ) &
!!$                                        * SixJSymbol_dinp( J  , T , K  , 0.5d0 , 0.5d0 , 0.5d0 )
!!$                         
!!$                         if(  abs(sum1-sum2) .gt. CGCONTRACTION_THRESHOLD ) then
!!$                            write(*,*) "-----"
!!$                            write(*,"(*(x,f8.1))") S, Sa, Sb, K, T, J
!!$                            write(*,"(*(x,d24.16))") sum1, sum2
!!$                         end if
!!$
!!$                      end do
!!$                   end do
!!$                end do
!!$
!!$             end do
!!$          end do
!!$       enddo
!!$    end do
!!$
!!$
!!$  end subroutine CheckClebschGordanContractions_1st1B
!!$  
!!$
!!$  !.. Check the summation 
!!$  !   \[
!!$  !   -\Pi^{-1}_{S_A}\sum_{\Sigma_A\Sigma_B\pi\theta\tau}(-1)^{\frac{1}{2}-\pi}
!!$  !     C_{S_A \Sigma_A,1/2 \pi}^{S \Sigma} C_{S_B \Sigma_B,1/2\theta}^{S  \Sigma  }
!!$  !     C_{1/2 \theta  ,1/2-\pi}^{T \tau  } C_{S_B \Sigma_B,T  \tau  }^{S_A\Sigma_A}=
!!$  !    =-(-1)^{S_B+1/2+S} \Pi_{T}\sjs{T}{1/2}{1/2}{S}{S_B}{S_A}
!!$  !   \]
!!$  !
!!$  
!!$  subroutine CheckClebschGordanContractions_3rd1B( lmax )
!!$
!!$    implicit none
!!$    integer, intent(in) :: lmax
!!$    real(kind(1d0)) :: dlmax
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: fun0,fun1,fun2,fun3,fun4,fun6,fun5,sum1,sum2,sum3
!!$    real(kind(1d0)) :: dl1,dl2,dl3,dm1,dm2,dm3
!!$    real(kind(1d0)) :: dl,dl12,dl23,dm,dm12,dm23
!!$    real(kind(1d0)) :: Sa,Sia,Sb,Sib,S,Sig,factor1,factor2
!!$    real(kind(1d0)) :: T,ta,K,ka,J,mu,te,pi,si
!!$
!!$    dlmax = dble( lmax ) + 0.1d0
!!$
!!$    do S = 0.d0, dlmax, 0.5d0
!!$       do Sig = -S, S + 0.1d0, 1.d0
!!$          do Sa = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$             do Sb = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$                !
!!$                do T = 0.d0, 1.1d0, 1.d0
!!$                   !
!!$                   sum1 = 0.d0
!!$                   do Sia = -Sa, Sa + 0.1d0, 1.d0
!!$                      do Sib = -Sb, Sb + 0.1d0, 1.d0
!!$                         do pi = -0.5d0, 0.51d0, 1.d0
!!$                            do te = -0.5d0, 0.51d0, 1.d0
!!$                               do ta = -T, T + 0.1d0, 1.d0
!!$                                  !
!!$                                  call WignerClebschGordan(Sa   ,0.5d0, S , Sia , pi , Sig , fun1)
!!$                                  call WignerClebschGordan(Sb   ,0.5d0, S , Sib , te , Sig , fun2)
!!$                                  call WignerClebschGordan(Sb   ,  T  , Sa, Sib , ta , Sia , fun3)
!!$                                  call WignerClebschGordan(0.5d0,0.5d0, T , te  ,-pi , ta  , fun4)
!!$
!!$                                  sum1 = sum1 + ((-1.d0)**(0.5d0-pi))*fun1*fun2*fun3*fun4
!!$
!!$                               enddo
!!$                            enddo
!!$                         enddo
!!$                      end do
!!$                   end do
!!$
!!$                   factor1 = -1.d0/sqrt(2.d0*Sa+1.d0)
!!$                   sum1 = factor1 * sum1
!!$
!!$                   factor2 = -((-1.d0)**(S + Sb + 0.5d0))*sqrt(2.d0*T + 1.d0)
!!$                   sum2 = factor2 * SixJSymbol_dinp( T  , 0.5d0 , 0.5d0, S , Sb ,Sa )
!!$                   if(  abs(sum1-sum2) .gt. CGCONTRACTION_THRESHOLD ) then
!!$                      write(*,*) "-----"
!!$                      write(*,"(*(x,f8.1))") S, Sa, Sb, K, T, J
!!$                      write(*,"(*(x,d24.16))") sum1, sum2
!!$                   end if
!!$                end do
!!$
!!$             end do
!!$          end do
!!$       enddo
!!$    end do
!!$
!!$
!!$  end subroutine CheckClebschGordanContractions_3rd1B
!!$
!!$
!!$
!!$  !.. Check the summation 
!!$  !   \[
!!$  !    \sum_{\Sigma_A\Sigma_B\pi\theta}
!!$  !    C_{S_A \Sigma_A,1/2\pi   }^{S\Sigma}
!!$  !    C_{S_B \Sigma_B,1/2\theta}^{S\Sigma} \times
!!$  !    \delta_{S_A,S_B}\delta_{\Sigma_A,\Sigma_B}\delta_{\theta,\pi} = \delta_{S_A,S_B}
!!$  !   \]
!!$  !
!!$  
!!$  subroutine CheckClebschGordanContractions_1st2B( lmax )
!!$
!!$    implicit none
!!$    integer, intent(in) :: lmax
!!$    real(kind(1d0)) :: dlmax
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: fun0,fun1,fun2,fun3,fun4,fun6,fun5,sum1,sum2,sum3
!!$    real(kind(1d0)) :: dl1,dl2,dl3,dm1,dm2,dm3
!!$    real(kind(1d0)) :: dl,dl12,dl23,dm,dm12,dm23
!!$    real(kind(1d0)) :: Sa,Sia,Sb,Sib,S,Sig,factor1,factor2
!!$    real(kind(1d0)) :: T,ta,K,ka,J,mu,te,pi,si
!!$
!!$    dlmax = dble( lmax ) + 0.1d0
!!$
!!$    do S = 0.d0, dlmax, 0.5d0
!!$       do Sig = -S, S + 0.1d0, 1.d0
!!$          do Sa = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$             do Sb = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$                !
!!$                   sum1 = 0.d0
!!$                   do Sia = -Sa, Sa + 0.1d0, 1.d0
!!$                      do Sib = -Sb, Sb + 0.1d0, 1.d0
!!$                         do pi = -0.5d0, 0.51d0, 1.d0
!!$                            do te = -0.5d0, 0.51d0, 1.d0
!!$                               !
!!$                               if(Sa .ne.Sb )cycle
!!$                               if(Sia.ne.Sib)cycle
!!$                               if(te .ne.pi )cycle
!!$                                  call WignerClebschGordan(Sa   ,0.5d0, S , Sia , pi , Sig , fun1)
!!$                                  call WignerClebschGordan(Sb   ,0.5d0, S , Sib , te , Sig , fun2)
!!$
!!$                                  sum1 = sum1 + fun1*fun2
!!$
!!$                            enddo
!!$                         enddo
!!$                      end do
!!$                   end do
!!$
!!$                   sum2 = 0.d0
!!$                   if(Sa.eq.Sb) sum2 = 1.d0
!!$                   
!!$                   if(  abs(sum1-sum2) .gt. CGCONTRACTION_THRESHOLD ) then
!!$                      write(*,*) "-----"
!!$                      write(*,"(*(x,f8.1))") S, Sa, Sb, K, T, J
!!$                      write(*,"(*(x,d24.16))") sum1, sum2
!!$                   end if
!!$
!!$             end do
!!$          end do
!!$       enddo
!!$    end do
!!$
!!$
!!$  end subroutine CheckClebschGordanContractions_1st2B
!!$
!!$
!!$    !.. Check the summation 
!!$  !   \[
!!$  !   \sum_{\Sigma_A\Sigma_B\pi\tau}\sum_{\rho}(-1)^{1/2-\rho}
!!$  !   C_{S_A \Sigma_A,1 /2 \pi  }^{S   \Sigma  }C_{S_B \Sigma_B,1/2  \pi }^{S \Sigma}
!!$  !   C_{S_B \Sigma_B, T   \tau }^{S_A \Sigma_A}C_{1/2 \rho    ,1/2 -\rho}^{T \tau  }
!!$  !   =\sqrt{2}\delta_{S_AS_B}\delta_{T,0}
!!$  !   \]
!!$  !
!!$  
!!$  subroutine CheckClebschGordanContractions_2nd2B( lmax )
!!$
!!$    implicit none
!!$    integer, intent(in) :: lmax
!!$    real(kind(1d0)) :: dlmax
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: fun0,fun1,fun2,fun3,fun4,fun6,fun5,sum1,sum2,sum3
!!$    real(kind(1d0)) :: dl1,dl2,dl3,dm1,dm2,dm3
!!$    real(kind(1d0)) :: dl,dl12,dl23,dm,dm12,dm23
!!$    real(kind(1d0)) :: Sa,Sia,Sb,Sib,S,Sig,factor1,factor2
!!$    real(kind(1d0)) :: T,ta,K,ka,J,mu,te,pi,si,ro
!!$
!!$    dlmax = dble( lmax ) + 0.1d0
!!$
!!$    do S = 0.d0, dlmax, 0.5d0
!!$       do Sig = -S, S + 0.1d0, 1.d0
!!$          do Sa = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$             do Sb = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$                !
!!$                do T = 0.d0, 1.1d0, 1.d0
!!$                   !
!!$                   sum1 = 0.d0
!!$                   do Sia = -Sa, Sa + 0.1d0, 1.d0
!!$                      do Sib = -Sb, Sb + 0.1d0, 1.d0
!!$                         do pi = -0.5d0, 0.51d0, 1.d0
!!$                            do ta = -T, T + 0.1d0, 1.d0
!!$                               !
!!$                               call WignerClebschGordan(Sa   ,0.5d0, S , Sia , pi , Sig , fun1)
!!$                               call WignerClebschGordan(Sb   ,0.5d0, S , Sib , pi , Sig , fun2)
!!$                               call WignerClebschGordan(Sb   ,  T  , Sa, Sib , ta , Sia , fun3)
!!$
!!$                               fun4 = 0.d0
!!$                               do ro = -0.5d0, 0.51d0, 1.d0
!!$                                  call WignerClebschGordan(0.5d0,0.5d0, T , ro  ,-ro , ta  , fun5)
!!$                                  fun4 = fun4 + ((-1.d0)**(0.5d0-ro))*fun5
!!$                               enddo
!!$
!!$                               sum1 = sum1 + fun1*fun2*fun3*fun4
!!$
!!$                            enddo
!!$                         enddo
!!$                      end do
!!$                   end do
!!$
!!$                   factor1 = 1.d0/sqrt(2.d0*Sa+1.d0)
!!$                   sum1 = factor1 * sum1
!!$
!!$                   sum2 = 0.d0
!!$                   if((Sa.eq.Sb).and.(T.eq.0.d0)) sum2 = sqrt(2.d0)/sqrt(2.d0*Sa+1.d0)
!!$                   !write(*,*) sum1,sum2
!!$                   if(  abs(sum1-sum2) .gt. CGCONTRACTION_THRESHOLD ) then
!!$                      write(*,*) "-----"
!!$                      write(*,"(*(x,f8.1))") S, Sa, Sb, K, T, J
!!$                      write(*,"(*(x,d24.16))") sum1, sum2
!!$                   end if
!!$                end do
!!$
!!$             end do
!!$          end do
!!$       enddo
!!$    end do
!!$
!!$  end subroutine CheckClebschGordanContractions_2nd2B
!!$
!!$
!!$
!!$  !.. Check the summation 
!!$  !   \[
!!$  !    -\Pi^{-1}_{S_A}\sum_{\Sigma_A\Sigma_B\pi\theta\tau}(-1)^{\frac{1}{2}-\pi}
!!$  !    C_{S_A \Sigma_A, 1/2  \pi }^{ S \Sigma} C_{S_B \Sigma_B, 1/2 \theta}^{S   \Sigma  }
!!$  !    C_{1/2 \theta  , 1/2 -\pi }^{ T \tau  } C_{S_B \Sigma_B, T   \tau  }^{S_A \Sigma_A}
!!$  !    =-(-1)^{S+S_B+1/2}\Pi_T\sjs{S_A}{T}{S_B}{1/2}{S}{1/2}
!!$  !   \]
!!$  !
!!$  
!!$  subroutine CheckClebschGordanContractions_3rd2B( lmax )
!!$
!!$    implicit none
!!$    integer, intent(in) :: lmax
!!$    real(kind(1d0)) :: dlmax
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: fun0,fun1,fun2,fun3,fun4,fun6,fun5,sum1,sum2,sum3
!!$    real(kind(1d0)) :: dl1,dl2,dl3,dm1,dm2,dm3
!!$    real(kind(1d0)) :: dl,dl12,dl23,dm,dm12,dm23
!!$    real(kind(1d0)) :: Sa,Sia,Sb,Sib,S,Sig,factor1,factor2
!!$    real(kind(1d0)) :: T,ta,K,ka,J,mu,te,pi,si,ro
!!$
!!$    dlmax = dble( lmax ) + 0.1d0
!!$
!!$    do S = 0.d0, dlmax, 0.5d0
!!$       do Sig = -S, S + 0.1d0, 1.d0
!!$          do Sa = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$             do Sb = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$                !
!!$                do T = 0.d0, 1.1d0, 1.d0
!!$                   !
!!$                   sum1 = 0.d0
!!$                   do Sia = -Sa, Sa + 0.1d0, 1.d0
!!$                      do Sib = -Sb, Sb + 0.1d0, 1.d0
!!$                         do pi = -0.5d0, 0.51d0, 1.d0
!!$                            do te = -0.5d0, 0.51d0, 1.d0
!!$                               do ta = -T, T + 0.1d0, 1.d0
!!$                                  !
!!$                                  call WignerClebschGordan(Sa   ,0.5d0, S , Sia , pi , Sig , fun1)
!!$                                  call WignerClebschGordan(Sb   ,0.5d0, S , Sib , te , Sig , fun2)
!!$                                  call WignerClebschGordan(0.5d0,0.5d0, T , te  ,-pi , ta  , fun3)
!!$                                  call WignerClebschGordan(Sb   ,  T  , Sa, Sib , ta , Sia , fun4)
!!$                                  
!!$
!!$                                  sum1 = sum1 + ((-1.d0)**(0.5d0-pi))*fun1*fun2*fun3*fun4
!!$
!!$                               enddo
!!$                            enddo
!!$                         enddo
!!$                      end do
!!$                   end do
!!$
!!$                   factor1 = -1.d0/sqrt(2.d0*Sa+1.d0)
!!$                   sum1 = factor1 * sum1
!!$
!!$                   factor2 = -((-1.d0)**(S + Sb + 0.5d0))*sqrt(2.d0*T + 1.d0)
!!$                   sum2 = factor2 * SixJSymbol_dinp( Sa, T, Sb , 0.5d0 , S , 0.5d0 )
!!$                   !write(*,*) sum1,sum2
!!$                   if(  abs(sum1-sum2) .gt. CGCONTRACTION_THRESHOLD ) then
!!$                      write(*,*) "-----"
!!$                      write(*,"(*(x,f8.1))") S, Sa, Sb, K, T, J
!!$                      write(*,"(*(x,d24.16))") sum1, sum2
!!$                   end if
!!$                end do
!!$
!!$             end do
!!$          end do
!!$       enddo
!!$    end do
!!$
!!$
!!$  end subroutine CheckClebschGordanContractions_3rd2B
!!$
!!$
!!$  !.. Check the summation 
!!$  !   \[
!!$  !   -\Pi_{S_A}^{-1}\sum_{\mu\tau\kappa\Sigma_A\Sigma_B\pi\theta\rho}(-1)^{-\rho-\pi} 
!!$  !    C_{ S_A \Sigma_A, 1/2 \pi   }^{ S   \Sigma   } C_{ S_B \Sigma_B,1/2  \theta}^{S \Sigma }
!!$  !    C_{ S_B \Sigma_B, K   \kappa}^{ S_A \Sigma_A } C_{ T   \tau    ,J   -\mu   }^{K \kappa }
!!$  !    C_{ 1/2 -\rho   ,1/2  -\pi  }^{ J   -\mu     } C_{ 1/2 \rho    ,1/2 \theta }^{T \tau   }\\
!!$  !    =-(-1)^{S_B+1/2-S} (-1)^{K+T}\Pi_{KJT}\sjs{ S_A }{ K }{ S_B }{ 1/2 }{  S  }{ 1/2 }
!!$  !                                          \sjs{  K  }{ J }{  T  }{ 1/2 }{ 1/2 }{ 1/2 }
!!$  !   \]
!!$  !
!!$  
!!$  subroutine CheckClebschGordanContractions_4th2B( lmax )
!!$
!!$    implicit none
!!$    integer, intent(in) :: lmax
!!$    real(kind(1d0)) :: dlmax
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: fun0,fun1,fun2,fun3,fun4,fun6,fun5,sum1,sum2,sum3
!!$    real(kind(1d0)) :: dl1,dl2,dl3,dm1,dm2,dm3
!!$    real(kind(1d0)) :: dl,dl12,dl23,dm,dm12,dm23
!!$    real(kind(1d0)) :: Sa,Sia,Sb,Sib,S,Sig,factor1,factor2
!!$    real(kind(1d0)) :: T,ta,K,ka,J,mu,te,pi,si,ro
!!$
!!$    dlmax = dble( lmax ) + 0.1d0
!!$
!!$    do S = 0.d0, dlmax, 0.5d0
!!$       do Sig = -S, S + 0.1d0, 1.d0
!!$          do Sa = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$             do Sb = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$                !
!!$                do K = 0.d0, 2.1d0, 1.d0
!!$                   do T = 0.d0, 1.1d0, 1.d0
!!$                      do J = 0.d0, 1.1d0, 1.d0
!!$                         !
!!$                         sum1 = 0.d0
!!$                         do Sia = -Sa, Sa + 0.1d0, 1.d0
!!$                            do Sib = -Sb, Sb + 0.1d0, 1.d0
!!$                               do pi = -0.5d0, 0.51d0, 1.d0
!!$                                  do te = -0.5d0, 0.51d0, 1.d0
!!$                                     do ro = -0.5d0, 0.51d0, 1.d0
!!$                                        !
!!$                                        do ka = -K, K + 0.1d0, 1.d0
!!$                                           do ta = -T, T + 0.1d0, 1.d0
!!$                                              do mu = -J, J + 0.1d0, 1.d0
!!$                                                 !
!!$                                                 call WignerClebschGordan(Sa   ,0.5d0, S , Sia , pi , Sig , fun1)
!!$                                                 call WignerClebschGordan(Sb   ,0.5d0, S , Sib , te , Sig , fun2)
!!$                                                 call WignerClebschGordan(Sb   ,  K  , Sa, Sib , ka , Sia , fun3)
!!$                                                 call WignerClebschGordan(  T  ,  J  , K , ta  ,-mu , ka  , fun4)
!!$                                                 call WignerClebschGordan(0.5d0,0.5d0, J ,-ro  ,-pi ,-mu  , fun5)
!!$                                                 call WignerClebschGordan(0.5d0,0.5d0, T , ro  , te , ta  , fun6)
!!$
!!$                                                 sum1 = sum1 + ((-1.d0)**(-ro-pi))*fun1*fun2*fun3*fun4*fun5*fun6
!!$                                                 
!!$                                              enddo
!!$                                           enddo
!!$                                        enddo
!!$                                     end do
!!$                                  end do
!!$                               enddo
!!$                            enddo
!!$                         enddo
!!$
!!$                         factor1 = -1.d0/sqrt(2.d0*Sa+1.d0)
!!$                         sum1 = factor1 * sum1
!!$
!!$                         factor2 = -((-1.d0)**(S + Sb + 0.5d0 + K + T))*sqrt((2.d0*K + 1.d0)*(2.d0*T + 1.d0)*(2.d0*J + 1.d0))
!!$                         !factor2 = -factor2
!!$                         sum2 = factor2 * SixJSymbol_dinp( Sa , K , Sb , 0.5d0 ,   S   , 0.5d0 ) &
!!$                                        * SixJSymbol_dinp( J  , T , K  , 0.5d0 , 0.5d0 , 0.5d0 )
!!$                         !write(*,*) sum1,sum2
!!$                         if(  abs(sum1-sum2) .gt. CGCONTRACTION_THRESHOLD ) then
!!$                            write(*,*) "-----"
!!$                            write(*,"(*(x,f8.1))") S, Sa, Sb, K, T, J
!!$                            write(*,"(*(x,d24.16))") sum1, sum2
!!$                         end if
!!$
!!$                      end do
!!$                   end do
!!$                end do
!!$
!!$             end do
!!$          end do
!!$       enddo
!!$    end do
!!$
!!$
!!$  end subroutine CheckClebschGordanContractions_4th2B
!!$
!!$
!!$  !.. Check the summation 
!!$  !   \[
!!$  !   ...
!!$  !   \]
!!$  !
!!$  
!!$  subroutine CheckClebschGordanContractions_6th2B( lmax )
!!$
!!$    implicit none
!!$    integer, intent(in) :: lmax
!!$    real(kind(1d0)) :: dlmax
!!$    integer         :: l1,l2,l3,m1,m2,m3
!!$    real(kind(1d0)) :: fun0,fun1,fun2,fun3,fun4,fun5,fun6,fun7,fun8,sum1,sum2,sum3
!!$    real(kind(1d0)) :: dl1,dl2,dl3,dm1,dm2,dm3
!!$    real(kind(1d0)) :: dl,dl12,dl23,dm,dm12,dm23
!!$    real(kind(1d0)) :: Sa,Sia,Sb,Sib,S,Sig,factor1,factor2
!!$    real(kind(1d0)) :: T,ta,K,ka,J,mu,C,ga,D,de,F,fi,te,pi,si,ro,al
!!$
!!$    dlmax = dble( lmax ) + 0.1d0
!!$
!!$    do S = 0.d0, dlmax, 0.5d0
!!$       do Sig = -S, S + 0.1d0, 1.d0
!!$          do Sa = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$             do Sb = abs( S - 0.5d0 ), S + 0.51d0, 1.d0
!!$                !
!!$                do F = 0.d0, 1.1d0, 1.d0
!!$                   do J = abs( F - 0.5d0 ),F + 0.51d0, 1.d0
!!$                      do C = 0.d0, 1.1d0, 1.d0
!!$                         do D = abs( C - 0.5d0 ),C + 0.51d0, 1.d0
!!$                            do K =  abs( J - D ), J + D + 0.1d0, 1.d0
!!$                               !
!!$                               sum1 = 0.d0
!!$                               do Sia = -Sa, Sa + 0.1d0, 1.d0
!!$                                  do Sib = -Sb, Sb + 0.1d0, 1.d0
!!$                                     do pi = -0.5d0, 0.51d0, 1.d0
!!$                                        do te = -0.5d0, 0.51d0, 1.d0
!!$                                           do ro = -0.5d0, 0.51d0, 1.d0
!!$                                              do al = -0.5d0, 0.51d0, 1.d0
!!$                                                 !
!!$
!!$                                                 do fi = -F, F + 0.1d0, 1.d0
!!$                                                    do mu = -J, J + 0.1d0, 1.d0
!!$                                                       do ga = -C, C + 0.1d0, 1.d0
!!$                                                          do de = -D, D + 0.1d0, 1.d0
!!$                                                             do ka = -K, K + 0.1d0, 1.d0
!!$                                                                !
!!$                                                                call WignerClebschGordan(Sa   ,0.5d0, S , Sia , pi , Sig , fun1)
!!$                                                                call WignerClebschGordan(Sb   ,0.5d0, S , Sib , te , Sig , fun2)
!!$                                                                call WignerClebschGordan(0.5d0,0.5d0, C , al  , ro , ga  , fun3)
!!$                                                                call WignerClebschGordan(C    ,0.5d0, D , ga  , te , de  , fun4)
!!$                                                                call WignerClebschGordan(0.5d0,0.5d0, F , -pi ,-ro , fi  , fun5)
!!$                                                                call WignerClebschGordan(F    ,0.5d0, J , fi  ,-al , mu  , fun6)
!!$                                                                call WignerClebschGordan(D    ,   J , K , de  , mu , ka  , fun7)
!!$                                                                call WignerClebschGordan(Sb   ,   K , Sa, Sib , ka , Sia , fun8)
!!$                                                                sum1 = sum1 + ((-1.d0)**(0.5d0-pi-ro-al)) &
!!$                                                                     *fun1*fun2*fun3*fun4*fun5*fun6*fun7*fun8
!!$
!!$                                                             end do
!!$                                                          end do
!!$                                                       end do
!!$                                                    enddo
!!$                                                 enddo
!!$                                              enddo
!!$                                           end do
!!$                                        end do
!!$                                     enddo
!!$                                  enddo
!!$                               enddo
!!$
!!$                               factor1 = 0.5d0/sqrt(2.d0*Sa+1.d0)
!!$                               sum1 = factor1 * sum1
!!$
!!$                               factor2 = 0.5d0*((-1.d0)**(S + Sb + 0.5d0 + D + C+ 0.5d0 + K))*sqrt(&
!!$                                     (2.d0*C + 1.d0)&
!!$                                    *(2.d0*D + 1.d0)&
!!$                                    *(2.d0*F + 1.d0)&
!!$                                    *(2.d0*J + 1.d0)&
!!$                                    *(2.d0*K + 1.d0))
!!$                               sum2 = factor2 * SixJSymbol_dinp( S     , 0.5d0 , Sa , K , Sb    , 0.5d0 ) &
!!$                                              * SixJSymbol_dinp( 0.5d0 , 0.5d0 , K  , J , D     , C     ) &
!!$                                              * SixJSymbol_dinp( C     , 0.5d0 , J  , F , 0.5d0 , 0.5d0 )
!!$                                              !write(*,*) sum1,sum2
!!$                               if(  abs(sum1-sum2) .gt. CGCONTRACTION_THRESHOLD ) then
!!$                                  write(*,*) "-----"
!!$                                  write(*,"(*(x,f8.1))") S, Sa, Sb, K, T, J
!!$                                  write(*,"(*(x,d24.16))") sum1, sum2
!!$                               end if
!!$
!!$                            end do
!!$                         end do
!!$                      end do
!!$                   end do
!!$                end do
!!$
!!$             end do
!!$          end do
!!$       enddo
!!$    end do
!!$
!!$
!!$  end subroutine CheckClebschGordanContractions_6th2B
!!$
!!$  !.. This routine compare the formulas of the 6j symbos for the particular input cases
!!$  !   SixJSymbol_a1, SixJSymbol_a2, SixJSymbol_b1, SixJSymbol_b2, SixJSymbol_c
!!$  !   with the general input routine SixJSymbol_dinp(
!!$  subroutine Check6jSpecialFormulas( lmax )
!!$    
!!$    implicit none
!!$    integer, intent(in) :: lmax
!!$    real(kind(1d0)) :: dlmax, dl
!!$    real(kind(1d0)) :: f0, f1, f2, f3, f4
!!$    real(kind(1d0)) :: sy0, sy1, sy2, sy3, sy4
!!$    real(kind(1d0)) :: dif0, dif1, dif2, dif3, dif4
!!$    real(kind(1d0)) :: factor
!!$
!!$    dlmax = dble( lmax ) + 0.1d0
!!$
!!$    do dl = 0.d0, dlmax, 0.5d0
!!$       !
!!$       f0 = SixJSymbol_a1( dl )
!!$       f1 = SixJSymbol_a2( dl )
!!$       f2 = SixJSymbol_b1( dl )
!!$       f3 = SixJSymbol_b2( dl )
!!$       f4 = SixJSymbol_c(  dl )
!!$
!!$       sy0 = SixJSymbol_dinp( 0.5d0, 0.5d0, 0.d0, dl       , dl, dl + 0.5d0)
!!$       sy1 = SixJSymbol_dinp( 0.5d0, 0.5d0, 0.d0, dl       , dl, dl - 0.5d0)
!!$       sy2 = SixJSymbol_dinp( 0.5d0, 0.5d0, 1.d0, dl       , dl, dl + 0.5d0)
!!$       sy3 = SixJSymbol_dinp( 0.5d0, 0.5d0, 1.d0, dl       , dl, dl - 0.5d0)
!!$       sy4 = SixJSymbol_dinp( 0.5d0, 0.5d0, 1.d0, dl + 1.d0, dl, dl + 0.5d0)
!!$
!!$       dif0 = f0 - sy0
!!$       dif1 = f1 - sy1
!!$       dif2 = f2 - sy2
!!$       dif3 = f3 - sy3
!!$       dif4 = f4 - sy4
!!$       factor = abs(dif0) + abs(dif1) + abs(dif2) + abs(dif3) + abs(dif4)
!!$       
!!$       if(factor .gt. CGCONTRACTION_THRESHOLD )then
!!$          !
!!$          write(*,*) "----------------------", factor, dl
!!$          write(*,*) f0, sy0, dif0
!!$          write(*,*) f1, sy1, dif1
!!$          write(*,*) f2, sy2, dif2
!!$          write(*,*) f3, sy3, dif3
!!$          write(*,*) f4, sy4, dif4
!!$       endif
!!$    enddo
!!$    
!!$       
!!$  end subroutine Check6jSpecialFormulas

end module ModuleAngularMomentum
