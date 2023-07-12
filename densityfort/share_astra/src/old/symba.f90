!Il file symba contiene qualche subroutine per il calcolo dei
!coefficienti e delle funzioni angolari
module symba

  !implicit none
  private
  logical, private :: NOT_INITIALIZED_YET=.TRUE.
  integer        , private, parameter :: DIMFACT=170
  real(kind(1d0)), private :: fact(0:DIMFACT)
  integer        , private, parameter :: DIMLOGFACT=500
  !PGF90 doesn't support quadruple precision!
  !real(kind(1q0)), private :: logfact(0:DIMLOGFACT)
  real(kind(1d0)), private :: logfact(0:DIMLOGFACT)

  private :: initsymba, delta

  public :: cgc,tjsy,sjsy,njsy,sjsyf,subsjsy,subnjsy,tjc,cgchq,cgcsc

  public :: sbangolm,fangol,MOTAF3,MOTSAF3,bispharm
  public :: RAssLegFun, RLegPol, LegPolF, LegPolS
  public :: Plm,Ylm,BSH,SABSH,WRME,Plm_B
  public :: mtssabsh,smmtssabsh,spsabsh,angcoup2,angcouppm

contains

  subroutine initsymba()
    integer :: i
    fact(0)=1.d0
    do i = 1,DIMFACT
       fact(i)=fact(i-1)*dble(i)
    enddo
    logfact(0:DIMFACT)=log(fact(0:DIMFACT))
    do i=DIMFACT+1,DIMLOGFACT
       !PGF90 doesn't support quadruple precision!
       !logfact(i)=logfact(i-1)+log(1.q0*i)
       logfact(i)=logfact(i-1)+log(1.d0*i)
    enddo

    NOT_INITIALIZED_YET=.FALSE.
    return
  end subroutine initsymba

  real(kind(1d0)) function CGC(J1,J2,J3,M1,M2)
    !================================================================
    !  calcolo coefficienti di clebsch-gordan
    !  Questa funzione e` stata testata esaustivamente per confronto
    !  con CLEBSCH
    IMPLICIT REAL*8(A-H,O-Z)

    if(Not_Initialized_Yet)call initsymba

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

  END FUNCTION CGC


  real(kind(1d0)) function CGCHQ(J1,J2,J3,M1,M2)
    !================================================================
    !  calcolo coefficienti di clebsch-gordan
    !  Questa funzione e` stata testata esaustivamente per confronto
    !  con CLEBSCH
    !PGF90 doesn't support quadruple precision!
    !IMPLICIT REAL(kind(1q0))(A-H,O-Z)
    IMPLICIT REAL(kind(1d0))(A-H,O-Z)
    integer, parameter :: MAXANG=50

    if(Not_Initialized_Yet)call initsymba

    cgchq = 0.d0

    if(J1>MAXANG.and.J2>MAXANG.and.J3>MAXANG)then
       cgchq=cgcsc(j1,j2,j3,m1,m2)
       return
    endif

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
    CC=logfact(L1)-logfact(L2)+logfact(IA1)+logfact(L3)-logfact(L4)-logfact(L5)+&
         logfact(IA2)-logfact(L6)+logfact(L7)-logfact(L8)
    IP1=J2+J3+M1-NI
    B1=logfact(IP1)
    IP2=J1-M1+NI
    B2=logfact(IP2)
    IP2=IP2+1
    D1=logfact(NI)
    IR1=NI+1
    IR2=J3-J1+J2-NI
    D2=logfact(IR2)
    IR3=J3+M3-NI
    D3=logfact(IR3)
    IR4=J1-J2-M3+NI
    D4=logfact(IR4)
    IR4=IR4+1
    FAC=1.d0
    IF(BTEST(NI+J2+M2,0)) FAC=-FAC
    !PGF90 doesn't support quadruple precision!
    !CC=FAC*sqrt(dble(L9))*exp(0.5q0*CC+B1+B2-D1-D2-D3-D4)
    CC=FAC*sqrt(dble(L9))*exp(0.5d0*CC+B1+B2-D1-D2-D3-D4)

    N=NM-NI

    S2=1.d0
    IF(N/=0)THEN
       FA=S2
       DO I=1,N
          FA=-FA*IP2*IR2/IP1*IR3/(IR1*IR4)
          S2=S2+FA
          IP1=IP1-1
          IP2=IP2+1
          IR1=IR1+1
          IR2=IR2-1
          IR3=IR3-1
          IR4=IR4+1
       ENDDO
    ENDIF
    CGCHQ=CC*S2
    
    !PGF90 doesn't support quadruple precision!
    !if(abs(CGCHQ)<1.q-15)CGCHQ=0.d0
    if(abs(CGCHQ)<1.d-15)CGCHQ=0.d0

    RETURN

  END FUNCTION CGCHQ

  real(kind(1d0)) function cgcsc(l1,l2,l3,m1_,m2_)
    !Implementa la formula approssimata semiclassica
    !per i coefficienti di Clebsch-Gordan
    ![formula n. 12, pag 265, Varshalovich]
    !(ristretta al caso di soli argomenti interi)
    implicit none
    integer, intent(in) :: l1,l2,l3,m1_,m2_
    
    real(kind(1d0)), parameter :: PIGRECO=3.14159265358979323844d0
    real(kind(1d0)) :: j1,j2,j3,m1,m2,m3,sgn,S,w1,w2,j1q,j2q,j3q
    real(kind(1d0)) :: t1,t2,t3,f1,f2,ct1,ct2,ct3,cf1,cf2

    cgcsc=0.d0
    if(l1<abs(m1).or.l2<abs(m2).or.l3<abs(m1+m2)) return
    if(l3>(l1+l2).or.l3<abs(l1-l2))return
    if((abs(m1)+abs(m2))==0.and.btest(l1+l2+l3,0)) return 
    
    j1=dble(l1)+0.5d0
    j2=dble(l2)+0.5d0
    j3=dble(l3)+0.5d0
    j1q=j1*j1
    j2q=j2*j2
    j3q=j3*j3
    m1=dble(m1_)
    m2=dble(m2_)
    m3=-(m1+m2)

    sgn=1.d0-2.d0*dble(mod(l3+m1_+m2_+1,2))

    S=sqrt(4.d0*(j1q*m2*m3+j2q*m1*m3+j3q*m1*m2)&
           +2.d0*(j1q*j2q+j1q*j3q+j2q*j3q)&
           -j1q*j1q-j2q*j2q-j3q*j3q)/4.d0

    w1=sqrt(dble(2*l3+1)/(2.d0*PIGRECO*S))

    ct1=(2.d0*j1q*m3+m1*(j1q-j2q+j3q))/sqrt((j1q-m1*m1)*(4.d0*j1q*j3q-(j1q-j2q+j3q)**2))
    ct2=(2.d0*j2q*m1+m2*(j2q-j3q+j1q))/sqrt((j2q-m2*m2)*(4.d0*j2q*j1q-(j2q-j3q+j1q)**2))
    ct3=(2.d0*j3q*m1+m2*(j3q-j1q+j2q))/sqrt((j3q-m3*m3)*(4.d0*j3q*j2q-(j3q-j1q+j2q)**2))

    cf1=0.5d0*(j1q-j2q-j3q-2.d0*m2*m3)/sqrt((j2q-m2*m2)*(j3q-m3*m3))
    cf2=0.5d0*(j2q-j3q-j1q-2.d0*m3*m1)/sqrt((j3q-m3*m3)*(j1q-m1*m1))

    t1=acos(ct1)
    t2=acos(ct2)
    t3=acos(ct3)
    f1=acos(cf1)
    f2=acos(cf2)

    w2=cos(j1*t1+j2*t2+j3*t3-m2*f1+m1*f2+PIGRECO/4.d0)

    cgcsc=sgn*w1*w2
    
    return
  end function cgcsc

  real(kind(1d0)) function tjc(j1,j2,j3,m1,m2)
    implicit none
    integer, intent(in) :: j1,j2,j3,m1,m2
    if(Not_Initialized_Yet)call initsymba
    tjc=cgc(j1,j2,j3,m1,m2)/sqrt(dble(2*j3+1))*(1-2*mod(abs(j1-j2+m1+m2),2))
    return
  end function tjc

  real(kind(1d0)) function tjsy(j1,j2,j3,m1,m2)
    implicit none
    integer, intent(in) :: j1,j2,j3,m1,m2
    if(Not_Initialized_Yet)call initsymba
    tjsy=cgc(j1,j2,j3,m1,m2)/sqrt(dble(2*j3+1))*(1-2*mod(abs(j1-j2+m1+m2),2))
    return
  end function tjsy

  real(kind(1d0)) function Racah(ja,jb,je,jd,jc,jf)
    integer, intent(in) :: ja,jb,je,jd,jc,jf
    Racah=sjsy(ja,jb,je,jd,jc,jf)
    if(mod(ja+jb+jc+jd,2)==1)Racah=-Racah
    return
  end function Racah
  
  logical function lsjsy(ja,jb,jc,jd,je,jf)
    !Controlla se sono soddisfatte le condizioni triangolari
    !per il coefficiente 6-j
    !  | JA JB JC |
    ! <            >
    !  | JD JE JF |
    implicit none
    integer, intent(in) :: ja,jb,jc,jd,je,jf
    lsjsy=.false.
    if(ja>jb+jc.or.ja<abs(jb-jc))return
    if(ja>je+jf.or.ja<abs(je-jf))return
    if(jb>jd+jf.or.jb<abs(jd-jf))return
    if(jc>jd+je.or.jc<abs(jd-je))return
    lsjsy=.true.
    return
  end function lsjsy

  real(kind(1d0)) function sjsy(ja,jb,je,jd,jc,jf)
    !SE INVOCATA COME SJSY(JA,JB,JC,JD,JE,JF)
    !CALCOLA PROPRIO IL COEFFICIENTE 6-J:
    !  | JA JB JC |
    ! <            >
    !  | JD JE JF |
    ! Ho controllato in modo esaustivo per confronto con sbwrac
    implicit real*8(a-h,o-z)
    integer, intent(in) :: ja,jb,je,jd,jc,jf

    if(Not_Initialized_Yet)call initsymba

    sjsy = 0.d0
    if((ja+jb-je)*(ja-jb+je)*(jb+je-ja).lt.0) return
    if((jd+jb-jf)*(jd-jb+jf)*(jb+jf-jd).lt.0) return
    if((jc+jd-je)*(jc-jd+je)*(jd+je-jc).lt.0) return
    if((ja+jc-jf)*(ja-jc+jf)*(jc+jf-ja).lt.0) return

    w = 0.d0

    iabe=ja+jb+je
    icde=jc+jd+je
    iacf=ja+jc+jf
    ibdf=jb+jd+jf

    jacdb=ja+jc+jb+jd
    jadef=ja+jd+je+jf
    jbcef=jb+jc+je+jf

    ia=max(iabe,icde,iacf,ibdf)
    ib=min(jacdb,jadef,jbcef)

    iin=ia
    ifin=ib

    sig=1.d0
    sigf=sig

    !if(btest(jacdb,0)) sigf=-sigf
    if(btest(iin,0)) sig=-sig

    aa=0.d0
    do i=iin,ifin
       bb=fact(i-iabe)*fact(i-icde)*fact(i-iacf)*&
            fact(i-ibdf)*fact(jacdb-i)*fact(jadef-i)*fact(jbcef-i)
       aa=aa+sig*fact(i+1)/bb
       sig=-sig
    enddo

    sjsy=aa*delta(ja,jb,je)*delta(jb,jd,jf)*&
         delta(jc,jd,je)*delta(ja,jc,jf)*sigf

    if(abs(sjsy)<1.d-14)sjsy=0.d0

    return

  end function sjsy

  subroutine subsjsy(ja,jb,je,jd,jc,jf,w,lv)
    !CALCOLA PROPRIO IL COEFFICIENTE 6-J:
    !  | JA JB JC |
    ! <            >
    !  | JD JE JF |
    ! Ho controllato in modo esaustivo per confronto con sbwrac
    implicit real*8(a-h,o-z)
    logical :: lv

    if(Not_Initialized_Yet)call initsymba

    lv=.FALSE.
    w = 0.d0
    if((ja+jb-je)*(ja-jb+je)*(jb+je-ja).lt.0) return
    if((jd+jb-jf)*(jd-jb+jf)*(jb+jf-jd).lt.0) return
    if((jc+jd-je)*(jc-jd+je)*(jd+je-jc).lt.0) return
    if((ja+jc-jf)*(ja-jc+jf)*(jc+jf-ja).lt.0) return

    iabe=ja+jb+je
    icde=jc+jd+je
    iacf=ja+jc+jf
    ibdf=jb+jd+jf

    jacdb=ja+jc+jb+jd
    jadef=ja+jd+je+jf
    jbcef=jb+jc+je+jf

    ia=max(iabe,icde,iacf,ibdf)
    ib=min(jacdb,jadef,jbcef)

    iin=ia
    ifin=ib

    sig=1.d0
    sigf=sig

    !if(btest(jacdb,0)) sigf=-sigf
    if(btest(iin,0))   sig =-sig

    aa=0.d0
    do i=iin,ifin
       bb=fact(i-iabe)*fact(i-icde)*fact(i-iacf)*&
            fact(i-ibdf)*fact(jacdb-i)*fact(jadef-i)*fact(jbcef-i)
       aa=aa+sig*fact(i+1)/bb
       sig=-sig
    enddo

    w=aa*delta(ja,jb,je)*delta(jb,jd,jf)*&
         delta(jc,jd,je)*delta(ja,jc,jf)*sigf

    if(abs(w)<1.d-14)then
       w=0.d0
    else
       lv=.TRUE.
    endif

    return

  end subroutine subsjsy

  function delta(j1,j2,j3)
    !CALCOLA LA FUNZIONE
    ![(j1+j2-j3)!*(j1-j2+j3)!*(-j1+j2+j3)!/(j1+j2+j3+1)!]**1/2
    implicit real*8 (a-h,o-z)
    if(Not_Initialized_Yet)call initsymba
    delta=dsqrt(fact(j1+j2-j3)*fact(j1-j2+j3)*fact(j2+j3-j1)/fact(j1+j2+j3+1))
    return
  end function delta

  real(kind(1d0)) function sjsyf(ja,jb,jc,jd,je,jf)
    !CALCOLA IL COEFFICIENTE 6-J "FISICO":
    !  | JA JB JE |
    ! <            >
    !  | JD JC JF |
    ! Ho controllato in modo esaustivo per confronto con sbwrac
    implicit real*8(a-h,o-z)
    if(Not_Initialized_Yet)call initsymba

    sjsyf = 0.d0
    if((ja+jb-je)*(ja-jb+je)*(jb+je-ja).lt.0) return
    if((jd+jb-jf)*(jd-jb+jf)*(jb+jf-jd).lt.0) return
    if((jc+jd-je)*(jc-jd+je)*(jd+je-jc).lt.0) return
    if((ja+jc-jf)*(ja-jc+jf)*(jc+jf-ja).lt.0) return

    w = 0.d0

    iabe=ja+jb+je
    icde=jc+jd+je
    iacf=ja+jc+jf
    ibdf=jb+jd+jf

    jacdb=ja+jc+jb+jd
    jadef=ja+jd+je+jf
    jbcef=jb+jc+je+jf

    ia=max(iabe,icde,iacf,ibdf)
    ib=min(jacdb,jadef,jbcef)

    iin=ia
    ifin=ib

    sig=1.d0
    sigf=sig

    !if(btest(jacdb,0)) sigf=-sigf
    if(btest(iin,0)) sig=-sig

    aa=0.d0
    do i=iin,ifin
       bb=fact(i-iabe)*fact(i-icde)*fact(i-iacf)*&
            fact(i-ibdf)*fact(jacdb-i)*fact(jadef-i)*fact(jbcef-i)
       aa=aa+sig*fact(i+1)/bb
       sig=-sig
    enddo

    sjsyf=aa*delta(ja,jb,je)*delta(jb,jd,jf)*&
         delta(jc,jd,je)*delta(ja,jc,jf)*sigf

    if(abs(sjsyf)<1.d-14)sjsyf=0.d0

    return
  end function sjsyf

  logical function ltjsy(j1,j2,j3)
    !Controlla se e` soddisfatta la condizione triangolare
    !per il coefficiente 3-j
    !    
    !     { j1 j2 j3 }
    !
    implicit none
    integer, intent(in) :: j1,j2,j3
    ltjsy=.FALSE.
    if(ntrico(j1,j2,j3))return
    ltjsy=.TRUE.
    return
  contains

    logical function ntrico(j1,j2,j3)
      integer :: j1,j2,j3
      ntrico = (j1<(abs(j2-j3))).or.(j1>(j2+j3))
      return
    end function ntrico

  end function ltjsy
    

  logical function lnjsy(j11,j12,j13,j21,j22,j23,j31,j32,j33)
    !Controlla se sono soddisfatte le condizioni triangolari
    !per il coefficiente 9-j
    !    
    !     | j11 j12 j13 |
    !    <  j21 j22 j23  >
    !     | j31 j32 j33 |
    !
    implicit none
    integer, intent(in) :: j11,j12,j13,j21,j22,j23,j31,j32,j33
    lnjsy=.false.
    if(ntrico(j11,j12,j13))return
    if(ntrico(j21,j22,j23))return
    if(ntrico(j31,j32,j33))return
    if(ntrico(j11,j21,j31))return
    if(ntrico(j12,j22,j32))return
    if(ntrico(j13,j23,j33))return
    lnjsy=.true.
    return
  contains

    logical function ntrico(j1,j2,j3)
      integer :: j1,j2,j3
      ntrico = (j1<(abs(j2-j3))).or.(j1>(j2+j3))
      return
    end function ntrico

  end function lnjsy


  real(kind(1d0)) function njsy(j11,j12,j13,j21,j22,j23,j31,j32,j33)
    !Calcola i simboli 9j o di Wigner o coefficienti di Fano:
    !    
    !     | j11 j12 j13 |
    !    <  j21 j22 j23  >
    !     | j31 j32 j33 |
    !
    !e sfrutta l'equazione n 5 a pagina 305 del Varshalovich 
    ! (Quantum Theory of angular momentum)
    ! Ho testato questa subroutine in qualche caso particolare
    implicit none
    integer, intent(in) :: j11,j12,j13,j21,j22,j23,j31,j32,j33
    integer :: ell, ellmi, ellma

    if(Not_Initialized_Yet)call initsymba

    njsy=0.d0

    if(ntrico(j11,j12,j13))return
    if(ntrico(j21,j22,j23))return
    if(ntrico(j31,j32,j33))return
    if(ntrico(j11,j21,j31))return
    if(ntrico(j12,j22,j32))return
    if(ntrico(j13,j23,j33))return

    ellmi = max(abs(j11-j33),abs(j32-j21),abs(j23-j12))
    ellma = min(    j11+j33 ,    j32+j21 ,    j23+j12 )
    do ell = ellmi, ellma
       njsy = njsy + dble(2*ell+1) * &
            sjsy(j11,j33,ell,j32,j21,j31) * &
            sjsy(j32,j21,ell,j23,j12,j22) * &
            sjsy(j23,j12,ell,j11,j33,j13)
    enddo

    return

  contains

    logical function ntrico(j1,j2,j3)
      integer :: j1,j2,j3
      ntrico = (j1<(abs(j2-j3))).or.(j1>(j2+j3))
      return
    end function ntrico

  end function njsy


  subroutine subnjsy(j11,j12,j13,j21,j22,j23,j31,j32,j33,w,lv)
    !Calcola i simboli 9j o di Wigner o coefficienti di Fano:
    !    
    !     | j11 j12 j13 |
    !    <  j21 j22 j23  >
    !     | j31 j32 j33 |
    !
    !e sfrutta l'equazione n 5 a pagina 305 del Varshalovich 
    ! (Quantum Theory of angular momentum)
    ! Ho testato questa subroutine in qualche caso particolare
    implicit none
    integer, intent(in) :: j11,j12,j13,j21,j22,j23,j31,j32,j33
    logical        , intent(out):: lv
    real(kind(1d0)), intent(out):: w
    integer :: ell, ellmi, ellma
    if(Not_Initialized_Yet)call initsymba

    w=0.d0
    lv=.FALSE.

    if(ntrico(j11,j12,j13))return
    if(ntrico(j21,j22,j23))return
    if(ntrico(j31,j32,j33))return
    if(ntrico(j11,j21,j31))return
    if(ntrico(j12,j22,j32))return
    if(ntrico(j13,j23,j33))return

    ellmi = max(abs(j11-j33),abs(j32-j21),abs(j23-j12))
    ellma = min(    j11+j33 ,    j32+j21 ,    j23+j12 )
    do ell = ellmi, ellma
       w = w + dble(2*ell+1) * &
            sjsy(j11,j33,ell,j32,j21,j31) * &
            sjsy(j32,j21,ell,j23,j12,j22) * &
            sjsy(j23,j12,ell,j11,j33,j13)
    enddo

    if(abs(w)<1.d-12)then
       w=0.d0
    else
       lv=.TRUE.
    endif

    return

  contains

    logical function ntrico(j1,j2,j3)
      integer :: j1,j2,j3
      ntrico = (j1<(abs(j2-j3))).or.(j1>(j2+j3))
      return
    end function ntrico

  end subroutine subnjsy


  subroutine sbangolm(la1,la2,la3,la12,lg,l1,l2,l3,l12,ll,acoe,lv)

    !CALCOLA I COEFFICIENTI C(L) DERIVANTI DA L'INTEGRAZIONE ANGOLARE
    !PER CIASCUNO DEI POSSIBILI MULTIPOLI L, DATE LE QUADRUPLETTE
    !(la1,la2,la12,la3) E (l1,l2,l12,l3) PER UN ASSEGNATO VALORE
    !DI Ltot. lg. PER I POSSIBILI 18 CASI I COEFFICIENTI SONO
    !CONTENUTI IN acoe(l,i) (i=1-18); IL VETTORE LOGICO*1 lv(i)
    !INDICHERA' SE IL GENERICO CONTRIBUTO iesimo E' DIVERSO DA 0 lv(i)=
    !TRUE
    !O SE E' NULLO lv(i)=FALSE; LA MATRICE incoe(2,18) CONTERRA' IN
    !incoe(1,i) IL VALORE INIZIALE DI l MENTRE incoe(2,i) IL VALORE
    !FINALE
    !ATTENZIONE! SI SUPPONE CHE GLI ACCOPPIAMENTI SIANO TALI DA SODDISFARE
    !SIA LE CONDIZIONI TRIANGOLARI D(la1,la2,la12), D(l1,l2,l12)
    !D(la12,la3,lg), E D(l12,l3,lg) SIA DI PARITA'
    !SPIEGAZIONE PIU' DETTAGLIATA:
    !SE LE FUNZIONI SPAZIALI PSI E PSI' SONO DEFINITE COME:
    !
    !PSI = somme(su m)[somme(su m1)[ C(la12,la3,lg;m,M-m) * C(la1,la2,la12;m1,m-m1) *
    !f[n1,la1,m1](1)*f[n2,la2,m-m1](2)*f[n3,la3,M-m](3)]]
    !
    !PSI' = somme(su m2)[somme(su m3)[ C(l12,l3,lg;m2,M-m2) * C(l1,l2,l12;m3,m2-m3) *
    !F[n1,l1,m3](1)*F[n2,l2,m2-m3](2)*F[n3,l3,M-m2](3)]]
    !
    !IL CALCOLO DEGLI ELEMENTI  ( 'A' E' L'ANTISIMMETRIZATORE )
    !< PSI * THETA | G | A (PSI' * THETA) > 
    !SI RIDUCE A SOMME CON OPPORTUNI COEFFICIENTI DI TERMINI DEL TIPO
    !< PSI | 1/rij | P ( PSI' ) >  dove P e' una permutazione
    !TALI TERMINI SONO 18 (I TRE 1/rij PER LE SEI PERMUTAZIONI).
    !ESPLICITANDO PSI E PSI' SALTAN FUORI SOMME SU QUATTRO INDICI DEL PRODOTTO
    !DI QUATTRO COEFF. DI CG PER UN MONOELETTRONICO PER UN BIELETTRONICO
    !SVILUPPANDO 1/rij NEL BIELETTRONICO IN MULTIPOLI ED ESEGUENDO L'INTEGRAZIONE 
    !ANGOLARE SALTANO FUORI UN TERMINE DELTIFORME (MONOELETTRONICO), ALTRI QUATTRO 
    !COEFF. DI CG, UN RADICALE E IL BIELETTRONICO MULTIPOLARE. 
    !IL TUTTO SI PUO' SCRIVERE COME SOMMA SUI MULTIPOLI DI UN FATTORE RADIALE 
    !MONOELETTRONICO PER UN BIELETTRONICO MULTIPOLARE PER UN COEFFICIENTE CON TUTTO 
    !IL RESTO (IL RADICALE PER LA TRIPLA SOMMATORIA DEL PRODOTTO DI 8 COEFF DI CG). 
    !QUEST'ULTIMO E' PROPRIO IL COEFFICIENTE FORNITO DALLA SBANGOL

    implicit none

    integer        , intent(in) :: la1,la2,la3,la12,lg,l1,l2,l3,l12,ll
    logical        , intent(out):: lv(1:18)
    real(kind(1d0)), intent(out):: acoe(1:18)

    integer :: ki,km,j
    real(kind(1d0)) :: sq12,w1,w2,w3,a1,a2,sig,sigd,aa

    if(Not_Initialized_Yet)call initsymba

    lv=.FALSE.
    acoe=0.d0

    sq12=sqrt(dble((2*la12+1)*(2*l12+1)))

    !CASO 1 (F3|f3) (F1f1|F2f2)(l)
    if(la3==l3.and.la12==l12)then
       if(delparf(la1,l1,la2,l2,ll).and.trif(l2,l1,la12))then
          lv(1)=.TRUE.
          sig=1.d0-2.d0*mod(l12+ll,2)
          a1=sbcsn(la1,l1,ll)
          a2=sbcsn(la2,l2,ll)
          w1=racah(la1,l1,ll,l2,la2,la12)
          acoe(1)=a1*a2*sig*w1
       endif
    endif

    !CASO 7 (F2|f2) (F1f1|F3f3)(l)
    if(la2==l2)then
       if(delparf(la1,l1,la3,l3,ll))then
          lv(7)=.TRUE.
          sig=1.d0-2.d0*mod(la1+la3+la12+l12,2)
          a1=sbcsn(la1,l1,ll)
          a2=sbcsn(la3,l3,ll)
          w1=racah(ll,la1,l1,l2,l12,la12) 
          w2=racah(l3,ll,la3,la12,lg,l12) 
          acoe(7)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 13 (F1|f1) (F2f2|F3f3)(l)
    if(la1==l1)then
       if(delparf(la2,l2,la3,l3,ll).and.tritparf(la2,l2,la3,l3,la12,l12,ll))then
          lv(13)=.TRUE.
          sig=1.d0-2.d0*mod(la12+l2+lg+la3+ll+l3,2)
          a1=sbcsn(la2,l2,ll)
          a2=sbcsn(la3,l3,ll)
          w1=racah(la12,la1,la2,l2,ll,l12)
          w2=racah(la3,l3,ll,l12,la12,lg)
          acoe(13)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 2 (F3|f1) (F1f3|F2f2)(l)
    if(la3==l1)then
       if(delparf(la1,l3,la2,l2,ll).and.trif(l3,la12,l2))then
          lv(2)=.TRUE.
          sig=1.d0-2.d0*mod(la1+la3+la2+l2+l12+ll,2)
          a1=sbcsn(la1,l3,ll)
          a2=sbcsn(la2,l2,ll)
          w1=racah(l3,la1,ll,la2,l2,la12)
          w2=racah(l3,la12,l2,la3,l12,lg)
          acoe(2)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 8 (F2|f2) (F1f3|F3f1)(l)
    if(la2==l2)then
       if(delparf(la1,l3,la3,l1,ll).and.tritf(l1,l3,la1,la3,la2,lg,ki,km))then
          lv(8)=.TRUE.
          sig=1.d0-2.d0*mod(l12+la12+l1+la3+ll,2)
          a1=sbcsn(la1,l3,ll)
          a2=sbcsn(la3,l1,ll)
          aa=0.d0
          do j=ki,km
             w1=racah(la2,l1,l12,l3,lg,j)
             w2=racah(la3,l1,ll,l3,la1,j)
             w3=racah(la2,la1,la12,la3,lg,j)
             aa=aa+w1*w2*w3*(2*j+1)
          enddo
          acoe(8)=a1*a2*sig*aa*sq12
       endif
    endif

    !CASO 14 (F1|f3) (F2f2|F3f1)(l)
    if(la1==l3)then
       if(delparf(la2,l2,la3,l1,ll).and.trif(la3,l12,la2))then
          lv(14)=.TRUE.
          sig=1.d0-2.d0*mod(la1+la12+ll+la3,2)
          a1=sbcsn(la2,l2,ll)
          a2=sbcsn(la3,l1,ll)
          w1=racah(la3,l12,la2,la1,la12,lg)
          w2=racah(la2,l2,ll,l1,la3,l12)
          acoe(14)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 3 (F3|f2) (F1f1|F2f3)(l)
    if(la3==l2)then
       if(delparf(la1,l1,la2,l3,ll).and.trif(la12,l1,l3))then
          lv(3)=.TRUE.
          sig=1.d0-2.d0*mod(la1+l12+lg+la2+l3+ll,2)
          a1=sbcsn(la1,l1,ll)
          a2=sbcsn(la2,l3,ll)
          w1=racah(l1,la1,ll,la2,l3,la12)
          w2=racah(l3,l1,la12,la3,lg,l12)
          acoe(3)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 9 (F2|f3) (F1f1|F3f2)(l)
    if(la2==l3)then
       if(delparf(la1,l1,la3,l2,ll).and.trif(l12,la1,la3))then
          lv(9)=.TRUE.
          sig=1.d0-2.d0*mod(l2+la3+l1+l12+la1+ll,2)
          a1=sbcsn(la1,l1,ll)
          a2=sbcsn(la3,l2,ll)
          w1=racah(la1,l1,ll,l2,la3,l12)
          w2=racah(la3,l12,la1,la2,la12,lg)
          acoe(9)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 15 (F1|f1) (F2f3|F3f2)(l)
    if(la1==l1)then
       if(delparf(la2,l3,la3,l2,ll).and.tritf(l2,l3,la1,lg,la2,la3,ki,km))then
          lv(15)=.TRUE.
          sig=1.d0-2.d0*mod(la2+la3+ll,2)
          a1=sbcsn(la2,l3,ll)
          a2=sbcsn(la3,l2,ll)
          aa=0.d0
          do j=ki,km
             w1=racah(la1,l2,l12,l3,lg,j)
             w2=racah(la3,l2,ll,l3,la2,j)
             w3=racah(la1,la2,la12,la3,lg,j)
             aa=aa+w1*w2*w3*(2*j+1)
          enddo
          acoe(15)=a1*a2*sig*aa*sq12
       endif
    endif

    !CASO 4 (F3|f3) (F1f2|F2f1)(l)
    if(la3==l3.and.la12==l12)then
       if(delparf(la1,l2,la2,l1,ll))then
          lv(4)=.TRUE.
          sig=1.d0-2.d0*mod(la2+la1+ll,2)
          a1=sbcsn(la1,l2,ll)
          a2=sbcsn(la2,l1,ll)
          w1=racah(la1,l2,ll,l1,la2,l12)
          acoe(4)=a1*a2*sig*w1
       endif
    endif

    !CASO 10 (F2|f1) (F1f2|F3f3)(l)
    if(la2==l1)then
       if(delparf(la1,l2,la3,l3,ll).and.tritparf(la1,l2,la3,l3,l12,la12,ll))then
          lv(10)=.TRUE.
          sig=1.d0-2.d0*mod(l1+lg+ll,2)
          a1=sbcsn(la1,l2,ll)
          a2=sbcsn(la3,l3,ll)
          w1=racah(la12,la2,la1,l2,ll,l12)
          w2=racah(la12,l12,ll,l3,la3,lg)
          acoe(10)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 16 (F1|f2) (F2f1|F3f3)(l)
    if(la1==l2)then
       if(delparf(la2,l1,la3,l3,ll).and.tritf(l1,l3,lg,la1,la2,la3,ki,km))then
          lv(16)=.TRUE.
          sig=1.d0-2.d0*mod(l2+lg+ll,2)
          a1=sbcsn(la2,l1,ll)
          a2=sbcsn(la3,l3,ll)
          aa=0.d0
          do j=ki,km
             w1=racah(la1,l12,l1,l3,j,lg)
             w2=racah(la2,l1,ll,l3,la3,j)
             w3=racah(la1,la2,la12,la3,lg,j)
             aa=aa+w1*w2*w3*(2*j+1)
          enddo
          acoe(16)=a1*a2*sig*aa*sq12
       endif
    endif

    !CASO 5 (F3|f1) (F1f2|F2f3)(l)
    if(la3==l1)then
       if(delparf(la1,l2,la2,l3,ll).and.trif(la12,l2,l3))then
          lv(5)=.TRUE.
          sig=1.d0-2.d0*mod(l1+lg+ll,2)
          a1=sbcsn(la1,l2,ll)
          a2=sbcsn(la2,l3,ll)
          w1=racah(l2,la1,ll,la2,l3,la12)
          w2=racah(la3,l2,l12,l3,lg,la12)
          acoe(5)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 11 (F2|f3) (F1f2|F3f1)(l)
    if(la2==l3)then
       if(delparf(la1,l2,la3,l1,ll).and.trif(l12,la1,la3))then
          lv(11)=.TRUE.
          sig=1.d0-2.d0*mod(la3+la1+ll,2)
          a1=sbcsn(la1,l2,ll)
          a2=sbcsn(la3,l1,ll)
          w1=racah(la1,l2,ll,l1,la3,l12)
          w2=racah(la3,l12,la1,la2,la12,lg)
          acoe(11)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 17 (F1|f2) (F2f3|F3f1)(l)
    if(la1==l2)then
       if(delparf(la2,l3,la3,l1,ll).and.tritf(l1,l3,lg,la1,la2,la3,ki,km))then
          lv(17)=.TRUE.
          sig=1.d0-2.d0*mod(la1+la2+la3+lg+ll,2)
          a1=sbcsn(la2,l3,ll)
          a2=sbcsn(la3,l1,ll)
          aa=0.d0
          sigd=1.d0-2.d0*mod(ki,2)
          do j=ki,km
             w1=racah(la1,l12,l1,l3,j,lg)
             w2=racah(la3,l1,ll,l3,la2,j)
             w3=racah(la1,la2,la12,la3,lg,j)
             aa=aa+w1*w2*w3*sigd*(2*j+1)
             sigd=-sigd
          enddo
          acoe(17)=a1*a2*sig*aa*sq12
       endif
    endif

    !CASO 6 (F3|f2) (F1f3|F2f1)(l)
    if(la3==l2)then
       if(delparf(la1,l3,la2,l1,ll).and.trif(la12,l3,l1))then
          lv(6)=.TRUE.
          sig=1.d0-2.d0*mod(la1+la2+ll,2)
          a1=sbcsn(la1,l3,ll)
          a2=sbcsn(la2,l1,ll)
          w1=racah(l3,la1,ll,la2,l1,la12)
          w2=racah(la12,l3,l1,l12,la3,lg)
          acoe(6)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    !CASO 12 (F2|f1) (F1f3|F3f2)(l)
    if(la2==l1)then
       if(delparf(la1,l3,la3,l2,ll).and.tritf(l2,l3,la1,la3,la2,lg,ki,km))then
          lv(12)=.TRUE.
          sig=1.d0-2.d0*mod(la2+la3+la12+ll,2)
          a1=sbcsn(la1,l3,ll)
          a2=sbcsn(la3,l2,ll)
          aa=0.d0
          do j=ki,km
             w1=racah(la2,l2,l12,l3,lg,j)
             w2=racah(la1,l3,ll,l2,la3,j)
             w3=racah(la2,la1,la12,la3,lg,j)
             aa=aa+w1*w2*w3*(2*j+1)
          enddo
          acoe(12)=a1*a2*sig*aa*sq12
       endif
    endif

    !CASO 18 (F1|f3) (F2f1|F3f2)(l)
    if(la1==l3)then
       if(delparf(la2,l1,la3,l2,ll).and.trif(l12,la2,la3))then
          lv(18)=.TRUE.
          sig=1.d0-2.d0*mod(l3+lg+ll,2)
          a1=sbcsn(la2,l1,ll)
          a2=sbcsn(la3,l2,ll)
          w1=racah(la1,la2,la12,la3,lg,l12)
          w2=racah(la2,l1,ll,l2,la3,l12)
          acoe(18)=a1*a2*sig*w1*w2*sq12
       endif
    endif

    return

  end subroutine sbangolm

  logical function fangol(la1,la2,la3,la12,lg,l1,l2,l3,l12,ll,acoe,ncase)
    !CALCOLA I COEFFICIENTI C(L) DERIVANTI DA L'INTEGRAZIONE ANGOLARE
    !PER CIASCUNO DEI POSSIBILI MULTIPOLI L, DATE LE QUADRUPLETTE
    !(la1,la2,la12,la3) E (l1,l2,l12,l3) PER UN ASSEGNATO VALORE
    !DI Ltot. lg. PER IL CASO ncase (1<=ncase<=18). I COEFFICIENTI SONO
    !CONTENUTI IN acoe(l); LA VARIABILE LOGICA*1 lv INDICHERA' SE IL 
    !CONTRIBUTO E' DIVERSO DA 0 lv=TRUE O SE E' NULLO lv=FALSE; 
    !IL VETTORE incoe(2) CONTERRA' IN incoe(1) IL VALORE INIZIALE DI l 
    !MENTRE incoe(2) IL VALORE FINALE
    !ATTENZIONE! SI SUPPONE CHE GLI ACCOPPIAMENTI SIANO TALI DA SODDISFARE
    !SIA LE CONDIZIONI TRIANGOLARI D(la1,la2,la12), D(l1,l2,l12)
    !D(la12,la3,lg), E D(l12,l3,lg) SIA DI PARITA'
    implicit none
    integer        , intent(in) :: la1,la2,la3,la12,lg,l1,l2,l3,l12,ll,ncase
    real(kind(1d0)), intent(out):: acoe

    integer :: ki,km,j
    real(kind(1d0)) :: sq12,w1,w2,w3,a1,a2,sig,sigd,aa

    if(Not_Initialized_Yet)call initsymba
    fangol=.FALSE.
    acoe=0.d0
    sq12=sqrt(dble((2*la12+1)*(2*l12+1)))

    select case(ncase)

    case(1 )!(F3|f3) (F1f1|F2f2)(l)
       if(la3==l3.and.la12==l12)then
          if(delparf(la1,l1,la2,l2,ll).and.trif(l2,l1,la12))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(l12+ll,2)
             a1=sbcsn(la1,l1,ll)
             a2=sbcsn(la2,l2,ll)
             w1=racah(la1,l1,ll,l2,la2,la12)
             acoe=a1*a2*sig*w1
          endif
       endif

    case(7 )!(F2|f2) (F1f1|F3f3)(l)
       if(la2==l2)then
          if(delparf(la1,l1,la3,l3,ll))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la1+la3+la12+l12,2)
             a1=sbcsn(la1,l1,ll)
             a2=sbcsn(la3,l3,ll)
             w1=racah(ll,la1,l1,l2,l12,la12) 
             w2=racah(l3,ll,la3,la12,lg,l12) 
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(13)!(F1|f1) (F2f2|F3f3)(l)
       if(la1==l1)then
          if(delparf(la2,l2,la3,l3,ll).and.tritparf(la2,l2,la3,l3,la12,l12,ll))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la12+l2+lg+la3+ll+l3,2)
             a1=sbcsn(la2,l2,ll)
             a2=sbcsn(la3,l3,ll)
             w1=racah(la12,la1,la2,l2,ll,l12)
             w2=racah(la3,l3,ll,l12,la12,lg)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(2 )!(F3|f1) (F1f3|F2f2)(l)
       if(la3==l1)then
          if(delparf(la1,l3,la2,l2,ll).and.trif(l3,la12,l2))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la1+la3+la2+l2+l12+ll,2)
             a1=sbcsn(la1,l3,ll)
             a2=sbcsn(la2,l2,ll)
             w1=racah(l3,la1,ll,la2,l2,la12)
             w2=racah(l3,la12,l2,la3,l12,lg)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(8 )!(F2|f2) (F1f3|F3f1)(l)
       if(la2==l2)then
          if(delparf(la1,l3,la3,l1,ll).and.tritf(l1,l3,la1,la3,la2,lg,ki,km))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(l12+la12+l1+la3+ll,2)
             a1=sbcsn(la1,l3,ll)
             a2=sbcsn(la3,l1,ll)
             aa=0.d0
             do j=ki,km
                w1=racah(la2,l1,l12,l3,lg,j)
                w2=racah(la3,l1,ll,l3,la1,j)
                w3=racah(la2,la1,la12,la3,lg,j)
                aa=aa+w1*w2*w3*(2*j+1)
             enddo
             acoe=a1*a2*sig*aa*sq12
          endif
       endif

    case(14)!(F1|f3) (F2f2|F3f1)(l)
       if(la1==l3)then
          if(delparf(la2,l2,la3,l1,ll).and.trif(la3,l12,la2))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la1+la12+ll+la3,2)
             a1=sbcsn(la2,l2,ll)
             a2=sbcsn(la3,l1,ll)
             w1=racah(la3,l12,la2,la1,la12,lg)
             w2=racah(la2,l2,ll,l1,la3,l12)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(3 )!(F3|f2) (F1f1|F2f3)(l)
       if(la3==l2)then
          if(delparf(la1,l1,la2,l3,ll).and.trif(la12,l1,l3))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la1+l12+lg+la2+l3+ll,2)
             a1=sbcsn(la1,l1,ll)
             a2=sbcsn(la2,l3,ll)
             w1=racah(l1,la1,ll,la2,l3,la12)
             w2=racah(l3,l1,la12,la3,lg,l12)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(9 )!(F2|f3) (F1f1|F3f2)(l)
       if(la2==l3)then
          if(delparf(la1,l1,la3,l2,ll).and.trif(l12,la1,la3))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(l2+la3+l1+l12+la1+ll,2)
             a1=sbcsn(la1,l1,ll)
             a2=sbcsn(la3,l2,ll)
             w1=racah(la1,l1,ll,l2,la3,l12)
             w2=racah(la3,l12,la1,la2,la12,lg)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(15)!(F1|f1) (F2f3|F3f2)(l)
       if(la1==l1)then
          if(delparf(la2,l3,la3,l2,ll).and.tritf(l2,l3,la1,lg,la2,la3,ki,km))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la2+la3+ll,2)
             a1=sbcsn(la2,l3,ll)
             a2=sbcsn(la3,l2,ll)
             aa=0.d0
             do j=ki,km
                w1=racah(la1,l2,l12,l3,lg,j)
                w2=racah(la3,l2,ll,l3,la2,j)
                w3=racah(la1,la2,la12,la3,lg,j)
                aa=aa+w1*w2*w3*(2*j+1)
             enddo
             acoe=a1*a2*sig*aa*sq12
          endif
       endif

    case(4 )!(F3|f3) (F1f2|F2f1)(l)
       if(la3==l3.and.la12==l12)then
          if(delparf(la1,l2,la2,l1,ll))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la2+la1+ll,2)
             a1=sbcsn(la1,l2,ll)
             a2=sbcsn(la2,l1,ll)
             w1=racah(la1,l2,ll,l1,la2,l12)
             acoe=a1*a2*sig*w1
          endif
       endif

    case(10)!(F2|f1) (F1f2|F3f3)(l)
       if(la2==l1)then
          if(delparf(la1,l2,la3,l3,ll).and.tritparf(la1,l2,la3,l3,l12,la12,ll))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(l1+lg+ll,2)
             a1=sbcsn(la1,l2,ll)
             a2=sbcsn(la3,l3,ll)
             w1=racah(la12,la2,la1,l2,ll,l12)
             w2=racah(la12,l12,ll,l3,la3,lg)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(16)!(F1|f2) (F2f1|F3f3)(l)
       if(la1==l2)then
          if(delparf(la2,l1,la3,l3,ll).and.tritf(l1,l3,lg,la1,la2,la3,ki,km))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(l2+lg+ll,2)
             a1=sbcsn(la2,l1,ll)
             a2=sbcsn(la3,l3,ll)
             aa=0.d0
             do j=ki,km
                w1=racah(la1,l12,l1,l3,j,lg)
                w2=racah(la2,l1,ll,l3,la3,j)
                w3=racah(la1,la2,la12,la3,lg,j)
                aa=aa+w1*w2*w3*(2*j+1)
             enddo
             acoe=a1*a2*sig*aa*sq12
          endif
       endif

    case(5 )!(F3|f1) (F1f2|F2f3)(l)
       if(la3==l1)then
          if(delparf(la1,l2,la2,l3,ll).and.trif(la12,l2,l3))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(l1+lg+ll,2)
             a1=sbcsn(la1,l2,ll)
             a2=sbcsn(la2,l3,ll)
             w1=racah(l2,la1,ll,la2,l3,la12)
             w2=racah(la3,l2,l12,l3,lg,la12)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(11)!(F2|f3) (F1f2|F3f1)(l)
       if(la2==l3)then
          if(delparf(la1,l2,la3,l1,ll).and.trif(l12,la1,la3))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la3+la1+ll,2)
             a1=sbcsn(la1,l2,ll)
             a2=sbcsn(la3,l1,ll)
             w1=racah(la1,l2,ll,l1,la3,l12)
             w2=racah(la3,l12,la1,la2,la12,lg)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(17)!(F1|f2) (F2f3|F3f1)(l)
       if(la1==l2)then
          if(delparf(la2,l3,la3,l1,ll).and.tritf(l1,l3,lg,la1,la2,la3,ki,km))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la1+la2+la3+lg+ll,2)
             a1=sbcsn(la2,l3,ll)
             a2=sbcsn(la3,l1,ll)
             aa=0.d0
             sigd=1.d0-2.d0*mod(ki,2)
             do j=ki,km
                w1=racah(la1,l12,l1,l3,j,lg)
                w2=racah(la3,l1,ll,l3,la2,j)
                w3=racah(la1,la2,la12,la3,lg,j)
                aa=aa+w1*w2*w3*sigd*(2*j+1)
                sigd=-sigd
             enddo
             acoe=a1*a2*sig*aa*sq12
          endif
       endif

    case(6 )!(F3|f2) (F1f3|F2f1)(l)
       if(la3==l2)then
          if(delparf(la1,l3,la2,l1,ll).and.trif(la12,l3,l1))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la1+la2+ll,2)
             a1=sbcsn(la1,l3,ll)
             a2=sbcsn(la2,l1,ll)
             w1=racah(l3,la1,ll,la2,l1,la12)
             w2=racah(la12,l3,l1,l12,la3,lg)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    case(12)!(F2|f1) (F1f3|F3f2)(l)
       if(la2==l1)then
          if(delparf(la1,l3,la3,l2,ll).and.tritf(l2,l3,la1,la3,la2,lg,ki,km))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(la2+la3+la12+ll,2)
             a1=sbcsn(la1,l3,ll)
             a2=sbcsn(la3,l2,ll)
             aa=0.d0
             do j=ki,km
                w1=racah(la2,l2,l12,l3,lg,j)
                w2=racah(la1,l3,ll,l2,la3,j)
                w3=racah(la2,la1,la12,la3,lg,j)
                aa=aa+w1*w2*w3*(2*j+1)
             enddo
             acoe=a1*a2*sig*aa*sq12
          endif
       endif

    case(18)!(F1|f3) (F2f1|F3f2)(l)
       if(la1==l3)then
          if(delparf(la2,l1,la3,l2,ll).and.trif(l12,la2,la3))then
             fangol=.TRUE.
             sig=1.d0-2.d0*mod(l3+lg+ll,2)
             a1=sbcsn(la2,l1,ll)
             a2=sbcsn(la3,l2,ll)
             w1=racah(la1,la2,la12,la3,lg,l12)
             w2=racah(la2,l1,ll,l2,la3,l12)
             acoe=a1*a2*sig*w1*w2*sq12
          endif
       endif

    end select

    return
  end function fangol
 
  subroutine MOTAF3(l1_,l2_,l3_,La,L_,l1,l2,l3,Lb,L,T,P,cv,lv)
    !Monoelectronic Tensor operator Angular Factor for 3 electrons
    !Le formule sono quelle nell'appendice della tesi per il 
    !calcolo degli elementi di matrice, con un'avvertenza: è
    !scorporato il termine \frac{C_{LM,Tt}^{L'M'}}{\sqrt{2L'+1}}.
    !In altri termini viene calcolato il coefficiente angolare
    !necessario a costruire l'elemento di matrice ridotto.
    !Prendiamo ad esempio il primo caso: o(1) E
    !La formula corrispondente, se non nulla, e`:
    !
    ! <[[f1'Xf2']Xf3']|o(1)E|[[f1Xf2]Xf3]>=
    ! =  <f1'||o||f1><f2'|f2><f3'|f3> X
    !  X (-1)^(Na+Nb+N) Pi_{La Lb L L'} X
    !  X \sjsy{La,Lb,T,l1,l1',l2} \sjsy{L,T,L',La,l3,Lb} X
    !  \frac{C_{LM,Tt}^{L'M'}}{\sqrt{2L'+1}}
    !
    ! Quindi, evidentemente:
    !
    ! <[[f1'Xf2']Xf3']||o(1)E||[[f1Xf2]Xf3]>=
    ! =  <f1'||o||f1><f2'|f2><f3'|f3> X
    !  X (-1)^(Na+Nb+N) Pi_{La Lb L L'} X
    !  X \sjsy{La,Lb,T,l1,l1',l2} \sjsy{L,T,L',La,l3,Lb}
    ! 
    ! Questa subroutine calcola solo la parte angolare, ossia
    ! in questo caso il coefficiente
    !
    ! (-1)^(Na+Nb+N) Pi_{La Lb L L'} \sjsy{La,Lb,T,l1,l1',l2} \sjsy{L,T,L',La,l3,Lb}
    !
    ! il flag corrispondente è vero solo se l2=l2' ed l3=l3'

    implicit none
    !Parametri angolari dello stato di bra
    integer        , intent(in) :: l1_,l2_,l3_,La,L_
    !Parametri angolari dello stato di ket
    integer        , intent(in) :: l1,l2,l3,Lb,L
    !Parametri tensoriali dell'operatore di transizione
    integer        , intent(in) :: T,P
    real(kind(1d0)), intent(out):: cv(18)
    logical        , intent(out):: lv(18)
    real(kind(1d0)) :: w1,w2,w3
    real(kind(1d0)) :: pifact
    integer :: Na,N_,Nb,N,No,i
    logical :: lo1,lo2,lo3

    cv = 0.d0
    lv =.false.

    !Conservazione momento angolare totale
    if(L_>L+T.or.L_<abs(L-T))return

    !Conservazione parità
    if(mod(l1_+l2_+l3_+l1+l2+l3+P,2)/=0)return

    Na =l1_+l2_+La
    N_ =l1_+l2_+l3_+L_
    Nb =l1+l2+Lb
    N  =l1+l2+l3+L
    No =T+P

    pifact=sqrt(dble((2*La+1)*(2*Lb+1)*(2*L+1)*(2*L_+1)))

    !caso  1: o(1) E
    i=1
    lo1 = l2_==l2 .and. l3_==l3
    lo2 = lsjsy(La,Lb,T,l1,l1_,l2)
    lo3 = lsjsy(L,T,L_,La,l3,Lb)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(Na+Nb+N)
       w2=sjsy(La,Lb,T,l1,l1_,l2)
       w3=sjsy(L,T,L_,La,l3,Lb)
       cv(i)=w1*w2*w3
    endif

    !caso  2: o(1) P13 
    i=2
    lo1 = l2_==l2 .and. l3_==l1
    lo2 = lsjsy(l1,l2,Lb,l1_,L_,La)
    lo3 = lsjsy(L,T,L_,l1_,Lb,l3)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(Nb+No+N_)
       w2=sjsy(l1,l2,Lb,l1_,L_,La)
       w3=sjsy(L,T,L_,l1_,Lb,l3)
       cv(i)=w1*w2*w3
    endif

    !caso  3: o(1) P23 
    i=3
    lo1 = l2_==l3 .and. l3_==l2
    lo2 = lnjsy(L,T,L_,l3,l1_,La,Lb,l1,l2)
    lv(i)=lo1.and.lo2
    if(lv(i))then
       w1=pifact*mypar(Na+Nb+N)
       w2=njsy(L,T,L_,l3,l1_,La,Lb,l1,l2)
       cv(i)=w1*w2
    endif

    !caso  4: o(1) P12
    i=4
    lo1 = l2_==l1 .and. l3_==l3
    lo2 = lsjsy(La,T,Lb,l2,l1,l1_)
    lo3 = lsjsy(L,T,L_,La,l3,Lb)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(Na+N)
       w2=sjsy(La,T,Lb,l2,l1,l1_)
       w3=sjsy(L,T,L_,La,l3,Lb)
       cv(i)=w1*w2*w3
    endif

    !caso  5: o(1) P123 
    i=5
    lo1 = l2_==l3 .and. l3_==l1
    lo2 = lnjsy(L,T,L_,l3,l1_,La,Lb,l2,l1)
    lv(i)=lo1.and.lo2
    if(lv(i))then
       w1=pifact*mypar(Na+N)
       w2=njsy(L,T,L_,l3,l1_,La,Lb,l2,l1)
       cv(i)=w1*w2
    endif

    !caso  6: o(1) P321 
    i=6
    lo1 = l2_==l1 .and. l3_==l2
    lo2 = lsjsy(La,l1,l1_,Lb,L_,l2)
    lo3 = lsjsy(L,T,L_,l1_,Lb,l3)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(No+N_)
       w2=sjsy(La,l1,l1_,Lb,L_,l2)
       w3=sjsy(L,T,L_,l1_,Lb,l3)
       cv(i)=w1*w2*w3
    endif

    !caso  7: o(2) E 
    i=7
    lo1 = l1_==l1 .and. l3_==l3
    lo2 = lsjsy(La,Lb,T,l2,l2_,l1)
    lo3 = lsjsy(L,T,L_,La,l3,Lb)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(N)
       w2=sjsy(La,Lb,T,l2,l2_,l1)
       w3=sjsy(L,T,L_,La,l3,Lb)
       cv(i)=w1*w2*w3
    endif

    !caso  8: o(2) P13 
    i=8
    lo1 = l1_==l3 .and. l3_==l1
    lo2 = lnjsy(L,T,L_,l3,l2_,La,Lb,l2,l1)
    lv(i)=lo1.and.lo2
    if(lv(i))then
       w1=pifact*mypar(N)
       w2=njsy(L,T,L_,l3,l2_,La,Lb,l2,l1)
       cv(i)=w1*w2
    endif

    !caso  9: o(2) P23 
    i=9
    lo1 = l1_==l1 .and. l3_==l2
    lo2 = lsjsy(l1,l2,Lb,L_,l2_,La)
    lo3 = lsjsy(L,T,L_,l2_,Lb,l3)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(N_+No+Na)
       w2=sjsy(l1,l2,Lb,L_,l2_,La)
       w3=sjsy(L,T,L_,l2_,Lb,l3)
       cv(i)=w1*w2*w3
    endif

    !caso 10: o(2) P12 
    i=10
    lo1 = l1_==l2 .and. l3_==l3
    lo2 = lsjsy(Lb,T,La,l2_,l2,l1)
    lo3 = lsjsy(L,T,L_,La,l3,Lb)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(N+Nb)
       w2=sjsy(Lb,T,La,l2_,l2,l1)
       w3=sjsy(L,T,L_,La,l3,Lb)
       cv(i)=w1*w2*w3
    endif

    !caso 11: o(2) P123 
    i=11
    lo1 = l1_==l2 .and. l3_==l1
    lo2 = lsjsy(l1,l2,Lb,l2_,L_,La)
    lo3 = lsjsy(L,T,L_,l2_,Lb,l3)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(Na+Nb+No+N_)
       w2=sjsy(l1,l2,Lb,l2_,L_,La)
       w3=sjsy(L,T,L_,l2_,Lb,l3)
       cv(i)=w1*w2*w3
    endif

    !caso 12: o(2) P321 
    i=12
    lo1 = l1_==l3 .and. l3_==l2
    lo2 = lnjsy(L,T,L_,l3,l2_,La,Lb,l1,l2)
    lv(i)=lo1.and.lo2
    if(lv(i))then
       w1=pifact*mypar(N+Nb)
       w2=njsy(L,T,L_,l3,l2_,La,Lb,l1,l2)
       cv(i)=w1*w2
    endif

    !caso 13: o(3) E 
    i=13
    lo1 = l1_==l1 .and. l2_==l2 .and. La==Lb
    lo2 = ltjsy(l1,l2,Lb)
    lo3 = lsjsy(L,T,L_,l3_,Lb,l3)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(N_+Na+No)/dble(2*La+1)
       w2=sjsy(L,T,L_,l3_,Lb,l3)
       cv(i)=w1*w2
    endif

    !caso 14: o(3) P13
    i=14
    lo1 = l1_==l3 .and. l2_==l2
    lo2 = lsjsy(l1,l2,Lb,l3,L,La)
    lo3 = lsjsy(L,T,L_,l3_,La,l1)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(Na+N_+No)
       w2=sjsy(l1,l2,Lb,l3,L,La)
       w3=sjsy(L,T,L_,l3_,La,l1)
       cv(i)=w1*w2*w3
    endif

    !caso 15: o(3) P23 
    i=15
    lo1 = l1_==l1 .and. l2_==l3
    lo2 = lsjsy(l1,l2,Lb,L,l3,La)
    lo3 = lsjsy(L,T,L_,l3_,La,l2)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(No+N_+Nb)
       w2=sjsy(l1,l2,Lb,L,l3,La)
       w3=sjsy(L,T,L_,l3_,La,l2)
       cv(i)=w1*w2*w3
    endif

    !caso 16: o(3) P12 
    i=16
    lo1 = l1_==l2 .and. l2_==l1 .and. La==Lb
    lo2 = ltjsy(l1,l2,Lb)
    lo3 = lsjsy(L,T,L_,l3_,Lb,l3)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(N_)/dble(2*La+1)
       w2=sjsy(L,T,L_,l3_,Lb,l3)
       cv(i)=w1*w2
    endif

    !caso 17: o(3) P123 
    i=17
    lo1 = l1_==l2 .and. l2_==l3
    lo2 = lsjsy(l1,l2,Lb,l3,L,La)
    lo3 = lsjsy(L,T,L_,l3_,La,l1)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(No+N_)
       w2=sjsy(l1,l2,Lb,l3,L,La)
       w3=sjsy(L,T,L_,l3_,La,l1)
       cv(i)=w1*w2*w3
    endif

    !caso 18: o(3) P321 
    i=18
    lo1 = l1_==l3 .and. l2_==l1
    lo2 = lsjsy(l1,l2,Lb,L,l3,La)
    lo3 = lsjsy(L,T,L_,l3_,La,l2)
    lv(i)=lo1.and.lo2.and.lo3
    if(lv(i))then
       w1=pifact*mypar(Na+Nb+N_+No)
       w2=sjsy(l1,l2,Lb,L,l3,La)
       w3=sjsy(L,T,L_,l3_,La,l2)
       cv(i)=w1*w2*w3
    endif

    return

  end subroutine MOTAF3


  subroutine MOTSAF3(l1_,l2_,l3_,La,L_,l1,l2,l3,Lb,L,T,P,cv,lv,i)
    !Monoelectronic Tensor operator Single Angular Factor for 3 electrons
    !Le formule sono quelle nell'appendice della tesi per il 
    !calcolo degli elementi di matrice, con un'avvertenza: è
    !scorporato il termine \frac{C_{LM,Tt}^{L'M'}}{\sqrt{2L'+1}}.
    !In altri termini viene calcolato il coefficiente angolare
    !necessario a costruire l'elemento di matrice ridotto.
    !Prendiamo ad esempio il primo caso: o(1) E
    !La formula corrispondente, se non nulla, e`:
    !
    ! <[[f1'Xf2']Xf3']|o(1)E|[[f1Xf2]Xf3]>=
    ! =  <f1'||o||f1><f2'|f2><f3'|f3> X
    !  X (-1)^(Na+Nb+N) Pi_{La Lb L L'} X
    !  X \sjsy{La,Lb,T,l1,l1',l2} \sjsy{L,T,L',La,l3,Lb} X
    !  \frac{C_{LM,Tt}^{L'M'}}{\sqrt{2L'+1}}
    !
    ! Quindi, evidentemente:
    !
    ! <[[f1'Xf2']Xf3']||o(1)E||[[f1Xf2]Xf3]>=
    ! =  <f1'||o||f1><f2'|f2><f3'|f3> X
    !  X (-1)^(Na+Nb+N) Pi_{La Lb L L'} X
    !  X \sjsy{La,Lb,T,l1,l1',l2} \sjsy{L,T,L',La,l3,Lb}
    ! 
    ! Questa subroutine calcola solo la parte angolare, ossia
    ! in questo caso il coefficiente
    !
    ! (-1)^(Na+Nb+N) Pi_{La Lb L L'} \sjsy{La,Lb,T,l1,l1',l2} \sjsy{L,T,L',La,l3,Lb}
    !
    ! il flag corrispondente è vero solo se l2=l2' ed l3=l3'

    implicit none
    !Parametri angolari dello stato di bra
    integer        , intent(in) :: l1_,l2_,l3_,La,L_
    !Parametri angolari dello stato di ket
    integer        , intent(in) :: l1,l2,l3,Lb,L
    !Parametri tensoriali dell'operatore di transizione
    integer        , intent(in) :: T,P
    real(kind(1d0)), intent(out):: cv
    logical        , intent(out):: lv
    integer        , intent(in) :: i
    real(kind(1d0)) :: w1,w2,w3
    real(kind(1d0)) :: pifact
    integer :: Na,N_,Nb,N,No
    logical :: lo1,lo2,lo3

    cv = 0.d0
    lv =.false.

    !Conservazione momento angolare totale
    if(L_>L+T.or.L_<abs(L-T))return

    !Conservazione parità
    if(mod(l1_+l2_+l3_+l1+l2+l3+P,2)/=0)return

    !Calcola la "Naturalità" dei vari sottosistemi
    Na =l1_+l2_+La
    N_ =l1_+l2_+l3_+L_
    Nb =l1+l2+Lb
    N  =l1+l2+l3+L
    No =T+P

    pifact=sqrt(dble((2*La+1)*(2*Lb+1)*(2*L+1)*(2*L_+1)))

    select case(i)

    case(1)! o(1) E
       lo1=l2_==l2.and.l3_==l3
       lo2=lsjsy(La,Lb,T,l1,l1_,l2)
       lo3=lsjsy(L,T,L_,La,l3,Lb)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(Na+Nb+N)
          w2=sjsy(La,Lb,T,l1,l1_,l2)
          w3=sjsy(L,T,L_,La,l3,Lb)
          cv=w1*w2*w3
       endif

    case(2)! o(1) P13 
       lo1=l2_==l2.and.l3_==l1
       lo2=lsjsy(l1,l2,Lb,l1_,L_,La)
       lo3=lsjsy(L,T,L_,l1_,Lb,l3)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(Nb+No+N_)
          w2=sjsy(l1,l2,Lb,l1_,L_,La)
          w3=sjsy(L,T,L_,l1_,Lb,l3)
          cv=w1*w2*w3
       endif

    case(3)! o(1) P23 
       lo1=l2_==l3.and.l3_==l2
       lo2=lnjsy(L,T,L_,l3,l1_,La,Lb,l1,l2)
       lv=lo1.and.lo2
       if(lv)then
          w1=pifact*mypar(Na+Nb+N)
          w2=njsy(L,T,L_,l3,l1_,La,Lb,l1,l2)
          cv=w1*w2
       endif

    case(4)! o(1) P12
       lo1=l2_==l1.and.l3_==l3
       lo2=lsjsy(La,T,Lb,l2,l1,l1_)
       lo3=lsjsy(L,T,L_,La,l3,Lb)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(Na+N)
          w2=sjsy(La,T,Lb,l2,l1,l1_)
          w3=sjsy(L,T,L_,La,l3,Lb)
          cv=w1*w2*w3
       endif

    case(5)! o(1) P123 
       lo1=l2_==l3.and.l3_==l1
       lo2=lnjsy(L,T,L_,l3,l1_,La,Lb,l2,l1)
       lv=lo1.and.lo2
       if(lv)then
          w1=pifact*mypar(Na+N)
          w2=njsy(L,T,L_,l3,l1_,La,Lb,l2,l1)
          cv=w1*w2
       endif

    case(6)! o(1) P321 
       lo1=l2_==l1.and.l3_==l2
       lo2=lsjsy(La,l1,l1_,Lb,L_,l2)
       lo3=lsjsy(L,T,L_,l1_,Lb,l3)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(No+N_)
          w2=sjsy(La,l1,l1_,Lb,L_,l2)
          w3=sjsy(L,T,L_,l1_,Lb,l3)
          cv=w1*w2*w3
       endif

    case(7)! o(2) E 
       lo1=l1_==l1.and.l3_==l3
       lo2=lsjsy(La,Lb,T,l2,l2_,l1)
       lo3=lsjsy(L,T,L_,La,l3,Lb)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(N)
          w2=sjsy(La,Lb,T,l2,l2_,l1)
          w3=sjsy(L,T,L_,La,l3,Lb)
          cv=w1*w2*w3
       endif

    case(8)! o(2) P13 
       lo1=l1_==l3.and.l3_==l1
       lo2=lnjsy(L,T,L_,l3,l2_,La,Lb,l2,l1)
       lv=lo1.and.lo2
       if(lv)then
          w1=pifact*mypar(N)
          w2=njsy(L,T,L_,l3,l2_,La,Lb,l2,l1)
          cv=w1*w2
       endif

    case(9)! o(2) P23 
       lo1=l1_==l1.and.l3_==l2
       lo2=lsjsy(l1,l2,Lb,L_,l2_,La)
       lo3=lsjsy(L,T,L_,l2_,Lb,l3)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(N_+No+Na)
          w2=sjsy(l1,l2,Lb,L_,l2_,La)
          w3=sjsy(L,T,L_,l2_,Lb,l3)
          cv=w1*w2*w3
       endif

    case(10)! o(2) P12 
       lo1=l1_==l2.and.l3_==l3
       lo2=lsjsy(Lb,T,La,l2_,l2,l1)
       lo3=lsjsy(L,T,L_,La,l3,Lb)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(N+Nb)
          w2=sjsy(Lb,T,La,l2_,l2,l1)
          w3=sjsy(L,T,L_,La,l3,Lb)
          cv=w1*w2*w3
       endif

    case(11)! o(2) P123 
       lo1=l1_==l2.and.l3_==l1
       lo2=lsjsy(l1,l2,Lb,l2_,L_,La)
       lo3=lsjsy(L,T,L_,l2_,Lb,l3)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(Na+Nb+No+N_)
          w2=sjsy(l1,l2,Lb,l2_,L_,La)
          w3=sjsy(L,T,L_,l2_,Lb,l3)
          cv=w1*w2*w3
       endif

    case(12)! o(2) P321 
       lo1=l1_==l3.and.l3_==l2
       lo2=lnjsy(L,T,L_,l3,l2_,La,Lb,l1,l2)
       lv=lo1.and.lo2
       if(lv)then
          w1=pifact*mypar(N+Nb)
          w2=njsy(L,T,L_,l3,l2_,La,Lb,l1,l2)
          cv=w1*w2
       endif

    case(13)! o(3) E 
       lo1=l1_==l1.and.l2_==l2.and.La==Lb
       lo2=ltjsy(l1,l2,Lb)
       lo3=lsjsy(L,T,L_,l3_,Lb,l3)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(N_+Na+No)/dble(2*La+1)
          w2=sjsy(L,T,L_,l3_,Lb,l3)
          cv=w1*w2
       endif

    case(14)! o(3) P13
       lo1=l1_==l3.and.l2_==l2
       lo2=lsjsy(l1,l2,Lb,l3,L,La)
       lo3=lsjsy(L,T,L_,l3_,La,l1)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(Na+N_+No)
          w2=sjsy(l1,l2,Lb,l3,L,La)
          w3=sjsy(L,T,L_,l3_,La,l1)
          cv=w1*w2*w3
       endif

    case(15)! o(3) P23 
       lo1=l1_==l1.and.l2_==l3
       lo2=lsjsy(l1,l2,Lb,L,l3,La)
       lo3=lsjsy(L,T,L_,l3_,La,l2)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(No+N_+Nb)
          w2=sjsy(l1,l2,Lb,L,l3,La)
          w3=sjsy(L,T,L_,l3_,La,l2)
          cv=w1*w2*w3
       endif

    case(16)! o(3) P12 
       lo1=l1_==l2.and.l2_==l1.and.La==Lb
       lo2=ltjsy(l1,l2,Lb)
       lo3=lsjsy(L,T,L_,l3_,Lb,l3)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(N_)/dble(2*La+1)
          w2=sjsy(L,T,L_,l3_,Lb,l3)
          cv=w1*w2
       endif

    case(17)! o(3) P123 
       lo1=l1_==l2.and.l2_==l3
       lo2=lsjsy(l1,l2,Lb,l3,L,La)
       lo3=lsjsy(L,T,L_,l3_,La,l1)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(No+N_)
          w2=sjsy(l1,l2,Lb,l3,L,La)
          w3=sjsy(L,T,L_,l3_,La,l1)
          cv=w1*w2*w3
       endif

    case(18)! o(3) P321 
       lo1=l1_==l3.and.l2_==l1
       lo2=lsjsy(l1,l2,Lb,L,l3,La)
       lo3=lsjsy(L,T,L_,l3_,La,l2)
       lv=lo1.and.lo2.and.lo3
       if(lv)then
          w1=pifact*mypar(Na+Nb+N_+No)
          w2=sjsy(l1,l2,Lb,L,l3,La)
          w3=sjsy(L,T,L_,l3_,La,l2)
          cv=w1*w2*w3
       endif

    end select

    return

  end subroutine MOTSAF3


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
    sbcsn=CGC(j1,j2,j3,0,0)*sqrt(dble((2*j1+1)*(2*j2+1))/dble(2*j3+1))
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
     w1=cgc(l1,l2,L,mu,M-mu)
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


!> Computes the associeted Legendre polynomials. 
!! Computes Pl-|m| from Pl|m|.
real(kind(1d0)) function Plm_B(x,l,m) result(res)
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
  Pmm=1.d0-2.d0*mod(abs(m),2)
  do j=1,abs(m)
     Pmm=Pmm*0.5d0*dble(abs(m)+j)*st
  enddo
  if(l==abs(m))then
     res=Pmm
     if ( m < 0 ) then
        res = (-1.d0)**m * Factorial_B(l-abs(m))/Factorial_B(l+abs(m)) * res
     end if
     return
  endif
  !Compute Plm with the recursive formula
  !(l-m)P_{l,m}=(2l-1)xP_{l-1,m}-(l-1+m)P_{l-2,m}
  c0=Pmm
  c1=Pmm*x*dble(2*abs(m)+1)
  do j=abs(m)+2,l
     c2=(dble(2*j-1)*x*c1-dble(j-1+abs(m))*c0)/dble(j-abs(m))
     c0=c1
     c1=c2
  enddo
  res=c1
  if ( m < 0 ) then
     res = (-1.d0)**m * Factorial_B(l-abs(m))/Factorial_B(l+abs(m)) * res
  end if
  return
end function Plm_B

!> Factorials function.
real(kind(1d0)) function Factorial_B( n ) result (res)
  integer, intent(in) :: n
  integer :: i
  res = 1.d0
  if ( n < 0 ) then
     write(*,*) 'The factorial is implemented for non-negative numbers.'
     stop
  end if
  do i = 1, n
    res = res * dble(i) 
  end do
end function Factorial_B


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
!!$  Plm_=Plm(x,l,m) ! does not work for m<0
  Plm_=Plm_B(x,l,m)
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

complex(kind(1d0)) function BSH(theta1,phi1,theta2,phi2,l1,l2,L,M) result(res)
  !Compute the bipolar spherical harmonics
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
     res=res+cgc(l1,l2,L,mu,M-mu)*z1*z2
  enddo
  return
end function BSH

complex(kind(1d0)) function SABSH(theta1,phi1,theta2,phi2,l1,l2,L,M,s) result(res)
  !Compute the symmetry adapted bipolar spherical harmonics
  !Yl1l2^LM+/-(Omega1,Omega2)=
  ! =[Y_{l1,l2}^{LM}(Omega1,Omega2)+/-Y_{l1,l2}^{LM}(Omega2,Omega1)]/sqrt(2)
  !If s is even => Y+, if s is odd => Y-
  implicit none
  real(kind(1d0)), intent(in) :: theta1,phi1,theta2,phi2
  integer        , intent(in) :: l1,l2,L,M,s
  complex(kind(1d0)) :: z1,z2
  z1=BSH(theta1,phi1,theta2,phi2,l1,l2,L,M)
  z2=BSH(theta2,phi2,theta1,phi1,l1,l2,L,M)*(1.d0-2.d0*mod(s,2))
  res=(z1+z2)/sqrt(2.d0)
  return
end function SABSH

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
          cgc(l1,l3,l,0,0)*cgc(l2,l4,l,0,0)*sjsy(l1,l2,J,l4,l3,l)+&
          cgc(l1,l4,l,0,0)*cgc(l2,l3,l,0,0)*sjsy(l1,l2,J,l3,l4,l)*&
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




end module symba
