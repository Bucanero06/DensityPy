!From Luca's hmk_standalone


SUBROUTINE DCGAMMA (X,Y,GR,GI)
  !
  !     COMPUTES THE GAMMA FUNCTION OF A COMPLEX ARGUMENT
  !     GR+iGI=GAMMA(X+iY)
  !
  IMPLICIT NONE
  COMPLEX(KIND(1d0)):: Z,ZZ,IZ,IZ2,IZ3,IZ4,FDIV,ZUN,Z1
  REAL(KIND(1d0)):: E,PI,DUEPI,X,Y,GR,GI,AA,XP,RO,UN
  INTEGER :: K,NA,NB,I
  E=2.718281828590452353D0
  PI=3.1415926535897932384D0
  DUEPI=6.2831853071795864769D0
  ZUN= DCMPLX(1.D0,0.D0)
  UN=1.D0
  Z=DCMPLX(X,Y)
  AA=DABS(DATAN(Y/X))
  IF(AA.GT.PI) THEN
     NA=IDINT(Y)+1
     XP=X+NA
     RO=DSQRT(XP*XP+Y*Y)
     IF(RO.LT.40) THEN
        K=IDINT(RO)
        NB=40-K
        XP=XP+NB
        Z1=DCMPLX(XP,Y)
        FDIV=ZUN
        DO I=1,NA+NB
           FDIV=FDIV*Z
           Z=Z+ZUN  
        ENDDO
        IZ=1/Z
        IZ2=IZ*IZ
        IZ3=IZ2*IZ
        IZ4=IZ3*IZ
        ZZ=(Z/E)**Z*CDSQRT(IZ*DUEPI)*(1+IZ/12+IZ2/288-139*IZ3/51840- &
             571*IZ4/2488320)
        ZZ=ZZ/FDIV
     ELSEIF(RO.GE.40) THEN
        Z1=DCMPLX(XP,Y)
        FDIV=ZUN
        DO I=1,NA
           FDIV=FDIV*Z
           Z=Z+ZUN
        ENDDO
        IZ=1/Z
        IZ2=IZ*IZ
        IZ3=IZ2*IZ
        IZ4=IZ3*IZ
        ZZ=(Z/E)**Z*CDSQRT(IZ*DUEPI)*(1+IZ/12+IZ2/288-139*IZ3/51840- &
             571*IZ4/2488320)
        ZZ=ZZ/FDIV
     ENDIF
  ELSEIF(AA.LE.PI)  THEN 
     Z=DCMPLX(X,Y)
     RO=DSQRT(X*X+Y*Y)
     IF(RO.LT.40)  THEN
        NA=40-IDINT(RO)+1
        XP=X+NA
        Z1=DCMPLX(XP,Y)
        FDIV=ZUN
        DO I=1,NA
           FDIV=FDIV*Z
           Z=Z+ZUN
        ENDDO
        IZ=1/Z
        IZ2=IZ*IZ
        IZ3=IZ2*IZ
        IZ4=IZ3*IZ
        ZZ=(Z/E)**Z*CDSQRT(IZ*DUEPI)*(1+IZ/12+IZ2/288-139*IZ3/51840- &
             571*IZ4/2488320)
        ZZ=ZZ/FDIV
     ELSEIF(RO.GE.40) THEN            
        Z=DCMPLX(X,Y)
        IZ=1/Z
        IZ2=IZ*IZ
        IZ3=IZ2*IZ
        IZ4=IZ3*IZ
        ZZ=(Z/E)**Z*CDSQRT(IZ*DUEPI)*(1+IZ/12+IZ2/288-139*IZ3/51840- &
             571*IZ4/2488320)
     ENDIF
  ENDIF
  GR=DREAL(ZZ)
  GI=DIMAG(ZZ)
  RETURN
END SUBROUTINE DCGAMMA

complex(kind(1d0)) function Gamma(z) result(Gam)
  implicit none
  complex(kind(1d0)), intent(in) :: z
  real(kind(1d0)) :: x,y,u,v
  external DCGAMMA
  x=dble(z)
  y=aimag(z)
  call DCGAMMA(x,y,u,v)
  Gam=(1.d0,0.d0)*u+(0.d0,1.d0)*v
  return
end function Gamma

!> Computes the Coulomb Phase of an electron (charge = -|e|)
!! with orbital angular momentum quantum number l (in units 
!! of $\hbar$) and energy above threshold E (au), in the field 
!! of a point charge Z|e|.
!! The Coulomb phase $\sigma_\ell(E)$ is defined as 
!! $$
!! e^{i\sigma_\ell(E)} = \arg[\Gamma( \ell + 1 - i Z/k )]
!! $$
!! with $k=\sqrt(2E)$
!! (see, e.g., Friedrich, Theoretical Atomic Physics, Eqns. 1.117 - 1.119)
real(kind(1d0)) function CoulombPhase( l, Z, E ) result( sigma )
   
   integer        , intent(in) :: l
   real(kind(1d0)), intent(in) :: Z
   real(kind(1d0)), intent(in) :: E
  
   complex(kind(1d0)) :: zk, zGam
   complex(kind(1d0)), external :: Gamma

   sigma = 0.d0
   if( l < 0 .or. E < 0.d0 )return

   zk    = (1.d0,0.d0) * (dble(l)+1.d0) - (0.d0,1.d0) * dble(Z) / sqrt(2.d0*E) 
   zGam  = Gamma( zk ) 
   sigma = atan2( aimag( zGam ), dble( zGam ) )

   return 
end function CoulombPhase 
