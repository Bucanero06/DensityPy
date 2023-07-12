! Function Irr_lm(l,m) result(i)
!   implicit none
!   integer, intent(in) :: l,m ! input
!   integer             :: i ! output
!   i=0
!   If((mod(l,2).eq.0).and.(mod(m,2).eq.0))then
!      If(m.ge.0)then
!         i=1 
!      Else
!         i=4 
!      End IF
!   ElseIf((mod(l,2).eq.0).and.(abs(mod(m,2)).eq.1))then
!      If(m.ge.0)then
!         i=6 
!      Else
!         i=7 
!      End IF
!   ElseIf((mod(l,2).eq.1).and.(mod(m,2).eq.0))then
!      If(m.ge.0)then
!         i=5 
!      Else
!         i=8 
!      End IF
!   ElseIf((mod(l,2).eq.1).and.(abs(mod(m,2)).eq.1))then
!      If(m.ge.0)then
!         i=2 
!      Else
!         i=3 
!      End IF
!   End IF
! end function irr_lm

! !This is the order found in ukrmol (same order than dalton):
! !#      Ag    Xlm       l even m even  N =  l/2+1      1  ->  1
! !#      B3u   Xlm       l odd  m odd   N = (l+1)/2     2  ->  8
! !#      B2u   Xl-m      l odd  m odd   N = (l+1)/2     3  ->  7
! !#      B1g   Xl-m      l even m even  N =  l/2        4  ->  2 
! !#      B1u   Xlm       l odd  m even  N = (l+1)/2     5  ->  6
! !#      B2g   Xlm       l even m odd   N =  l/2        6  ->  3
! !#      B3g   Xl-m      l even m odd   N =  l/2        7  ->  4
! !#      Au    Xl-m      l odd  m even  N = (l-1)/2     8  ->  5
! !This is the order found in CloseCouplingBasis file
! !#      Ag    Xlm       l even m even  N =  l/2+1
! !#      B1g   Xl-m      l even m even  N =  l/2
! !#      B2g   Xlm       l even m odd   N =  l/2
! !#      B3g   Xl-m      l even m odd   N =  l/2
! !#      Au    Xl-m      l odd  m even  N = (l-1)/2
! !#      B1u   Xlm       l odd  m even  N = (l+1)/2
! !#      B2u   Xl-m      l odd  m odd   N = (l+1)/2
! !#      B3u   Xlm       l odd  m odd   N = (l+1)/2

!This is the integral of three spherical harmonics (none of them conjugate)
function ThreeYlmIntegral(l1,m1,l2,m2,l3,m3) result(inte)
  use ModuleAngularMomentum
  implicit none
  integer, intent(in) :: l1,m1,l2,m2,l3,m3 ! input
  real(kind(1d0))     :: inte ! output
  real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
  inte=0.d0
  If((m1+m2).eq.(-m3))then
     inte=((-1.d0)**(m3))*sqrt(((2*l1+1)*(2*l2+1))/((4*PI)*(2*l3+1)))
     inte=inte*ClebschGordanCoefficient(l1,l2,l3,m1,m2)*ClebschGordanCoefficient(l1,l2,l3,0,0)
  End IF
end function ThreeYlmIntegral

!This is the integral of Xl1m1  Xl2m2 Yl3m3 
function TwoXlmOneYlmIntegral(l1,m1,l2,m2,l3,m3) result(inte)
  implicit none
  integer, intent(in) :: l1,m1,l2,m2,l3,m3 ! input
  real(kind(1d0))     :: inte,amp1,fas1,amp2,fas2 ! output
  real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
  real   (kind(1d0)), external :: ThreeYlmIntegral
  write(*,*) "TwoXlmOneYlmIntegral must be complex"
  stop
  inte=0.d0
  If(m1.eq.0)then
     amp1=1.d0
     fas1=0.d0
  ElseIf(m1.lt.0)then
     amp1=1.d0/sqrt(2.d0)
     fas1=amp1*((-1.d0)**abs(m1))
  Else
     amp1=1.d0/sqrt(2.d0)
     fas1=-amp1*((-1.d0)**abs(m1))
  End IF
  If(m2.eq.0)then
     amp2=1.d0
     fas2=0.d0
  ElseIf(m2.lt.0)then
     amp2=1.d0/sqrt(2.d0)
     fas2=amp2*((-1.d0)**abs(m2))
  Else
     amp2=1.d0/sqrt(2.d0)
     fas2=-amp2*((-1.d0)**abs(m2))
  End IF
  inte=0.d0
  inte=     amp1*amp2*ThreeYlmIntegral(l1,abs(m1),l2, abs(m2),l3,m3)
  inte=inte+amp1*fas2*ThreeYlmIntegral(l1,abs(m1),l2,-abs(m2),l3,m3)
  inte=inte+fas1*amp2*ThreeYlmIntegral(l1,-abs(m1),l2, abs(m2),l3,m3)
  inte=inte+fas1*fas2*ThreeYlmIntegral(l1,-abs(m1),l2,-abs(m2),l3,m3)    
end function TwoXlmOneYlmIntegral


!This is the integral of Xl1m1  Xl2m2 Yl3m3 
function ThreeXlmIntegral(l1,m1,l2,m2,l3,m3) result(inteR)
  implicit none
  integer, intent(in) :: l1,m1,l2,m2,l3,m3 ! input
  ! real(kind(1d0))     :: inte,amp1,fas1,amp2,fas2,amp3,fas3 ! output
  real(kind(1d0))     :: inteR
  complex(kind(1d0))     :: inte,amp1,fas1,amp2,fas2,amp3,fas3 ! output
  complex(kind(1d0)), parameter  :: ci=(0.d0,1.d0),c1=(1.d0,0.d0),c0=(0.d0,0.d0)
  real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
  real   (kind(1d0)), external :: ThreeYlmIntegral
  inte=c0
  If(m1.eq.0)then
     amp1=c1*1.d0
     fas1=c1*0.d0
  ElseIf(m1.lt.0)then
     amp1=-ci*1.d0/sqrt(2.d0)
     fas1=-amp1*((-1.d0)**abs(m1))
  Else
     amp1=c1*1.d0/sqrt(2.d0)
     fas1=amp1*((-1.d0)**abs(m1))
  End IF
  If(m2.eq.0)then
     amp2=c1*1.d0
     fas2=c1*0.d0
  ElseIf(m2.lt.0)then
     amp2=-ci*1.d0/sqrt(2.d0)
     fas2=-amp2*((-1.d0)**abs(m2))
  Else
     amp2=c1*1.d0/sqrt(2.d0)
     fas2=amp2*((-1.d0)**abs(m2))
  End IF
  If(m3.eq.0)then
     amp3=c1*1.d0
     fas3=c1*0.d0
  ElseIf(m3.lt.0)then
     amp3=-ci*1.d0/sqrt(2.d0)
     fas3=-amp3*((-1.d0)**abs(m3))
  Else
     amp3=c1*1.d0/sqrt(2.d0)
     fas3=amp3*((-1.d0)**abs(m3))
  End IF
  inte=c0
  ! m1=m2=m3=0
  inte=     amp1*amp2*amp3*ThreeYlmIntegral(l1, abs(m1),l2, abs(m2),l3, abs(m3))
  inte=inte+fas1*fas2*fas3*ThreeYlmIntegral(l1,-abs(m1),l2,-abs(m2),l3,-abs(m3)) !this is exactly 0 because fas1,2,3=0

  !    If m1=0 this is multiplied by i^{-1} because either m2 or m3 are negative
  inte=inte+amp1*amp2*fas3*ThreeYlmIntegral(l1, abs(m1),l2, abs(m2),l3,-abs(m3)) 
  inte=inte+amp1*fas2*amp3*ThreeYlmIntegral(l1, abs(m1),l2,-abs(m2),l3, abs(m3))

  inte=inte+amp1*fas2*fas3*ThreeYlmIntegral(l1, abs(m1),l2,-abs(m2),l3,-abs(m3))

  inte=inte+fas1*amp2*amp3*ThreeYlmIntegral(l1,-abs(m1),l2, abs(m2),l3, abs(m3)) 
  inte=inte+fas1*amp2*fas3*ThreeYlmIntegral(l1,-abs(m1),l2, abs(m2),l3,-abs(m3))

  inte=inte+fas1*fas2*amp3*ThreeYlmIntegral(l1,-abs(m1),l2,-abs(m2),l3, abs(m3))

  inteR=real(inte)
  If(abs(aimag(inte)).gt.(10.d0**(-10)))then

     write(*,*) "matrix element is complex, something is wrong",inte
     write(*,*) l1,m1,l2,m2,l3,m3
     stop
  End If


end function ThreeXlmIntegral

!This is the integral of Xl1m1  Xl2m2 Yl3m3 
function Xlm(theta,phi,l,m) result(inteR)
  use ModuleAngularMomentum
  implicit none
  real   (kind(1d0)), intent(in) :: theta,phi
  integer, intent(in) :: l,m ! input
  ! real(kind(1d0))     :: inte,amp1,fas1,amp2,fas2,amp3,fas3 
  real(kind(1d0))     :: inteR! output
  complex(kind(1d0))     :: inte,amp1,fas1,amp2,fas2,amp3,fas3 
  complex(kind(1d0)), parameter  :: ci=(0.d0,1.d0),c1=(1.d0,0.d0),c0=(0.d0,0.d0)
  real   (kind(1d0)), parameter  :: PI=3.14159265358979323844d0
  inte=c0
  If(m.eq.0)then
     amp1=c1*1.d0
     fas1=c1*0.d0
  ElseIf(m.lt.0)then
     amp1=-ci*1.d0/sqrt(2.d0)
     fas1=amp1*((-1.d0)**abs(m))
  Else
     amp1=c1*1.d0/sqrt(2.d0)
     fas1=-amp1*((-1.d0)**abs(m))
  End IF

  inte=c0
  ! m1=m2=m3=0
  inte=amp1*Ylm(theta,phi,l,abs(m))+fas1*Ylm(theta,phi,l,-abs(m))

  inteR=real(inte)
  If(abs(aimag(inte)).gt.0.d0)then

     write(*,*) "Xlm is complex, something is wrong"
     write(*,*) theta,phi,l,m
     stop
  End If


end function Xlm

function kron(i,j) result(k)
  implicit none
  integer, intent(in) :: i,j ! input
  integer    :: k ! output
  k=0
  If(i.eq.j)then
     k=1
  End IF
end function kron


