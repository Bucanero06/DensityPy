 real(kind(1d0)) function NCD_Phi(t,t_0,tau) result(Phi)
    !
    !.. Return the Normal Cumulative Distribution Phi(x),
    !   which is defined as 
    !   $$
    !   \Phi(x)\,=\,\frac{1}{2\pi}\int_{-\infty}^x\,e^{-t^2/2} dt,
    !   $$
    !   for $x=(t-t_0)/\tau$.
    !
    !   To compute Phi, the program uses the relation 
    !   with the error distribution function
    !   $$
    !   \Phi(x)\,=\,\frac{1}{2}\left[\,1\,+\,\erf(x/\sqrt{2})\,\right]
    !   $$
    !..
    !
    implicit none
    !
    real(kind(1d0)), intent(in) :: t,t_0,tau
    !
    real(kind(1d0)), parameter :: TAU_THRESHOLD=1.d-20
    !real(kind(1d0)), external :: DERF
    real(kind(1d0)) :: x
    !
    !.. If tau is very small, reduces the definition
    !   of Phi to that of the Theta function
    !..
    if(abs(tau)<TAU_THRESHOLD)then
       if(t<t_0)then
          Phi=0.d0
       else
          Phi=1.d0
       endif
       return
    endif
    !
    x=(t-t_0)/tau
    !
    Phi=0.5d0*(1.d0+DERF(x))
    !
    return
    !
  end function NCD_Phi


  real(kind(1d0)) function dNCD_Phi(t,t_0,tau) result(dPhi)
    !.. Return the derivative of NCD_Phi
    !..
    implicit none
    !
    real(kind(1d0)), intent(in) :: t,t_0,tau
    !
    real(kind(1d0)), parameter :: TAU_THRESHOLD=1.d-20
    real(kind(1d0)), parameter :: One_Over_Sqrt_2Pi = 0.39894228040143267794d0
    real(kind(1d0)) :: x
    !
    dPhi=0.d0
    if(abs(tau)<TAU_THRESHOLD)return
    !
    x=(t-t_0)/tau
    !
    dPhi = One_Over_Sqrt_2Pi / tau * Exp( -0.5d0 * x * x )
    !
    return
    !
  end function DNCD_Phi
