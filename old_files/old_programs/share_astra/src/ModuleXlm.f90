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
!> Symmetry adapted spherical harmonics (SASH) \f$ X_{\ell m} \f$ are defined as
!! 
!! For the ASTRA code, we will use the convention in MEST [Olsen2000] ch.6,
!! page 210, Eqns. 6.4.19-6.4.21: 
!! for $ m > 0 $,  
!! \f{eqnarray}
!!    X_{\ell  0} &=& Y_{\ell 0}, \\
!!    X_{\ell  m} &=& (-1)^m \sqrt{2} \Re[ Y_{\ell m} ] \\
!!    X_{\ell -m} &=& (-1)^m \sqrt{2} \Im[ Y_{\ell m} ] \\
!! \f}
!! or
!! \f{eqnarray}
!!    X_{\ell  m} &=& (-1)^m \frac{1}{ \sqrt{2}} [ Y_{\ell m} + (-1)^m Y_{\ell -m}] \\
!!    X_{\ell -m} &=& (-1)^m \frac{1}{i\sqrt{2}} [ Y_{\ell m} - (-1)^m Y_{\ell -m}] \\
!! \f}
!!    
!! - **** - OLD XCHEM CONVENTION - WARNING - 
!! \f{eqnarray}
!!    X_{\ell  0} &=& Y_{\ell 0}, \\
!!    X_{\ell  m} &=& \frac{1}{\sqrt{2}   }[ Y_{\ell m} + ( -1 )^{m} Y_{\ell -m} ], m > 0 \\
!!    X_{\ell -m} &=& \frac{1}{\sqrt{2}\,i}[ Y_{\ell m} - ( -1 )^{m} Y_{\ell -m} ], m > 0.
!! \f}
!! Notice that the phase convention on the SASH differes in Olsen and XCHEM
!! - **** - 
!!
!! where $Y_{\ell m}$ are ordinary spherical harmonics, \f$ |m| \leq \ell \f$.
!!
!! Each SASH span a single irreducible representation of the D2h group or any of its subgroup,
!! if defined with respect to the same set of coordinate cartesian axes. 
!!
!!
!!  C1    
!!      A     Xlm Xl-m  l any  m any   N = (2l+1) 
!!
!!  Cs 
!!      A'    Xlm Xl-m  l+m  even      N = l+1
!!      A''    Xlm Xl-m  l+m  odd       N = l
!!
!!  C2 
!!      A     Xlm Xl-m  l any  m even  N = 2[l/2]+1
!!      B     Xlm Xl-m  l any  m odd   N = 2[(l+1)/2]
!!
!!  Ci
!!      Ag    Xlm Xl-m  l even         N = 2l+1
!!      Au    Xlm Xl-m  l odd          N = 2l+1
!!
!!  C2v
!!      A1    Xlm       l any  m even  N = [l/2] + 1
!!      A2    Xl-m      l any  m even  N = [l/2]
!!      B1    Xlm       l any  m odd   N = [(l+1)/2]
!!      B2    Xl-m      l any  m odd   N = [(l+1)/2]
!!
!!  C2h
!!      Ag    Xlm Xl-m  l even m even  N = l+1
!!      Bg    Xlm Xl-m  l even m odd   N = l
!!      Au    Xlm Xl-m  l odd  m even  N = l
!!      Bu    Xlm Xl-m  l odd  m odd   N = l+1
!!
!!  D2
!!      A     Xlm       l even m even  N =  l/2+1  ||  
!!            Xl-m      l odd  m even  N = (l-1)/2
!!      B1    Xlm       l odd  m even  N = (l+1)/2 ||  
!!            Xl-m      l even m even  N =  l/2
!!      B2    Xlm       l even m odd   N =  l/2    ||  
!!            Xl-m      l odd  m odd   N = (l+1)/2
!!      B3    Xlm       l odd  m odd   N = (l+1)/2 ||  
!!            Xl-m      l even m odd   N =  l/2
!!
!!  D2h
!!      Ag    Xlm       l even m even  N =  l/2+1
!!      B1g   Xl-m      l even m even  N =  l/2
!!      B2g   Xlm       l even m odd   N =  l/2
!!      B3g   Xl-m      l even m odd   N =  l/2
!!      Au    Xl-m      l odd  m even  N = (l-1)/2
!!      B1u   Xlm       l odd  m even  N = (l+1)/2
!!      B2u   Xl-m      l odd  m odd   N = (l+1)/2
!!      B3u   Xlm       l odd  m odd   N = (l+1)/2
!!
module ModuleXlm
  
  use, intrinsic :: ISO_FORTRAN_ENV
  
  !.. Level 0.0
  use ModuleErrorHandling
  use ModuleString
  use ModuleConstants
  use ModuleAngularMomentum
  !.. Level 0.1
  use ModuleGroups

  implicit none

  private

  !> Defines the class of a definite symmetry adapted spherical harmonic.
  type, public :: ClassXlm
     ! {{{ 

     !> Angular momentum.
     integer :: l
     !> Angular momentum projection.
     integer :: m

     ! }}}
   contains
     !> Initializes the Xlm class.
     procedure :: Init => ClassXlmInit
     generic   :: Getl => ClassXlmGetl
     generic   :: Getm => ClassXlmGetm
     generic   :: Show => ClassXlmShow
     generic   :: Save => ClassXlmSave
     generic   :: Load => ClassXlmLoad
     generic   :: Is   => ClassXlmIs
     generic   :: GetLabel => ClassXlmGetLabel
     generic   :: GetIrrep => ClassXlmGetIrrep
     generic, public :: assignment(=)  => ClassXlmCopyXlm
     ! {{{ 

     procedure, private :: ClassXlmGetl
     procedure, private :: ClassXlmGetm
     procedure, public  :: ClassXlmShow
     procedure, public  :: ClassXlmSave
     procedure, public  :: ClassXlmLoad
     procedure, private :: ClassXlmIs
     procedure, private :: ClassXlmGetLabel
     procedure, private :: ClassXlmGetIrrep
     procedure, private :: ClassXlmCopyXlm

     ! }}}
  end type ClassXlm


  !> Defines the class of the symmetry adapted spherical harmonics set with a definite symmetry.
  type, public :: ClassXlmSymmetricSet
     ! {{{ 

     private
     !> Maximum angular momentum
     integer                     :: Lmax
     !> Number of Symmetric Xlm for each angular momentum
     integer, allocatable        :: dim(:)
     !> Absolute index of the last Xlm with angular momentum 
     !! lower than a given l
     integer, allocatable        :: idim(:)
     !> Irrep class associated to the Xml set.
     type(ClassIrrep), pointer   :: Irrep
     !> Total number of Xlm functions for a given symmetry.
     integer                     :: NumXlm
     !> List of the Xlm functions belonging to the set of definite symmetry.
     type(ClassXlm), allocatable :: XlmList(:)

     ! }}}
   contains
     !
     generic, public :: init      =>  ClassXlmSymmetricSetInit
     generic, public :: GetIrrep  =>  ClassXlmSymmetricSetGetIrrep
     generic, public :: Show      =>  ClassXlmSymmetricSetShow
     generic, public :: GetNXlm   =>  ClassXlmSymmetricSetGetNXlm
     generic, public :: GetXlm    =>  ClassXlmSymmetricSetGetXlm
     generic, public :: GetMList  =>  ClassXlmSymmetricSetGetMList
     generic, public :: GetNm     =>  ClassXlmSymmetricSetGetNm
     generic, public :: ValidXlm  =>  ClassXlmSymmetricSetValidXlm
!!$     generic, public :: GetCartesianExpan => ClassXlmSymmetricSetGetCartesianExpan
     procedure, public :: Free => ClassXlmSymmetricSetFree

     !> Retrieves if the ClassXlmSymmetricSet has been initialized
     procedure         :: Initialized => ClassXlmSymmetricSetInitialized
     ! {{{ 

     procedure, private :: ClassXlmSymmetricSetInit
     procedure, private :: ClassXlmSymmetricSetGetIrrep
     procedure, private :: ClassXlmSymmetricSetShow
     procedure, private :: ClassXlmSymmetricSetGetNXlm
     procedure, private :: ClassXlmSymmetricSetGetXlm
     procedure, private :: ClassXlmSymmetricSetGetMList
     procedure, private :: ClassXlmSymmetricSetGetNm
     procedure, private :: ClassXlmSymmetricSetValidXlm
!!$     procedure, private :: ClassXlmSymmetricSetGetCartesianExpan
     procedure, private :: Set_DimVectors => ClassXlmSymmetricSetDimension
     !
     procedure, private :: Set_C1_A
     !
     procedure, private :: Set_Cs_Ap
     procedure, private :: Set_Cs_App
     !
     procedure, private :: Set_C2_A
     procedure, private :: Set_C2_B
     !
     procedure, private :: Set_Ci_Ag
     procedure, private :: Set_Ci_Au
     !
     procedure, private :: Set_C2v_A1
     procedure, private :: Set_C2v_A2
     procedure, private :: Set_C2v_B1
     procedure, private :: Set_C2v_B2
     !
     procedure, private :: Set_C2h_Ag
     procedure, private :: Set_C2h_Bg
     procedure, private :: Set_C2h_Au
     procedure, private :: Set_C2h_Bu
     !
     procedure, private :: Set_D2_A
     procedure, private :: Set_D2_B1
     procedure, private :: Set_D2_B2
     procedure, private :: Set_D2_B3
     !
     procedure, private :: Set_D2h_Ag
     procedure, private :: Set_D2h_B1g
     procedure, private :: Set_D2h_B2g
     procedure, private :: Set_D2h_B3g
     procedure, private :: Set_D2h_Au
     procedure, private :: Set_D2h_B1u
     procedure, private :: Set_D2h_B2u
     procedure, private :: Set_D2h_B3u
     !
     final :: ClassXlmSymmetricSetFinal
     ! }}}
  end type ClassXlmSymmetricSet


  type, public :: ClassXlmSet
     ! {{{ 

     private
     logical                             :: Initialized = .false.
     type(ClassGroup)          , pointer :: Group
     !> Maximum angular momentum that can be represented.
     integer                             :: Lmax
     integer                             :: NIrreps
     type(ClassXlmSymmetricSet), pointer :: XlmSymSet(:)

     ! }}}
   contains

     generic,   public :: GetSymSet  => ClassXlmSetGetSymSet, ClassXlmSetGetSymSetFromPtr

     !> Given the global group and the maximum L, initializes 
     !! the Xlm symmetric sets that correspond to all the group's irreps
     procedure, public :: init       => ClassXlmSetInit
     !> Given an irrep of the group, returns a pointer to the
     !! Xlm symmetric set corresponding to that irrep
     procedure, public :: ClassXlmSetGetSymSet
     procedure, private:: ClassXlmSetGetSymSetFromPtr
     procedure, public :: Show       => ClassXlmSetShow
     !> For a given irreducible representation and an angular momentum retrieve for all te possible Xlm ( r^l*Xlm ), its expansion in terms of the cartessian monomials. 
!!$     procedure, public :: GetCartesianExpan
  end type ClassXlmSet

  type(ClassXlmSet), public :: GlobalXlmSet




  public :: TripleXlmIntegral, EvalXlm
  private :: delta
  private :: signm
  private :: Ylm_Xlm_Overlap_Factor1
  
contains

  !> Initializes the Xlm class.
  subroutine ClassXlmInit( Xlm, l, m )
    !
    !> Class of Xlm with definite l and m
    class (ClassXlm),intent(inout) :: Xlm
    !> Angular momentum.
    integer,         intent(in)    :: l
    !> Angular momentum projection.
    integer,         intent(in)    :: m
    !
    Xlm%l = l
    Xlm%m = m
    !
  end subroutine ClassXlmInit


  !> Evaluates the Xlm function for specific angles.
  real(kind(1d0)) function EvalXlm( l, m, theta, phi ) result(res)
    !
    integer        , intent(in) :: l, m
    real(kind(1d0)), intent(in) :: theta, phi
    !
    !! \f{eqnarray}
    !!    X_{\ell  0} &=& Y_{\ell 0}, \\
    !!    X_{\ell  m} &=& (-1)^m \sqrt{2} \Re[ Y_{\ell m} ] \\
    !!    X_{\ell -m} &=& (-1)^m \sqrt{2} \Im[ Y_{\ell m} ] \\
    !! \f}
    !
    if(m==0)then
       res = dble(Ylm(theta,phi,l,0))
    elseif(m>0)then
       res = sqrt(2.d0) * dble( Ylm(theta,phi,l,m))
       if(mod(m,2)==1) res = -res
    else
       res = sqrt(2.d0) * aimag(Ylm(theta,phi,l,-m))
       if(mod(-m,2)==1) res = -res
    endif
    !
  end function EvalXlm
  

  !> Tell to which irrep the Xlm belong
  function ClassXlmGetIrrep( Xlm, Group ) result( irrep )
    !
    !> Class of Xlm with definite l and m
    class(ClassXlm)  , intent(in) :: Xlm
    !> input group
    class(ClassGroup), intent(in) :: Group
    class(ClassIrrep), pointer    :: irrep
    !
    type(ClassXlmSymmetricSet)   :: XlmSymSet
    type(ClassIrrep), pointer    :: irrepv(:)
    integer :: iIrrep
    !
    irrepv => Group%GetIrrepList()
    do iIrrep = 1, size(irrepv)
       call XlmSymSet%init( Xlm%GetL(), irrepv( iIrrep ))
       if( XlmSymSet%ValidXlm( Xlm ))then
          irrep => irrepv( iIrrep )
          return
       endif
    enddo
    !
  end function ClassXlmGetIrrep


  !> Copies the information of Xlm class to another.
  subroutine ClassXlmCopyXlm( XlmOut, XlmIn )
    !> Class of Xlm to be assigned.
    class(ClassXlm)  , intent(inout) :: XlmOut
    class(ClassXlm)  , intent(in)    :: XlmIn
    !
    XlmOut%l = XlmIn%l
    XlmOut%m = XlmIn%m
    !
  end subroutine ClassXlmCopyXlm


  subroutine ClassXlmSetInit( XlmSet, Group, Lmax )
    !
    class(ClassXlmSet),         intent(out) :: XlmSet
    class(ClassGroup) , target, intent(in)  :: Group
    integer           ,         intent(in)  :: Lmax
    integer :: iIrrep
    type(ClassIrrep), pointer, dimension(:) :: irrepv
    !
    XlmSet%Group => Group
    XlmSet%Lmax  =  Lmax

    irrepv => Group%GetIrrepList()
    XlmSet%NIrreps = size( irrepv )
    allocate( XlmSet%XlmSymSet( XlmSet%NIrreps ) )
    do iIrrep = 1, size( irrepv )
       call XlmSet%XlmSymSet( iIrrep )%init( Lmax, irrepv( iIrrep ) )
    enddo
    XlmSet%Initialized = .true.

  end subroutine ClassXlmSetInit



  subroutine ClassXlmSymmetricSetInit( XlmSym, Lmax, IrrepPtr )
    class(ClassXlmSymmetricSet), intent(out):: XlmSym
    integer                    , intent(in) :: Lmax
    type(ClassIrrep)           , intent(in) :: irrepPtr
    !
    character(len=16) :: IdStrn
    
    IdStrn = IrrepPtr%GetIdStrn()

    XlmSym%lmax = lmax

    if( IdStrn .is. "C1.A"    ) call XlmSym%Set_C1_A   ( LMax, irrepPtr )
    if( IdStrn .is. "Cs.Ap"   ) call XlmSym%Set_Cs_Ap  ( LMax, irrepPtr )
    if( IdStrn .is. "Cs.App"  ) call XlmSym%Set_Cs_App ( LMax, irrepPtr )
    if( IdStrn .is. "C2.A"    ) call XlmSym%Set_C2_A   ( LMax, irrepPtr )
    if( IdStrn .is. "C2.B"    ) call XlmSym%Set_C2_B   ( LMax, irrepPtr )
    if( IdStrn .is. "Ci.Ag"   ) call XlmSym%Set_Ci_Ag  ( LMax, irrepPtr )
    if( IdStrn .is. "Ci.Au"   ) call XlmSym%Set_Ci_Au  ( LMax, irrepPtr )
    if( IdStrn .is. "C2v.A1"  ) call XlmSym%Set_C2v_A1 ( LMax, irrepPtr )
    if( IdStrn .is. "C2v.B1"  ) call XlmSym%Set_C2v_B1 ( LMax, irrepPtr )
    if( IdStrn .is. "C2v.A2"  ) call XlmSym%Set_C2v_A2 ( LMax, irrepPtr )
    if( IdStrn .is. "C2v.B2"  ) call XlmSym%Set_C2v_B2 ( LMax, irrepPtr )
    if( IdStrn .is. "C2h.Ag"  ) call XlmSym%Set_C2h_Ag ( LMax, irrepPtr )
    if( IdStrn .is. "C2h.Bg"  ) call XlmSym%Set_C2h_Bg ( LMax, irrepPtr )
    if( IdStrn .is. "C2h.Au"  ) call XlmSym%Set_C2h_Au ( LMax, irrepPtr )
    if( IdStrn .is. "C2h.Bu"  ) call XlmSym%Set_C2h_Bu ( LMax, irrepPtr )
    if( IdStrn .is. "D2.A"    ) call XlmSym%Set_D2_A   ( LMax, irrepPtr )
    if( IdStrn .is. "D2.B1"   ) call XlmSym%Set_D2_B1  ( LMax, irrepPtr )
    if( IdStrn .is. "D2.B2"   ) call XlmSym%Set_D2_B2  ( LMax, irrepPtr )
    if( IdStrn .is. "D2.B3"   ) call XlmSym%Set_D2_B3  ( LMax, irrepPtr )
    if( IdStrn .is. "D2h.Ag"  ) call XlmSym%Set_D2h_Ag ( LMax, irrepPtr )
    if( IdStrn .is. "D2h.B1g" ) call XlmSym%Set_D2h_B1g( LMax, irrepPtr )
    if( IdStrn .is. "D2h.B2g" ) call XlmSym%Set_D2h_B2g( LMax, irrepPtr )
    if( IdStrn .is. "D2h.B3g" ) call XlmSym%Set_D2h_B3g( LMax, irrepPtr )
    if( IdStrn .is. "D2h.Au"  ) call XlmSym%Set_D2h_Au ( LMax, irrepPtr )
    if( IdStrn .is. "D2h.B1u" ) call XlmSym%Set_D2h_B1u( LMax, irrepPtr )
    if( IdStrn .is. "D2h.B2u" ) call XlmSym%Set_D2h_B2u( LMax, irrepPtr )
    if( IdStrn .is. "D2h.B3u" ) call XlmSym%Set_D2h_B3u( LMax, irrepPtr )

    call XlmSym%Set_DimVectors()

  end subroutine ClassXlmSymmetricSetInit



  subroutine Set_C1_A( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + 2*l+1
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m= -l, l
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_C1_A


  subroutine Set_Cs_Ap( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    XlmSym%Lmax = LMax
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + l+1
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m= -l, l
          if ( mod(l+m,2) == 0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_Cs_Ap


  subroutine Set_Cs_App( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + l
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m= -l, l
          if ( mod(l+m+1,2) == 0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_Cs_App




  subroutine Set_C2_A( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + 2*floor(dble(l/2))+1
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m= -l, l
          if ( mod(m,2) == 0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_C2_A



  subroutine Set_C2_B( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + 2*floor(dble((l+1)/2))
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m= -l, l
          if ( mod(m+1,2) == 0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_C2_B




  subroutine Set_Ci_Ag( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax, 2
       nXlm = nXlm + 2*l+1
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax, 2
       do m = -l, l
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_Ci_Ag



  subroutine Set_Ci_Au( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 1, lmax, 2
       nXlm = nXlm + 2*l+1
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 1, lmax, 2
       do m = -l, l
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_Ci_Au



  subroutine Set_C2v_A1( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + floor(dble(l/2))+1
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m = 0, l, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_C2v_A1



  subroutine Set_C2v_A2( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + floor(dble(l/2))
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m = -l, -2
          if ( mod(-m,2) == 0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_C2v_A2



  subroutine Set_C2v_B1( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + floor(dble((l+1)/2))
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m = 1, l, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_C2v_B1



  subroutine Set_C2v_B2( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax
       nXlm = nXlm + floor(dble((l+1)/2))
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m = -l, -1
          if ( mod(-m+1,2) == 0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_C2v_B2



  subroutine Set_C2h_Ag( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax, 2
       nXlm = nXlm + l+1
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax, 2
       do m = -l, l, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_C2h_Ag



  subroutine Set_C2h_Bg( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax, 2
       nXlm = nXlm + l
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax, 2
       do m = -l+1, l-1, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_C2h_Bg



  subroutine Set_C2h_Au( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 1, lmax, 2
       nXlm = nXlm + l
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 1, lmax, 2
       do m = -l+1, l-1, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_C2h_Au



  subroutine Set_C2h_Bu( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 1, lmax, 2
       nXlm = nXlm + l+1
    enddo 
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 1, lmax, 2
       do m = -l, l, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_C2h_Bu



  subroutine Set_D2_A( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    !
    do l = 0, lmax, 2
       nXlm = nXlm + l/2 + 1
    enddo 
    do l = 1, lmax, 2
       nXlm = nXlm + (l-1)/2
    enddo 
    !
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m = -l, l
          if ( mod(l,2)==0 .and. mod(m,2)==0 .and. m>=0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          elseif ( mod(l+1,2)==0 .and. mod(m,2)==0 .and. m<0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_D2_A



  subroutine Set_D2_B1( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    !
    do l = 0, lmax, 2
       nXlm = nXlm + l/2
    enddo 
    do l = 1, lmax, 2
       nXlm = nXlm + (l+1)/2
    enddo 
    !
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m = -l, l
          if ( mod(l,2)==0 .and. mod(m,2)==0 .and. m<0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          elseif ( mod(l+1,2)==0 .and. mod(m,2)==0 .and. m>=0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_D2_B1



  subroutine Set_D2_B2( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    !
    do l = 0, lmax, 2
       nXlm = nXlm + l/2
    enddo 
    do l = 1, lmax, 2
       nXlm = nXlm + (l+1)/2
    enddo 
    !
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m = -l, l
          if ( mod(l,2)==0 .and. mod(m+1,2)==0 .and. m>=0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          elseif ( mod(l+1,2)==0 .and. mod(m+1,2)==0 .and. m<0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_D2_B2



  subroutine Set_D2_B3( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    !
    XlmSym%Irrep => Irrep
    nXlm = 0
    !
    do l = 0, lmax, 2
       nXlm = nXlm + l/2
    enddo 
    do l = 1, lmax, 2
       nXlm = nXlm + (l+1)/2
    enddo 
    !
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) ) 
    nXlm = 0
    do l = 0, lmax
       do m = -l, l
          if ( mod(l+1,2)==0 .and. mod(m+1,2)==0 .and. m>=0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          elseif ( mod(l,2)==0 .and. mod(m+1,2)==0 .and. m<0 ) then
             nXlm = nXlm + 1
             call XlmSym%XlmList( nXlm )%init( l, m )
          end if
       enddo
    enddo
  end subroutine Set_D2_B3



  subroutine Set_D2h_Ag( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax, 2
       nXlm = nXlm + l/2 + 1
    enddo
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) )
    nXlm = 0
    do l = 0, lmax, 2
       do m = 0, l, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_D2h_Ag



  subroutine Set_D2h_B1g( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax, 2
       nXlm = nXlm + l/2
    enddo
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) )
    nXlm = 0
    do l = 0, lmax, 2
       do m = -l, -2, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_D2h_B1g



  subroutine Set_D2h_B2g( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax, 2
       nXlm = nXlm + l/2
    enddo
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) )
    nXlm = 0
    do l = 0, lmax, 2
       do m = 1, l, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_D2h_B2g



  subroutine Set_D2h_B3g( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 0, lmax, 2
       nXlm = nXlm + l/2
    enddo
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) )
    nXlm = 0
    do l = 0, lmax, 2
       do m = -l+1, -1, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_D2h_B3g



  subroutine Set_D2h_Au( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 1, lmax, 2
       nXlm = nXlm + (l-1)/2
    enddo
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) )
    nXlm = 0
    do l = 1, lmax, 2
       do m = -l+1, -2, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_D2h_Au



  subroutine Set_D2h_B1u( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 1, lmax, 2
       nXlm = nXlm + ( l + 1 ) / 2 
    enddo
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) )
    nXlm = 0
    do l = 1, lmax, 2
       do m = 0, l, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_D2h_B1u



  subroutine Set_D2h_B2u( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 1, lmax, 2
       nXlm = nXlm + ( l + 1 ) / 2 
    enddo
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) )
    nXlm = 0
    do l = 1, lmax, 2
       do m = -l, -1, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_D2h_B2u



  subroutine Set_D2h_B3u( XlmSym, Lmax, Irrep )
    class( ClassXlmSymmetricSet ), intent(inout) :: XlmSym
    integer                      , intent(in)    :: Lmax
    class( ClassIrrep ), target  , intent(in)    :: Irrep
    !
    integer :: l,m,nXlm
    XlmSym%Irrep => Irrep
    nXlm = 0
    do l = 1, lmax, 2
       nXlm = nXlm + ( l + 1 ) / 2 
    enddo
    XlmSym%NumXlm = nXlm
    allocate( XlmSym%XlmList( nXlm ) )
    nXlm = 0
    do l = 1, lmax, 2
       do m = 1, l, 2
          nXlm = nXlm + 1
          call XlmSym%XlmList( nXlm )%init( l, m )
       enddo
    enddo
  end subroutine Set_D2h_B3u


  function ClassXlmGetl( Xlm ) result( l )
    class(ClassXlm), intent(in) :: Xlm
    integer :: l
    l = Xlm%l
  end function ClassXlmGetl


  function ClassXlmGetm( Xlm ) result( m )
    class(ClassXlm), intent(in) :: Xlm
    integer :: m
    m = Xlm%m
  end function ClassXlmGetm


  subroutine ClassXlmSymmetricSetDimension( XlmSym )
    class(ClassXlmSymmetricSet), intent(inout) :: XlmSym
    integer :: iXlm, lmax, il
    !
    lmax = XlmSym%lmax

    if(allocated(XlmSym%dim))deallocate(XlmSym%dim)
    allocate(XlmSym%dim(0:lmax))
    XlmSym%dim = 0
    do iXlm = 1, XlmSym%NumXlm
       il = XlmSym%XlmList( iXlm )%Getl()
       XlmSym%dim( il ) = XlmSym%dim( il ) + 1
    enddo
    if(allocated(XlmSym%idim))deallocate(XlmSym%idim)
    allocate(XlmSym%idim(0:lmax+1))
    XlmSym%idim = 0
    do il = 0, lmax
       XlmSym%idim(il+1)=XlmSym%idim(il)+XlmSym%dim(il)
    enddo
  end subroutine ClassXlmSymmetricSetDimension

     !> Retrieves if the ClassXlmSymmetricSet has been initialized
  logical function  ClassXlmSymmetricSetInitialized( XlmSet ) result( Initialized )
    !
    !> Class of the Xlm set with a definite symmetry.
    class(ClassXlmSymmetricSet), intent(in) :: XlmSet
    !
    if ( allocated( XlmSet%XlmList ) ) then
       !
       Initialized = .true.
       !
    else
       !
       Initialized = .false.
       !
    end if
    !
  end function ClassXlmSymmetricSetInitialized



  !> Computes the angular integral 
  !! \f[
  !!     \langle X_{\ell_{Bra}m_{Bra}}|X_{\ell_{Op}m_{Op}}|X_{\ell_{Ket}m_{Ket}}\rangle
  !! \f]
  real(kind(1d0)) function TripleXlmIntegral(LBra,MBra,LOp,MOp,LKet,MKet) result( Res )
    !
    integer, intent(in) :: LBra,MBra,LOp,MOp,LKet,MKet
    !
    complex(kind(1d0)) :: FactMult
    real(kind(1d0))    :: SignFactorMBra, SignFactorMOp, SignFactorMKet, SignFactor, FactMult2
    !
    Res = 0.d0
    !
    FactMult = Ylm_Xlm_Overlap_Factor1( MKet ) * &
         Ylm_Xlm_Overlap_Factor1( MOp ) * &
         conjg(Ylm_Xlm_Overlap_Factor1( MBra ))
    !
    FactMult2 = sqrt(dble((2*LKet+1)*(2*LOp+1))/(4.d0*PI*dble(2*LBra+1))) * &
         ClebschGordanCoefficient(LKet,LOp,LBra,0,0)
    !
    SignFactorMBra = signm(MBra)*dble((-1)**MBra)
    SignFactorMOp  = signm(MOp) *dble((-1)**MOp)
    SignFactorMKet = signm(MKet)*dble((-1)**MKet)
    !
    SignFactor = dble((-1)**(LKet+LOp-LBra))
    !
    Res = ( ClebschGordanCoefficient(LKet,LOp,LBra,abs(MKet),abs(MOp))*&
         delta(abs(MBra),abs(MKet)+abs(MOp)) * &
         (1.d0 + SignFactor*SignFactorMBra*SignFactorMOp*SignFactorMKet) + &
         ClebschGordanCoefficient(LKet,LOp,LBra,abs(MKet),abs(MOp))*&
         delta(-abs(MBra),abs(MKet)+abs(MOp)) * &
         (SignFactorMBra + SignFactor*SignFactorMKet*SignFactorMOp) + &
         ClebschGordanCoefficient(LKet,LOp,LBra,abs(MKet),-abs(MOp))*&
         delta(abs(MBra),abs(MKet)-abs(MOp)) * &
         (SignFactorMOp + SignFactor*SignFactorMBra*SignFactorMKet) + &
         ClebschGordanCoefficient(LKet,LOp,LBra,-abs(MKet),abs(MOp))*&
         delta(abs(MBra),-abs(MKet)+abs(MOp)) * &
         (SignFactorMKet + SignFactor*SignFactorMBra*SignFactorMOp) )
    !
    if ( abs(aimag(FactMult))>1.d-12 .and. abs(Res)>1.d-12 ) then
       write(output_unit,*) "Complex factor: ", FactMult
       write(output_unit,*) "Clebsh Gordan part: ", Res
       call Assert( "The triple integral of Xlm has to be real." )
    else
       Res = Res * dble(FactMult) * FactMult2
!!$write(output_unit,*) "triple",lbra,mbra,lket,mket,lop,mop,FactMult2,Res
    end if
    !
  end function TripleXlmIntegral



!!$  real(kind(1d0)) function TripleXlmYlmXlmIntegral(LBra,MBra,LOp,MOp,LKet,MKet) result( Res )
!!$
!!$
!!$  end function TripleXlmYlmXlmIntegral


  complex(kind(1d0)) function Ylm_Xlm_Overlap_Factor1( m ) result( zres )
    integer, intent(in) :: m
    if( m < 0 )then
       zres = -Zi / sqrt(2.d0)
    else if( m == 0 )then
       zres = Z1
    else !(m>0)
       zres = Z1 / sqrt(2.d0)
    endif
  end function Ylm_Xlm_Overlap_Factor1


  real(kind(1d0)) function Inv_YXlm_Overlap_Factor1( m ) result( res )
    integer, intent(in) :: m
    if( m < 0 )then
       res = dble((-1)**m) / sqrt(2.d0)
    else if( m == 0 )then
       res = 1.d0
    else !(m>0)
       res = 1.d0 / sqrt(2.d0)
    endif
  end function Inv_YXlm_Overlap_Factor1


  real(kind(1d0)) function delta( x, y ) result(d)
    !
    integer, intent(in) :: x, y
    !
    if ( x == y ) then
       d = 1.d0
    else
       d = 0.d0
    end if
    !
  end function delta



  real(kind(1d0)) function signm( m )
    integer, intent(in) :: m
    if ( m == 0 ) then
       signm = 0.d0
    elseif ( m > 0 ) then
       signm = 1.d0
    else
       signm = -1.d0
    end if
  end function signm


  !> Builds the transformation matrix that expresses the spherical gaussians in terms of the Xlm symmetry adapted functions. Supposes that the Xlm functions have the same ordening than the harmonics, according to l and m.
  subroutine SpherXlmConversionMatrix( SpherHarmIndexes, TransMat )
    !
    !> Array of the spherical harmonics indexes n,l,m in r^{n} Y_{lm} in the rows, and in the columns the different functions.
    integer,                         intent(in)  :: SpherHarmIndexes(:,:)
    !> If this array is a matrix B, and de vectors S and X represent the spherical harmonics and Xlm functions respectively then S = B x X. So considering the comment and the others above: M = A x B. 
    complex(kind(1d0)), allocatable, intent(out) :: TransMat(:,:)
    !
    complex(kind(1d0)), allocatable :: ArrayTransMat(:,:)
    integer :: NumSpherHarm, i, n, l, m, j
    !
    NumSpherHarm = size(SpherHarmIndexes,2)
    !
    allocate( ArrayTransMat(NumSpherHarm,NumSpherHarm) )
    ArrayTransMat = Z0
    !
    do i = 1, NumSpherHarm
       !
       n = SpherHarmIndexes(1,i)
       l = SpherHarmIndexes(2,i)
       m = SpherHarmIndexes(3,i)
       !
       do j = 1, NumSpherHarm
          !
          if ( (SpherHarmIndexes(1,j)==n) .and. (SpherHarmIndexes(2,j)==l) ) then
             !
             if ( m == 0 ) then
                !
                if ( SpherHarmIndexes(3,j) == 0 ) then
                   ArrayTransMat(i,j) = Z1
                end if
                !
             elseif ( m > 0 ) then
                !
                if ( SpherHarmIndexes(3,j) == m ) then
                   ArrayTransMat(i,j) = Z1/sqrt(2.d0)
                end if
                if ( SpherHarmIndexes(3,j) == -m ) then
                   ArrayTransMat(i,j) = Zi/sqrt(2.d0)
                end if
                !
             elseif( m < 0 ) then
                !
                if ( SpherHarmIndexes(3,j) == m ) then
!!$                   ArrayTransMat(i,j) = -Zi/sqrt(2.d0)
                   ArrayTransMat(i,j) = -Zi/sqrt(2.d0)*(-1.d0)**(-m)
                end if
                if ( SpherHarmIndexes(3,j) == -m ) then
!!$                   ArrayTransMat(i,j) = Z1/sqrt(2.d0)
                   ArrayTransMat(i,j) = Z1/sqrt(2.d0)*(-1.d0)**(-m)
                end if
                !
             end if
             !
          end if
          !
       end do
       !
    end do
    !
    allocate( TransMat, source = ArrayTransMat )
    deallocate( ArrayTransMat )
    !
  end subroutine SpherXlmConversionMatrix



  !> Factorial function
  real(kind(1d0)) function fact(n) 
    !
    implicit none
    !> Factorial function argument.
    integer, intent(in) :: n
    !
    real(kind(1d0)) :: factorial
    integer :: i
    !
    factorial = 1.d0
    if (n==0) then
       fact = factorial
    else
       do i = 1, n
          factorial = factorial*dble(i) 
       end do
       fact = factorial
    end if
    !
  end function fact


  !> Combinatorics.
  real(kind(1d0)) function Comb(a,b)
    !
    !> a <= b
    integer, intent(in) :: a
    !> a <= b
    integer, intent(in) :: b
    !
    Comb = fact(b)/(fact(a)*fact(b-a))
    !
  end function Comb




  !> Prints on a unit the relevant info of ClassXlmSet%
  subroutine ClassXlmSetShow( XlmSet, unit )
    !
    !> Class ClassXlmSet
    class(ClassXlmSet), intent(in) :: XlmSet
    !> Unit in which the info will be printed.
    integer, optional,  intent(in) :: unit
    !
    integer :: OutputUnit, i
    !
    if ( .not.  XlmSet%Initialized ) then
       call Assert( "The class has not been initialized, so cannot be shown." )
    end if
    !
    if ( present(unit) ) then
       OutputUnit = unit
    else
       OutputUnit = OUTPUT_UNIT
    end if
    !
    write(OutputUnit,fmt='(a)') "ClassXlmSet info:"
    write(OutputUnit,fmt='(a)') 
    write(OutputUnit,fmt='(a,i2)') 'Lmax =', XlmSet%Lmax
    write(OutputUnit,fmt='(a,i2)') 'Num Irreps =', XlmSet%NIrreps
    write(OutputUnit,fmt='(a,2x,a)') 'Group =>', XlmSet%Group%GetName()
    !
    do i = 1, XlmSet%NIrreps
       call XlmSet%XlmSymSet(i)%Show( OutputUnit )
    end do
    !
    write(OutputUnit,fmt='(a)') 
    !
  end subroutine ClassXlmSetShow

  !> Prints on a unit the relevant info of ClassXlmSymmetricSet%
  subroutine ClassXlmSymmetricSetShow( XlmSymSet, unit )
    !
    class(ClassXlmSymmetricSet), intent(in) :: XlmSymSet
    integer, optional,           intent(in) :: unit
    !
    integer :: i, outunit
    !
    if ( present(unit) ) then
       outunit = unit
    else
       outunit = OUTPUT_UNIT
    end if
    !
    call XlmSymSet%Irrep%Show( unit )
    !
    write(outunit,fmt='(a,2x,i5)') "Number of Xlm", XlmSymSet%NumXlm
    !
    do i = 1, XlmSymSet%NumXlm
       write(outunit,fmt='(2(a,x),x,2i2)') "l", "m", XlmSymSet%XlmList(i)%l, XlmSymSet%XlmList(i)%m
    end do
    !
  end subroutine ClassXlmSymmetricSetShow


  !> Given an irrep of the group, returns a pointer to the
  !! Xlm symmetric set corresponding to that irrep
  function ClassXlmSetGetSymSetFromPtr( XlmSet, IrrepPtr ) result( XlmSymSetPtr )
    !
    class(ClassXlmSet), target , intent(in) :: XlmSet
    type(ClassIrrep)  , pointer, intent(in) :: IrrepPtr
    type(ClassXlmSymmetricSet) , pointer    :: XlmSymSetPtr

    XlmSymSetPtr => XlmSet%GetSymSet( IrrepPtr%GetName() )

  end function ClassXlmSetGetSymSetFromPtr

  !> Given an irrep of the group, returns a pointer to the
  !! Xlm symmetric set corresponding to that irrep
  function ClassXlmSetGetSymSet( XlmSet, IrrepName ) result( XlmSymSetPtr )
    !
    class(ClassXlmSet), target, intent(in) :: XlmSet
    character(len=*)          , intent(in) :: IrrepName
    type(ClassXlmSymmetricSet), pointer    :: XlmSymSetPtr
    !
    integer :: i
    logical :: Found
    !
    Found = .false.
    !
    XlmSymSetPtr => NULL()
    !
    do i = 1, XlmSet%NIrreps
       if ( XlmSet%XlmSymSet(i)%Irrep%NameIs(IrrepName) ) then
          XlmSymSetPtr => XlmSet%XlmSymSet(i)
          Found = .true.
       end if
    end do
    !
    if ( .not. Found ) then
       call Assert( "The requested irreducible representation "//trim(IrrepName)//" is not present in the Xlm set." )
    end if
    !
  end function ClassXlmSetGetSymSet


  function ClassXlmSymmetricSetGetIrrep( XlmSymSet ) result( IrrepPtr )
    !
    class(ClassXlmSymmetricSet), intent(in) :: XlmSymSet
    type(ClassIrrep), pointer :: IrrepPtr
    !
    IrrepPtr => XlmSymSet%Irrep
    !
  end function ClassXlmSymmetricSetGetIrrep


  function ClassXlmSymmetricSetGetNXlm( XlmSymSet ) result( NXlm )
    class(ClassXlmSymmetricSet), intent(in) :: XlmSymSet
    integer :: NXlm
    NXlm = XlmSymSet%NumXlm
  end function ClassXlmSymmetricSetGetNXlm


  function ClassXlmSymmetricSetGetXlm( XlmSymSet, iXlm ) result( Xlm )
    class(ClassXlmSymmetricSet), intent(in) :: XlmSymSet
    integer                    , intent(in) :: iXlm
    type(ClassXlm) :: Xlm
    Xlm = XlmSymSet%XlmList(iXlm)
  end function ClassXlmSymmetricSetGetXlm


  subroutine ClassXlmSymmetricSetGetMList( XlmSymSet, l, Mlist ) 
    class(ClassXlmSymmetricSet), intent(in)  :: XlmSymSet
    integer                    , intent(in)  :: l
    integer, allocatable       , intent(out) :: Mlist(:)
    integer :: i, icount
    !
    if(l>XlmSymSet%Lmax)call Assert( &
         " l "//AlphabeticNumber(l)//&
         " > Lmax "//AlphabeticNumber(XlmSymSet%Lmax)//&
         "in ClassXlmSymmetricSetGetMList")
    !
    if(allocated(MList))deallocate(Mlist)
    if(XlmSymSet%dim(l)<=0)return
    allocate(Mlist(XlmSymSet%dim(l)))
    icount=0
    do i = XlmSymSet%idim(l)+1, XlmSymSet%idim(l)+XlmSymSet%dim(l)
       icount=icount+1
       Mlist(icount) = XlmSymSet%XlmList(i)%Getm()
    enddo
    !
  end subroutine ClassXlmSymmetricSetGetMList


  integer function ClassXlmSymmetricSetGetNm( XlmSymSet, l ) result( nm )
    class(ClassXlmSymmetricSet), intent(in)  :: XlmSymSet
    integer                    , intent(in)  :: l
    !
    if(l>XlmSymSet%Lmax)call Assert( &
         " l "//AlphabeticNumber(l)//&
         " > Lmax "//AlphabeticNumber(XlmSymSet%Lmax)//&
         "in ClassXlmSymmetricSetGetMList")
    nm = XlmSymSet%dim(l)
    !
  end function ClassXlmSymmetricSetGetNm

  !> Builds the transformation matrix M, such that a vector of  of Xlm symmetry adapted spherical harmonics X is expressed in terms of a vector C of spherical monomials through: X = M x C.
  subroutine BuildTransMatfromSpherToXlm( Mat, MatInv )
    !
    complex(kind(1d0)),              intent(in)  :: Mat(:,:)
    complex(kind(1d0)), allocatable, intent(out) :: MatInv(:,:)
    !
    complex(kind(1d0)), allocatable :: work(:), MatCopy(:,:)
    integer, allocatable :: ipiv(:)
    integer :: lwork, info, NumRows, i, j
    real(kind(1d0)), parameter :: Threshold = 1.d-12
    !
    allocate( MatCopy, source = Mat )
    !
    NumRows = size(MatCopy,1)
    !
    allocate( ipiv(NumRows) )
    !
    call zgetrf( NumRows, NumRows, MatCopy, NumRows, ipiv, info )
    !
    lwork = NumRows
    allocate( work(lwork) )
    !
    call zgetri( NumRows, MatCopy, NumRows, ipiv, work, lwork, info )
    !
    !
    ! Removes posibly very close to zero elements arising from inversion proces.
    do j = 1, size(MatCopy,2)
       do i = 1, size(MatCopy,1)
          if ( abs(MatCopy(i,j)) < Threshold ) then
             MatCopy(i,j) = Z0
          end if
       end do
    end do
    !
    allocate( MatInv, source = MatCopy )
    !
  end subroutine BuildTransMatfromSpherToXlm

  
  !> Given a monomial indexes matrix, containing in the first, the second and the third rows the x, y and z exponents respectively ( ix, iy, iz ), and in each column storing a diferent monomial; this subroutines reorders the columns in a way that the columns (monomials) are sorted in an increasing semi-monotonic fashion of ix+iy+iz = l. 
  subroutine Reorder( MonInd )
    !
    integer, intent(inout) :: MonInd(:,:)
    !
    integer              :: i, j, k, Counter, NumCart, IER
    integer, allocatable :: Sum(:), iperm(:), OldSum(:), ForbiddenIndexes(:)
    integer, allocatable :: NewMonInd(:,:)
    logical              :: NewMon
    integer, parameter   :: SORT_INCR_ORD = 2
    !
    NumCart = size(MonInd,2)
    !
    allocate( Sum(NumCart) )
    !
    do i = 1, NumCart
       !
       iperm(i)=i
       Sum(i) = MonInd(1,i) + MonInd(2,i) + MonInd(3,i)
       !
    end do
    !
    allocate( OldSum, source = Sum )
    !
    allocate( ForbiddenIndexes(NumCart) )
    ForbiddenIndexes = 0
    !
    allocate( NewMonInd(3,NumCart) )
    !
    !.. More standard routine to sort a vector of integers ..
    call IPSORT(Sum,NumCart,iperm,SORT_INCR_ORD,IER)
    if(IER/=0)then
       write(*,*) "ERROR IN IPSORT IN Reorder"
       stop
    endif
    deallocate(iperm)
    !
    Counter = 0
    !
    first : do i = 1, NumCart
       do j = 1, NumCart
          !
          if ( Sum(i) == OldSum(j) ) then
             !
             NewMon = .false.
             !
             do k = 1, NumCart
                if ( j == ForbiddenIndexes(k) ) exit
                if ( k == NumCart ) NewMon = .true.
             end do
             !
             if ( NewMon ) then
                NewMonInd(:,i) = MonInd(:,j) 
                Counter = Counter + 1
                ForbiddenIndexes(Counter) = j
                cycle first
             end if
             !
          end if
          !
       end do
    end do first
    !
    MonInd = NewMonInd
    deallocate(Sum,OldSum,ForbiddenIndexes,NewMonInd)
    !
  end subroutine Reorder



  !> Given a monomial indexes matrix, containing in the first, the second and the third rows the n,l and m values of r^nYlm  respectively, and in each column storing a diferent monomial; this subroutines reorders the columns in a way that the columns (monomials) are sorted in an increasing order of n, l and m with that priority (i.e.: 000,11-1,110,111,200,22-2,...)% If l is even n has to be even too, the same occurs when it are odd. 
  subroutine SpherReorder( MonInd )
    !
    integer, intent(inout) :: MonInd(:,:)
    !
    integer :: n, l, m, MaxN, Counter, i
    integer :: NumSpher, PosN, PosL, PosM
    integer, allocatable ::  CopyMonInd(:,:), NewMonInd(:,:)
    character(len=32) :: Cstrn, Nstrn
    !
    NumSpher = size(MonInd,2)
    !
    allocate( CopyMonInd, source = MonInd )
    allocate( NewMonInd(3,NumSpher) )
    !
    MaxN = maxval( CopyMonInd(1,:) )
    !
    Counter = 0
    !
    ! Assumes the monomial are not repeated.
    do n = 0, MaxN
       do l = 0, n
          if (mod(n-l,2) == 0 ) then
             do m = -l, l
                !
                do i = 1, NumSpher
                   if ( (n == CopyMonInd(1,i)) .and. &
                        (l == CopyMonInd(2,i)) .and. &
                        (m == CopyMonInd(3,i)) )  then
                      Counter = Counter + 1
                      NewMonInd (1,Counter) = n
                      NewMonInd (2,Counter) = l
                      NewMonInd (3,Counter) = m
                      exit
                   end if
                end do
                !
             end do
          end if
       end do
    end do
    !
    write(Cstrn,*) Counter
    write(Nstrn,*) NumSpher
    if ( Counter /= NumSpher ) then
       call Assert( "The number of spherical monomials after the ordering "//&
               trim(adjustl(Cstrn))//" is not the same than before "//&
               trim(adjustl(Nstrn))//"." )
    end if
    !
    !
    MonInd = NewMonInd
    !
  end subroutine SpherReorder




  subroutine PrintsComplexMatrixOnScreen( Mat )
    !
    complex(kind(1d0)), intent(in) :: Mat(:,:)
    !
    integer :: i, j, unit
    real(kind(1d0)) :: RealPart, ImagPart
    !
    unit = OUTPUT_UNIT
    !
    do j = 1, size(Mat,2)
       do i = 1, size(Mat,1)
          !
          RealPart = dble(Mat(i,j))
          ImagPart = aimag(Mat(i,j))
          write(unit,fmt='(2i2,2f24.16)') i, j, RealPart, ImagPart
          !
       end do
    end do
    !
  end subroutine PrintsComplexMatrixOnScreen



  subroutine PrintsIntegerMatrixOnScreen( Mat )
    !
    integer, intent(in) :: Mat(:,:)
    !
    integer :: i, j, unit
    !
    unit = OUTPUT_UNIT
    !
    do j = 1, size(Mat,2)
       do i = 1, size(Mat,1)
          !
          write(unit,fmt='(2i2,i6)') i, j, Mat(i,j)
          !
       end do
    end do
    !
  end subroutine PrintsIntegerMatrixOnScreen




  subroutine ClassXlmShow( Xlm, unit )
    !
    class(ClassXlm),   intent(in) :: Xlm
    integer, optional, intent(in) :: unit
    !
    integer :: outunit
    !
    if ( present(unit) ) then
       outunit = unit
    else
       outunit = OUTPUT_UNIT
    end if
    !
    write(outunit,fmt="(a,a,i4,4x,a,i4)") "Xlm info : ", "l", Xlm%Getl(), "m", Xlm%Getm()
    !
  end subroutine ClassXlmShow



  subroutine ClassXlmSave( Xlm, unit )
    !
    class(ClassXlm),   intent(in) :: Xlm
    integer,           intent(in) :: unit
    !
    write(unit,*) Xlm%Getl(), Xlm%Getm()
    !
  end subroutine ClassXlmSave


  subroutine ClassXlmLoad( Xlm, unit )
    !
    class(ClassXlm),   intent(inout) :: Xlm
    integer,           intent(in) :: unit
    !
    read(unit,*) Xlm%l, Xlm%m
    !
  end subroutine ClassXlmLoad



  !> Provided an Xlm Class retrieves whether it appears in the Xlm symmetric Set (true) or not (false)% 
  logical function ClassXlmSymmetricSetValidXlm( XlmSymSet, Xlm ) result( Valid )
    !
    class(ClassXlmSymmetricSet), intent(in) :: XlmSymSet
    class(ClassXlm),             intent(in) :: Xlm
    !
    integer :: i
    type(ClassXlm) :: XlmTest
    Valid = .false.
    do i = 1, XlmSymSet%GetNXlm()
       XlmTest = XlmSymSet%GetXlm(i)
       if ( Xlm%Is(XlmTest) ) then
          Valid = .true.
          return
       end if
    end do
    !
  end function ClassXlmSymmetricSetValidXlm


  !> Retrieves True if the two Xlm classes coincide and False if do not.
  logical function ClassXlmIs( Xlm, XlmTest ) result( Same )
    !
    class(ClassXlm), intent(in) :: Xlm
    class(ClassXlm), intent(in) :: XlmTest
    !
    Same = .false.
    if ( (Xlm%GetL() == XlmTest%GetL()) .and. &
         (Xlm%GetM() == XlmTest%GetM()) ) Same = .true. 
    !
  end function ClassXlmIs


  !> Gets the label of the Xlm functions in the form "Xl.m".
  function ClassXlmGetLabel( Xlm ) result( Label )
    !
    class(ClassXlm), intent(in) :: Xlm
    character(len=:), allocatable :: Label
    !
    character(len=32) :: lstrn, mstrn
    !
    write(lstrn,*) Xlm%GetL()
    write(mstrn,*) Xlm%GetM()
    !
    allocate( Label, source = trim(adjustl(lstrn))//"."//trim(adjustl(mstrn)) )
    !
  end function ClassXlmGetLabel



  subroutine ClassXlmSymmetricSetFree( Self )
    class(ClassXlmSymmetricSet), intent(inout) :: Self
    Self%Lmax = -1
    if ( allocated(Self%dim) ) deallocate(Self%dim)
    if ( allocated(Self%idim) ) deallocate(Self%idim)
    Self%Irrep => NULL()
    Self%NumXlm = -1
    if ( allocated(Self%XlmList) ) deallocate(Self%XlmList)
  end subroutine ClassXlmSymmetricSetFree

  subroutine ClassXlmSymmetricSetFinal( Self )
    type(ClassXlmSymmetricSet) :: Self
    call Self%Free()
  end subroutine ClassXlmSymmetricSetFinal

  integer function FindValIndex( IntArray, IntVal )
    !
    integer, intent(in) :: IntArray(:)
    integer, intent(in) :: IntVal
    !
    integer :: i
    !
    ! Meaning it did not find the value
    FindValIndex = -1
    !
    do i = 1, size(IntArray)
       if ( IntVal == IntArray(i) ) then
          FindValIndex = i
          return
       end if
    end do
    !
  end function FindValIndex


  
end module ModuleXlm
