 Module ModuleDensityMatricesJUAN

  use ModuleBasisJUAN
  implicit none

  !Multiplication matrix for the irr
  integer, dimension(8,8)   :: irrmult=(/ 1, 2, 3, 4, 5, 6, 7, 8, &
                                          2, 1, 4, 3, 6, 5, 8, 7, &
                                          3, 4, 1, 2, 7, 8, 5, 6, &
                                          4, 3, 2, 1, 8, 7, 6, 5, &
                                          5, 6, 7, 8, 1, 2, 3, 4, &
                                          6, 5, 8, 7, 2, 1, 4, 3, &
                                          7, 8, 5, 6, 3, 4, 1, 2, &
                                          8, 7, 6, 5, 4, 3, 2, 1 /)

  type, public :: ClassDensityMatrices

     !character(len=120)          :: filename
     !> Directory in which the results will be stored
     character(len=:)             , allocatable :: StorageDir
     Class(BasisElementInfo), allocatable :: Orbitals(:)
     integer                              :: Nbasis,irr_sizes(8),imin_irr(8),imax_irr(8)
     real(kind(1d0)), allocatable         :: RhoDM(:,:),PiDM(:,:,:,:)

  ! D2h   Ag  B3u B2u B1g B1u B2g B3g Au
  !
  ! Ag    Ag  B3u B2u B1g B1u B2g B3g Au
  ! B3u   B3u Ag  B1g B2u B2g B1u Au  B3g
  ! B2u   B2u B1g Ag  B3u B3g Au  B1u B2g
  ! B1g   B1g B2u B3u Ag  Au  B3g B2g B1u 
  ! B1u   B1u B2g B3g Au  Ag  B3u B2u B1g  
  ! B2g   B2g B1u Au  B3g B3u Ag  B1g B2u
  ! B3g   B3g Au  B1u B2g B2u B1g Ag  B3u
  ! Au    Au  B3g B2g B1u B1g B2u B3u Ag
     
   contains

    procedure :: InitOrb                  => InitializeOrbitalsData
    procedure :: rho_stex                 => RhoDensityMatrixSTEX
    procedure :: pi_stex                  => PiDensityMatrixSTEX
    procedure :: Ga_stex                  => GammaDensityMatrixSTEX
     

     procedure, private :: InitializeOrbitalsData
  end type ClassDensityMatrices

contains


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !         Procedures   ClassDensityMatrices
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


  subroutine InitializeOrbitalsData( Basis,Orbitals,Nbasis )
    class(ClassDensityMatrices), intent(inout) :: Basis
    integer,                     intent(in)    :: Nbasis
    Class(BasisElementInfo),     intent(in) :: Orbitals(1:Ndim)
    Basis.Nbasis=Nbasis
    Basis.irr_sizes=irr_total*2   !spin
    Basis.imin_irr=imin_irr
    Basis.imax_irr=imax_irr
    allocate(Basis.Orbitals, source = Orbitals)
  end subroutine InitializeOrbitalsData

  real(kind(1.d0)) function RhoDensityMatrixSTEX(DM,iB,iA,iQ,iP) result(Element)
    Class(ClassDensityMatrices), intent(in) :: DM
    integer           , intent(in) :: iA, iB, iQ, iP

    Element =          delta(iA,iB)*delta(iP,iQ)*(1.d0-delta(iA,iP))
    Element = Element -delta(iA,iQ)*delta(iP,iB)*(1.d0-delta(iA,iB))*(1.d0-delta(iP,iQ))
    ! Element=0.d0
    ! If((iA.eq.iB).and.(i.eq.j).and.(iA.ne.j))Then
    !    Element=1.d0
    ! End IF
    ! If((iA.ne.iB).and.(i.ne.j))Then
    !    If((iA.eq.i).and.(iB.eq.j))Then
    !       Element=Element-1.d0
    !    End IF
    ! End IF
  end function RhoDensityMatrixSTEX

  real(kind(1.d0)) function PiDensityMatrixSTEX(DM,iB,iA,iR,iS,iP,iQ) result(Element)
    Class(ClassDensityMatrices), intent(in) :: DM
    integer           , intent(in) :: iA,iB,iR,iS,iP,iQ
    
    Element = delta(iA,iB)*(1.d0-delta(iA,iP))*(1.d0-delta(iA,iQ))*(1.d0-delta(iP,iQ))* &
         (delta(iP,iR)*delta(iQ,iS)-delta(iQ,iR)*delta(iP,iS))
    Element = Element + (1.d0 - delta(iA,iB))*delta(iQ,iB)*(1.d0-delta(iA,iP))*(1.d0-delta(iP,iB)) * &
         (delta(iA,iR)*delta(iP,iS)-delta(iP,iR)*delta(iA,iS))
    Element = Element + delta(iP,iB)*(1.d0 - delta(iA,iB))*(1.d0-delta(iA,iQ))*(1.d0-delta(iB,iQ)) * &
         (delta(iQ,iR)*delta(iA,iS)-delta(iA,iR)*delta(iQ,iS))


  end function PiDensityMatrixSTEX

  

  real(kind(1.d0)) function GammaDensityMatrixSTEX(DM,iA,iB,iR,iS,iT,iO,iP,iQ) result(Element)
    Class(ClassDensityMatrices), intent(in) :: DM
    integer           , intent(in) :: iA,iB,iR,iS,iT,iO,iP,iQ

    Element = (1.d0-delta(iA,iO))*(1.d0-delta(iA,iP))*(1.d0-delta(iA,iQ))*                 &
         (1.d0-delta(iO,iP))*(1.d0-delta(iO,iQ))*(1.d0-delta(iO,iP))*(1.d0-delta(iP,iQ))*  &
         (1.d0-delta(iB,iR))*(1.d0-delta(iB,iS))*(1.d0-delta(iB,iT))*(1.d0-delta(iR,iS))*  &
         (1.d0-delta(iR,iT))*(1.d0-delta(iR,iT))*(1.d0-delta(iS,iT))
    if(nint(abs(Element)).eq.1)then
       
       Element =           delta(iA,iB)*delta(iO,iR)*(delta(iP,iS)*delta(iQ,iT)-delta(iQ,iS)*delta(iP,iT))
       Element = Element + delta(iA,iB)*delta(iP,iR)*(delta(iQ,iS)*delta(iO,iT)-delta(iO,iS)*delta(iQ,iT))
       Element = Element + delta(iA,iB)*delta(iQ,iR)*(delta(iO,iS)*delta(iP,iT)-delta(iP,iS)*delta(iO,iT))
       
       Element = Element + delta(iO,iB)*delta(iP,iR)*(delta(iA,iS)*delta(iQ,iT)-delta(iQ,iS)*delta(iA,iT))
       Element = Element + delta(iO,iB)*delta(iQ,iR)*(delta(iP,iS)*delta(iA,iT)-delta(iA,iS)*delta(iP,iT))
       Element = Element + delta(iO,iB)*delta(iA,iR)*(delta(iQ,iS)*delta(iP,iT)-delta(iP,iS)*delta(iQ,iT))
       
       Element = Element + delta(iP,iB)*delta(iQ,iR)*(delta(iA,iS)*delta(iO,iT)-delta(iO,iS)*delta(iA,iT))
       Element = Element + delta(iP,iB)*delta(iA,iR)*(delta(iO,iS)*delta(iQ,iT)-delta(iQ,iS)*delta(iO,iT))
       Element = Element + delta(iP,iB)*delta(iO,iR)*(delta(iQ,iS)*delta(iA,iT)-delta(iA,iS)*delta(iQ,iT))
       
       Element = Element + delta(iQ,iB)*delta(iA,iR)*(delta(iP,iS)*delta(iO,iT)-delta(iO,iS)*delta(iP,iT))
       Element = Element + delta(iQ,iB)*delta(iP,iR)*(delta(iO,iS)*delta(iA,iT)-delta(iA,iS)*delta(iO,iT))
       Element = Element + delta(iQ,iB)*delta(iO,iR)*(delta(iA,iS)*delta(iP,iT)-delta(iP,iS)*delta(iA,iT))
      
    endif


    
  end function GammaDensityMatrixSTEX

    real(kind(1d0)) function delta(i,j) result(res)
      integer, intent(in) :: i, j
      res = 1.d0
      if(i.ne.j) res = 0.d0
    end function delta




  
end module ModuleDensityMatricesJUAN

!This is the order found in ukrmol (same order than dalton):
!#      Ag    Xlm       l even m even  N =  l/2+1      1  ->  1
!#      B3u   Xlm       l odd  m odd   N = (l+1)/2     2  ->  8
!#      B2u   Xl-m      l odd  m odd   N = (l+1)/2     3  ->  7
!#      B1g   Xl-m      l even m even  N =  l/2        4  ->  2 
!#      B1u   Xlm       l odd  m even  N = (l+1)/2     5  ->  6
!#      B2g   Xlm       l even m odd   N =  l/2        6  ->  3
!#      B3g   Xl-m      l even m odd   N =  l/2        7  ->  4
!#      Au    Xl-m      l odd  m even  N = (l-1)/2     8  ->  5
!This is the order found in CloseCouplingBasis file
!#      Ag    Xlm       l even m even  N =  l/2+1
!#      B1g   Xl-m      l even m even  N =  l/2
!#      B2g   Xlm       l even m odd   N =  l/2
!#      B3g   Xl-m      l even m odd   N =  l/2
!#      Au    Xl-m      l odd  m even  N = (l-1)/2
!#      B1u   Xlm       l odd  m even  N = (l+1)/2
!#      B2u   Xl-m      l odd  m odd   N = (l+1)/2
!#      B3u   Xlm       l odd  m odd   N = (l+1)/2

!This is from   XCHEM_ModuleGroups.f90
  ! D2h   Ag  B1g B2g B3g   Au  B1u B2u B3u
  !
  ! Ag    Ag  B1g B2g B3g   Au  B1u B2u B3u
  ! B1g   B1g Ag  B3g B2g   B1u Au  B3u B2u
  ! B2g   B2g B3g Ag  B1g   B2u B3u Au  B1u
  ! B3g   B3g B2g B1g Ag    B3u B2u B1u Au
  !
  ! Au    Au  B1u B2u B3u   Ag  B1g B2g B3g
  ! B1u   B1u Au  B3u B2u   B1g Ag  B3g B2g
  ! B2u   B2u B3u Au  B1u   B2g B3g Ag  B1g
  ! B3u   B3u B2u B1u Au    B3g B2g B1g Ag

!The table according to dalton and ukrmol would be

  ! D2h   Ag  B3u B2u B1g B1u B2g B3g Au
  !
  ! Ag    Ag  B3u B2u B1g B1u B2g B3g Au
  ! B3u   B3u Ag  B1g B2u B2g B1u Au  B3g
  ! B2u   B2u B1g Ag  B3u B3g Au  B1u B2g
  ! B1g   B1g B2u B3u Ag  Au  B3g B2g B1u 
  ! B1u   B1u B2g B3g Au  Ag  B3u B2u B1g  
  ! B2g   B2g B1u Au  B3g B3u Ag  B1g B2u
  ! B3g   B3g Au  B1u B2g B2u B1g Ag  B3u
  ! Au    Au  B3g B2g B1u B1g B2u B3u Ag

