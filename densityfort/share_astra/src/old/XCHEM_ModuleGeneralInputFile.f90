
!> \file
!!
!! Contains the variables and routines devoted to read the [GeneralInputFile](@ref GeneralInputFile).

module moduleGeneralInputFile

   use ModuleParameterList

   implicit none

   private

   type, public :: ClassSAEGeneralInputData
      private
    contains
      generic, public :: ParseFile => ClassSAEGeneralInputDataParseFile
      generic, public :: GetStore => ClassSAEGeneralInputDataGetStore
      generic, public :: GetBasisFile => ClassSAEGeneralInputDataGetBasisFile
      ! {{{ Private Procedures
      procedure, private :: ClassSAEGeneralInputDataParseFile
      procedure, private :: ClassSAEGeneralInputDataGetStore 
      procedure, private :: ClassSAEGeneralInputDataGetBasisFile
      ! }}}
   end type ClassSAEGeneralInputData 
   type(ClassSAEGeneralInputData), public :: SAEGeneralInputData

   
   !> Maximum length of a name to be used.
   integer, parameter :: NAME_LENGTH=500

   !> Directory where the results of the programs will be stored in.
   character(len=NAME_LENGTH), public :: storedir
   !> Radial basis configuration file ([Basis](@ref Basis)).
   character(len=NAME_LENGTH), public :: BasisFile
   !> Absorption potential configuration file ([AbsorptionPotential](@ref AbsorptionPotential)).
   character(len=NAME_LENGTH), public :: AbsorptionPotentialFile
   !> [Time propagation program](@ref Propagation1) configuration file.
   character(len=NAME_LENGTH), public :: propf
   !> Configuration file for the analysis of time propagation results program.
   character(len=NAME_LENGTH), public :: aresf


   public ParseGeneralInputFile

contains


  subroutine ClassSAEGeneralInputDataParseFile( this, ConfigurationFile )
    Class(ClassSAEGeneralInputData), intent(in) :: this
    character(len=*), intent(in) :: ConfigurationFile
    call ParseGeneralInputFile( ConfigurationFile )
  end subroutine ClassSAEGeneralInputDataParseFile


  subroutine ClassSAEGeneralInputDataGetBasisFile( this, SAEBasisFile )
    Class(ClassSAEGeneralInputData), intent(in) :: this
    character(len=:), allocatable, intent(out) :: SAEBasisFile
    allocate(SAEBasisFile, source = BasisFile)
  end subroutine ClassSAEGeneralInputDataGetBasisFile


  subroutine ClassSAEGeneralInputDataGetStore( this, store )
    Class(ClassSAEGeneralInputData), intent(in) :: this
    character(len=:), allocatable, intent(out) :: store
    allocate(store,source=storedir)
  end subroutine ClassSAEGeneralInputDataGetStore


  !> Parses the [GeneralInputFile](@ref GeneralInputFile).
  subroutine ParseGeneralInputFile( ConfigurationFile )
    !
    use, intrinsic :: ISO_FORTRAN_ENV
    !> [GeneralInputFile](@ref GeneralInputFile).
    character(len=*) , intent(in) :: ConfigurationFile
    !
    type(ClassParameterList) :: List
    !
    call List%Add( "StoreDirectory" ,"A"          , "required" )
    call List%Add( "BasisFile"      ,"Basis"      , "optional" )
    call List%Add( "AbsorptionFile" ,"Absorption" , "optional" )
    call List%Add( "PropagationFile","Propagation", "optional" )
    call List%Add( "AnalysisFile"   ,"Analysis"   , "optional" )
    !
    call List%Parse( ConfigurationFile )
    !
    call List%Get( "StoreDirectory" , storedir  )
    call List%Get( "BasisFile"      , BasisFile )
    call List%Get( "AbsorptionFile" , AbsorptionPotentialFile     )
    call List%Get( "PropagationFile", propf     )
    call List%Get( "AnalysisFile"   , aresf     )
    !
  end subroutine ParseGeneralInputFile

end module
