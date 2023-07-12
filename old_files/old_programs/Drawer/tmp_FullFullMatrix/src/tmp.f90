
use intrinsic :: iso_c_binding

enum, bind(c) :: orbital_type
    enumerator :: ORB_TYPE_MO
    enumerator :: ORB_TYPE_HYBRID
    enumerator :: ORB_TYPE_BSPLINE_INT
    enumerator :: ORB_TYPE_BSPLINE_EXT 
end enum

integer, parameter :: ORB_TYPE_MO     = 1
integer, parameter :: ORB_TYPE_HYBRID = 2


integer :: myorbital_type

logical :: isMO
logical :: is..



if( myorbital_type == ORB_TYPE_MO )then

  .....

else(...


select case (myorbital_type)

case( ORB_TYPE_MO )
  ....
case( ORB_TYPE_BSPLINE_INT ) 
  ....



