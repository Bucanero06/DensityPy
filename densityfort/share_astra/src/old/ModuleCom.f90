module ModuleCom

  logical, parameter :: LOG_INFO = .FALSE.
  logical, parameter :: TER_INFO = .TRUE.
  
  integer :: inp_com_uid
  integer :: log_com_uid
  
  character(len=:), allocatable :: CLIST_INP 
  character(len=:), allocatable :: CLIST_OUT 
  character(len=:), allocatable :: CLIST_NEW 

  character(len=:), allocatable :: CLIST_LOG 

end module ModuleCom
