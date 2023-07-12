program Proc2

  implicit none

  character(len=*), parameter :: DATA_FIFO ="data.fifo"
  integer            :: i, data_uid
  integer, parameter :: N = 100000000
  real(kind(1d0))    :: mat(N)

  open(newunit=data_uid,file=DATA_FIFO,access="stream")
  read(data_uid) (mat(i),i=1,N)
  write(*,*)"data received", sum(mat)
  close(data_uid)

end program Proc2
