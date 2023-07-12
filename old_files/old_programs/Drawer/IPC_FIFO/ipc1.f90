!.. Fifo's (named pipes) are created with mkfifo name
!   The max pipe buffer size is around 1MB
!   to extend it, it is necessary to invoke the
!   following command as root
!      sudo sysctl fs.pipe-max-size=4194304
!   That said, the default pipe buffer size is likely
!   smaller than the maximum. To read and, if necessary,
!   chage the buffer size for a specific named pipe, 
!   it is necessary to make a call to fcntl with
!   F_GETPIPE_SZ and F_SETPIPE_SZ, respectively
!..


!.. Writing to unformatted file takes about 50% more time
!   than writing to a named pipe. Since the process of
!   writing and reading is sequential, however, it takes
!   twice that amount of time to transfer the data from
!   one process to another, whereas in a FIFO the transfer
!   is synchronous. This means that communication through
!   binary files is three times slower than through pipes.
!
!   Writing to formatted units takes about TWENTY times
!   longer compared to unformatted units. Communication
!   through formatted files, therefore, is about sixty
!   times slower than through (unformatted) named pipes.
!
!   Finally, FIFOs take around 40 nanoseconds to transfer
!   a double across programs (5 nanoseconds per byte).
program Proc1

  implicit none

  character(len=*), parameter :: DATA_FIFO ="data.fifo"
  integer            :: i, data_uid
  integer, parameter :: N = 100000000
  real(kind(1d0))    :: mat(N)

  open(newunit=data_uid,file=DATA_FIFO,access="stream")
  mat=1.d0
  write(*,*) "sending data",sum(mat)
  write(data_uid) (mat(i),i=1,N)
  close(data_uid)

end program Proc1
