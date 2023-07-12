!> Determines, by bisections, the interval in 
!! which x falls
!! =>   0,   if to the left of all points
!! =>   i,   if comprised between x_i and x_{i+1}, i < NDAT
!! =>  NDAT, if to the right of all points
!<
function WhichInterval(x,NDAT,XDAT) result (interval)
  implicit none
  real(kind(1d0)), intent(in) :: x
  integer        , intent(in) :: NDAT
  real(kind(1d0)), intent(in) :: XDAT(:)
  integer :: interval
  !
  integer :: imi, ima, ime
  !
  if(x<XDAT(1))then
     imi=0
  elseif(x>XDAT(NDAT))then
     imi=NDAT
  else
     imi=1
     ima=NDAT
     do while(ima>imi+1)
        ime=(imi+ima+1)/2
        if(x>=XDAT(ime))then
           imi=ime
        else
           ima=ime
        endif
     enddo
  endif
  interval=imi
end function WhichInterval


!> Evaluates the linear interpolant of the data set
!! at a given value of the independent variable
!<
complex(kind(1d0)) function LinEval(x,XDAT,ZDAT,NDAT) result (res)
  implicit none
  real(kind(1d0)), intent(in)  :: x
  real(kind(1d0)), intent(in)  :: XDAT(:)
  complex(kind(1d0)), intent(in)  :: ZDAT(:)
  integer        , intent(in)  :: NDAT
  interface
     function WhichInterval(x, NDAT, XDAT) result( interval )
       real(kind(1d0)), intent(in) :: x
       integer        , intent(in) :: NDAT
       real(kind(1d0)), intent(in) :: XDAT(:)
     end function WhichInterval
  end interface
  real(kind(1d0)), parameter :: eps=1.d-24
  integer :: imi

  !.. Extrapolates below the first and above the last point,
  !   Interpolates between points
  !..
  imi = min(max(WhichInterval( x, NDAT, XDAT ),1),NDAT-1)

  res=ZDAT(imi)+(x-XDAT(imi))/(XDAT(imi+1)-XDAT(imi))*(ZDAT(imi+1)-ZDAT(imi))

end function LinEval


subroutine InterpolateAmpl( xv1, za1, xv2, za2 )
  
  implicit none
  
  real   (kind(1d0)), intent(in) :: xv1(:)
  complex(kind(1d0)), intent(in) :: za1(:,:)
  real   (kind(1d0)), intent(in) :: xv2(:)
  complex(kind(1d0)), intent(out):: za2(:,:)
  
  interface
     complex(kind(1d0)) function LinEval(x,XDAT,YDAT,NDAT) result (res)
       real(kind(1d0)), intent(in)  :: x
       real(kind(1d0)), intent(in)  :: XDAT(:)
       complex(kind(1d0)), intent(in) :: YDAT(:)
       integer        , intent(in)  :: NDAT
     end function LinEval
  end interface

  integer :: irow, icol, ncol, nrow
  
  if(size(za1,2).ne.size(za2,2))then
     write(*,*) "Inconsistent number of columns in InterpolateAmpl"
     stop
  endif
  if(size(xv1,1).ne.size(za1,1))then
     write(*,*) "Inconsistent number of rows in #1 entries of InterpolateAmpl"
     stop
  endif
  if(size(xv2,1).ne.size(za2,1))then
     write(*,*) "Inconsistent number of rows in #2 entries of InterpolateAmpl"
     stop
  endif

  ncol = size(za1,2)
  nrow = size(xv1,1)
  do icol = 1, ncol
     do irow = 1, size(xv2,1)
        za2(irow,icol) = LinEval(xv2(irow),xv1,za1(:,icol),nrow)
     enddo
  enddo

end subroutine InterpolateAmpl




