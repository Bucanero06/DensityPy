!From Luca's hmk_standalone

module chebmod
  implicit none
  integer, parameter :: MAXCHEBORDER=511

  type cheb_series
     integer :: order
     real(kind(1d0)) :: a,b
     real(kind(1d0)) :: c(0:MAXCHEBORDER)
  end type cheb_series

contains

  subroutine set_cheb_zero(cs)
    type(cheb_series) :: cs
    cs%order=0
    cs%a=0.d0
    cs%b=0.d0
    return
  end subroutine set_cheb_zero

  subroutine free_cheb_series(cs)
    type(cheb_series) :: cs
    cs%order=0
    cs%a=0.d0
    cs%b=0.d0
    return
  end subroutine free_cheb_series

  integer function save_cheb_series(c,STREAM)
    type(cheb_series), intent(in) :: c
    integer          , intent(in) :: STREAM
    integer :: stat,i
    logical :: alloc
    save_cheb_series=0
    write(STREAM,IOSTAT=stat)c%order,c%a,c%b
    if(stat/=0)then
       save_cheb_series=1
       return
    endif
       write(STREAM,IOSTAT=stat)(c%c(i),i=0,c%order)
       if(stat/=0)save_cheb_series=2
    return
  end function save_cheb_series

  integer function load_cheb_series(c,STREAM)
    type(cheb_series), intent(out):: c
    integer          , intent(in) :: STREAM
    integer :: stat,i
    logical :: alloc
    load_cheb_series=0
    read(STREAM,IOSTAT=stat)c%order,c%a,c%b
    if(stat/=0)then
       load_cheb_series=1
       return
    endif
       read(STREAM,IOSTAT=stat)(c%c(i),i=0,c%order)
       if(stat/=0)load_cheb_series=2
    return
  end function load_cheb_series

  integer function cheb_copy(c1,c2)
    type(cheb_series), intent(in) :: c1
    type(cheb_series), intent(out):: c2
    cheb_copy=0
    c2%a=c1%a
    c2%b=c1%b
    c2%order=c1%order
    c2%c(0:c2%order)=c1%c(0:c1%order)
    return
  end function cheb_copy

  real(kind(1d0)) function cheb_eval_e(x,cs)
    !Input variables
    real(kind(1d0))  , intent(in) :: x
    type(cheb_series), intent(in) :: cs
    !Local variables
    integer j
    real(kind(1d0)) :: d,dd,y,y2,temp!,e
    d =0.d0
    dd=0.d0
    y =(2.d0*x-cs%a-cs%b)/(cs%b-cs%a)
    y2= 2.d0*y
    do j=cs%order,1,-1
       temp=d
       d=y2*d-dd+cs%c(j)
       dd=temp
    enddo
    d=y*d-dd+0.5d0*cs%c(0)
    cheb_eval_e=d
    return
  end function cheb_eval_e

  real(kind(1d0)) function cheb_eval(x,cs)
    real(kind(1d0))  , intent(in) :: x
    type(cheb_series), intent(in) :: cs
    cheb_eval=cheb_eval_e(x,cs)
    return
  end function cheb_eval
  
  subroutine cheb_init_grid(a,b,order,xgrid)
    real(kind(1d0)), parameter  :: PI=3.1415926535897932d0
    real(kind(1d0)), intent(in) :: a,b
    integer        , intent(in) :: order
    real(kind(1d0)), intent(out):: xgrid(0:order)
    integer :: i
    real(kind(1d0)) :: wp,wm,ropu
    wm=0.5d0*(b-a)
    wp=0.5d0*(b+a)
    ropu=1.d0/dble(order+1)
    do i=0,order
       xgrid(i)=wp+wm*cos(PI*(dble(i)+0.5d0)*ropu)
    enddo
    return
  end subroutine cheb_init_grid
  
  subroutine cheb_init(a,b,order,yg,cs)
    real(kind(1d0))  , parameter  :: PI=3.1415926535897932d0
    real(kind(1d0))  , intent(in) :: a,b
    integer          , intent(in) :: order
    real(kind(1d0))  , intent(in) :: yg(0:order)
    type(cheb_series), intent(out):: cs
    integer :: j, k
    real(kind(1d0)) :: sum,fac,ropu,jr
    if(order>MAXCHEBORDER)call crash
    cs%a=a
    cs%b=b
    cs%order=order
    cs%c=0.d0
    ropu=1.d0/dble(order+1)
    fac =2.d0*ropu
    ropu=PI*ropu
    do j=0,order
       jr=dble(j)
       sum=0.d0
       do k=0,order
          sum=sum+yg(k)*cos(jr*(dble(k)+0.5d0)*ropu)
       enddo
       cs%c(j)=fac*sum
    enddo
    return
  end subroutine cheb_init

end module chebmod
