! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)                               
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!
module bspline_grid_mod
   use precisn
   use utils, only: xermsg
   use mpi_mod
   use const, only: stdout

   private
   
   public bspline_grid_obj, read_BTO, write_BTO, print_BTO

   !> Parameters of the B-spline orbital. The intention is that the user declares this type in the main progam and specifies:
   !> A, B, order, l,m and ind, no_bps. This object can then be passed to bspl_function%init, which checks the parameters
   !> and initializes the parameters in the particular bspl_function. The default values of some of the variables below
   !> have been chosen so that an error will be triggered by the member procedure check if they have not been set to meaningful
   !> values.
   type :: bspline_grid_obj
      !> The range of the radial basis of which this function is a member.
      real(kind=cfp) :: A, B = -1
      !> Number of break-points in the basis.
      integer :: no_bps = -1
      !> Order of the radial B-spline basis.
      integer :: order = -1
      !> Number of B-splines.
      integer :: n
      !> Index of the first radial B-spline whose first derivative at point r=A is 0.
      integer :: ind_0_der = 0
      !> Array of coefficients in the B-spline basis that can be used to evaluate B-spline basis functions.
      real(kind=cfp), allocatable :: bcoef(:)
      !> Number of knots.
      integer :: no_knots = 0
      !> Array of knots.
      real(kind=cfp), allocatable :: knots(:)
      !> Work array used for evaluation of this spline.
      real(kind=cfp), allocatable :: work(:)
      !> Helper variable used for evaluation of B-splines.
      integer :: inbv = 0
      !> Relative precision (tolerance) for numerical calculation of integrals involving this B-spline grid.
      real(kind=cfp) :: tol = 0.0_cfp
   contains
      !> Checks that the parameters of this function are acceptable. Returns 0 on success, positive integer on fail;
      !> the integer is the number of the error message described in the routine. The parameters that are checked are:
      !> A, B, order, l,m, ind and no_bps. All the other data components are set by bspline_grid_obj%init_grid and
      !> bspline_grid_obj%normalize. The only value that is not checked is 'norm'.
      !> \memberof bspline_grid_obj
      procedure :: check
      !> Returns the norm of the radial B-spline with a given index.
      procedure :: normalize
      !> Reads-in the BTO data from the stream file and position given.
      procedure :: read
      !> Writes the BTO data into the stream file and position given.
      procedure :: write
      !> Prints the basic grid parameters to stdout.
      procedure :: print
      !> Prints all grid parameters to stdout.
      procedure :: print_grid
      !> This procedure sets up the grid of radial B-splines given the values: A, B, order, no_bps. The variables initialized
      !> by this routine are: knots, n, no_knots, bcoef, work. Note that it is here where different breakpoint sequences
      !> can be implemented. At the moment only linear sequence is implemented. This routine does not have to be used, i.e.
      !> the user can set-up the grid parameters by hand to whatever grid needed.
      procedure :: init_grid
      !> For a given index of the radial B-spline it returns the values r1, r2 corresponding to its radial range.
      procedure :: bspline_range
      !> Calculates reduced boundary amplitude of a B-spline function (i.e. without the l,m-dependent factors
      !> delta_{m,mp}*delta_{l,lp}). The last argument der is the order of the derivative requried.
      !> For a simple evaluation of the B-spline value set der = 0.
      procedure :: bspline_amplitude
   end type bspline_grid_obj

contains

!ROUTINES for bspline_grid_obj
   function check(this)
      implicit none
      class(bspline_grid_obj) :: this
      integer :: check

         check = 0

         if (this%no_bps .le. 0) then
            check = 1
            return
         endif

         if (this%order .le. 0) then
            check = 2
            return
         endif

         if ( (this%no_bps + this%order - 2) .le. 0) then !total number of B-splines
            check = 3
            return
         endif

         if ((this%A < 0) .or. (this%B .le. 0)) then
            check = 4
            return
         endif

         if (.not.(allocated(this%knots)) .or. .not.(allocated(this%bcoef)) .or. .not.(allocated(this%work))) then
            check = 5
            return
         endif

         if (size(this%knots) .ne. this%no_knots) then
            check = 6
            return
         endif

         if (size(this%bcoef) .ne. this%n) then
            check = 7
            return
         endif

         if (size(this%work) .ne. 3*this%order) then
            check = 8
            return
         endif

   end function check

   subroutine init_grid(this,A,B,order,no_bsplines)
      use bspline_base, only: bvalu
      implicit none
      class(bspline_grid_obj) :: this
      real(kind=cfp), intent(in) :: A,B
      integer, intent(in) :: order,no_bsplines

      integer :: err, i
      real(kind=cfp) :: cfp_dummy, val

         this%A = A
         this%B = B
         this%order = order
         this%n = no_bsplines
         this%no_bps = no_bsplines - order + 2 !number of break-points: we assume maximal possible continuity at the break-points, i.e. order-2

         !calculate the parameters for the grid of radial B-splines:
         this%no_knots = this%n + this%order   !number of knots

         if (allocated(this%knots)) deallocate(this%knots)
         if (allocated(this%bcoef)) deallocate(this%bcoef)
         if (allocated(this%work)) deallocate(this%work)

         allocate(this%knots(1:this%no_knots),this%bcoef(1:this%n),this%work(1:3*this%order),stat=err)
         if (err .ne. 0) call xermsg ('bto_function', 'init_grid', 'Memory allocation failed', err, 1)

         this%work = 0.0_cfp

         !set-up knots so that the B-splines have order-2 continuuous derivatives
         !Important: note that this choice influences the behaviour of the first and last '#order' B-splines. 
         !In particular the first B-spline is discontinuous at r=A and the second has a discont. first derivative at r=A.
         !Note that these properties are used in bspl_function to test whether this B-spline
         !function satisfies the boundary conditions at r=A. For the B-splines that start at r>A we also test that they are smooth and continuous.
         do i = 1,this%order-1
            this%knots(i) = this%A
            this%knots(this%no_knots+1-i) = this%B
         enddo
         
         !uniform spacing of knots
         do i = this%order, this%no_knots-(this%order-1)
            this%knots(i) = min(this%A + (this%B-this%A)/(this%no_bps-1)*(i-this%order),this%B) !min ensures that the knot value does not overflow this%B even by a tiny amount.
         enddo
         !todo temp
         !this%knots(10)=0.255207514048475_cfp
         !this%knots(11)=2.281938789998041_cfp

         !obtain tolerance for numerical quadratures
         this%tol = 10*f1mach(4,cfp_dummy) !max(1.0d-18,f1mach(4,cfp_dummy))

         !Find the first radial B-spline which has zero first derivative at r=A.
         !This is needed to ensure that the Bloch operator is only needed at r=B.
         this%ind_0_der = 0
         this%bcoef = 0.0_cfp
         i = 0
         do
            i=i+1
            if (i > this%n) exit
            this%bcoef(i) = 1.0_cfp
            val = bvalu(this%knots,this%bcoef,this%n,this%order,1,this%A,this%inbv,this%work)
            this%bcoef(i) = 0.0_cfp
            if (val .eq. 0.0_cfp) then
               this%ind_0_der = i
               exit
            endif
         enddo

         if (this % ind_0_der == 0) then
            call xermsg ('bto_function', 'init_grid', &
                         'Bad B-spline basis: there is no radial B-spline in the basis with zero first derivative at r=A.', 1, 1)
         end if

         this%bcoef = 0.0_cfp;

         !finally check that everything has been set-up as required
         err = this%check()
         if (err /= 0) then
            call xermsg ('bspline_grid_obj', 'init_grid', &
                         'The bspline parameters have been set-up incorrectly; see bspline_grid_obj%check().', err, 1)
         end if

         call this%print_grid

   end subroutine init_grid

   function normalize(this,ind)
      use function_integration
      use quadrature_module, only: cfp_bsqad
      implicit none
      class(bspline_grid_obj) :: this
      integer, intent(in) :: ind
      real(kind=cfp) :: normalize

      integer :: i, err
      real(kind=cfp) :: quad, r1, r2
      type(power_function) :: f_poly

         !we assume that the grid parameters have meaningful values so we don't call this%check() but only check here for correctness of the ind value.
         if (ind <= 0 .or. ind > this % n) then
            call xermsg ('bspline_grid_obj', 'normalize', 'On input the value of ind was out of range.', 1, 1)
         end if
     
         !set up the power function for the numerical quadrature
         f_poly%l = 0 !L=0; the normalization is independent of the BTO L.
         i = 0 !no derivative of the B-spline
         !calculate the normalization factor for this orbital; quad = \int_{r1}^{r2} dr (B_{i}(r))^{2}
         this%bcoef = 0.0_cfp
         this%bcoef(ind) = 1.0_cfp
         call this%bspline_range(ind,r1,r2)
         call cfp_bsqad(f_poly, this%knots, this%bcoef, this%n, this%order, i, r1, r2, this%tol, quad, err, this%work)
         if (err /= 1) then
            call xermsg ('bto_function', 'normalize', &
                         'Calculation of normalization of the B-spline orbital does not meet the requested precision.', err, 1)
         end if

         !normalization
         normalize = quad**(-0.5_cfp)

         this%inbv = 1

   end function normalize

   subroutine read(this,lunit,posit,pos_after_rw)
      implicit none
      class(bspline_grid_obj) :: this
      integer, intent(in) :: lunit, posit
      integer, intent(out) :: pos_after_rw

      integer :: err
      real(kind=cfp) :: bcast(1:3)

         if (myrank .eq. master) then
            read(lunit,pos=posit,err=10) this%A, this%B
            read(lunit,err=10) this%no_bps
            read(lunit,err=10) this%order
            read(lunit,err=10) this%n, this%no_knots
            read(lunit,err=10) this%tol
            this%inbv = 0
   
            if (allocated(this%bcoef)) deallocate(this%bcoef)
            if (allocated(this%knots)) deallocate(this%knots)
            if (allocated(this%work)) deallocate(this%work)
            allocate(this%bcoef(1:this%n),this%knots(1:this%no_knots),this%work(1:this%order*3),stat=err)
            if (err .ne. 0) call xermsg ('bspline_grid_obj', 'read', 'Memory allocation has failed.', err, 1)
   
            read(lunit,err=10) this%knots(1:this%no_knots)
            inquire(lunit,pos=pos_after_rw)
         endif

         !master broadcasts all its data to the other processes
         !todo replace by one 9-element array broadcast
         call mpi_mod_bcast(this%no_bps,master)
         call mpi_mod_bcast(this%order,master)
         call mpi_mod_bcast(this%n,master)
         call mpi_mod_bcast(this%no_knots,master)
         call mpi_mod_bcast(pos_after_rw,master)
         this%inbv = 0

         bcast(1:3) = (/this%A,this%B,this%tol/)
         call mpi_mod_bcast(bcast,master)
         this%A = bcast(1); this%B = bcast(2); this%tol = bcast(3)

         if (myrank .ne. master) then
            if (allocated(this%bcoef)) deallocate(this%bcoef)
            if (allocated(this%knots)) deallocate(this%knots)
            if (allocated(this%work)) deallocate(this%work)
            allocate(this%bcoef(1:this%n),this%knots(1:this%no_knots),this%work(1:this%order*3),stat=err)
            if (err .ne. 0) call xermsg ('bspline_grid_obj', 'read', 'Memory allocation 2 has failed.', err, 1)
         endif

         this%bcoef = 0.0_cfp

         call mpi_mod_bcast(this%knots,master)
         
         err = this%check()
         if (err /= 0) then
            call xermsg ('bspline_grid_obj', 'read', &
                         'BTO data read-in but bspline_grid_obj%check() has failed. &
                         &See bspline_grid_obj%check for details.', err, 1)
         end if

         return

 10      call xermsg ('bspline_grid_obj', 'read', 'Error reading the B-spline grid data from the file and position given.', 1, 1)

   end subroutine read

   subroutine write(this,lunit,posit,pos_after_rw)
      implicit none
      class(bspline_grid_obj) :: this
      integer, intent(in) :: lunit, posit
      integer, intent(out) :: pos_after_rw

      integer :: err

         err = this%check()
         if (err /= 0) then
            call xermsg ('bspline_grid_obj', 'write', &
                         'bspline_grid_obj%check() has failed. Erroneous data will not be written. &
                         &See bspline_grid_obj%check for details.', err, 1)
         end if

         if (myrank .eq. master) then
            write(lunit,pos=posit,err=10) this%A, this%B
            write(lunit,err=10) this%no_bps
            write(lunit,err=10) this%order
            write(lunit,err=10) this%n, this%no_knots
            write(lunit,err=10) this%tol
            write(lunit,err=10) this%knots(1:this%no_knots)
            inquire(lunit,pos=pos_after_rw)
         endif

         !master ensures all processes know where the record ends
         call mpi_mod_bcast(pos_after_rw,master)

         return

 10      call xermsg ('bspline_grid_obj', 'write', 'Error writing the B-spline grid data into the file and position given.', 1, 1)

   end subroutine write

   subroutine print(this)
      implicit none
      class(bspline_grid_obj) :: this

      integer :: err

         err = this%check()
         if (err /= 0) then
            call xermsg ('bspline_grid_obj', 'print', &
                         'bspline_grid_obj%check() has failed. Erroneous data will not be written. &
                         &See bspline_grid_obj%check for details.', err, 1)
         end if

         write(stdout,'("Parameters of the B-spline grid:")')
         write(stdout,'("A, B: ",2e20.10)') this%A, this%B
         write(stdout,'("Order of the radial B-splines: ",i0)') this%order
         write(stdout,'("Number of break-points: ",i0)') this%no_bps
         write(stdout,'("Number of B-splines in the basis: ",i0)') this%n

   end subroutine print

   subroutine print_grid(this)
      implicit none
      class(bspline_grid_obj) :: this

      integer :: err, i

         err = this%check()
         if (err /= 0) then
            call xermsg ('bspline_grid_obj', 'print', &
                         'bspline_grid_obj%check() has failed. Erroneous data will not be written. &
                         &See bspline_grid_obj%check for details.', err, 1)
         end if

         write(stdout,'("Parameters of the B-spline grid:")')
         write(stdout,'("A, B: ",2e20.10)') this%A, this%B
         write(stdout,'("Order of the radial B-splines: ",i0)') this%order
         write(stdout,'("Number of break-points: ",i0)') this%no_bps
         write(stdout,'("Number of B-splines in the basis: ",i0)') this%n

         write(stdout,'(/,"Array of knots: ")')
         do i=1,size(this%knots)
            write(stdout,'(i5,e25.15)') i,this%knots(i)
         enddo !i

   end subroutine print_grid

   subroutine bspline_range(this,ind,r1,r2)
      implicit none
      class(bspline_grid_obj) :: this
      integer, intent(in) :: ind
      real(kind=cfp), intent(out) :: r1, r2

         if (ind <= 0 .or. ind > this % n) then
            call xermsg ('bspline_grid_obj', 'write', 'On input the value of ind was out of range.', 1, 1)
         end if

         if (.not. allocated(this % knots)) then
            call xermsg ('bspline_grid_obj', 'write', 'The array of knots has not been allocated.', 2, 1)
         end if

         r1 = this%knots(ind)
         r2 = this%knots(ind+this%order)

   end subroutine bspline_range

   subroutine read_BTO(bspline_grid,l,bspline_index,number_of_functions,norm,non_zero_at_boundary,lunit,posit,pos_after_rw)
      implicit none
      class(bspline_grid_obj) :: bspline_grid
      integer, intent(in) :: lunit, posit
      integer, intent(out) :: l, bspline_index, number_of_functions, pos_after_rw
      real(kind=cfp), intent(out) :: norm
      logical, intent(out) :: non_zero_at_boundary

      integer :: start

         if (myrank .eq. master) then
            read(lunit,pos=posit,err=10) l
            read(lunit,err=10) bspline_index
            read(lunit,err=10) number_of_functions
            read(lunit,err=10) norm
            read(lunit,err=10) non_zero_at_boundary
            inquire(lunit,pos=pos_after_rw)
         endif

         !master broadcasts all its data to the other processes
         call mpi_mod_bcast(l,master)
         call mpi_mod_bcast(bspline_index,master)
         call mpi_mod_bcast(number_of_functions,master)
         call mpi_mod_bcast(norm,master)
         call mpi_mod_bcast(non_zero_at_boundary,master)
         call mpi_mod_bcast(pos_after_rw,master)

         start = pos_after_rw
         call bspline_grid%read(lunit,start,pos_after_rw)

         return

 10      call xermsg ('bspline_grid_obj', 'read_BTO', 'Error reading the BTO data from the file and position given.', 1, 1)

   end subroutine read_BTO

   subroutine write_BTO(bspline_grid,l,bspline_index,number_of_functions,norm,non_zero_at_boundary,lunit,posit,pos_after_rw)
      implicit none
      class(bspline_grid_obj) :: bspline_grid
      integer, intent(in) :: l, bspline_index, number_of_functions, lunit, posit
      integer, intent(out) :: pos_after_rw
      real(kind=cfp), intent(in) :: norm
      logical, intent(in) :: non_zero_at_boundary

      integer :: start

         if (myrank .eq. master) then
            write(lunit,pos=posit,err=10) l
            write(lunit,err=10) bspline_index
            write(lunit,err=10) number_of_functions
            write(lunit,err=10) norm
            write(lunit,err=10) non_zero_at_boundary
            inquire(lunit,pos=pos_after_rw)
         endif

         call mpi_mod_bcast(pos_after_rw,master)

         start = pos_after_rw
         call bspline_grid%write(lunit,start,pos_after_rw)

         return

 10      call xermsg ('bspline_grid_obj', 'write_BTO', 'Error writing the BTO data from the file and position given.', 1, 1)

   end subroutine write_BTO

   subroutine print_BTO(bspline_grid,l,bspline_index,number_of_functions,norm,non_zero_at_boundary)
      implicit none
      class(bspline_grid_obj) :: bspline_grid
      integer, intent(in) :: l, bspline_index, number_of_functions
      real(kind=cfp), intent(in) :: norm
      logical, intent(in) :: non_zero_at_boundary

         write(stdout,'("L: ",i0)') l
         write(stdout,'("Index of the radial B-spline, total number of radial B-splines &
                      &generated by the present grid: ",i0,1x,i0)') bspline_index, bspline_grid%n
         write(stdout,'("Normalization factor: ",e20.10)') norm
         write(stdout,'("Is non-zero at the boundary: ",l)') non_zero_at_boundary

         call bspline_grid%print

   end subroutine print_BTO

   function bspline_amplitude(this,r,norm,ind,der)
      use bspline_base, only: bvalu
      implicit none
      class(bspline_grid_obj) :: this
      real(kind=cfp), intent(in) :: r, norm
      real(kind=cfp) :: bspline_amplitude
      integer, intent(in) :: ind
      integer, intent(in) :: der

      real(kind=cfp) :: r1, r2

         bspline_amplitude = 0.0_cfp

         call this%bspline_range(ind,r1,r2)

         if (r .ge. r1 .and. r .le. r2) then
            this%bcoef = 0.0_cfp
            this%bcoef(ind) = 1.0_cfp
            bspline_amplitude = norm*bvalu(this%knots,this%bcoef,this%n,this%order,der,r,this%inbv,this%work)
         endif

   end function bspline_amplitude

end module bspline_grid_mod
