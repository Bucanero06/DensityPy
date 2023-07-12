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
module free_scattering_mod
   use precisn
   use mpi_mod, only: mpi_xermsg, mpi_mod_barrier
   use molecular_basis_mod
   implicit none

   private

   public free_scattering

contains

   !> Performs the free potential scattering given the location of the transformed overlap and kinetic energy integrals, orbital basis set and R-matrix radius.
   subroutine free_scattering(mo_integral_storage,molecular_orbital_basis,a,min_energy,max_energy,nE)
      use const, only: stdout, kinetic_ints, one_p_sym_ints
      use phys_const, only: to_ev
      use integral_storage_mod
      use molecular_basis_mod
      use blas_lapack, only: blasint, syev
      use parallel_arrays
      implicit none
      type(integral_storage_obj), intent(inout) :: mo_integral_storage
      type(molecular_orbital_basis_obj), intent(inout) :: molecular_orbital_basis
      real(kind=cfp), intent(in) :: a, min_energy, max_energy
      integer, intent(in) :: nE

      integer :: err, i, j, k, channel, sym, ind(1:1), sym_start, sym_end
      type(p2d_array_obj), pointer :: mo_integrals
      type(p2d_array_obj), target :: ao_integrals
      type(integral_storage_obj), target :: mo_integral_disk
      type(data_header_obj) :: header
      type(integral_options_obj) :: options
      real(kind=cfp), allocatable :: ham(:,:), bamps_all(:,:)
      integer, allocatable :: chan_lm_all(:,:), l_list(:)
      logical, parameter :: tgt_is_local = .true.

      real(kind=cfp), allocatable :: work(:), e(:), w(:,:), rmatrix(:,:), kmatrix(:,:), eigphase(:)
      real(kind=cfp) :: energy, p
      integer(blasint) :: lwork, info, n, n_chan
      integer :: tot_chan, col_kei
      logical, parameter :: normalize_to_a = .true.

         write(stdout,'(/,"--------->free_scattering_mod:free_scattering")')

         if (min_energy < 0.0_cfp .or. max_energy < min_energy .or. nE < 1) then
            call mpi_xermsg ('free_scattering_mod', 'free_scattering', &
                             'Incorrect energy grid parameters on input (min_energy,max_energy,nE).', 1, 1)
         endif

         call mpi_mod_barrier(err)

         if (mo_integral_storage%in_memory()) then
            mo_integrals => mo_integral_storage%integral
         else !on disk
            err = mo_integral_disk%init(memory=ao_integrals)
            if (err .ne. 0) call mpi_xermsg('free_scattering_mod','free_scattering','Memory allocation 1 failed',err,1)

            !we need to get the full header corresponding to the transformed 1-electron integrals.
            header%name = mo_integral_storage%contruct_header_string(molecular_orbital_basis%get_basis_name(),one_p_sym_ints)

            !read the integrals from the disk into mo_integrals
            call mo_integral_disk%read(mo_integral_storage,header%name,options,tgt_is_local)

            mo_integrals => ao_integrals
         endif

         !Get the column number for the mo_integrals%a array where the KEI are stored.
         col_kei = mo_integrals%find_column_matching_name(kinetic_ints)

         !FROM THIS POINT ONWARDS mo_integrals POINTS TO THE STORAGE IN MEMORY OF THE KINETIC ENERGY INTEGRALS

         !EVALUATE ALL BOUNDARY AMPLITUDES
         call molecular_orbital_basis%calculate_amplitudes(a,normalize_to_a,bamps_all,chan_lm_all)

         tot_chan = size(chan_lm_all,2) !total number of channels

         !LOOP OVER ALL SYMMETRIES

         do sym=1,molecular_orbital_basis%no_irr

            !BUILD THE HAMILTONIAN = KINETIC ENERGY MATRIX

            write(stdout,'(/,"Symmetry: ",i1)') sym

            n = molecular_orbital_basis%get_number_of_orbitals(sym)

            if (n .eq. 0) then
               write(stdout,'("No orbitals in this symmetry: skipping")')
               cycle
            endif

            if (allocated(ham)) deallocate(ham)
            allocate(ham(n,n),stat=err)
            if (err .ne. 0) call mpi_xermsg('free_scattering_mod','free_scattering','Memory allocation 2 failed',err,1)

            ham = 0.0_cfp

            sym_start = 1
            do i=1,sym-1
               sym_start = sym_start+molecular_orbital_basis%get_number_of_orbitals(i)
            enddo
            sym_end = sym_start-1 + n

            do i=sym_start,sym_end
               do j=sym_start,i
                  ind(1:1) = molecular_orbital_basis%integral_index(kinetic_ints,reshape(source=(/j,i/),shape=(/2,1/)))
                  ham(j-sym_start+1,i-sym_start+1) = mo_integrals%a(ind(1),col_kei)
                  ham(i-sym_start+1,j-sym_start+1) = ham(j-sym_start+1,i-sym_start+1)
               enddo
            enddo

            !DIAGONALIZE THE HAMILTONIAN

            if (allocated(work)) deallocate(work)
            if (allocated(e)) deallocate(e)
            allocate(work(1:1),e(1:n),stat=err)
            if (err .ne. 0) call mpi_xermsg('free_scattering_mod','free_scattering','Memory allocation 3 failed',err,1)
            !determine the size of the auxiliary array WORK needed for diagonalization of H+L
            lwork = -1
            call syev('V','U',n,ham,n,e,work,lwork,info)
            if (info /= 0) then
                call mpi_xermsg ('free_scattering_mod', 'free_scattering', &
                                 'SYEV call 1 failed with an error or you are linking against 32bit integer MKL &
                                 &while using 64bit default integer.', int(info), 1)
            end if

            lwork=work(1)
            deallocate(work); allocate(work(1:lwork),stat=err)
            if (err .ne. 0) call mpi_xermsg('free_scattering_mod','free_scattering','Memory allocation 4 failed',err,1)
   
            !eigenvectors of H+L are normalized; eigenvectors are held in columns of ham
            call syev('V','U',n,ham,n,e,work,lwork,info)
            if (info /= 0) then
                call mpi_xermsg ('free_scattering_mod', 'free_scattering', &
                                 'SYEV call 2 failed with an error or you are linking against 32bit integer MKL &
                                 &while using 64bit default integer.', int(info), 1)
            end if

            write(stdout,'(/,"Eigenvalues of the Kinetic energy matrix:")')
            do i=1,n
               write(stdout,'(i0,e25.15)') i,e(i)
               if (e(i) < 0.0_cfp) then
                  call mpi_xermsg ('free_scattering_mod', 'free_scattering', &
                                   'At least one eigenvalue of the kinetic energy matrix is negative: &
                                   &severe numerical instability or R-matrix radius is too small.', 2, 1)
               end if
            enddo

            !GET CHANNEL DATA AND CONSTRUCT THE TRUE BOUNDARY AMPLITUDES

            if (allocated(w)) deallocate(w)
            if (allocated(l_list)) deallocate(l_list)

            write(stdout,'(/,"Partial waves (m,l) in the current symmetry:")')
            n_chan = 0
            do i=1,tot_chan
               if (chan_lm_all(3,i) .eq. sym) then 
                  n_chan = n_chan + 1
                  write(stdout,'(2i5)') chan_lm_all(1:2,i)
               endif
            enddo

            write(stdout,'(/,"Number of channels in the current symmetry: ",i0)') n_chan

            if (n_chan .eq. 0) cycle

            allocate(w(n_chan,n),l_list(n_chan),stat=err)
            if (err .ne. 0) call mpi_xermsg('free_scattering_mod','free_scattering','Memory allocation 5 failed',err,1)

            w(:,:) = 0.0_cfp
            channel = 0
            do i=1,tot_chan
               if (chan_lm_all(3,i) .eq. sym) then 
                  channel = channel + 1
                  l_list(channel) = chan_lm_all(2,i) !save the L-values for each channel
                  !calculate the R-matrix boundary amplitudes in each channel
                  do k=1,n
                     do j=1,n !over all coefficients of the R-matrix basis functions
                        w(channel,k) = w(channel,k) + ham(j,k)*bamps_all(i,sym_start-1+j)
                        !write(stdout,'("cf",i10,2e25.15)') j,ham(j,k),bamps_all(i,sym_start-1+j)
                     enddo
                     !write(stdout,'("amplitude",2i10,e25.15)') channel,k,w(channel,k)
                  enddo
               endif
            enddo

            if (allocated(rmatrix)) deallocate(rmatrix)
            if (allocated(kmatrix)) deallocate(kmatrix)
            if (allocated(eigphase)) deallocate(eigphase)

            allocate(rmatrix(1:n_chan,1:n_chan),kmatrix(1:n_chan,1:n_chan),eigphase(1:n_chan),stat=err)
            if (err .ne. 0) call mpi_xermsg('free_scattering_mod','free_scattering','Memory allocation 6 failed',err,1)

            if (allocated(work)) deallocate(work)
            lwork=3*n_chan-1
            allocate(work(1:lwork),stat=err)
            if (err .ne. 0) call mpi_xermsg('free_scattering_mod','free_scattering','Memory allocation 7 failed',err,1)

            write(stdout,'(/,"EIGENPHASE SUMS")')

            do i=1,nE+1
               energy = min_energy + (i-1)*(max_energy-min_energy)/real(nE,kind=cfp) !from min_energy [H] to max_energy [H]: nE+1 points
               p = sqrt(2.0_cfp*energy)
               call calc_rmat_e(rmatrix,energy,a,w,e,int(n),int(n_chan)) !get rmatrix
               call calc_kmat_e(kmatrix,a,energy,int(n_chan),l_list,rmatrix) !get kmatrix
  
               call syev('N','U',n_chan,kmatrix,n_chan,eigphase,work,lwork,info) !diagonalize the k-matrix
               if (info /= 0) then
                  call mpi_xermsg ('free_scattering_mod', 'free_scattering', &
                                   'SYEV call 3 failed with an error or you are linking against 32bit integer MKL &
                                   &while using 64bit default integer.', int(info), 1)
               end if

               write(stdout,'(e25.15,e25.15)') energy*to_ev,sum(atan(eigphase))
            enddo

         enddo !sym

         call mpi_mod_barrier(err)

         write(stdout,'(/,"<---------done:free_scattering_mod:free_scattering")')

   end subroutine free_scattering

   subroutine calc_rmat_e(rmatrix,e,a0,u,ek,n,n_chan)
      implicit none
      integer, intent(in) :: n_chan, n
      real(kind=cfp), intent(in) :: u(:,:), ek(:), a0, e
      real(kind=cfp), intent(out) :: rmatrix(:,:)
 
      integer :: i, j, k
 
         rmatrix = 0.0_cfp
         do i=1,n_chan
            do j=1,n_chan
               do k=1,n
                  rmatrix(j,i) = rmatrix(j,i) + u(j,k)*u(i,k)/(ek(k)-e)
               enddo
            enddo
         enddo
         rmatrix = rmatrix/a0*0.5_cfp
 
   end subroutine calc_rmat_e
 
   subroutine calc_kmat_e(kmatrix,a0,energy,n_chan,l_list,rmatrix)
      use phys_const, only: pi
      use blas_lapack, only: blasint, trsm
      implicit none
      integer, intent(in) :: n_chan, l_list(:)
      real(kind=cfp), intent(in) :: a0, energy, rmatrix(:,:)
      real(kind=cfp), intent(out) :: kmatrix(:,:)
 
      integer :: i,j,k
      integer(kind=4) :: l
      integer(blasint) :: n
      real(kind=cfp) :: bes, dbes, bes2, dbes2, p, rb(1:n_chan), rn(1:n_chan), drb(1:n_chan), drn(1:n_chan), &
                        A_mat(1:n_chan,1:n_chan), B_mat(1:n_chan,1:n_chan), f(1:n_chan,1:n_chan,1:2)
      real(kind=cfp) :: alpha, df, fp(1:n_chan,1:n_chan,1:2), mv, xv
 
 
         bes=0.0_cfp; dbes=0.0_cfp; bes2=0.0_cfp; dbes2=0.0_cfp
         
         rb=0.0_cfp !riccatti-bessel function
         rn=0.0_cfp !riccatti-neumann function
         !and their first derivatives:
         drb=0.0_cfp 
         drn=0.0_cfp
    
         p=sqrt(2.0_cfp*energy)
         
         if (p.ne.0.0_cfp) then
            do i=1,n_chan
               bes=0.0_cfp; dbes=0.0_cfp; bes2=0.0_cfp; dbes2=0.0_cfp
               l = l_list(i)
               !returns bessel functions y(pr) of the first (order=l+1/2) and second kind (order=-l-1/2) and their first derivatives d/dr[y(pr)]
               call bessel(a0,p,l,0.5_cfp,bes,dbes,bes2,dbes2)
            
               rb(i)=sqrt(pi*p*a0/2.0_cfp)*bes   !riccati-bessel function
               rn(i)=-sqrt(pi*p*a0/2.0_cfp)*bes2 !riccati-neumann function 
            
               drb(i)=pi*p/4.0_cfp*((pi*p*a0/2.0_cfp)**(-0.5_cfp))*bes+sqrt(pi*p*a0/2.0_cfp)*dbes     !derivative of the r-b function
               drn(i)=-pi*p/4.0_cfp*((pi*p*a0/2.0_cfp)**(-0.5_cfp))*bes2-sqrt(pi*p*a0/2.0_cfp)*dbes2  !derivative of the r-n function
            enddo
 
            f = 0.0_cfp
            fp = 0.0_cfp
            do i=1,n_chan
               f(i,i,1) = rb(i)
               fp(i,i,1) = drb(i)
               do j=1,n_chan
                  f(j,i,2) = rn(i)
                  fp(j,i,2) = drn(i)
               enddo
            enddo
 
            A_mat = 0.0_cfp
            B_mat = 0.0_cfp
            do j=1,n_chan
               do i=1,n_chan
                  A_mat(i,j) = f(i,j,2)
               enddo
            enddo
            
            do j=1,n_chan
               do k=1,n_chan
                  df = fp(k,j,2) !-bsto(k)*f(k,j,2)
                  do i=1,n_chan
                     A_mat(i,j) = A_mat(i,j)-a0*rmatrix(i,k)*df
                  enddo
               enddo
            enddo
            
            do j=1,n_chan
               do i=1,n_chan
                  B_mat(i,j)=-f(i,j,1)
               enddo
            enddo
            
            do j=1,n_chan
               do k=1,n_chan
                  df = fp(k,j,1) !-bsto(k)*f(k,j,1)
                  do i=1,n_chan
                     B_mat(i,j)=B_mat(i,j)+a0*rmatrix(i,k)*df
                  enddo
               enddo
            enddo
 
            alpha = 1.0_cfp
            n = n_chan
            call trsm('L','U','N','N',n,n,alpha,A_mat,n,B_mat,n)           
            kmatrix = B_mat
            xv = maxval(kmatrix)
            mv = minval(kmatrix)
            xv = max(abs(xv),abs(mv))
         else
            kmatrix = 0.0_cfp
         endif
 
   end subroutine calc_kmat_e

   subroutine bessel(r,p,nu,nufrac,j,dj,y,dy)
   !returns the values of the bessel functions of the first and second kind of order 'nu+nufrac' at 'r*p' and its derivatives over 'r' there,
   !where 'nu' is the integer part of the order and 'nufrac' is the fractional part of the order.
      implicit none
      real(kind=cfp) :: j, dj, r, nufrac, y, dy,p,z
      integer(kind=4) :: nu, nb,ncalc
      real(kind=cfp),allocatable :: b(:)
      
      ncalc=0
      
      nb=nu+2
      allocate(b(1:nb))
      z=r*p
      
      call rjbesl(z,nufrac,nb,b,ncalc)     !bessel functions of the first kind of various orders
      
      if ((ncalc.ne.nb).and.(r.ne.0.0_cfp)) then
         print *, 'ncalc =',ncalc
         call mpi_xermsg('free_scattering_mod','bessel','Error in calling rjbesl',1,1)
      endif
      
      j=b(nb-1)                           !bessel function of the first kind of the order we want
      dj=-p*b(nb)+(nu+nufrac)/r*j         !derivative of the bessel function (from forward recurrence relations):
      
      call rybesl(z,nufrac,nb,b,ncalc) !bessel functions of the second kind of various orders
      
      if ((ncalc.ne.nb).and.(r.ne.0.0_cfp)) then
         print *, 'ncalc =',ncalc
         call mpi_xermsg('free_scattering_mod','bessel','Error in calling rybesl',2,1)
      endif
      
      y=b(nb-1)                           !bessel function of the second kind of the order we want
      dy=-p*b(nb)+(nu+nufrac)/r*y         !derivative of the bessel function (from forward recurrence relations):
      
      deallocate(b)
   
   end subroutine bessel


   subroutine rybesl ( x, alpha, nb, by, ncalc )
   
   !*****************************************************************************80
   !
   !! RYBESL calculates Y Bessel function with non-integer orders.
   !
   !  Discussion:
   !
   !    This routine calculates Bessel functions Y SUB(N+ALPHA) (X)
   !    for non-negative argument X, and non-negative order N+ALPHA.
   !
   !    This program draws heavily on Temme's Algol program for Y(a,x)
   !    and Y(a+1,x) and on Campbell's programs for Y_nu(x).  Temme's
   !    scheme is used for x < THRESH, and Campbell's scheme is used
   !    in the asymptotic region.  Segments of code from both sources
   !    have been translated into Fortran77, merged, and heavily modified.
   !    Modifications include parameterization of machine dependencies,
   !    use of a new approximation for ln(gamma(x)), and built-in
   !    protection against over/underflow.
   !
   !    In case of an error, NCALC .NE. NB, and not all Y's are
   !    calculated to the desired accuracy.
   !
   !    NCALC < -1:  An argument is out of range. For example,
   !    NB <= 0, IZE is not 1 or 2, or IZE=1 and ABS(X) .GE.
   !    XMAX.  In this case, BY(1) = 0.0, the remainder of the
   !    BY-vector is not calculated, and NCALC is set to
   !    MIN0(NB,0)-2  so that NCALC .NE. NB.
   !
   !    NCALC = -1:  Y(ALPHA,X) .GE. XINF.  The requested function
   !    values are set to 0.0.
   !
   !    1 < NCALC < NB: Not all requested function values could
   !    be calculated accurately.  BY(I) contains correct function
   !    values for I <= NCALC, and and the remaining NB-NCALC
   !    array elements contain 0.0.
   !
   !  Modified:
   !
   !    03 April 2007
   !
   !  Author:
   !
   !    William Cody
   !
   !  Reference:
   !
   !    JB Campbell,
   !    Bessel functions J_nu(x) and Y_nu(x) of real order and real argument,
   !    Computational Physics Communications,
   !    Volume 18, 1979, pages 133-142.
   !
   !    NM Temme,
   !    On the numerical evaluation of the ordinary Bessel function
   !    of the second kind,
   !    Journal of Computational Physics,
   !    Volume 21, 1976, pages 343-350.
   !
   !  Parameters:
   !
   !    Input, real ( kind = 8 ) X, the argument.  0 <= X.
   !
   !    Input, real ( kind = 8 ) ALPHA, the fractional part of the order
   !    for which the Y's are to be calculated.  0 <= ALPHA < 1.0.
   !
   !    Input, integer ( kind = 4 ) NB, the number of functions to be calculated, NB .GT. 0.
   !    The first function calculated is of order ALPHA, and the
   !    last is of order (NB - 1 + ALPHA).
   !
   !    Output, real ( kind = 8 ) BY(NB).  If the routine terminates normally
   !    (NCALC=NB), the vector BY contains the functions Y(ALPHA,X) through
   !    Y(NB-1+ALPHA,X),  If (0 < NCALC < NB), BY(I) contains correct
   !    function values for I <= NCALC, and contains the ratios
   !    Y(ALPHA+I-1,X)/Y(ALPHA+I-2,X) for the rest of the array.
   !
   !    Output, integer ( kind = 4 ) NCALC, error flag.  Before using the vector BY, the
   !    user should check that NCALC=NB, i.e., all orders have been calculated
   !    to the desired accuracy.
   !
     implicit none
   
     integer ( kind = 4 ) nb
   
     real(kind=cfp) alfa
     real(kind=cfp) alpha
     real(kind=cfp) aye
     real(kind=cfp) b
     real(kind=cfp) by(nb)
     real(kind=cfp) c
     real(kind=cfp) ch(21)
     real(kind=cfp) cosmu
     real(kind=cfp) d
     real(kind=cfp) del
     real(kind=cfp) den
     real(kind=cfp) ddiv
     real(kind=cfp) div
     real(kind=cfp) dmu
     real(kind=cfp) d1
     real(kind=cfp) d2
     real(kind=cfp) e
     real(kind=cfp) eight
     real(kind=cfp) en
     real(kind=cfp) enu
     real(kind=cfp) en1
     real(kind=cfp) eps
     real(kind=cfp) even
     real(kind=cfp) ex
     real(kind=cfp) f
     real(kind=cfp) fivpi
     real(kind=cfp) g
     real(kind=cfp) gamma
     real(kind=cfp) h
     real(kind=cfp) half
     integer ( kind = 4 ) i
     integer ( kind = 4 ) k
     integer ( kind = 4 ) na
     integer ( kind = 4 ) ncalc
     real(kind=cfp) odd
     real(kind=cfp) onbpi
     real(kind=cfp) one
     real(kind=cfp) one5
     real(kind=cfp) p
     real(kind=cfp) pa
     real(kind=cfp) pa1
     real(kind=cfp) pi
     real(kind=cfp) piby2
     real(kind=cfp) pim5
     real(kind=cfp) q
     real(kind=cfp) qa
     real(kind=cfp) qa1
     real(kind=cfp) q0
     real(kind=cfp) r
     real(kind=cfp) s
     real(kind=cfp) sinmu
     real(kind=cfp) sq2bpi
     real(kind=cfp) ten9
     real(kind=cfp) term
     real(kind=cfp) three
     real(kind=cfp) thresh
     real(kind=cfp) two
     real(kind=cfp) twobyx
     real(kind=cfp) x
     real(kind=cfp) xinf
     real(kind=cfp) xlarge
     real(kind=cfp) xmin
     real(kind=cfp) xna
     real(kind=cfp) x2
     real(kind=cfp) ya
     real(kind=cfp) ya1
     real(kind=cfp) zero
     logical, save :: first = .true.
   !
   !  Mathematical constants
   !    FIVPI = 5*PI
   !    PIM5 = 5*PI - 15
   !    ONBPI = 1/PI
   !    PIBY2 = PI/2
   !    SQ2BPI = SQUARE ROOT OF 2/PI
   !
     data zero / 0.0d0 /
     data half / 0.5d0 /
     data one / 1.0d0 /
     data two / 2.0d0 /
     data three / 3.0d0 /
     data eight /8.0d0 /
     data one5 / 15.0d0 /
     data ten9 / 1.9d1/
     data fivpi /1.5707963267948966192d1 /
     data piby2 / 1.5707963267948966192d0/
     data pi /3.1415926535897932385d0 /
     data sq2bpi / 7.9788456080286535588d-1/
     data pim5 /7.0796326794896619231d-1/
     data onbpi / 3.1830988618379067154d-1/
   !
   !  Machine-dependent constants
   !
     data del / 1.0d-8 /
     data xmin / 4.46d-308 /
     data xinf / 1.79d308 /
     data eps / 1.11d-16 /
     data thresh / 16.0d0 /
     data xlarge / 1.0d8 /
   !
   !  Coefficients for Chebyshev polynomial expansion of
   !  1/gamma(1-x), abs(x) <= .5
   !
     data ch/-0.67735241822398840964d-23,-0.61455180116049879894d-22, &
              0.29017595056104745456d-20, 0.13639417919073099464d-18, &
              0.23826220476859635824d-17,-0.90642907957550702534d-17, &
             -0.14943667065169001769d-14,-0.33919078305362211264d-13, &
             -0.17023776642512729175d-12, 0.91609750938768647911d-11, &
              0.24230957900482704055d-09, 0.17451364971382984243d-08, &
             -0.33126119768180852711d-07,-0.86592079961391259661d-06, &
             -0.49717367041957398581d-05, 0.76309597585908126618d-04, &
              0.12719271366545622927d-02, 0.17063050710955562222d-02, &
             -0.76852840844786673690d-01,-0.28387654227602353814d+00, &
              0.92187029365045265648d+00/

     if (cfp .eq. ep1 .and. first) then
        call mpi_xermsg ('free_scattering', 'rybesl', &
                         'Quad precision version not implemented: precision degradation into double. &
                         &This warning message will be supressed for subsequent calls to rybesl.', 1, 0)
        first = .false.
     endif
   
     ex = x
     enu = alpha
   
     if ( 0 < nb .and. &
       xmin <= x .and. &
       ex < xlarge .and. &
       zero <= enu .and. &
       enu < one )  then
   
       xna = aint ( enu + half )
       na = int ( xna )
   
       if ( na == 1 ) then
         enu = enu - xna
       end if
   
       if ( enu == -half ) then
   
         p = sq2bpi / sqrt ( ex )
         ya = p * sin ( ex )
         ya1 = -p * cos ( ex )
   !
   !  Use Temme's scheme for small X.
   !
       else if ( ex < three ) then
   
         b = ex * half
         d = - log ( b )
         f = enu * d
         e = b**( -enu )
   
         if ( abs ( enu ) < del ) then
           c = onbpi
         else
           c = enu / sin ( enu * pi )
         end if
   !
   !  Computation of sinh(f)/f.
   !
         if ( abs ( f ) < one ) then
           x2 = f * f
           en = ten9
           s = one
           do i = 1, 9
             s = s * x2 / en / ( en - one ) + one
             en = en - two
           end do
         else
           s = ( e - one / e ) * half / f
         end if
   !
   !  Computation of 1/gamma(1-a) using Chebyshev polynomials.
   !
         x2 = enu * enu * eight
         aye = ch(1)
         even = zero
         alfa = ch(2)
         odd = zero
   
         do i = 3, 19, 2
           even = -( aye + aye + even )
           aye = -even * x2 - aye + ch(i)
           odd = -( alfa + alfa + odd )
           alfa = -odd * x2 - alfa + ch(i+1)
         end do
   
         even = ( even * half + aye ) * x2 - aye + ch(21)
         odd = ( odd + alfa ) * two
         gamma = odd * enu + even
   !
   !  End of computation of 1/gamma(1-a).
   !
         g = e * gamma
         e = ( e + one / e ) * half
         f = two * c * ( odd * e + even * s * d )
         e = enu * enu
         p = g * c
         q = onbpi / g
         c = enu * piby2
   
         if ( abs ( c ) < del ) then
           r = one
         else
           r = sin ( c ) / c
         end if
   
         r = pi * c * r * r
         c = one
         d = - b * b
         h = zero
         ya = f + r * q
         ya1 = p
         en = zero
   
     100     continue
   
         en = en + one
   
         if ( eps < abs ( g / ( one + abs ( ya ) ) ) &
           + abs ( h / ( one + abs ( ya1 ) ) ) ) then
           f = ( f * en + p + q ) / ( en * en - e )
           c = c * d / en
           p = p / ( en - enu )
           q = q / ( en + enu )
           g = c * ( f + r * q )
           h = c * p - en * g
           ya = ya + g
           ya1 = ya1 + h
           go to 100
         end if
   
         ya = -ya
         ya1 = -ya1 / b
   
       else if ( ex < thresh ) then
   !
   !  Use Temme's scheme for moderate X.
   !
         c = ( half - enu ) * ( half + enu )
         b = ex + ex
         e = ( ex * onbpi * cos ( enu * pi ) / eps )
         e = e * e
         p = one
         q = -ex
         r = one + ex * ex
         s = r
         en = two
   
     200     continue
   
         if ( r * en * en < e ) then
           en1 = en + one
           d = ( en - one + c / en ) / s
           p = ( en + en - p * d ) / en1
           q = ( -b + q * d ) / en1
           s = p * p + q * q
           r = r * s
           en = en1
           go to 200
         end if
   
         f = p / s
         p = f
         g = -q / s
         q = g
   
     220     continue
   
         en = en - one
   
         if ( zero < en ) then
           r = en1 * ( two - p ) - two
           s = b + en1 * q
           d = ( en - one + c / en ) / ( r * r + s * s )
           p = d * r
           q = d * s
           e = f + one
           f = p * e - g * q
           g = q * e + p * g
           en1 = en
           go to 220
         end if
   
         f = one + f
         d = f * f + g * g
         pa = f / d
         qa = -g / d
         d = enu + half -p
         q = q + ex
         pa1 = ( pa * q - qa * d ) / ex
         qa1 = ( qa * q + pa * d ) / ex
         b = ex - piby2 * ( enu + half )
         c = cos ( b )
         s = sin ( b )
         d = sq2bpi / sqrt ( ex )
         ya = d * ( pa * s + qa * c )
         ya1 = d * ( qa1 * s - pa1 * c )
       else
   !
   !  Use Campbell's asymptotic scheme.
   !
         na = 0
         d1 = aint ( ex / fivpi )
         i = int ( d1 )
         dmu = (( ex - one5 * d1 ) - d1 * pim5 ) &
           - ( alpha + half ) * piby2
   
         if ( i - 2 * ( i / 2 ) == 0 ) then
           cosmu = cos ( dmu )
           sinmu = sin ( dmu )
         else
           cosmu = -cos ( dmu )
           sinmu = -sin ( dmu )
         end if
   
         ddiv = eight * ex
         dmu = alpha
         den = sqrt ( ex )
   
         do k = 1, 2
   
           p = cosmu
           cosmu = sinmu
           sinmu = -p
           d1 = ( two * dmu - one ) * ( two * dmu + one )
           d2 = zero
           div = ddiv
           p = zero
           q = zero
           q0 = d1 / div
           term = q0
   
           do i = 2, 20
   
             d2 = d2 + eight
             d1 = d1 - d2
             div = div + ddiv
             term = -term * d1 / div
             p = p + term
             d2 = d2 + eight
             d1 = d1 - d2
             div = div + ddiv
             term = term * d1 / div
             q = q + term
   
             if ( abs ( term ) <= eps ) then
               exit
             end if
   
           end do
   
           p = p + one
           q = q + q0
   
           if ( k == 1 ) then
             ya = sq2bpi * ( p * cosmu - q * sinmu ) / den
           else
             ya1 = sq2bpi * ( p * cosmu - q * sinmu ) / den
           end if
   
           dmu = dmu + one
   
         end do
   
       end if
   
       if ( na == 1 ) then
         h = two * ( enu + one ) / ex
         if ( one < h ) then
           if ( xinf / h < abs ( ya1 ) ) then
             h = zero
             ya = zero
           end if
         end if
         h = h * ya1 - ya
         ya = ya1
         ya1 = h
       end if
   !
   !  Now have first one or two Y's.
   !
       by(1) = ya
       by(2) = ya1
   
       if ( ya1 == zero ) then
   
         ncalc = 1
   
       else
   
         aye = one + alpha
         twobyx = two / ex
         ncalc = 2
   
         do i = 3, nb
   
           if ( twobyx < one ) then
             if ( xinf / aye <= abs ( by(i-1) ) * twobyx ) then
               go to 450
             end if
           else
             if ( xinf / aye / twobyx <= abs ( by(i-1) ) ) then
               go to 450
             end if
           end if
   
           by(i) = twobyx * aye * by(i-1) - by(i-2)
           aye = aye + one
           ncalc = ncalc + 1
   
         end do
   
       end if
   
     450   continue
   
       do i = ncalc + 1, nb
         by(i) = zero
       end do
   
     else
   
       by(1) = zero
       ncalc = min ( nb, 0_4 ) - 1
   
     end if
   
     return
   end subroutine
   
   subroutine rjbesl ( x, alpha, nb, b, ncalc )
   
   !*****************************************************************************80
   !
   !! RJBESL calculates J Bessel function with non-integer orders.
   !
   !  Discussion:
   !
   !    This routine calculates Bessel functions J sub(N+ALPHA) (X)
   !    for non-negative argument X, and non-negative order N+ALPHA.
   !
   !    This program is based on a program written by David Sookne
   !    that computes values of the Bessel functions J or I of real
   !    argument and integer order.  Modifications include the restriction
   !    of the computation to the J Bessel function of non-negative real
   !    argument, the extension of the computation to arbitrary positive
   !    order, and the elimination of most underflow.
   !
   !    In case of an error,  NCALC .NE. NB, and not all J's are
   !    calculated to the desired accuracy.
   !
   !    NCALC < 0:  An argument is out of range. For example,
   !    NBES <= 0, ALPHA < 0 or .GT. 1, or X is too large.
   !    In this case, B(1) is set to zero, the remainder of the
   !    B-vector is not calculated, and NCALC is set to
   !    MIN(NB,0)-1 so that NCALC .NE. NB.
   !
   !    NB .GT. NCALC .GT. 0: Not all requested function values could
   !    be calculated accurately.  This usually occurs because NB is
   !    much larger than ABS(X).  In this case, B(N) is calculated
   !    to the desired accuracy for N <= NCALC, but precision
   !    is lost for NCALC < N <= NB.  If B(N) does not vanish
   !    for N .GT. NCALC (because it is too small to be represented),
   !    and B(N)/B(NCALC) = 10**(-K), then only the first NSIG-K
   !    significant figures of B(N) can be trusted.
   !
   !  Modified:
   !
   !    03 April 2007
   !
   !  Author:
   !
   !    William Cody
   !
   !  Reference:
   !
   !    Frank Olver, David Sookne,
   !    A Note on Backward Recurrence Algorithms,
   !    Mathematics of Computation,
   !    Volume 26, 1972, pages 941-947.
   !
   !    David Sookne,
   !    Bessel Functions of Real Argument and Integer Order,
   !    NBS Journal of Res. B,
   !    Volume 77B, 1973, pages 125-132.
   !
   !  Parameters:
   !
   !    Input, real ( kind = 8 ) X, the argument for which the
   !    J's are to be calculated.
   !
   !    Input, real ( kind = 8 ) ALPHA, the fractional part of order for which
   !    the J's or exponentially scaled J's (J*exp(X)) are to be calculated.
   !    0 <= ALPHA < 1.0.
   !
   !    Input, integer ( kind = 4 ) NB, the number of functions to be calculated.
   !    0 < NB.  The first function calculated is of order ALPHA, and the
   !    last is of order (NB - 1 + ALPHA).
   !
   !    Output, real ( kind = 8 ) B(NB).  If RJBESL terminates normally, with
   !    NCALC = NB, then B contains the functions J/ALPHA/(X) through
   !    J/NB-1+ALPHA/(X), or the corresponding exponentially scaled functions.
   !
   !    Output, integer ( kind = 4 ) NCALC, error indicator.  If NCALC = NB, then all the
   !    requested values were calculated to the desired accuracy.
   !
   !  Local Parameters:
   !
   !    IT, the number of bits in the mantissa of a working precision
   !    variable.
   !
   !    NSIG, the decimal significance desired.  Should be set to
   !    INT(LOG10(2)*IT+1).  Setting NSIG lower will result
   !    in decreased accuracy while setting NSIG higher will
   !    increase CPU time without increasing accuracy.  The
   !    truncation error is limited to a relative error of
   !    T=.5*10**(-NSIG).
   !
   !    ENTEN = 10.0**K, where K is the largest integer such that
   !    ENTEN is machine-representable in working precision
   !
   !    ENSIG = 10.0**NSIG
   !
   !    RTNSIG = 10.0**(-K) for the smallest integer K such that
   !    K .GE. NSIG/4
   !
   !    ENMTEN, the smallest ABS(X) such that X/4 does not underflow
   !
   !    XLARGE, the upper limit on the magnitude of X.  If ABS(X)=N,
   !    then at least N iterations of the backward recursion
   !    will be executed.  The value of 10000.0 is used on
   !    every machine.
   !
     implicit none
   
     integer ( kind = 4 ) nb
   
     real(kind=cfp) alpha
     real(kind=cfp) alpem
     real(kind=cfp) alp2em
     real(kind=cfp) b(nb)
     real(kind=cfp) capp
     real(kind=cfp) capq
     real(kind=cfp) eighth
     real(kind=cfp) em
     real(kind=cfp) en
     real(kind=cfp) enmten
     real(kind=cfp) ensig
     real(kind=cfp) enten
     real(kind=cfp) fact(25)
     real(kind=cfp) four
     real(kind=cfp) gnu
     real(kind=cfp) half
     real(kind=cfp) halfx
     integer ( kind = 4 ) i
     integer ( kind = 4 ) j
     integer ( kind = 4 ) k
     integer ( kind = 4 ) l
     integer ( kind = 4 ) m
     integer ( kind = 4 ) magx
     integer ( kind = 4 ) n
     integer ( kind = 4 ) nbmx
     integer ( kind = 4 ) ncalc
     integer ( kind = 4 ) nend
     integer ( kind = 4 ) nstart
     real(kind=cfp) one
     real(kind=cfp) one30
     real(kind=cfp) p
     real(kind=cfp) pi2
     real(kind=cfp) plast
     real(kind=cfp) pold
     real(kind=cfp) psave
     real(kind=cfp) psavel
     real(kind=cfp) rtnsig
     real(kind=cfp) s
     real(kind=cfp) sum
     real(kind=cfp) t
     real(kind=cfp) t1
     real(kind=cfp) tempa
     real(kind=cfp) tempb
     real(kind=cfp) tempc
     real(kind=cfp) test
     real(kind=cfp) three
     real(kind=cfp) three5
     real(kind=cfp) tover
     real(kind=cfp) two
     real(kind=cfp) twofiv
     real(kind=cfp) twopi1
     real(kind=cfp) twopi2
     real(kind=cfp) x
     real(kind=cfp) xc
     real(kind=cfp) xin
     real(kind=cfp) xk
     real(kind=cfp) xlarge
     real(kind=cfp) xm
     real(kind=cfp) vcos
     real(kind=cfp) vsin
     real(kind=cfp) z
     real(kind=cfp) zero
     logical, save :: first = .true.
   !
   !  Mathematical constants
   !
   !   PI2    - 2 / PI
   !   TWOPI1 - first few significant digits of 2 * PI
   !   TWOPI2 - (2*PI - TWOPI) to working precision, i.e.,
   !            TWOPI1 + TWOPI2 = 2 * PI to extra precision.
   !
     data pi2 / 0.636619772367581343075535d0 /
     data twopi1 / 6.28125d0 /
     data twopi2 / 1.935307179586476925286767d-3 /
     data zero /0.0d0 /
     data eighth / 0.125d0 /
     data half / 0.5d0 /
     data one / 1.0d0/
     data two /2.0d0 /
     data three / 3.0d0 /
     data four / 4.0d0 /
     data twofiv /25.0d0/
     data one30 /130.0d0 /
     data three5 / 35.0d0/
   !
   !  Machine-dependent parameters
   !
     data enten /1.0d38 /
     data ensig / 1.0d17 /
     data rtnsig / 1.0d-4/
     data enmten /1.2d-37 /
     data xlarge / 1.0d4/
   !
   !  Factorial(N)
   !
     data fact / &
      1.0d0, &
      1.0d0, &
      2.0d0, &
      6.0d0, &
      24.0d0, &
      1.2d2, &
      7.2d2, &
      5.04d3, &
      4.032d4, &
      3.6288d5,3.6288d6,3.99168d7,4.790016d8,6.2270208d9, &
      8.71782912d10,1.307674368d12,2.0922789888d13,3.55687428096d14, &
      6.402373705728d15,1.21645100408832d17,2.43290200817664d18, &
      5.109094217170944d19,1.12400072777760768d21, &
      2.585201673888497664d22, &
      6.2044840173323943936d23/

     if (cfp .eq. ep1 .and. first) then
        call mpi_xermsg ('free_scattering', 'rjbesl', &
                         'Quad precision version not implemented: precision degradation into double. &
                         &This warning message will be supressed for subsequent calls to rjbesl.', 1, 0)
        first = .false.
     endif
   !
   !  Check for out of range arguments.
   !
     magx = int ( x )
   
     if ( &
       0 < nb .and. &
       zero <= x .and. &
       x <= xlarge .and. &
       zero <= alpha .and. &
       alpha < one ) then
   !
   !  Initialize result array to zero.
   !
       ncalc = nb
       do i = 1, nb
         b(i) = zero
       end do
   !
   !  Branch to use 2-term ascending series for small X and asymptotic
   !  form for large X when NB is not too large.
   !
       if ( x < rtnsig ) then
   !
   !  Two-term ascending series for small X.
   !
         tempa = one
         alpem = one + alpha
         halfx = zero
   
         if ( enmten < x ) then
           halfx = half * x
         end if
   
         if ( alpha /= zero ) then
           tempa = halfx**alpha / ( alpha * gamma ( alpha ) )
         end if
   
         tempb = zero
   
         if ( one < x + one ) then
           tempb = -halfx * halfx
         end if
   
         b(1) = tempa + tempa * tempb / alpem
   
         if ( x /= zero .and. b(1) == zero ) then
           ncalc = 0
         end if
   
         if ( nb /= 1 ) then
   
           if ( x <= zero ) then
   
             do n = 2, nb
               b(n) = zero
             end do
   !
   !  Calculate higher order functions.
   !
           else
   
             tempc = halfx
             tover = ( enmten + enmten ) / x
   
             if ( tempb /= zero ) then
               tover = enmten / tempb
             end if
   
             do n = 2, nb
   
               tempa = tempa / alpem
               alpem = alpem + one
               tempa = tempa * tempc
   
               if ( tempa <= tover * alpem ) then
                 tempa = zero
               end if
   
               b(n) = tempa + tempa * tempb / alpem
   
               if ( b(n) == zero .and. n < ncalc ) then
                 ncalc = n - 1
               end if
   
             end do
   
           end if
         end if
   !
   !  Asymptotic series for 21 < X.
   !
       else if ( twofiv < x .and. nb <= magx + 1 ) then
   
         xc = sqrt ( pi2 / x )
         xin = ( eighth / x )**2
         m = 11
   
         if ( three5 <= x ) then
           m = 8
         end if
   
         if ( one30 <= x ) then
           m = 4
         end if
   
         xm = four * real ( m, kind = 8 )
   !
   !  Argument reduction for SIN and COS routines.
   !
         t = aint ( x / ( twopi1 + twopi2 ) + half )
         z = ( ( x - t * twopi1 ) - t * twopi2 ) &
           - ( alpha + half ) / pi2
         vsin = sin ( z )
         vcos = cos ( z )
         gnu = alpha + alpha
   
         do i = 1, 2
   
           s = ( ( xm - one ) - gnu ) * ( ( xm - one ) + gnu ) &
             * xin * half
           t = ( gnu - ( xm - three ) ) * ( gnu + ( xm - three ) )
           capp = s * t / fact(2*m+1)
           t1 = ( gnu - ( xm + one ) ) * ( gnu + ( xm + one ) )
           capq = s * t1 / fact(2*m+2)
           xk = xm
           k = m + m
           t1 = t
   
           do j = 2, m
             xk = xk - four
             s = ( ( xk - one ) - gnu ) * ( ( xk - one ) + gnu )
             t = ( gnu - ( xk - three ) ) * ( gnu + ( xk - three ) )
             capp = ( capp + one / fact(k-1) ) * s * t * xin
             capq = ( capq + one / fact(k) ) * s * t1 * xin
             k = k - 2
             t1 = t
           end do
   
           capp = capp + one
           capq = ( capq + one ) * ( gnu * gnu - one ) * ( eighth / x )
           b(i) = xc * ( capp * vcos - capq * vsin )
   
           if ( nb == 1 ) then
             return
           end if
   
           t = vsin
           vsin = -vcos
           vcos = t
           gnu = gnu + two
   
         end do
   !
   !  If 2 < NB, compute J(X,ORDER+I)  I = 2, NB-1.
   !
         if ( 2 < nb ) then
           gnu = alpha + alpha + two
           do j = 3, nb
             b(j) = gnu * b(j-1) / x - b(j-2)
             gnu = gnu + two
           end do
         end if
   !
   !  Use recurrence to generate results.  First initialize the
   !  calculation of P's.
   !
       else
   
         nbmx = nb - magx
         n = magx + 1
         en = real ( n + n, kind = 8 ) + ( alpha + alpha )
         plast = one
         p = en / x
   !
   !  Calculate general significance test.
   !
         test = ensig + ensig
   !
   !  Calculate P's until N = NB-1.  Check for possible overflow.
   !
         if ( 3 <= nbmx ) then
   
           tover = enten / ensig
           nstart = magx + 2
           nend = nb - 1
           en = real ( nstart + nstart, kind = 8 ) - two + ( alpha + alpha )
   
           do k = nstart, nend
   
             n = k
             en = en + two
             pold = plast
             plast = p
             p = en * plast / x - pold
   !
   !  To avoid overflow, divide P's by TOVER.  Calculate P's until
   !  1 < ABS(P).
   !
             if ( tover < p ) then
   
               tover = enten
               p = p / tover
               plast = plast / tover
               psave = p
               psavel = plast
               nstart = n + 1
   
     100           continue
   
               n = n + 1
               en = en + two
               pold = plast
               plast = p
               p = en * plast / x - pold
   
               if ( p <= one ) then
                 go to 100
               end if
   
               tempb = en / x
   !
   !  Calculate backward test and find NCALC, the highest N such that
   !  the test is passed.
   !
               test = pold * plast &
                 * ( half - half / ( tempb * tempb ) )
               test = test / ensig
               p = plast * tover
               n = n - 1
               en = en - two
               nend = min ( nb, n )
   
               do l = nstart, nend
                 pold = psavel
                 psavel = psave
                 psave = en * psavel / x - pold
                 if ( test < psave * psavel ) then
                   ncalc = l - 1
                   go to 190
                 end if
               end do
   
               ncalc = nend
               go to 190
   
             end if
   
           end do
   
           n = nend
           en = real ( n + n, kind = 8 ) + ( alpha + alpha )
   !
   !  Calculate special significance test for 2 < NBMX.
   !
           test = max ( test, sqrt ( plast * ensig ) * sqrt ( p + p ) )
   
         end if
   !
   !  Calculate P's until significance test passes.
   !
     140     continue
   
         n = n + 1
         en = en + two
         pold = plast
         plast = p
         p = en * plast / x - pold
         if ( p < test ) then
           go to 140
         end if
   !
   !  Initialize the backward recursion and the normalization sum.
   !
     190     continue
   
         n = n + 1
         en = en + two
         tempb = zero
         tempa = one / p
         m = 2 * n - 4 * ( n / 2 )
         sum = zero
         em = real ( n / 2, kind = 8 )
         alpem = ( em - one ) + alpha
         alp2em = ( em + em ) + alpha
   
         if ( m /= 0 ) then
           sum = tempa * alpem * alp2em / em
         end if
   
         nend = n - nb
   !
   !  Recur backward via difference equation, calculating (but not
   !  storing) B(N), until N = NB.
   !
         if ( 0 < nend ) then
   
           do l = 1, nend
   
             n = n - 1
             en = en - two
             tempc = tempb
             tempb = tempa
             tempa = ( en * tempb ) / x - tempc
             m = 2 - m
   
             if ( m /= 0 ) then
               em = em - one
               alp2em = ( em + em ) + alpha
               if ( n == 1 ) then
                 go to 210
               end if
               alpem = ( em - one ) + alpha
               if ( alpem == zero ) then
                 alpem = one
               end if
               sum = ( sum + tempa * alp2em ) * alpem / em
             end if
   
           end do
   
         end if
   !
   !  Store B(NB).
   !
     210     continue
   
         b(n) = tempa
   
         if ( 0 <= nend ) then
   
           if ( nb <= 1 ) then
   
             alp2em = alpha
             if ( alpha + one == one ) then
               alp2em = one
             end if
             sum = sum + b(1) * alp2em
             go to 250
   
           else
   !
   !  Calculate and store B(NB-1).
   !
             n = n - 1
             en = en - two
             b(n) = ( en * tempa ) / x - tempb
   
             if ( n == 1 ) then
               go to 240
             end if
   
             m = 2 - m
   
             if ( m /= 0 ) then
               em = em - one
               alp2em = ( em + em ) + alpha
               alpem = ( em - one ) + alpha
               if ( alpem == zero ) then
                 alpem = one
               end if
               sum = ( sum + b(n) * alp2em ) * alpem / em
             end if
   
           end if
   
         end if
   
         nend = n - 2
   !
   !  Calculate via difference equation and store B(N), until N = 2.
   !
         if ( nend /= 0 ) then
   
           do l = 1, nend
             n = n - 1
             en = en - two
             b(n) = ( en * b(n+1) ) / x - b(n+2)
             m = 2 - m
             if ( m /= 0 ) then
               em = em - one
               alp2em = ( em + em ) + alpha
               alpem = ( em - one ) + alpha
               if ( alpem == zero ) then
                 alpem = one
               end if
               sum = ( sum + b(n) * alp2em ) * alpem / em
             end if
           end do
   
         end if
   !
   !  Calculate B(1).
   !
         b(1) = two * ( alpha + one ) * b(2) / x - b(3)
   
     240     continue
   
         em = em - one
         alp2em = ( em + em ) + alpha
   
         if ( alp2em == zero ) then
           alp2em = one
         end if
   
         sum = sum + b(1) * alp2em
   !
   !  Normalize.  Divide all B(N) by sum.
   !
     250     continue
   
         if ( alpha + one /= one ) then
           sum = sum * gamma ( alpha ) * ( x * half )**( -alpha )
         end if
   
         tempa = enmten
         if ( one < sum ) then
           tempa = tempa * sum
         end if
   
         do n = 1, nb
           if ( abs ( b(n) ) < tempa ) then
             b(n) = zero
           end if
           b(n) = b(n) / sum
         end do
   
       end if
   !
   !  Error return: X, NB, or ALPHA is out of range.
   !
     else
       b(1) = zero
       ncalc = min ( nb, 0_4 ) - 1
     end if
   
     return
   end subroutine

end module free_scattering_mod
