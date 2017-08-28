!======================================================================!
!     Title    : ho_basis.f90                                          !
!     Coded by : Yusa Shusaku                                          !
!     Date     : 2009-4-30-Fri ~ 2009-5-6-Wed                          !
!     Last modification : 2009-7-26-Sun                                !
!                                                                      !
!     using libraries : lapack                                         !
!                                                                      !
!     This program computes binding energy and rms radius of deuteron  !
!     by expanding the wave function with eigen functions of harmonic  !
!     oscillator. The Schroedinger equation reduces to the seqular     !
!     equation of the Hamiltonian matrix. The eigenvalues of the       !
!     Hamiltonian matrix give the energy and eigenvectors give the     !
!     expansion coefficients.                                          !
!                                                                      !
!     ** Structure of this program **                                  !
!                                                                      !
!     module      global_const                                         !
!     module      pot_prof                                             !
!     module      basis                                                !
!     program     main                                                 !
!     subroutine  read_input                                           !
!     subroutine  read_int_prof                                        !
!     subroutine  energy_surfice                                       !
!     subroutine  base_expand                                          !
!     subroutine  optimize_omega                                       !
!       -> contains  subroutine  xp_min                                !
!     subroutine  Hamiltonian                                          !
!     subroutine  mat_elem                                             !
!     subroutine  Integral_V                                           !
!     function    V                                                    !
!     subroutine  eps                                                  !
!     subroutine  wave_fct                                             !
!     subroutine  rms_only                                             !
!     function    PSI                                                  !
!     function    phi_ho                                               !
!     function    Vint                                                 !
!     function    NH                                                   !
!======================================================================!
      module global_const
!----------------------------------------------------------------------!
!     Definition of global constants.                                  !
!     thread    :  Number of threads in the parallization.             !
!     Nbase_max :  Maximum number of basis.                            !
!     hw_min    :  Minimum value of frequency of the h.o..             !
!     hw_max    :  Maximum value of frequency of the h.o..             !
!     b         :  Upper limit of the integration in the calculation   !
!                  of the matrix element of the Hamiltonian.           !
!----------------------------------------------------------------------!
      implicit none
!$    integer :: thread
      integer, parameter :: Nbase_max = 250
      real(kind=8), parameter :: hbar = 197.327d0
      real(kind=8), parameter, private :: mass_p = 938.272d0
      real(kind=8), parameter, private :: mass_n = 939.565d0
      real(kind=8), parameter :: rmass = mass_p*mass_n/(mass_p+mass_n)
      real(kind=8), parameter :: hw_min = 10.0d0
      real(kind=8), parameter :: hw_max = 40.0d0
      real(kind=8), parameter :: b = 27.0d0
      end module
!======================================================================!
      module pot_prof
      implicit none
      type prof
        integer :: int_form
        integer :: Nterm
        real(kind=8), allocatable :: V(:)
        real(kind=8), allocatable :: mu(:)
      end type
      end module
!======================================================================!
      module basis
      interface
        subroutine base_expand(p, eps, Nbase, step, Nopt, hw, Egs, Cn)
        use pot_prof, only : prof
        integer, intent(in) :: Nopt, step
        integer, intent(inout) :: Nbase
        real(kind=8), intent(in) :: eps
        real(kind=8), intent(inout) :: hw
        real(kind=8), intent(out) :: Egs
        real(kind=8), intent(out), allocatable :: Cn(:)
        type(prof), intent(in) :: p
        end subroutine
      end interface
      end module
!======================================================================!
      program main
      use global_const, only : rmass, hbar
!$    use global_const, only : thread
      use pot_prof, only : prof
      use basis, only : base_expand
      implicit none
      integer :: Nbase, Nopt, step, int_type, ios
      real(kind=8), allocatable, dimension(:) :: Cn
      real(kind=8) :: rms, hw, ak, Egs, eps
      character(len=20), parameter :: FM = '(1x,a,f8.4,a)'
      character(len=20), parameter :: FM2 = '(1x,a,i4)'
      type(prof) :: p

      write(6,*)
      write(6,*) '1 : MT-V'
      write(6,*) '2 : Volkov'
      write(6,*) '3 : ATS3'
      write(6,*) '4 : Minnesota'
      write(6,*)
      write(6,'(a)',advance='no') ' Choose the interaction : '
      read(5,*) int_type
      write(6,*)

      call read_input(Nbase, step, Nopt, hw, eps, ios)
      call read_int_prof(int_type, p)
!     call energy_surfice(10, step, 'energy_surfice', p)
      call base_expand(p, eps, Nbase, step, Nopt, hw, Egs, Cn)
      call wf_rms(Nbase, hw, Cn, 'wave-fct.dat', 70.0d0, rms)
      ak = sqrt(rmass * hw / (hbar * hbar))

      write(6,*)
      write(6,'(a)',advance='no') ' Interaction : '
      select case(int_type)
        case(1) ; write(6,*) 'MT-V'
        case(2) ; write(6,*) 'Volkov'
        case(3) ; write(6,*) 'ATS3'
        case(4) ; write(6,*) 'Minnesota'
      end select
      write(6,*)
      write(6,FM) 'Omega =', hw, ' MeV'
      write(6,FM) 'l =',1.0d0/ak,' fm (characteristic length of h.o.)'
      write(6,FM2) 'Number of basis =', Nbase
      write(6,*)
      write(6,FM) 'E   =', Egs, ' MeV'
      write(6,FM) 'rms =', rms, ' fm'
      write(6,*)

      stop
      end program
!======================================================================!
      subroutine read_input(Nbase, step, Nopt, hw, eps, ios)
!$    use global_const, only : thread
      implicit none
      integer, intent(out) :: Nbase, step, Nopt, ios
      real(kind=8), intent(out) :: hw, eps

      open(7, file='input', status='old', action='read', iostat=ios)
      read(7,*) Nbase, step, Nopt
      read(7,*) hw, eps
!$    read(7,*) thread
      close(7)

      return
      end subroutine
!======================================================================!
      subroutine read_int_prof(int_type, p)
      use pot_prof, only : prof
      implicit none
      integer, intent(in) :: int_type
      type(prof), intent(out) :: p

      select case(int_type)
      case(1)   ! MT-V
        p%int_form = 0
        p%Nterm = 2
        allocate(p%V(p%Nterm), p%mu(p%Nterm))
        p%V(1) =  1458.05d0 ; p%mu(1) = 3.11d0
        p%V(2) = - 578.09d0 ; p%mu(2) = 1.55d0
      case(2)   ! Volkov
        p%int_form = 1
        p%Nterm = 2
        allocate(p%V(p%Nterm), p%mu(p%Nterm))
        p%V(1) =  144.86d0 ; p%mu(1) = 0.82d0 ** (-2)
        p%V(2) = - 83.34d0 ; p%mu(2) = 1.60d0 ** (-2)
      case(3)   ! ATS3
        p%int_form = 1
        p%Nterm = 3
        allocate(p%V(p%Nterm), p%mu(p%Nterm))
        p%V(1) =  1000.0d0 ; p%mu(1) = 3.0d0  
        p%V(2) = - 326.7d0 ; p%mu(2) = 1.05d0
        p%V(3) = -  43.0d0 ; p%mu(3) = 0.60d0
      case(4)   ! Minnesota
        p%int_form = 1
        p%Nterm = 2
        allocate(p%V(p%Nterm), p%mu(p%Nterm))
        p%V(1) =   200.0d0 ; p%mu(1) = 1.487d0
        p%V(2) = - 178.0d0 ; p%mu(2) = 0.639d0
      case default
        stop 'specify the interaction correctly.'
      end select

      return
      end subroutine
!======================================================================!
      subroutine energy_surfice(Nbase, step, Fname, p)
      use global_const, only : rmass, hbar
      use pot_prof, only : prof
      use f95_lapack
      implicit none
      integer, intent(in) :: Nbase, step
      real(kind=8), allocatable, dimension(:,:) :: H, Cn
      real(kind=8), parameter :: epsr = 1.0d-15
      real(kind=8) :: c, hw, ak, E(Nbase)
      character(len=*), intent(in) :: Fname
      type(prof), intent(in) :: p

      allocate(H(Nbase,Nbase), Cn(Nbase,Nbase))
      c = rmass / (hbar * hbar)
      hw = 0.01d0
      ak = sqrt(c * hw)
      open(7,file=Fname)
      do while(hw <= 50.0d0)
        call Hamiltonian(Nbase, step, hw, ak, H, p)
        call dsyev_f95(H, E, 'N', 'L')
        write(7,*) hw, E(1)
        hw = hw + 0.1d0
        ak = sqrt(c * hw)
      end do
      close(7)
      deallocate(H, Cn)

      return
      end subroutine
!======================================================================!
      subroutine base_expand(p, eps, Nbase, step, Nopt, hw, Egs, Cn)
      use global_const, only : Nbase_max, rmass, hbar
      use pot_prof, only : prof
      use f95_lapack
      implicit none
      logical :: not_optimized
      integer, intent(in) :: Nopt, step
      integer, intent(inout) :: Nbase
      real(kind=8), intent(in) :: eps
      real(kind=8), intent(inout) :: hw
      real(kind=8), intent(out) :: Egs
      real(kind=8), intent(out), allocatable :: Cn(:)
      real(kind=8), allocatable, dimension(:,:) :: H
      real(kind=8), allocatable :: E(:)
      real(kind=8), parameter :: epsr = 1.0d-15
      real(kind=8) :: E0, dE, ak
      character(len=30), parameter :: FM1 = '(1x,a,i4,a,f10.6)'
      character(len=30), parameter :: FM2 = '(1x,a,f9.4,a)'
      type(prof), intent(in) :: p

      E0 = 0.0d0
      not_optimized = .true.
      ak = sqrt(rmass * hw / (hbar * hbar))

      do Nbase=Nbase, Nbase_max, step
        allocate(H(Nbase,Nbase), E(Nbase))

        if (not_optimized .and. Nbase >= Nopt) then
          write(6,'(a)',advance='no') '   Optimizing omega ...'
          call optimize_omega(Nopt, step, hw, ak, eps*10.0d0, p)
          write(6,FM2) '=> omega =', hw, ' MeV'
          not_optimized = .false.
        end if
        
        call Hamiltonian(Nbase, step, hw, ak, H, p)
        call dsyev_f95(H, E, 'N', 'L')
        write(6,FM1) 'Nbase =', Nbase, ',  E =', E(1)
        
        dE = abs(E(1) - E0) / abs(E(1))
        if (dE < eps .or. Nbase + step > Nbase_max) exit
        E0 = E(1)
        deallocate(H, E)
      end do
      call Hamiltonian(Nbase, step, hw, ak, H, p)
      call dsyev_f95(H, E, 'V', 'L')
      allocate(Cn(Nbase))
      Egs = E(1)
      Cn = H(:,1)
      deallocate(H, E)

      return
      end subroutine
!======================================================================!
      subroutine optimize_omega(Nbase, step, hw, ak, epsr_E, p)
      use global_const, only : hbar, rmass, hw_min, hw_max
      use pot_prof, only : prof
      use f95_lapack
      implicit none
      integer, intent(in) :: Nbase, step
      real(kind=8), intent(in) :: epsr_E
      real(kind=8), intent(out) :: hw, ak
      real(kind=8), parameter :: epsr = 1.0d-15
      real(kind=8), allocatable, dimension(:,:) :: H, Cn
      real(kind=8), dimension(Nbase) :: E0, E1, E2, E
      real(kind=8) :: c, hw0, hw1, hw2, ak0, ak1, ak2
      type(prof), intent(in) :: p

      allocate(H(Nbase,Nbase), Cn(Nbase,Nbase))
      c = rmass / (hbar * hbar)
      hw0 = hw_min
      hw2 = hw_max
      hw1 = 0.5d0 * (hw0 + hw2)
      ak0 = sqrt(c * hw0)
      ak1 = sqrt(c * hw1)
      ak2 = sqrt(c * hw2)
      call Hamiltonian(Nbase, step, hw0, ak0, H, p)
      call dsyev_f95(H, E0, 'N', 'L')
      call Hamiltonian(Nbase, step, hw1, ak1, H, p)
      call dsyev_f95(H, E1, 'N', 'L')
      call Hamiltonian(Nbase, step, hw2, ak2, H, p)
      call dsyev_f95(H, E2, 'N', 'L')
      if (E1(1) > E0(1) .or. E1(1) > E2(1)) then
        if (E0(1) > E2(1)) then
          hw = hw2 ; ak = ak2
        else
          hw = hw0 ; ak = ak0
        end if
        return
      end if
      do
        call xp_min(hw0, hw1, hw2, E0(1), E1(1), E2(1), hw)
        ak = sqrt(c * hw)
        call Hamiltonian(Nbase, step, hw, ak, H, p)
        call dsyev_f95(H, E, 'N', 'L')
        if (abs(E(1) - E1(1)) < abs(E(1)) * epsr_E) return
        if (hw < hw1) then
          if (E(1) < E1(1)) then
            hw2 = hw1 ; ak2 = sqrt(c * hw2)
            hw1 = hw  ; ak1 = sqrt(c * hw1)
          else
            hw0 = hw  ; ak0 = sqrt(c * hw0)
          end if
        else
          if (E(1) < E1(1)) then
            hw0 = hw1 ; ak0 = sqrt(c * hw0)
            hw1 = hw  ; ak1 = sqrt(c * hw1)
          else
            hw2 = hw  ; ak2 = sqrt(c * hw2)
          end if
        end if
        call Hamiltonian(Nbase, step, hw0, ak0, H, p)
        call dsyev_f95(H, E0, 'N', 'L')
        call Hamiltonian(Nbase, step, hw1, ak1, H, p)
        call dsyev_f95(H, E1, 'N', 'L')
        call Hamiltonian(Nbase, step, hw2, ak2, H, p)
        call dsyev_f95(H, E2, 'N', 'L')
      end do
      deallocate(H, Cn)

      contains
!**********************************************************************!
      pure subroutine xp_min(a, b, c, fa, fb, fc, x)
      implicit none
      real(kind=8), intent(in) :: a, b, c, fa, fb, fc
      real(kind=8), intent(out) :: x

      x = (b - a) ** 2 * (fb - fc) - (b - c) ** 2 * (fb - fa)
      x = b - 0.5d0 * x / ((b - a) * (fb - fc) - (b - c) * (fb - fa)) 

      return
      end subroutine
!**********************************************************************!
      end subroutine
!======================================================================!
      subroutine Hamiltonian(Nbase, step, hw, ak, H, p)
      use pot_prof, only : prof
      implicit none
      integer, intent(in) :: Nbase, step
      integer :: i
      type(prof), intent(in) :: p
      real(kind=8), intent(in) :: hw, ak
      real(kind=8), intent(out) :: H(Nbase,Nbase)
      real(kind=8) :: e(Nbase)

      call mat_elem(Nbase, step, hw, ak, H, p)
      call eps(Nbase, step, hw, e)
      forall (i=1:Nbase) H(i,i) = H(i,i) + e(i)

      return
      end subroutine
!======================================================================!
      subroutine mat_elem(Nbase, step, hw, ak, V, p)
      use pot_prof, only : prof
      use global_const, only : Nbase_max
!$    use global_const, only : thread
      implicit none
      integer, intent(in) :: Nbase, step
      integer :: n, m, n2, m2
      real(kind=8), intent(in) :: hw, ak
      real(kind=8), intent(out) :: V(Nbase,Nbase)
      real(kind=8), save :: Vp(Nbase_max,Nbase_max), hw0 = 0.0d0
      real(kind=8) :: S
      type(prof), intent(in) :: p

!$    call omp_set_num_threads(thread)
      if (hw == hw0) then
        V(1:Nbase-1,1:Nbase-1) = Vp(1:Nbase-1,1:Nbase-1)
!$OMP parallel shared(hw,ak,p,Nbase,V) private(m,n,m2,n2,S)
!$OMP do
        do m=0, Nbase-1
          m2 = 2 * m + 1
          do n=Nbase-step, Nbase-1
            n2 = 2 * n + 1
            call Integral_V(n2, m2, hw, ak, S, p)
            V(n+1,m+1) = (- 1.0d0) ** (n + m) * S
            Vp(n+1,m+1) = V(n+1,m+1)
          end do
        end do
!$OMP end do
!$OMP end parallel
      else
!$OMP parallel shared(hw,ak,p,Nbase,V,Vp) private(m,n,m2,n2,S)
!$OMP do
        do m=0, Nbase-1
          m2 = 2 * m + 1
          do n=m, Nbase-1
            n2 = 2 * n + 1
            call Integral_V(n2, m2, hw, ak, S, p)
            V(n+1,m+1) = (- 1.0d0) ** (n + m) * S
            Vp(n+1,m+1) = V(n+1,m+1)
          end do
        end do
!$OMP end do
!$OMP end parallel
        hw0 = hw
      end if

      return
      end subroutine
!======================================================================!
      subroutine Integral_V(n, m, hw, ak, S, p)
      use pot_prof, only : prof
      use global_const, only : b
      implicit none
      integer, intent(in) :: n, m
      integer :: j, ngrid
      real(kind=8), intent(in) :: ak, hw
      real(kind=8), intent(out) :: S
      real(kind=8), external :: V, NH
      real(kind=8), parameter :: a = 1.0d-100
      real(kind=8), parameter :: hh = 0.1d0
      real(kind=8) :: x
      type(prof), intent(in) :: p

      ngrid = nint((b - a) / hh)
      S = 0.5d0 * (exp(- a * a) * V(a,hw,ak,p) * NH(n,a) * NH(m,a)   &
     &           + exp(- b * b) * V(b,hw,ak,p) * NH(n,b) * NH(m,b))
      x = a
      do j=1, ngrid-1
        x = x + hh
        S = S + exp(- x * x) * V(x,hw,ak,p) * NH(n,x) * NH(m,x)
      end do
      S = S * hh

      return
      end subroutine
!======================================================================!
      function V(r, hw, ak, p) result(f)
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(in) :: r, hw, ak
      real(kind=8), external :: Vint
      real(kind=8) :: f
      type(prof), intent(in) :: p

      f = 2.0d0 * Vint(r/ak, p) - r * r * hw

      return
      end function
!======================================================================!
      subroutine eps(Nbase, step, hw, e)
      use global_const, only : Nbase_max
      implicit none
      integer, intent(in) :: Nbase, step
      integer :: n
      real(kind=8), intent(in) :: hw
      real(kind=8), intent(out) :: e(Nbase)
      real(kind=8), save :: ep(Nbase_max), hw0 = 0.0d0

      if (hw == hw0) then
        e(1:Nbase-1) = ep(1:Nbase-1)
        do n=Nbase-step, Nbase-1
          e(n+1) = 0.5d0 * hw * (3.0d0 + 4.0d0 * dble(n))
          ep(n+1) = e(n+1)
        end do
      else
        do n=0, Nbase-1
          e(n+1) = 0.5d0 * hw * (3.0d0 + 4.0d0 * dble(n))
          ep(n+1) = e(n+1)
        end do
        hw0 = hw
      end if

      return
      end subroutine
!======================================================================!
      subroutine wf_rms(Nbase, hw, Cn, Fname, rmax, rms)
      use global_const, only : rmass, hbar
      implicit none
      integer, intent(in) :: Nbase
      real(kind=8), intent(in) :: hw, Cn(Nbase), rmax
      real(kind=8), intent(out) :: rms
      real(kind=8), external :: PSI
      real(kind=8), parameter :: dr = 0.05d0
      real(kind=8) :: r, u, ak
      character(len=*), intent(in) :: Fname

      ak = sqrt(rmass / (hbar * hbar) * hw)
      open(7, file=Fname)
      r = 1.0d-200
      rms = 0.0d0
      do while(r <= rmax)
        u = PSI(Nbase, ak, Cn, r) * r
        write(7,*) r, u
        rms = rms + dr * (u * r) ** 2
        r = r + dr
      end do
      close(7)
      rms = (rms + 0.5d0 * dr * (u * r) ** 2)
      rms = 0.5d0 * sqrt(rms)


      return
      end subroutine
!======================================================================!
      function PSI(Nbase, ak, Cn, r) result(f)
      implicit none
      integer, intent(in) :: Nbase
      integer :: i
      real(kind=8), intent(in) :: r, ak, Cn(Nbase)
      real(kind=8), external :: phi_ho
      real(kind=8) :: f

      f = 0.0d0
      do i=0, Nbase-1
        f = f + Cn(i+1) * phi_ho(i, ak, r)
      end do

      return
      end function
!======================================================================!
      function phi_ho(n, ak, r) result(f)
      implicit none
      integer, intent(in) :: n
      real(kind=8), intent(in) :: r, ak
      real(kind=8), external :: NH
      real(kind=8) :: f, x

      x = ak * r
      f = (- 1.0d0) ** n * sqrt(2.0d0 * ak)            &
     &    * exp(- 0.5d0 * x * x) * NH(2*n+1,x) / r 

      return
      end function
!======================================================================!
      pure function Vint(r, p) result(V)
      use pot_prof, only : prof
      implicit none
      integer :: i
      real(kind=8), intent(in) :: r
      real(kind=8) :: V
      type(prof), intent(in) :: p

      select case(p%int_form)
        case(0)       ! Yukawa
          V = p%V(1) * exp(- p%mu(1) * r) / r
          do i=2, p%Nterm
            V = V + p%V(i) * exp(- p%mu(i) * r) / r
          end do
        case(1)       ! Gauss
          V = p%V(1) * exp(- p%mu(1) * r * r)
          do i=2, p%Nterm
            V = V + p%V(i) * exp(- p%mu(i) * r * r)
          end do
      end select

      return
      end function
!======================================================================!
      pure function NH(n, x) result(f)
      implicit none
      integer, intent(in) :: n
      integer :: i
      real(kind=8), intent(in) :: x
      real(kind=8), parameter :: c = 0.7511255444649425d0 ! 1 / PI^(0.25)
      real(kind=8) :: f, f0, f1

      if (n == 0) then
        f = c
        return
      else if (n == 1) then
        f = sqrt(2.0d0) * c * x
        return
      end if

      f0 = c
      f1 = sqrt(2.0d0) * c * x

      do i=1, n-1
        f = (x * f1 - sqrt(dble(i)*0.5d0) * f0) * sqrt(2.0d0/dble(i+1))
        f0 = f1
        f1 = f
      end do

      return
      end function
!======================================================================!
