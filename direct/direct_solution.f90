!======================================================================!
!     Title  : direct_solution.f90                                     !
!     Author : Yusa Shusaku                                            !
!     Date   : 2009-5-3-Sun                                            !
!     Last modification : 2009-5-4-Mon                                 !
!                                                                      !
!     ** Structure of this program **                                  !
!                                                                      !
!     module      global_const                                         !
!     module      pot_prof                                             !
!     program     main                                                 !
!     subroutine  read_int_prof                                        !
!     subroutine  bind_energy                                          !
!     subroutine  solve_eq                                             !
!     subroutine  wave_fct                                             !
!     subroutine  Init                                                 !
!     subroutine  Runge_Kutta                                          !
!     subroutine  fct                                                  !
!     function    Vint                                                 !
!     subroutine  Normal_const                                         !
!     subroutine  Normalized_wf                                        !
!======================================================================!
      module global_const
      implicit none
      real(kind=8), parameter :: hbar = 197.327d0
      real(kind=8), parameter, private :: mass_p = 938.272d0
      real(kind=8), parameter, private :: mass_n = 939.565d0
      real(kind=8), parameter :: rmass = mass_p*mass_n/(mass_p+mass_n)
      real(kind=8), parameter :: rmax = 70.0d0
      real(kind=8), parameter :: dr = 0.05d0
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
      program main
      use pot_prof, only : prof
      implicit none
      integer :: int_type
      real(kind=8), parameter :: epsr = 1.0d-10
      real(kind=8) :: E, rms, C
      type(prof) :: p
      character(len=20), parameter :: FM = '(1x,a,f8.4,a)'

      write(6,*)
      write(6,*) '1 : MT-V'
      write(6,*) '2 : Volkov'
      write(6,*) '3 : ATS3'
      write(6,*) '4 : Minnesota'
      write(6,*)
      write(6,'(a)',advance='no') ' Choose the interaction : '
      read(5,*) int_type

      call read_int_prof(int_type, p)
      call bind_energy(-10.0d0, epsr, E, p)
!     call wave_fct(E, 'wf.dat', rms, p)
      call Normal_const(E, C, p)
      call Normalized_wf(E, C, 'wave_fct', rms, p)

      write(6,*)
      write(6,'(a)',advance='no') ' Interaction : '
      select case(int_type)
        case(1) ; write(6,*) 'MT-V'
        case(2) ; write(6,*) 'Volkov'
        case(3) ; write(6,*) 'ATS3'
        case(4) ; write(6,*) 'Minnesota'
      end select
      write(6,*)
      write(6,FM) 'E   =', E, ' MeV'
      write(6,FM) 'rms =', rms, ' fm'
      write(6,*)

      stop
      end program
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
      subroutine bind_energy(Est, epsr, E, p)
      use pot_prof, only : prof
      implicit none
      integer :: Nnode, Nstate, Nstate_max
      real(kind=8), intent(in) :: Est, epsr
      real(kind=8), intent(out) :: E
      real(kind=8) :: E0, Emin, Emax
      type(prof), intent(in) :: p

      Emin = Est
      Nstate = 0
!     call solve_eq(0.0d0, Nstate_max, p)
!     do 
        E = 0.5d0 * Emin
        E0 = 0.0d0
        do 
          call solve_eq(E, Nnode, p)
          if (Nnode >= Nstate + 1) then
            Emax = E
          else if (Nnode == Nstate) then
            Emin = E
          end if
          E = 0.5d0 * (Emin + Emax)
          if (abs(E0 - E) < abs(E) * epsr) exit
          E0 = E
        end do
!       Nstate = Nstate + 1
!       if (Nstate == Nstate_max) exit
!       Emin = E
!     end do

      return
      end subroutine
!======================================================================!
      subroutine solve_eq(E, Nnode, p)
      use global_const, only : rmax, dr
      use pot_prof, only : prof
      implicit none
      logical :: flag_p, flag_n
      integer, intent(out) :: Nnode
      real(kind=8), intent(in) :: E
      real(kind=8) :: r, u(2)
      type(prof), intent(in) :: p

      call Init(r, u, p)
      Nnode = 0
      flag_p = .true.
      flag_n = .false.
      do while (r <= rmax)
        call Runge_Kutta(E, r, u, dr, p)
        if (flag_n .and. u(1) > 0.0d0) then
          flag_p = .true.
          flag_n = .false.
          Nnode = Nnode + 1
        else if (flag_p .and. u(1) < 0.0d0) then
          flag_p = .false.
          flag_n = .true.
          Nnode = Nnode + 1
        end if
        r = r + dr
      end do

      return
      end subroutine
!======================================================================!
      subroutine wave_fct(E, Fname, rms, p)
      use global_const, only : rmax, dr
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(in) :: E
      real(kind=8), intent(out) :: rms
      real(kind=8) :: r, u(2), S
      character(len=*), intent(in) :: Fname
      type(prof), intent(in) :: p

      call Init(r, u, p)
      S = 0.5d0 * u(1) * u(1)
      rms = 0.0d0
      open(7,file=Fname)
      do while(r <= rmax)
        write(7,*) r, u(1)
        call runge_kutta(E, r, u, dr, p)
        r = r + dr
        S = S + u(1) * u(1)
        rms = rms + (u(1) * r) ** 2
      end do
      close(7)
      call runge_kutta(E, r, u, dr, p)
      r = r + dr
      S = S + 0.5d0 * u(1) * u(1)
      rms = rms + 0.5d0 * (u(1) * r) ** 2 
      rms = 0.5d0 * sqrt(rms / S)

      return
      end subroutine
!======================================================================!
      pure subroutine Init(x, y, p)
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(out) :: x
      real(kind=8), intent(out) :: y(2)
      type(prof), intent(in) :: p

      select case(p%int_form)
        case(0)
          x = 1.0d-300
          y(1) = x
          y(2) = 1.0d0
        case(1)   
          x = 0.0d0
          y(1) = 0.0d0
          y(2) = 1.0d0
      end select

      return
      end subroutine
!======================================================================!
      subroutine Runge_Kutta(E, x, y, h, p)
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(in) :: x, h, E
      real(kind=8), intent(inout) :: y(2)
      real(kind=8), dimension(2) :: k1, k2, k3, k4
      real(kind=8) :: h2
      type(prof), intent(in) :: p

      h2 = 0.5d0 * h
      call fct(E, x, y, k1, p)
      call fct(E, x+h2, y+h2*k1, k2, p)
      call fct(E, x+h2, y+h2*k2, k3, p)
      call fct(E, x+h,  y+h*k3,  k4, p)
      y = y + h * (k1 + 2.0d0 * (k2 + k3) + k4) / 6.0d0

      contains
!**********************************************************************!
      subroutine fct(E, x, y, f, p)
      use global_const, only : rmass, hbar
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(in) :: x, E, y(2)
      real(kind=8), intent(out) :: f(2)
      real(kind=8), external :: Vint
      type(prof), intent(in) :: p

      f(1) = y(2)
      f(2) = - 2.0d0 * rmass / (hbar * hbar) * (E - Vint(x,p)) * y(1)

      return
      end subroutine
!**********************************************************************!
      end subroutine
!======================================================================!
      pure function Vint(r, p) result(V)
      use pot_prof, only : prof
      implicit none
      integer :: i
      real(kind=8), intent(in) :: r
      real(kind=8) :: V
      type(prof), intent(in) :: p

      select case(p%int_form)
        case(0)     ! Yukawa
          V = p%V(1) * exp(- p%mu(1) * r) / r
          do i=2, p%Nterm
            V = V + p%V(i) * exp(- p%mu(i) * r) / r
          end do
        case(1)     ! Gauss
          V = p%V(1) * exp(- p%mu(1) * r * r)
          do i=2, p%Nterm
            V = V + p%V(i) * exp(- p%mu(i) * r * r)
          end do
      end select

      return
      end function
!======================================================================!
      subroutine Normal_const(E, C, p)
      use global_const, only : rmax, dr
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(in) :: E
      real(kind=8), intent(out) :: C
      real(kind=8) :: r, u(2)
      type(prof), intent(in) :: p

      call Init(r, u, p)
      C = 0.5d0 * u(1) * u(1)
      do while (r <= rmax)
        call runge_kutta(E, r, u, dr, p)
        r = r + dr
        C = C + u(1) * u(1)
      end do
      call runge_kutta(E, r, u, dr, p)
      r = r + dr
      C = (C + 0.5d0 * u(1) * u(1)) * dr

      return
      end subroutine
!======================================================================!
      subroutine Normalized_wf(E, C, Fname, rms, p)
      use global_const, only : rmax, dr
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(in) :: E, C
      real(kind=8), intent(out) :: rms
      real(kind=8) :: r, u(2)
      character(len=*), intent(in) :: Fname
      type(prof), intent(in) :: p

      call Init(r, u, p)
      u = u / sqrt(C)    ! Normalization
      rms = 0.0d0
      open(7,file=Fname)
      do while(r <= rmax)
        write(7,*) r, u(1)
        call runge_kutta(E, r, u, dr, p)
        r = r + dr
        rms = rms + (u(1) * r) ** 2
      end do
      close(7)
      call runge_kutta(E, r, u, dr, p)
      r = r + dr
      rms = (rms + 0.5d0 * (u(1) * r) ** 2) * dr
      rms = 0.5d0 * sqrt(rms)

      return
      end subroutine
!======================================================================!
      subroutine Init_Numerov(E, x, y0, y, h, p)
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(in) :: E, h
      real(kind=8), intent(out) :: x, y0, y
      real(kind=8) :: yt(2)
      type(prof), intent(in) :: p

      x = 0.0d0
      yt(1) = 0.0d0
      yt(2) = 1.0d0
      call runge_kutta(E, x, yt, h, p)
      y0 = 0.0d0
      y = yt(1)

      return
      end subroutine
!======================================================================!
      subroutine Numerov(E, x, y0, y, h, p)
      use global_const, only : hbar, rmass
      use pot_prof, only : prof
      implicit none
      real(kind=8), intent(in) :: E, x, h
      real(kind=8), intent(inout) :: y0, y
      real(kind=8), external :: Vint
      real(kind=8) :: w1, w2, c, y1
      type(prof), intent(in) :: p

      c = 2.0d0 * rmass / (hbar * hbar) * h * h
      w1 = 1.0d0 + c * (E - Vint(x,p)) / 12.0d0
      w2 = 2.0d0 * (1.0d0 - 5.0d0/12.0d0 * c * (E - Vint(x+h,p)))
      y1 = y
      y = (w2 * y - w1 * y0) / (1.0d0 + c * (E - Vint(x+h+h,p))/12.0d0)
      y0 = y1

      return
      end subroutine
!======================================================================!
      subroutine solve_eq_Numerov(E, Nnode, p)
      use global_const, only : rmax, dr
      use pot_prof, only : prof
      implicit none
      logical :: flag_p, flag_n
      integer, intent(out) :: Nnode
      real(kind=8), intent(in) :: E
      real(kind=8) :: r, u(2), y0, y
      type(prof), intent(in) :: p

      call Init_Numerov(E, r, y0, y, dr, p)
      Nnode = 0
      flag_p = .true.
      flag_n = .false.
      do while (r <= rmax)
        call Numerov(E, r, y0, y, dr, p)
        if (y0 > 0.0d0 .and. flag_n) then
          flag_p = .true.
          flag_n = .false.
          Nnode = Nnode + 1
        else if (y0 < 0.0d0 .and. flag_p) then
          flag_p = .false.
          flag_n = .true.
          Nnode = Nnode + 1
        end if
        r = r + dr
      end do

      return
      end subroutine
