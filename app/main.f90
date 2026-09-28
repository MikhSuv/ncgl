program main
  use precision_mod
  use function_interfaces
  use test_functions
  use nclg
  implicit none

  ! Starting point.
  real(dp) :: x0(2) = [3.0_dp, 4.0_dp]
  ! Required accuracy: the numerical gradient is only good to about 1e-8.
  real(dp), parameter :: eps = 1.0e-9_dp
  ! Restart period and iteration limit of the method.
  integer, parameter :: m = 20
  integer, parameter :: max_iter = 1000
  real(dp) :: x_min(2)

  call nclg_use(parabola_c37, x0, eps, m, max_iter)

  contains 
    subroutine nclg_use(f, x0, eps, m, max_iter)
      procedure(multivariable_func) :: f
      real(dp), intent(in) :: x0(:)
      real(dp), intent(in) :: eps
      integer, intent(in), optional :: m
      integer, intent(in), optional :: max_iter
      real(dp) :: x_min(size(x0))

      integer :: restart = 10
      integer :: mi = 1000

      if (.not. present(m)) restart = max(1, m)
      if (.not. present(max_iter)) mi = max_iter

      x_min = cg_min(parabola_c37, x0, eps, restart,mi)
      print *, "x_min =", x_min
      print *, "f(x_min) =", parabola_c37(x_min)
      print *, "norm of grad =", norm2(grad(x_min, parabola_c37))

    end subroutine nclg_use
end program main
