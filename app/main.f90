program main
  use precision_mod
  use function_interfaces
  use test_functions
  use nclg
  implicit none

  ! Starting point.
  real(dp) :: x0(2) = [3.0_dp, 4.0_dp]
  ! Required accuracy: the numerical gradient is only good to about 1e-8.
  real(dp), parameter :: eps = 1.0e-6_dp
  ! Restart period and iteration limit of the method.
  integer, parameter :: m = 5
  integer, parameter :: max_iter = 100
  real(dp) :: x_min(2)

  x_min = cg_min(quadratic, x0, eps, m, max_iter)
  print *, "x_min =", x_min
  print *, "f(x_min) =", quadratic(x_min)
  print *, "norm of grad =", norm2(grad(x_min, quadratic))

end program main
