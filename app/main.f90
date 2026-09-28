! Demonstration program for the nclg module.
! Runs the line search nclg::find_alpha along the steepest descent direction
! for a function that has a minimum along the ray and for one that has none.
program main
  use precision_mod
  use function_interfaces
  use test_functions
  use nclg
  implicit none

  ! Starting point of both demonstrations.
  real(dp) :: p(2)

  ! One step of the line search along the steepest descent direction.
  p = [3.0_dp, 0.0_dp]
  call step(quadratic, p)

  ! The same for a function without a minimum along the ray.
  p = [0.0_dp, 0.0_dp]
  call step(exponential, p)

  contains

  ! Take a single steepest descent step and print the gradient at the starting
  ! point together with the step found by the line search.
  ! Arguments:
  !   func — multivariable function conforming to the multivariable_func interface
  !   x0(:)— point in R^n at which the step starts
  subroutine step(func, x0)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp) :: alpha
    ! Steepest descent direction p = -grad f(x0).
    real(dp) :: d(size(x0))
    d = grad(x0, func)
    alpha = find_alpha(x0, -d, 1.0e-14_dp, func)
    print *, "x =", x0, " grad =", d
    print *, "alpha =", alpha, " f(x + alpha*p) =", func(x0 - alpha * d)
  end subroutine step

end program main
