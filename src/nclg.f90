module nclg
  use precision_mod
  use function_interfaces
  implicit none
  private

  public :: grad
contains
   ! Compute the gradient of function func at point x0.
   ! Uses forward finite differences with step size delta_x = sqrt(machine epsilon).
   ! Arguments:
   !   x0(:) — point in R^n at which to evaluate the gradient
   !   F    — multivariable function conforming to the multivariable_func interface
   ! Returns:
   !   g(n) — numerical approximation of the gradient
  function grad(x0, func) result(g)
    real(dp), intent(in) :: x0(:)
    procedure(multivariable_func) :: func
    real(dp) :: g(size(x0))
    real(dp) :: x(size(x0))
    real(dp) :: f0
    real(dp), parameter :: delta_x = sqrt(epsilon(1.0_dp))
    integer i, n
    
    n = size(x0)
    f0 = func(x0)
    do i = 1, n
     x = x0 
     x(i) = x(i) + delta_x
     g(i) = func(x) - f0
    end do 

    g = g / delta_x
    
  end function grad
end module nclg
