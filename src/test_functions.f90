! Module providing a collection of test functions F: R^n -> R together with
! their analytic gradients. These are used to validate the numerical gradient
! routine by comparing numeric and analytic results.
module test_functions
  use precision_mod
  
  implicit none

  contains
    ! Constant function f(x) = 1.
    ! Analytic gradient: grad f(x) = 0.
    function const(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = 1.0_dp ! const
    end function const

    function const_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = 0.0_dp ! gradient of const
    end function const_grad

    ! Quadratic function f(x) = sum_i x_i^2.
    ! Analytic gradient: grad f(x) = 2*x.
    function quadratic(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = sum(x*x)
    end function quadratic

    function quadratic_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = 2*x
    end function quadratic_grad

    ! Linear function f(x) = a·x + 17 with constant coefficients a_i = 18.
    ! Analytic gradient: grad f(x) = 18 everywhere.
    function linear(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      real(dp) :: a(size(x))
      a = 18.0_dp
      y = dot_product(a, x) + 17.0_dp
    end function linear

    function linear_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = 18.0_dp
    end function linear_grad

    ! Exponential function f(x) = sum_i exp(x_i).
    ! Analytic gradient: grad f(x) = exp(x).
    function exponential(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = sum(exp(x))
    end function exponential

    function exponential_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = exp(x)
    end function exponential_grad

    ! Logarithm function f(x) = sum_i ln(x_i). Only defined for x_i > 0.
    ! Analytic gradient: grad f(x) = 1/x.
    function logarithm(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = sum(log(x))
    end function logarithm
    
    function logarithm_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = 1.0_dp/x
    end function logarithm_grad

    ! Sinus function f(x) = sum_i sin(x_i).
    ! Analytic gradient: grad f(x) = cos(x).
    function sinus(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = sum(sin(x))
    end function sinus

    function sinus_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = cos(x)
    end function sinus_grad
 
end module test_functions
