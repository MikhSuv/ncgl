! Module providing a collection of test functions F: R^n -> R, some of them
! together with their analytic gradients. The functions are used to check the
! numerical gradient of nclg::grad and the search for a minimum performed by
! nclg::cg_min.
! A function comes with a gradient only when a test compares the two; the
! analytic gradient is given in the comment of the function itself.
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

    ! Analytic gradient of const.
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

    ! Analytic gradient of quadratic.
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

    ! Analytic gradient of linear.
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

    ! Analytic gradient of exponential.
    function exponential_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = exp(x)
    end function exponential_grad

    ! Sinus function f(x) = sum_i sin(x_i).
    ! Analytic gradient: grad f(x) = cos(x).
    function sinus(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = sum(sin(x))
    end function sinus

    ! Analytic gradient of sinus.
    function sinus_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = cos(x)
    end function sinus_grad
 
    ! 1D shifted parabola f(x) = (x(1) - c)^2 / 2 with the minimum at x(1) = c.
    ! Used to test the search for a minimum in one dimension, where the method
    ! follows the only direction available.
    ! c = 0.3: minimum near the starting point x(1) = 0.
    function parabola_c03(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = 0.5_dp * (x(1) - 0.3_dp)**2
    end function parabola_c03

    ! c = 3.7: minimum far from the starting point x(1) = 0, so the line search
    ! has to grow the step several times to reach it.
    function parabola_c37(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = 0.5_dp * (x(1) - 3.7_dp)**2 + 20.0_dp 
    end function parabola_c37

    ! Stretched quadratic f(x) = x1^2 + 10*x2^2 with the minimum at the origin.
    ! Its condition number is 10, so the iterates of the conjugate method and of
    ! steep descent differ by orders of magnitude within ten iterations.
    ! Analytic gradient: grad f(x) = (2*x1, 20*x2).
    function quadratic_k10(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = x(1)**2 + 10.0_dp * x(2)**2
    end function quadratic_k10

    ! Analytic gradient of quadratic_k10. Components beyond the second are zero
    ! because f does not depend on them.
    function quadratic_k10_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = 0.0_dp
      y(1) = 2.0_dp * x(1)
      y(2) = 20.0_dp * x(2)
    end function quadratic_k10_grad

    ! Stretched quadratic f(x) = x1^2 + 10*x2^2 + 100*x3^2, condition number 100.
    ! Analytic gradient: grad f(x) = (2*x1, 20*x2, 200*x3).
    function quadratic_k100(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = x(1)**2 + 10.0_dp * x(2)**2 + 100.0_dp * x(3)**2
    end function quadratic_k100

    ! Analytic gradient of quadratic_k100.
    function quadratic_k100_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = 0.0_dp
      y(1) = 2.0_dp * x(1)
      y(2) = 20.0_dp * x(2)
      y(3) = 200.0_dp * x(3)
    end function quadratic_k100_grad

    ! Bounded below but nonconvex f(x) = x1^2 + 10*x2^2 + 50*sin(10*x1) + 20*sin(5*x3).
    ! The oscillations make the conjugate direction leave the cone of descent, so
    ! the restart on a bad direction is taken repeatedly here.
    ! Analytic gradient: grad f(x) = (2*x1 + 500*cos(10*x1), 20*x2, 100*cos(5*x3)).
    function wavy(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      y = x(1)**2 + 10.0_dp * x(2)**2 + 50.0_dp * sin(10.0_dp * x(1)) &
          + 20.0_dp * sin(5.0_dp * x(3))
    end function wavy

    ! Analytic gradient of wavy.
    function wavy_grad(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y(size(x))
      y = 0.0_dp
      y(1) = 2.0_dp * x(1) + 500.0_dp * cos(10.0_dp * x(1))
      y(2) = 20.0_dp * x(2)
      y(3) = 100.0_dp * cos(5.0_dp * x(3))
    end function wavy_grad

    function h(x) result(y)
      real(dp), intent(in) :: x(:)
      real(dp) :: y
      
      y = (1.0_dp - x(1))**2 + 100.0_dp * (x(2) - x(1)**2)**2 
    end function h

end module test_functions
