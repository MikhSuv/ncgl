! Test program for nclg::find_alpha.
! Verifies that find_alpha(x, p, eps, f) returns a step alpha >= 0 that satisfies
! the Armijo and Wolfe conditions for g(t) = f(x + t*p) and, where the line
! function is a parabola, that it lies in the interval where those conditions
! hold. Also checks that the search terminates on rays without a minimum.
program test_find_alpha
  use precision_mod
  use function_interfaces
  use test_functions
  use nclg
  implicit none

  ! Accuracy requested from the line search.
  real(dp), parameter :: eps = 1.0e-14_dp
  integer :: total = 0
  integer :: failed = 0
  integer :: passed = 0
  logical :: status

  ! The curvature condition |g'(alpha)| <= c2*|g'(0)| is satisfied by every step
  ! of a line function g(t) = g(0) + g'(0)*t + q*t^2/2 with g'(0) < 0 exactly for
  ! alpha in [(1 - c2)*|g'(0)|/q, (1 + c2)*|g'(0)|/q]. The intervals below are
  ! this range for the test functions.

  ! Test 1: quadratic f(x)=sum(x^2), x=(3,0), p=-grad=(-6,0).
  !   g(t) = (3-6t)^2 => g'(0) = -36, q = 72 => alpha in [0.45, 0.55].
  status = test_step(quadratic, [3.0_dp, 0.0_dp], [-6.0_dp, 0.0_dp], eps, 0.45_dp, 0.55_dp)
  call report(status, total, passed, failed, 1)

  ! Test 2: quadratic at the minimum x=(0,0): p is not a descent direction, so
  ! the search must return the empty step.
  status = test_step(quadratic, [0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], eps, 0.0_dp, 0.0_dp)
  call report(status, total, passed, failed, 2)

  ! Test 3: parabola with min at c=0.3 inside the initial step.
  !   g(t) = (t-0.3)^2/2 => g'(0) = -0.3, q = 1 => alpha in [0.27, 0.33].
  status = test_step(parabola_c03, [0.0_dp], [1.0_dp], eps, 0.27_dp, 0.33_dp)
  call report(status, total, passed, failed, 3)

  ! Test 4: parabola with min at c=1.0, exactly one step away.
  !   g(t) = (t-1)^2/2 => alpha in [0.9, 1.1].
  status = test_step(parabola_c10, [0.0_dp], [1.0_dp], eps, 0.9_dp, 1.1_dp)
  call report(status, total, passed, failed, 4)

  ! Test 5: parabola with min at c=3.7, far beyond the initial step: the
  ! bracketing stage has to grow the step to reach the minimum.
  !   g(t) = (t-3.7)^2/2 => alpha in [3.33, 4.07].
  status = test_step(parabola_c37, [0.0_dp], [1.0_dp], eps, 3.33_dp, 4.07_dp)
  call report(status, total, passed, failed, 5)

  ! Test 6: constant function: no descent direction, the empty step is correct.
  status = test_step(const, [0.0_dp], [1.0_dp], eps, 0.0_dp, 0.0_dp)
  call report(status, total, passed, failed, 6)

  ! Test 7: increasing linear function: minimum on the ray is at alpha = 0.
  status = test_step(linear, [0.0_dp], [1.0_dp], eps, 0.0_dp, 0.0_dp)
  call report(status, total, passed, failed, 7)

  ! Test 8: strictly decreasing linear function: no minimum on the ray.
  ! Regression test: the search must terminate and return a finite alpha >= 0.
  status = test_no_garbage(linear, [2.0_dp], [-1.0_dp], eps)
  call report(status, total, passed, failed, 8)

  ! Test 9: exponential function along a steepest descent direction: g decays
  ! without ever turning up, and the curvature condition is met long before the
  ! exponential underflows. Regression test for the search getting stuck while
  ! narrowing a bracket far from the origin, where the step became too small to
  ! move the point in floating point arithmetic.
  status = test_no_garbage(exponential, [0.0_dp, 0.0_dp], [-1.0_dp, -1.0_dp], eps)
  call report(status, total, passed, failed, 9)

  ! Test 10: quadratic with the minimum at x=(1000,0) and p = -grad: the line
  ! function decreases over a range of steps far wider than the first one.
  !   g(t) = (1000-2000t)^2 => g'(0) = -4e6, q = 8e6 => alpha in [0.45, 0.55].
  status = test_step(quadratic, [1000.0_dp, 0.0_dp], [-2000.0_dp, 0.0_dp], eps, 0.45_dp, 0.55_dp)
  call report(status, total, passed, failed, 10)

  ! Test 11: repeated calls with identical arguments return the same result.
  call report(test_deterministic(parabola_c37), total, passed, failed, 11)
  
  ! Test 12: coarse eps stops the zoom stage earlier but must not break the
  ! Wolfe conditions.
  status = test_step(parabola_c37, [0.0_dp], [1.0_dp], 1.0e-4_dp, 3.33_dp, 4.07_dp)
  call report(status, total, passed, failed, 12)

  print *, passed, "/", total, " tests passed"

contains

  ! Run find_alpha and check that the returned step is finite, non-negative,
  ! lies in [lo, hi] and satisfies the Armijo and Wolfe conditions.
  ! Returns .true. if the test passed.
  function test_step(func, x, p, eps, lo, hi) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: p(size(x))
    real(dp), intent(in) :: eps
    real(dp), intent(in) :: lo, hi
    real(dp), parameter :: c1 = 1.0e-4_dp
    real(dp), parameter :: c2 = 0.1_dp
    ! Slack absorbing the error of the numerical gradient.
    real(dp), parameter :: noise = 1.0e-6_dp
    real(dp) :: alpha, f0, ft, d0, dt
    logical :: flag

    alpha = find_alpha(x, p, eps, func)
    f0 = func(x)
    ft = func(x + alpha * p)
    d0 = dot_product(p, grad(x, func))
    dt = dot_product(p, grad(x + alpha * p, func))

    flag = (alpha >= lo) .and. (alpha <= hi) .and. (alpha <= huge(alpha))
    ! The conditions only constrain a step along a descent direction: otherwise
    ! the search returns alpha = 0 and there is nothing left to check.
    if (d0 < 0.0_dp .and. alpha > 0.0_dp) then
      flag = flag .and. (ft <= f0 + c1 * alpha * d0 + noise * abs(f0))
      flag = flag .and. (abs(dt) <= -c2 * d0 + noise * abs(d0))
    end if
    if (.not. flag) print *, "  expected alpha in [", lo, ",", hi, "] got:", alpha
  end function test_step

  ! Check that find_alpha terminates and does not return garbage when the
  ! function has no minimum along the ray (strictly decreasing): the result
  ! must be finite and non-negative.
  function test_no_garbage(func, x, p, eps) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: p(size(x))
    real(dp), intent(in) :: eps
    real(dp) :: alpha
    logical :: flag
    alpha = find_alpha(x, p, eps, func)
    flag = (alpha >= 0.0_dp) .and. (alpha <= huge(alpha))
    if (.not. flag) print *, "  got:", alpha
  end function test_no_garbage

  ! Two identical consecutive calls must return byte-identical results.
  function test_deterministic(func) result(flag)
    procedure(multivariable_func) :: func
    real(dp) :: a1, a2
    logical :: flag
    a1 = find_alpha([0.0_dp], [1.0_dp], eps, func)
    a2 = find_alpha([0.0_dp], [1.0_dp], eps, func)
    flag = a1 == a2
    if (.not. flag) print *, "  call1:", a1, " call2:", a2
  end function test_deterministic

  ! Count the result of a test and print the id of the one that failed.
  subroutine report(status, total, passed, failed, id)
    logical, intent(in) :: status
    integer, intent(inout) :: total, passed, failed
    integer, intent(in) :: id
    total = total + 1
    if (status) then
      passed = passed + 1
    else
      failed = failed + 1
      print *, "Test ", id, " failed"
    end if
  end subroutine report

end program test_find_alpha
