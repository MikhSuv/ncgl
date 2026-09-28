! Test program for nclg::cg_min.
! Verifies that the conjugate gradient method finds the minimum of a test
! function, returns the starting point untouched when the required accuracy is
! already met, and terminates on functions that have no minimum.
program test_cg_min
  use precision_mod
  use function_interfaces
  use test_functions
  use nclg
  implicit none

  ! Accuracy requested from the method. The numerical gradient of grad is
  ! accurate to about 1e-11, so this is a limit set by the tests and not by the
  ! precision of the gradient: the norm printed for every test is much smaller.
  real(dp), parameter :: eps = 1.0e-6_dp
  ! Tolerance on the point found.
  real(dp), parameter :: tol = 1.0e-6_dp
  integer :: total = 0
  integer :: failed = 0
  integer :: passed = 0
  logical :: status
  ! Deviation of the result of a test from the expected one.
  real(dp) :: residual
  ! Point found by the last test, reported together with the residual.
  real(dp), allocatable :: x_found(:)

  ! Test 1: quadratic f(x)=sum(x^2), whose minimum is the origin.
  status = test_min(quadratic, [3.0_dp, 4.0_dp], eps, [0.0_dp, 0.0_dp], tol, 5, 100, residual)
  call report(status, total, passed, failed, 1, residual, quadratic)

  ! Test 2: the same from a point close to the minimum, and with the default
  ! restart period and iteration limit of cg_min.
  status = test_min(quadratic, [0.1_dp, 0.1_dp], eps, [0.0_dp, 0.0_dp], tol, 0, 0, residual)
  call report(status, total, passed, failed, 2, residual, quadratic)

  ! Test 3: the same from a point with a negative component.
  status = test_min(quadratic, [-5.0_dp, 2.0_dp], eps, [0.0_dp, 0.0_dp], tol, 5, 100, residual)
  call report(status, total, passed, failed, 3, residual, quadratic)

  ! Test 4: the same in five dimensions, with the default parameters.
  status = test_min(quadratic, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], eps, &
                    [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], tol, 0, 0, residual)
  call report(status, total, passed, failed, 4, residual, quadratic)

  ! Test 5: the memory of the search direction is dropped at every iteration,
  ! which turns the method into steepest descent.
  status = test_min(quadratic, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], eps, &
                    [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], tol, 1, 100, residual)
  call report(status, total, passed, failed, 5, residual, quadratic)

  ! Test 6: the memory is never dropped, since m is larger than the number of
  ! iterations allowed.
  status = test_min(quadratic, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], eps, &
                    [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], tol, 1000, 100, residual)
  call report(status, total, passed, failed, 6, residual, quadratic)

  ! Test 7: starting from the minimum the gradient already satisfies the
  ! criterion, so the starting point must be returned unchanged.
  status = test_same(quadratic, [0.0_dp, 0.0_dp], eps, residual)
  call report(status, total, passed, failed, 7, residual, quadratic)

  ! Test 8: constant function, the gradient is zero everywhere.
  status = test_same(const, [1.0_dp, 2.0_dp], eps, residual)
  call report(status, total, passed, failed, 8, residual, const)

  ! Test 9: an accuracy that is already met at the starting point must be
  ! reported without taking a single step.
  status = test_same(quadratic, [3.0_dp, 4.0_dp], 1.0e+3_dp, residual)
  call report(status, total, passed, failed, 9, residual, quadratic)

  ! Test 10: 1D parabola with the minimum at x(1) = 0.3.
  status = test_min(parabola_c03, [0.0_dp], eps, [0.3_dp], tol, 5, 100, residual)
  call report(status, total, passed, failed, 10, residual, parabola_c03)

  ! Test 11: 1D parabola with the minimum at x(1) = 3.7, far from the start.
  status = test_min(parabola_c37, [0.0_dp], eps, [3.7_dp], tol, 5, 100, residual)
  call report(status, total, passed, failed, 11, residual, parabola_c37)

  ! Test 12: f(x) = sum(sin(x_i)) reaches the global minimum -2 in two
  ! dimensions. Which minimum is found depends on the starting point, so only
  ! the value of the function is checked.
  status = test_value(sinus, [1.0_dp, 1.0_dp], eps, -2.0_dp, tol, residual)
  call report(status, total, passed, failed, 12, residual, sinus)

  ! Test 13: the decreasing linear function has no minimum. The iteration limit
  ! must stop the method and the point returned must still decrease f.
  status = test_downhill(linear, [0.0_dp], eps, 5, 3, residual)
  call report(status, total, passed, failed, 13, residual, linear)

  ! Test 14: the same for the exponential function.
  status = test_downhill(exponential, [1.0_dp, 1.0_dp], eps, 5, 3, residual)
  call report(status, total, passed, failed, 14, residual, exponential)

  ! Test 15: two calls with identical arguments must return the same point.
  call report(test_deterministic(quadratic, [3.0_dp, 4.0_dp], eps, residual), &
              total, passed, failed, 15, residual, quadratic)

  ! The tests below use the stretched quadratics quadratic_k10 and quadratic_k100
  ! of test_functions. Unlike the plain quadratic their condition number is 10 and
  ! 100, so the conjugate directions make a difference: the method reaches the
  ! minimum in a few iterations, while steep descent needs far more.

  ! Test 16: the conjugate directions reach the minimum well within ten
  ! iterations, which a broken beta would not manage.
  status = test_min(quadratic_k10, [1.0_dp, 1.0_dp], eps, [0.0_dp, 0.0_dp], 1.0e-6_dp, 0, 10, residual)
  call report(status, total, passed, failed, 16, residual, quadratic_k10)

  ! Test 17: the same without any periodic restart, since m exceeds the
  ! iteration limit.
  status = test_min(quadratic_k10, [1.0_dp, 1.0_dp], eps, [0.0_dp, 0.0_dp], 1.0e-6_dp, 1000, 10, residual)
  call report(status, total, passed, failed, 17, residual, quadratic_k10)

  ! Test 18: cg_min raises the restart period m to at least ten, so asking for a
  ! shorter one must give exactly the result of m = 10. This pins down the
  ! smallest restart period the method accepts.
  status = test_same_run(quadratic_k10, [1.0_dp, 1.0_dp], eps, 1, 10, 20, residual)
  call report(status, total, passed, failed, 18, residual, quadratic_k10)

  ! Test 19: the same in three dimensions with the condition number 100, where
  ! the conjugate method still wins by several orders of magnitude.
  status = test_min(quadratic_k100, [1.0_dp, 1.0_dp, 1.0_dp], eps, &
                    [0.0_dp, 0.0_dp, 0.0_dp], 1.0e-5_dp, 0, 10, residual)
  call report(status, total, passed, failed, 19, residual, quadratic_k100)

  ! Test 20: a nonconvex function, on which the conjugate direction stops being
  ! a descent direction and cg_min has to drop its memory. The point returned
  ! must still be finite and much better than the start.
  status = test_downhill(wavy, [1.0_dp, 2.0_dp, 0.5_dp], eps, 4, 50, residual)
  call report(status, total, passed, failed, 20, residual, wavy)

  print *, passed, "/", total, " tests passed"

contains

  ! Run cg_min and check that it found the expected minimizer and decreased f.
  ! A restart period m <= 0 and an iteration limit max_iter <= 0 ask for the
  ! defaults of cg_min. The residual is the distance to the expected minimizer.
  ! Returns .true. if the test passed.
  function test_min(func, x0, eps, x_ref, tol, m, max_iter, residual) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    real(dp), intent(in) :: x_ref(size(x0))
    real(dp), intent(in) :: tol
    integer, intent(in) :: m, max_iter
    real(dp), intent(out) :: residual
    real(dp) :: x(size(x0))
    logical :: flag
    x = minimize(func, x0, eps, m, max_iter)
    x_found = x
    residual = norm2(x - x_ref)
    flag = (residual <= tol) .and. (func(x) <= func(x0))
    if (.not. flag) print *, "  expected:", x_ref, " got:", x
  end function test_min

  ! Run cg_min with the default parameters and check that the value of the
  ! function at the point found is the expected one and that f decreased.
  ! The residual is the difference between the value found and the expected one.
  function test_value(func, x0, eps, f_ref, tol, residual) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    real(dp), intent(in) :: f_ref, tol
    real(dp), intent(out) :: residual
    real(dp) :: x(size(x0))
    logical :: flag
    x = cg_min(func, x0, eps)
    x_found = x
    residual = abs(func(x) - f_ref)
    flag = (residual <= tol) .and. (func(x) <= func(x0))
    if (.not. flag) print *, "  expected f =", f_ref, " got:", func(x)
  end function test_value

  ! Check that the restart period m is raised to at least ten: a run with
  ! m = 1 must return exactly what a run with m = 10 returns. The residual is the
  ! distance between the two points found.
  function test_same_run(func, x0, eps, m, m_ref, max_iter, residual) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    integer, intent(in) :: m, m_ref, max_iter
    real(dp), intent(out) :: residual
    real(dp) :: x(size(x0)), x_ref(size(x0))
    logical :: flag
    x = cg_min(func, x0, eps, m, max_iter)
    x_ref = cg_min(func, x0, eps, m_ref, max_iter)
    x_found = x
    residual = maxval(abs(x - x_ref))
    flag = residual == 0.0_dp
    if (.not. flag) print *, "  m =", m, " gave", x, " m =", m_ref, " gave", x_ref
  end function test_same_run

  ! Check that cg_min returns the starting point unchanged when the required
  ! accuracy is already met there. The residual is the largest change made.
  function test_same(func, x0, eps, residual) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    real(dp), intent(out) :: residual
    real(dp) :: x(size(x0))
    logical :: flag
    x = cg_min(func, x0, eps)
    x_found = x
    residual = maxval(abs(x - x0))
    flag = residual == 0.0_dp
    if (.not. flag) print *, "  expected:", x0, " got:", x
  end function test_same

  ! Check that cg_min terminates and decreases f by more than one: used for
  ! functions without a minimum and for nonconvex ones, where only the fact that
  ! the method keeps making progress is known in advance. The residual is the
  ! decrease of f that was achieved.
  function test_downhill(func, x0, eps, m, max_iter, residual) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    integer, intent(in) :: m, max_iter
    real(dp), intent(out) :: residual
    real(dp) :: x(size(x0))
    logical :: flag
    x = cg_min(func, x0, eps, m, max_iter)
    x_found = x
    residual = func(x0) - func(x)
    flag = (norm2(x) <= huge(x)) .and. (residual >= 1.0_dp)
    if (.not. flag) print *, "  got:", x, " f =", func(x), " expected at most:", func(x0) - 1.0_dp
  end function test_downhill

  ! Two identical consecutive calls must return the same point. The residual is
  ! the largest difference between the two points returned.
  function test_deterministic(func, x0, eps, residual) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    real(dp), intent(out) :: residual
    real(dp) :: x1(size(x0)), x2(size(x0))
    logical :: flag
    x1 = cg_min(func, x0, eps)
    x2 = cg_min(func, x0, eps)
    x_found = x1
    residual = maxval(abs(x1 - x2))
    flag = residual == 0.0_dp
    if (.not. flag) print *, "  call1:", x1, " call2:", x2
  end function test_deterministic

  ! Call cg_min, passing the restart period and the iteration limit only when
  ! they are positive.
  function minimize(func, x0, eps, m, max_iter) result(x)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    integer, intent(in) :: m, max_iter
    real(dp) :: x(size(x0))
    if (m > 0 .and. max_iter > 0) then
      x = cg_min(func, x0, eps, m, max_iter)
    else
      x = cg_min(func, x0, eps)
    end if
  end function minimize

  ! Count the result of a test and print the residual it reports together with
  ! the point found and the norm of the gradient there. The norm shows at once
  ! whether the accuracy required by the method was reached.
  subroutine report(status, total, passed, failed, id, residual, func)
    logical, intent(in) :: status
    integer, intent(inout) :: total, passed, failed
    integer, intent(in) :: id
    real(dp), intent(in) :: residual
    procedure(multivariable_func) :: func
    real(dp) :: g(size(x_found))
    total = total + 1
    g = grad(x_found, func)
    print *, "Test ", id, " residual =", residual,&
           " |grad| =", norm2(g)
    if (status) then
      passed = passed + 1
    else
      failed = failed + 1
      print *, "  failed"
    end if
  end subroutine report

end program test_cg_min
