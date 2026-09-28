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

  ! Accuracy requested from the method. The numerical gradient of grad is only
  ! good to about 1e-8, so a smaller value is never reached.
  real(dp), parameter :: eps = 1.0e-6_dp
  ! Tolerance on the point found.
  real(dp), parameter :: tol = 1.0e-6_dp
  integer :: total = 0
  integer :: failed = 0
  integer :: passed = 0
  logical :: status

  ! Test 1: quadratic f(x)=sum(x^2), whose minimum is the origin.
  status = test_min(quadratic, [3.0_dp, 4.0_dp], eps, [0.0_dp, 0.0_dp], tol, 5, 100)
  call report(status, total, passed, failed, 1)

  ! Test 2: the same from a point close to the minimum, and with the default
  ! restart period and iteration limit of cg_min.
  status = test_min(quadratic, [0.1_dp, 0.1_dp], eps, [0.0_dp, 0.0_dp], tol, 0, 0)
  call report(status, total, passed, failed, 2)

  ! Test 3: the same from a point with a negative component.
  status = test_min(quadratic, [-5.0_dp, 2.0_dp], eps, [0.0_dp, 0.0_dp], tol, 5, 100)
  call report(status, total, passed, failed, 3)

  ! Test 4: the same in five dimensions, with the default parameters.
  status = test_min(quadratic, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], eps, &
                    [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], tol, 0, 0)
  call report(status, total, passed, failed, 4)

  ! Test 5: the memory of the search direction is dropped at every iteration,
  ! which turns the method into steepest descent.
  status = test_min(quadratic, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], eps, &
                    [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], tol, 1, 100)
  call report(status, total, passed, failed, 5)

  ! Test 6: the memory is never dropped, since m is larger than the number of
  ! iterations allowed.
  status = test_min(quadratic, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], eps, &
                    [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], tol, 1000, 100)
  call report(status, total, passed, failed, 6)

  ! Test 7: starting from the minimum the gradient already satisfies the
  ! criterion, so the starting point must be returned unchanged.
  status = test_same(quadratic, [0.0_dp, 0.0_dp], eps)
  call report(status, total, passed, failed, 7)

  ! Test 8: constant function, the gradient is zero everywhere.
  status = test_same(const, [1.0_dp, 2.0_dp], eps)
  call report(status, total, passed, failed, 8)

  ! Test 9: an accuracy that is already met at the starting point must be
  ! reported without taking a single step.
  status = test_same(quadratic, [3.0_dp, 4.0_dp], 1.0e+3_dp)
  call report(status, total, passed, failed, 9)

  ! Test 10: 1D parabola with the minimum at x(1) = 0.3.
  status = test_min(parabola_c03, [0.0_dp], eps, [0.3_dp], tol, 5, 100)
  call report(status, total, passed, failed, 10)

  ! Test 11: 1D parabola with the minimum at x(1) = 3.7, far from the start.
  status = test_min(parabola_c37, [0.0_dp], eps, [3.7_dp], tol, 5, 100)
  call report(status, total, passed, failed, 11)

  ! Test 12: f(x) = sum(sin(x_i)) reaches the global minimum -2 in two
  ! dimensions. Which minimum is found depends on the starting point, so only
  ! the value of the function is checked.
  status = test_value(sinus, [1.0_dp, 1.0_dp], eps, -2.0_dp, tol)
  call report(status, total, passed, failed, 12)

  ! Test 13: the decreasing linear function has no minimum. The iteration limit
  ! must stop the method and the point returned must still decrease f.
  status = test_downhill(linear, [0.0_dp], eps, 5, 3)
  call report(status, total, passed, failed, 13)

  ! Test 14: the same for the exponential function.
  status = test_downhill(exponential, [1.0_dp, 1.0_dp], eps, 5, 3)
  call report(status, total, passed, failed, 14)

  ! Test 15: two calls with identical arguments must return the same point.
  call report(test_deterministic(quadratic, [3.0_dp, 4.0_dp], eps), &
              total, passed, failed, 15)

  ! The tests below use a stretched quadratic f(x) = x1^2 + 10*x2^2 declared at
  ! the end of this file. Unlike the quadratic of test_functions its condition
  ! number is 10, so the conjugate directions make a difference: the method
  ! reaches the minimum in two iterations, while steep descent needs far more.

  ! Test 16: the conjugate directions reach the minimum well within ten
  ! iterations, which a broken beta would not manage.
  status = test_min(stretched, [1.0_dp, 1.0_dp], eps, [0.0_dp, 0.0_dp], 1.0e-6_dp, 0, 10)
  call report(status, total, passed, failed, 16)

  ! Test 17: the same without any periodic restart, since m exceeds the
  ! iteration limit.
  status = test_min(stretched, [1.0_dp, 1.0_dp], eps, [0.0_dp, 0.0_dp], 1.0e-6_dp, 1000, 10)
  call report(status, total, passed, failed, 17)

  ! Test 18: dropping the memory of the direction at every iteration turns the
  ! method into steep descent, which is not accurate yet after ten iterations.
  ! This pins down what the restart period m does.
  status = test_slow(stretched, [1.0_dp, 1.0_dp], eps, [0.0_dp, 0.0_dp], 1.0e-3_dp, 1, 10)
  call report(status, total, passed, failed, 18)

  ! Test 19: the same in three dimensions with the condition number 100, where
  ! the conjugate method still wins by several orders of magnitude.
  status = test_min(stretched3, [1.0_dp, 1.0_dp, 1.0_dp], eps, &
                    [0.0_dp, 0.0_dp, 0.0_dp], 1.0e-5_dp, 0, 10)
  call report(status, total, passed, failed, 19)

  ! Test 20: a nonconvex function, on which the conjugate direction stops being
  ! a descent direction and cg_min has to drop its memory. The point returned
  ! must still be finite and much better than the start.
  status = test_downhill(wavy, [1.0_dp, 2.0_dp, 0.5_dp], eps, 4, 50)
  call report(status, total, passed, failed, 20)

  print *, passed, "/", total, " tests passed"

contains

  ! Run cg_min and check that it found the expected minimizer and decreased f.
  ! A restart period m <= 0 and an iteration limit max_iter <= 0 ask for the
  ! defaults of cg_min.
  ! Returns .true. if the test passed.
  function test_min(func, x0, eps, x_ref, tol, m, max_iter) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    real(dp), intent(in) :: x_ref(size(x0))
    real(dp), intent(in) :: tol
    integer, intent(in) :: m, max_iter
    real(dp) :: x(size(x0))
    logical :: flag
    x = minimize(func, x0, eps, m, max_iter)
    flag = (norm2(x - x_ref) <= tol) .and. (func(x) <= func(x0))
    if (.not. flag) print *, "  expected:", x_ref, " got:", x
  end function test_min

  ! Run cg_min with the default parameters and check that the value of the
  ! function at the point found is the expected one and that f decreased.
  function test_value(func, x0, eps, f_ref, tol) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    real(dp), intent(in) :: f_ref, tol
    real(dp) :: x(size(x0))
    logical :: flag
    x = cg_min(func, x0, eps)
    flag = (abs(func(x) - f_ref) <= tol) .and. (func(x) <= func(x0))
    if (.not. flag) print *, "  expected f =", f_ref, " got:", func(x)
  end function test_value

  ! Check that cg_min returns a point that is still far from the expected
  ! minimizer after the given number of iterations.
  function test_slow(func, x0, eps, x_ref, tol, m, max_iter) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    real(dp), intent(in) :: x_ref(size(x0))
    real(dp), intent(in) :: tol
    integer, intent(in) :: m, max_iter
    real(dp) :: x(size(x0))
    logical :: flag
    x = minimize(func, x0, eps, m, max_iter)
    flag = norm2(x - x_ref) > tol
    if (.not. flag) print *, "  expected a distance above", tol, " got:", norm2(x - x_ref)
  end function test_slow

  ! Check that cg_min returns the starting point unchanged when the required
  ! accuracy is already met there.
  function test_same(func, x0, eps) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    real(dp) :: x(size(x0))
    logical :: flag
    x = cg_min(func, x0, eps)
    flag = all(x == x0)
    if (.not. flag) print *, "  expected:", x0, " got:", x
  end function test_same

  ! Check that cg_min terminates and decreases f by more than one: used for
  ! functions without a minimum and for nonconvex ones, where only the fact that
  ! the method keeps making progress is known in advance.
  function test_downhill(func, x0, eps, m, max_iter) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    integer, intent(in) :: m, max_iter
    real(dp) :: x(size(x0))
    logical :: flag
    x = cg_min(func, x0, eps, m, max_iter)
    flag = (norm2(x) <= huge(x)) .and. (func(x) <= func(x0) - 1.0_dp)
    if (.not. flag) print *, "  got:", x, " f =", func(x), " expected at most:", func(x0) - 1.0_dp
  end function test_downhill

  ! Two identical consecutive calls must return the same point.
  function test_deterministic(func, x0, eps) result(flag)
    procedure(multivariable_func) :: func
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    logical :: flag
    flag = all(cg_min(func, x0, eps) == cg_min(func, x0, eps))
    if (.not. flag) print *, "  the two calls disagree"
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

  ! Stretched quadratic f(x) = x1^2 + 10*x2^2 with the minimum at the origin.
  ! Its condition number is 10, so the iterates of the conjugate method and of
  ! steep descent differ by orders of magnitude within ten iterations.
  function stretched(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y
    y = x(1)**2 + 10.0_dp * x(2)**2
  end function stretched

  ! Stretched quadratic f(x) = x1^2 + 10*x2^2 + 100*x3^2, condition number 100.
  function stretched3(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y
    y = x(1)**2 + 10.0_dp * x(2)**2 + 100.0_dp * x(3)**2
  end function stretched3

  ! Bounded below but nonconvex f(x) = x1^2 + 10*x2^2 + 50*sin(10*x1) + 20*sin(5*x3).
  ! The oscillations make the conjugate direction leave the cone of descent, so
  ! the restart on a bad direction is taken repeatedly here.
  function wavy(x) result(y)
    real(dp), intent(in) :: x(:)
    real(dp) :: y
    y = x(1)**2 + 10.0_dp * x(2)**2 + 50.0_dp * sin(10.0_dp * x(1)) &
        + 20.0_dp * sin(5.0_dp * x(3))
  end function wavy

end program test_cg_min
