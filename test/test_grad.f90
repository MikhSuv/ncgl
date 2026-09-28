! Test program for nclg.
! Verifies the numerical gradient against the analytic gradient of each
! test function and prints the number of passed tests out of the total
! together with the residual of every one of them.
program check
  use precision_mod
  use function_interfaces
  use test_functions
  use nclg
implicit none

  ! Tolerance for the maximum absolute residual between the numeric and
  ! analytic gradients.
  real(dp), parameter :: tol = 1.0e-5_dp
  ! Test point shared by all tests (overridden per test below). The assignment
  ! sets its size, so a test may use a point of any dimension.
  real(dp), allocatable :: p(:)
  integer :: total = 0
  integer :: failed = 0
  integer :: passed = 0
  logical :: status
  real(dp) :: residual

  ! Test 1: constant function. Gradient should be identically zero.
  p = [1.0_dp]
  status = test_gradient(const, const_grad, p, residual)
  call report(status, total, passed, failed, 1, residual)

  ! Test 2: quadratic function f(x) = sum(x^2).
  p = [2.0_dp]
  status = test_gradient(quadratic, quadratic_grad, p, residual)
  call report(status, total, passed, failed, 2, residual)

  ! Test 3: linear function.
  p = [3.0_dp]
  status = test_gradient(linear, linear_grad, p, residual)
  call report(status, total, passed, failed, 3, residual)

  ! Test 4: exponential function.
  p = [4.0_dp]
  status = test_gradient(exponential, exponential_grad, p, residual)
  call report(status, total, passed, failed, 4, residual)

  ! Test 5: sine function.
  p = [6.0_dp]
  status = test_gradient(sinus, sinus_grad, p, residual)
  call report(status, total, passed, failed, 5, residual)

  ! Test 6: stretched quadratic in two dimensions.
  p = [2.0_dp, 3.0_dp]
  status = test_gradient(quadratic_k10, quadratic_k10_grad, p, residual)
  call report(status, total, passed, failed, 6, residual)

  ! Test 7: stretched quadratic in three dimensions.
  p = [1.0_dp, 2.0_dp, 3.0_dp]
  status = test_gradient(quadratic_k100, quadratic_k100_grad, p, residual)
  call report(status, total, passed, failed, 7, residual)

  ! Test 8: nonconvex function. Its oscillations make the residual an order of
  ! magnitude larger than for the smooth functions, but it stays below tol.
  p = [1.0_dp, 2.0_dp, 3.0_dp]
  status = test_gradient(wavy, wavy_grad, p, residual)
  call report(status, total, passed, failed, 8, residual)

  print *, passed, "/", total, " tests passed" 
  ! A nonzero exit status tells the build system that a test failed.
  if (failed > 0) error stop

  contains
    ! Compare the numerical gradient nclg::grad(func, point) with the
    ! analytic gradient func_grad(point). The residual is the maximum
    ! absolute component-wise difference. Returns .true. if it stays below tol.
    function test_gradient(func, func_grad, point, residual) result(flag)
      procedure(multivariable_func) :: func
      procedure(gradient) :: func_grad
      real(dp), intent(in) :: point(:)
      real(dp), intent(out) :: residual
      logical :: flag
      residual = maxval(abs(grad(point, func) - func_grad(point)))
      flag = residual < tol

    end function test_gradient

    ! Count the result of a test and print the residual it reports.
    subroutine report(status, total, passed, failed, id, residual)
      logical, intent(in) :: status
      integer, intent(inout) :: total, passed, failed
      integer, intent(in) :: id
      real(dp), intent(in) :: residual
      total = total + 1
      print *, "Test ", id, " residual =", residual
      if (status) then
        passed = passed + 1
      else
        failed = failed + 1
        print *, "  failed"
      end if
    end subroutine report
end program check
