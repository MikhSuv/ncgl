! Test program for nclg.
! Verifies the numerical gradient against the analytic gradient of each
! test function and prints the number of passed tests out of the total.
program check
  use precision_mod
  use function_interfaces
  use test_functions
  use nclg
implicit none

  ! Tolerance for the maximum absolute residual between the numeric and
  ! analytic gradients.
  real(dp), parameter :: tol = 1.0e-5_dp
  ! Test point shared by all tests (overridden per test below).
  real(dp), dimension(10) :: p 
  integer :: total = 0
  integer :: failed = 0
  integer :: passed = 0
  logical :: status

  ! Test 1: constant function. Gradient should be identically zero.
  p = 1.0_dp
  status = test_gradient(const, const_grad, p)
  total = total + 1
  if (status) then 
    passed = passed + 1
  else 
    failed = failed + 1
    print *, "Test ", total, " failed"
  end if

  ! Test 2: quadratic function f(x) = sum(x^2).
  p = 2.0_dp
  status = test_gradient(quadratic, quadratic_grad, p )
  total = total + 1
  if (status) then 
    passed = passed + 1
  else 
    failed = failed + 1
    print *, "Test ", total, " failed"
  end if


  ! Test 3: linear function.
  p = 3.0_dp
  status = test_gradient(linear, linear_grad, p )
  total = total + 1
  if (status) then 
    passed = passed + 1
  else 
    failed = failed + 1
    print *, "Test ", total, " failed"
  end if

  ! Test 4: exponential function.
  p = 4.0_dp
  status = test_gradient(exponential, exponential_grad, p )
  total = total + 1
  if (status) then 
    passed = passed + 1
  else 
    failed = failed + 1
    print *, "Test ", total, " failed"
  end if

  ! Test 5: logarithm function (test point must stay positive).
  p = 5.0_dp
  status = test_gradient(logarithm, logarithm_grad, p )
  total = total + 1
  if (status) then 
    passed = passed + 1
  else 
    failed = failed + 1
    print *, "Test ", total, " failed"
  end if

  ! Test 6: sine function.
  p = 6.0_dp
  status = test_gradient(sinus, sinus_grad, p )
  total = total + 1
  if (status) then 
    passed = passed + 1
  else 
    failed = failed + 1
    print *, "Test ", total, " failed"
  end if

  print *, passed, "/", total, " tests passed" 

  contains
    ! Compare the numerical gradient nclg::grad(func, point) with the
    ! analytic gradient func_grad(point). Returns .true. if the maximum
    ! absolute component-wise difference stays below tol.
    function test_gradient(func, func_grad, point) result(flag)
      procedure(multivariable_func) :: func
      procedure(gradient) :: func_grad
      real(dp), intent(in) :: point(:)
      real(dp) :: residual
      logical :: flag
        residual = maxval(abs(grad(point, func) - func_grad(point)))
      if (residual < tol ) then 
        flag = .true.
      else 
        flag = .false.
        print *, "residual:", residual
      end if

    end function test_gradient
end program check
