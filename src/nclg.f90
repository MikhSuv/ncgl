module nclg
  use precision_mod
  use function_interfaces
  implicit none
  private

  public :: grad, find_alpha, cg_min
contains

  ! Compute the gradient of function func at point x0.
  ! Uses central finite differences with step size delta_x = eps^(1/3), which
  ! balances the truncation error against the rounding error of f.
  ! Arguments:
  !   x0(:) — point in R^n at which to evaluate the gradient
  !   f     — multivariable function conforming to the multivariable_func interface
  ! Returns:
  !   g(n)  — numerical approximation of the gradient
  function grad(x0, f) result(g)
    real(dp), intent(in) :: x0(:)
    procedure(multivariable_func) :: f
    real(dp) :: g(size(x0))
    real(dp) :: x(size(x0)), x_(size(x0))
    real(dp) :: f0
    real(dp), parameter :: delta_x = epsilon(1.0_dp)**(1.0_dp/3.0_dp)
    integer i, n
    n = size(x0)
    do i = 1,n
      x = x0
      x_ = x0
      x(i) = x(i) + delta_x
      x_(i) = x_(i) - delta_x
      g(i) = f(x) - f(x_)
    end do
    g = 0.5_dp * g/delta_x
  end function grad

  ! Find alpha >= 0 minimizing g(t) = f(x + t*p) over the ray t >= 0 with the
  ! Armijo and Wolfe conditions. Two stages: the step is doubled until the
  ! minimum is bracketed, then the bracket is narrowed by safeguarded quadratic
  ! interpolation. Both stages are bounded, so the search always terminates.
  ! Arguments:
  !   x(:)  — point in R^n at which the search starts
  !   p(:)  — search direction, a descent direction of f for x to be useful
  !   eps   — narrowest bracket accepted by the zoom stage
  !   f     — multivariable function conforming to the multivariable_func interface
  ! Returns:
  !   alpha — step length >= 0 satisfying the Armijo and Wolfe conditions, or the
  !           point with the smallest g among those examined if there is no such step
  function find_alpha(x, p, eps, f) result(alpha)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: p(size(x))
    real(dp), intent(in) :: eps
    procedure(multivariable_func) :: f

    real(dp) :: alpha
    ! Sufficient decrease (Armijo) and curvature (Wolfe) parameters.
    real(dp), parameter :: c1 = 1.0e-4_dp
    real(dp), parameter :: c2 = 0.1_dp
    ! Factor the step grows by while the minimum is being bracketed.
    real(dp), parameter :: grow = 2.0_dp
    ! Longest step considered along the ray.
    real(dp), parameter :: alpha_max = 1.0e+3_dp
    integer, parameter :: max_iter = 50
    real(dp) :: d0, d_t, d_lo
    real(dp) :: f0, f_t, f_lo, f_hi
    real(dp) :: t, t_lo, t_hi
    logical  :: bracketed, found
    integer  :: k

    f0 = f(x)
    d0 = dot_product(p, grad(x, f))
    ! A direction that does not decrease f at x has no minimum to look for:
    ! the empty step is the best admissible one.
    if (d0 >= 0.0_dp) then
      alpha = 0.0_dp
      return
    end if

    ! A step of unit length in x-space keeps the search independent of the scale of p.
    t = min(1.0_dp, 1.0_dp / sqrt(dot_product(p, p)))
    ! Every accepted point is a new minimum of g, so t_lo always keeps the best
    ! point seen so far together with the slope of g there.
    t_lo = 0.0_dp
    f_lo = f0
    d_lo = d0
    f_hi = f0
    bracketed = .false.
    found = .false.

    ! Bracket the minimum by doubling the step while g keeps decreasing.
    do k = 1, max_iter
      f_t = f(x + t * p)
      if (f_t > f0 + c1 * t * d0 .or. f_t >= f_lo) then
        t_hi = t
        f_hi = f_t
        bracketed = .true.
        exit
      end if
      d_t = dot_product(p, grad(x + t * p, f))
      t_lo = t
      f_lo = f_t
      d_lo = d_t
      if (abs(d_t) <= -c2 * d0) then
        alpha = t
        found = .true.
        exit
      end if
      ! The ray has no minimum within the steps allowed: t_lo is the best of them.
      if (t >= alpha_max) exit
      t = min(grow * t, alpha_max)
    end do

    ! Zoom: every trial point is compared with the best one, which is never
    ! dropped from the bracket.
    if (bracketed .and. .not. found) then
      do k = 1, max_iter
        t = min_quadratic(t_lo, f_lo, d_lo, t_hi, f_hi)
        f_t = f(x + t * p)
        if (f_t > f0 + c1 * t * d0 .or. f_t >= f_lo) then
          t_hi = t
        else
          d_t = dot_product(p, grad(x + t * p, f))
          if (abs(d_t) <= -c2 * d0) then
            alpha = t
            found = .true.
            exit
          end if
          ! g'(t) turned positive: the minimum is on the other side of t_lo now.
          if (d_t * (t_hi - t_lo) >= 0.0_dp) t_hi = t_lo
          t_lo = t
          f_lo = f_t
          d_lo = d_t
        end if
        if (abs(t_hi - t_lo) <= eps) exit
      end do
    end if

    ! No step satisfies the Wolfe conditions: keep the point with the smallest g.
    if (.not. found) alpha = t_lo

  end function find_alpha

  ! Fletcher-Reeves coefficient beta = (g·g)/(g_prev·g_prev).
  ! Arguments:
  !   g_prev(:) — gradient at the previous iterate
  !   g(:)     — gradient at the current iterate
  ! Returns:
  !   beta — coefficient of the previous search direction
  function find_beta_fr(g_prev, g) result(beta)
   real(dp), intent(in) :: g_prev(:)
   real(dp), intent(in) :: g(:)
   real(dp) :: beta

   beta = dot_product(g, g) / dot_product(g_prev, g_prev)

  end function find_beta_fr

  ! Polak-Ribiere coefficient beta = (g·(g - g_prev))/(g_prev·g_prev).
  ! Arguments:
  !   g_prev(:) — gradient at the previous iterate
  !   g(:)     — gradient at the current iterate
  ! Returns:
  !   beta — coefficient of the previous search direction
  function find_beta_pr(g_prev, g) result(beta)
   real(dp), intent(in) :: g_prev(:)
   real(dp), intent(in) :: g(:)
   real(dp) :: beta

   beta = dot_product(g, g - g_prev) / dot_product(g_prev, g_prev)

  end function find_beta_pr

  ! Minimize f with the conjugate gradient method.
  ! Starting from x0 the search direction is p = -g, then it is updated with
  ! p = -g + beta*p using the Polak-Ribiere coefficient. The memory of the
  ! previous direction is dropped every m iterations and whenever the new
  ! gradient no longer decreases f enough, that is when g·g_prev > gamma*|g|^2.
  ! Arguments:
  !   f           — multivariable function conforming to the multivariable_func interface
  !   x0(:)       — starting point in R^n
  !   eps         — required accuracy, the iteration stops at |grad f(x)| <= eps
  !   m           — restart period, the default is n + 1 iterations
  !   max_iter    — iteration limit, the default is 1000
  ! Returns:
  !   x_min — point with the smallest f reached; the last one if the iteration
  !           limit is reached before the required accuracy
  function cg_min(f, x0, eps, m, max_iter) result(x_min)
    procedure(multivariable_func) :: f
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: eps
    integer, intent(in), optional :: m
    integer, intent(in), optional :: max_iter

    real(dp) :: x_min(size(x0))
    real(dp) :: g(size(x0)), g_prev(size(x0)), p(size(x0))
    real(dp) :: alpha, beta
    ! Restart on a direction that stops decreasing f enough.
    real(dp), parameter :: gamma = 0.2_dp
    integer :: i, iter_max, restart

    iter_max = 1000000
    restart = size(x0) + 1
    if (present(max_iter)) iter_max = max_iter
    if (present(m)) restart = max(m, 10)

    x_min = x0
    g = grad(x0, f)
    p = -1.0_dp * g

    do i = 1, iter_max
      if (norm2(g) < eps) return
      g_prev = g
      alpha = find_alpha(x_min, p, eps, f)
      x_min = x_min + alpha * p
      g = grad(x_min, f)
      if (norm2(g) <= eps) return
      beta = find_beta_pr(g_prev, g)
      p = -1.0_dp * g + beta * p
      if (mod(i, restart) == 0) p = -1.0_dp * g
      if (dot_product(g_prev, g) > gamma * norm2(g)**2) p = -1.0_dp * g
    end do

  end function cg_min

  ! Argument of the minimum of the quadratic passing through (t_a, f_a) with
  ! slope d_a at t_a and through (t_b, f_b). The midpoint of the two points is
  ! returned when the fit is not convex or its minimum falls outside the
  ! safeguarded part of the interval, which keeps the zoom stage bracketed.
  function min_quadratic(t_a, f_a, d_a, t_b, f_b) result(t)
    real(dp), intent(in) :: t_a, f_a, d_a, t_b, f_b
    real(dp) :: t
    real(dp), parameter :: safeguard = 0.1_dp
    real(dp) :: dt, c, lo, hi, t_try

    dt = t_b - t_a
    t = 0.5_dp * (t_a + t_b)
    c = 0.0_dp
    if (dt /= 0.0_dp) c = (f_b - f_a - d_a * dt) / (dt * dt)
    lo = min(t_a, t_b) + safeguard * abs(dt)
    hi = max(t_a, t_b) - safeguard * abs(dt)
    if (c > 0.0_dp) then
      t_try = t_a - 0.5_dp * d_a / c
      if (t_try > lo .and. t_try < hi) t = t_try
    end if
  end function min_quadratic
end module nclg
