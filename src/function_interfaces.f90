! Module defining abstract interfaces for multivariable functions.
! Used to declare functions F: R^n -> R for use in numerical methods.
module function_interfaces
   use precision_mod
   implicit none

   abstract interface
      ! Interface for a function F mapping R^n -> R.
      ! Input:  x(:)   — input vector of length n
      ! Output: y(:)   — real value result
      function multivariable_func(x) result(y)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp) :: y
      end function multivariable_func
   end interface

   abstract interface
      ! Interface for a function gradient.
      ! Input:  x(:)   — input vector of length n
      ! Output: g(:)   — output vector of length n
      function gradient(x) result(g)
        import :: dp
        real(dp), intent(in) :: x(:)
        real(dp) :: g(size(x))
      end function gradient
   end interface

end module function_interfaces
