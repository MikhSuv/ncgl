module nclg
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, nclg!"
  end subroutine say_hello
end module nclg
