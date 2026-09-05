module sbl_gate_module
  use, intrinsic :: iso_c_binding
  implicit none

  type, bind(C) :: gate_column_state
    real(c_double) :: prev_v_f(2)   ! Fast eigenvector continuation
    integer(c_int) :: is_active     ! Discrete gate state G_n
  end type gate_column_state

contains

  subroutine sbl_fast_gate(e, shear2, n2, state, lambda_f, c_v, gate_active) &
             bind(C, name="sbl_fast_gate")
    real(c_double), value, intent(in)      :: e, shear2, n2
    type(gate_column_state), intent(inout) :: state
    real(c_double), intent(out)            :: lambda_f, c_v
    integer(c_int), intent(out)            :: gate_active
    ! Computes fast eigenvalue \lambda_f, mode overlap c_v, and updates G_n
  end subroutine sbl_fast_gate
end module sbl_gate_module