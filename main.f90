

Subroutine err_finalize(err_message)
  use Global_Variables
  use communication
  implicit none
  character(*),intent(in) :: err_message
  if (comm_is_root()) then
    write(*,*) err_message
  endif
  call comm_finalize

  stop
End Subroutine Err_finalize


Program main
  use Global_Variables, only: cfunction
  use arted_sc, main_sc => main
  use arted_ms, main_ms => main
  implicit none
  
  read(*,*) cfunction
  select case(cfunction)
  case ("arted_sc")
    call main_sc
  case ("arted_ms")
    call main_ms
  case default
    write(*,*) "Invalid cfunction parameter"
  end select
  
  stop
End Program main
