! precision module that comes before all other compiled fortran scripts

! kind 4=single 4 byte float
! kind 8=double 8 byte float
! kind 16=quadruple 16 byte float
module prec
	integer, parameter :: sp = kind(1.0)
	integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
	integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))
end module prec
