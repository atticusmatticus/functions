! precision module that comes before all other compiled fortran scripts

! kind 4=single 4 byte float
! kind 8=double 8 byte float
! kind 16=quadruple 16 byte float
module prec
	integer, parameter	:: sp = 4	!kind(1.0)
	integer, parameter	:: dp = 8	!kind(1.0d0)
	integer, parameter	:: qp = 16
end module prec
