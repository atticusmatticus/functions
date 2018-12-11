! ideal solvent energy and force calculation


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!    MODS    !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! data from the config file.
module idealCfgCL3
	use prec
	integer			:: nAtoms=5
	real(kind=dp)	:: rMin, cosThStepSize, phiStepSize
	real(kind=dp)	:: rEq_CH = 1.1
	real(kind=dp)	:: rEq_CCl = 1.758
	real(kind=dp)	:: aEq_HCl = 1.87937134
	real(kind=dp)	:: AljH = 3.08283e+5
	real(kind=dp)	:: BljH = 2.45437e+2
	real(kind=dp)	:: AljC = 6.93952e+6
	real(kind=dp)	:: BljC = 1.89195e+3
	real(kind=dp)	:: AljCl = 1.53255e+7
	real(kind=dp)	:: BljCl = 3.69121e+3
end module idealCfgCL3


! angle arrays
module idealArray
	use prec
	real(kind=dp),allocatable	:: rAxis(:), pos_i(:,:), pos_r(:,:), pos_t(:,:), Alj(:), Blj(:), theta(:), phi(:)
end module idealArray


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!    SUBS    !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module idealSolv
	use prec
	implicit none

contains

	! run the whole ideal_CL3 code in order.
	subroutine ideal_CL3(rBins,rStepSize,cosThBins,cosTh_min,cosTh_max,phiBins,phi_min,phi_max, idealHist)
		use idealCfgCL3
		use idealArray
		implicit none
		integer			:: rBins, cosThBins, phiBins
		real(kind=dp)	:: rStepSize, cosTh_min, cosTh_max, phi_min, phi_max
		real(kind=dp)	:: idealHist(4,rBins,cosThBins,phiBins)

		! calculate initial geometry and initialize arrays
		call initial_geom

		! compute average force integral.
		call compute_ideal_hist(rBins,rStepSize,cosThBins,cosTh_min,cosTh_max,phiBins,phi_min,phi_max, idealHist)

	end subroutine ideal_CL3


	! generate starting configuration of CL3 @ theta=0, phi=0 and LJ coefficient arrays
	subroutine initial_geom
		use constants
		use idealCfgCL3
		use idealArray
		implicit none
		integer				:: i, j
		real(kind=dp)	:: psi

		allocate( pos_i(nAtoms,3), pos_r(nAtoms,3), pos_t(nAtoms,3), Alj(nAtoms), Blj(nAtoms) )

		pos_i = 0_dp
		do i = 1, nAtoms
			pos_i(i,3) = 1_dp
		end do

		! 1 is H
		pos_i(1,:) = [ 0_dp, -1_dp, 0_dp ]
		pos_i(1,:) = pos_i(1,:) * rEq_CH
		Alj(1) = AljH
		Blj(1) = BljH

		! 2 is C
		pos_i(2,:) = [ 0_dp, 0_dp, 0_dp ]
		Alj(2) = AljC
		Blj(2) = BljC

		! 3-5 are Cl
		! FIXME, BUG: the Cl-C-Cl angle isn't correct.
		psi = 0_dp
		do i = 3, nAtoms
			pos_i(i,1) = -dsin(psi)*dsin(aEq_HCl)
			pos_i(i,2) = -dcos(aEq_HCl)
			pos_i(i,3) = -dsin(aEq_HCl)*dcos(psi)

			pos_i(i,:) = pos_i(i,:) * rEq_CCl

			psi = psi + 2_dp*pi/3_dp

			Alj(i) = AljCl
			Blj(i) = BljCl
		end do

	end subroutine initial_geom


	! compute idealized force histogram
	subroutine compute_ideal_hist(rBins,rStepSize,cosThBins,cosTh_min,cosTh_max,phiBins,phi_min,phi_max, idealHist)
		use constants
		use functions
		use idealCfgCL3
		use idealArray
		implicit none
		integer							:: i, ir, ith, iphi, rBins, cosThBins, phiBins
		real(kind=dp),dimension(3)	:: r_vec, s_vec, t_vec, force_vec
		real(kind=dp)					:: rStepSize, cosTh_min, cosTh_max, phi_min, phi_max, potent_en, &
			& idealHist(4,rBins,cosThBins,phiBins)

		! Distance Axis
		allocate( rAxis(rBins) , theta(cosThBins), phi(phiBins) )
		rAxis = 0_dp

		rMin = 0_dp

		do i = 1, rBins
			rAxis(i) = rMin + ((i-1) * rStepSize + rStepSize/2_dp)
		end do

		! ANGLES
		! Theta
		! tilt off of z
		!cosTh_max = 1_dp
		!cosTh_min = -1_dp
		cosThStepSize = (cosTh_max - cosTh_min) / real(cosThBins, dp)
		do ith = 1, cosThBins
			theta(ith) = acos(((ith-1)*cosThStepSize + cosThStepSize/2_dp) + cosTh_min)
		end do
		!write(*,*) "Config Cos(Theta) Step Size: ", cosThStepSize

		! Phi
		! twist about z
		!phi_max = 2_dp*pi/3_dp ! note: this is sufficient for a molecule with C3 symmetry.
		!phi_min = 0_dp
		phiStepSize = (phi_max - phi_min) / real(phiBins, dp)
		do iphi = 1, phiBins
			phi(iphi) = phi_min + ((iphi-1)*phiStepSize + phiStepSize/2_dp)
		end do
		!write(*,*) "Config Phi Step Size: ", phiStepSize

		!flush(6)

		! Calculate the average force integral for top half of bisecting plane of cylinder
		do ith = 1, cosThBins
			do iphi = 1, phiBins
				! rotate initial geom
				do i = 1, nAtoms
					pos_r(i,:) = rotate_x( theta(ith), rotate_y( phi(iphi), pos_i(i,:) ) )
				end do
				do ir = 1, rBins
					! set the vector to be translated as the rotated vector at the origin, then translate in the y dimension
					pos_t = pos_r
					! translate each solvent atom along y
					do i = 1, nAtoms
						pos_t(i,2) = pos_r(i,2) + rAxis(ir)
					end do !i atoms
					r_vec = [0_dp,0_dp,0_dp] - pos_t(2,:)
					r_vec = r_vec / norm2(r_vec)
					t_vec = cross_product(r_vec, (pos_t(1,:)-pos_t(2,:)))
					t_vec = t_vec / norm2(t_vec)
					s_vec = cross_product(t_vec, r_vec)
					s_vec = s_vec / norm2(s_vec)
					! force calculation
					force_vec = 0_dp
					potent_en = 0_dp
					do i = 1, nAtoms
						potent_en = potent_en + lj_energy(i, pos_t(i,:))
						force_vec = force_vec - lj_force(i, pos_t(i,:))
					end do
					idealHist(1,ir,ith,iphi) = potent_en							! -T*ln(g(r)) == u_pmf or u_dir??? xxx
					idealHist(2,ir,ith,iphi) = dot_product(force_vec, r_vec)	! <f.r>
					idealHist(3,ir,ith,iphi) = dot_product(force_vec, s_vec)	! <f.s>
					idealHist(4,ir,ith,iphi) = dot_product(force_vec, t_vec)	! <f.t>
				end do !r
			end do !phi
		end do !theta

	end subroutine compute_ideal_hist


	! calculate LJ force, (-gradient of PE), for a particular atom type on the solute from a distance 'pos'
	function lj_force(i, pos)
		use idealArray
		implicit none
		integer,intent(in)			:: i											! inputs not to be changed
		real(kind=dp),intent(in)	:: pos(3)									! inputs not to be changed
		real(kind=dp)					:: posDist2, posDist6, lj_force(3)	! output

		posDist2 = pos(1)**2 + pos(2)**2 + pos(3)**2
		posDist6 = posDist2**(-3)

		lj_force = ( posDist6*( 12_dp*Alj(i)*posDist6 - 6_dp*Blj(i) )/posDist2 ) * pos
	end function lj_force


	! calculate LJ potential energy, for a particular atom type on the solute from a distance 'pos'
	function lj_energy(i, pos)
		use idealArray
		implicit none
		integer,intent(in)			:: i											! inputs not to be changed
		real(kind=dp),intent(in)	:: pos(3)									! inputs not to be changed
		real(kind=dp)					:: posDist2, posDist6, lj_energy		! output

		posDist2 = pos(1)**2 + pos(2)**2 + pos(3)**2
		posDist6 = posDist2**(-3)

		lj_energy = posDist6*( Alj(i)*posDist6 - Blj(i) )
	end function lj_energy

end module idealSolv
