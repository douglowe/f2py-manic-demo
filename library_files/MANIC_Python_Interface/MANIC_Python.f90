subroutine manic_python(N_P_SPEC, C_P, CNP_P, CFACTOR_P, &
						traj_length_P, &
						time_traj_P, temp_traj_P, pres_traj_P, &
						zmbl_traj_P, tide_traj_P, sun_traj_P, &
						DT_P, TSTEPMAX_P, PDFITE_P, &
						M_P, ND_P, DU_P, UEDGE_P, R_R_P, &
						tstart_P, tend_P, tout_P, start_P )

!f2py integer, intent(in) :: N_P_SPEC
!f2py real(kind=dp_alt), intent(in, out), depend(NP_SPEC) :: C_P(N_P_SPEC)
!f2py real(kind=dp_alt), intent(in, out), depend(NP_SPEC) :: CNP_P(N_P_SPEC)
!f2py real(kind=dp_alt), intent(in) :: CFACTOR_P
!f2py integer, intent(in) :: traj_length_P
!f2py real(kind=dp_alt), intent(in), depend(traj_length_P) :: time_traj_P(traj_length_P)
!f2py real(kind=dp_alt), intent(in), depend(traj_length_P) :: temp_traj_P(traj_length_P)
!f2py real(kind=dp_alt), intent(in), depend(traj_length_P) :: pres_traj_P(traj_length_P)
!f2py real(kind=dp_alt), intent(in), depend(traj_length_P) :: zmbl_traj_P(traj_length_P)
!f2py real(kind=dp_alt), intent(in), depend(traj_length_P) :: tide_traj_P(traj_length_P)
!f2py real(kind=dp_alt), intent(in), depend(traj_length_P) :: sun_traj_P(traj_length_P)
!f2py real(kind=dp_alt), intent(in) :: DT_P
!f2py real(kind=dp_alt), intent(in) :: TSTEPMAX_P
!f2py logical, intent(in) :: PDFITE_P
!f2py integer, intent(in) :: M_P
!f2py integer, intent(in) :: ND_P
!f2py real(kind=dp_alt), intent(in), depend(M_P, ND_P) :: DU_P(M_P, ND_P)
!f2py real(kind=dp_alt), intent(in), depend(M_P, ND_P) :: UEDGE_P(M_P+1, ND_P)
!f2py real(kind=dp_alt), intent(in), depend(ND_P) :: R_R_P(ND_P)
!f2py real(kind=dp_alt), intent(in) :: tstart_P
!f2py real(kind=dp_alt), intent(in) :: tend_P
!f2py real(kind=dp_alt), intent(out) :: tout_P
!f2py logical, intent(in) :: start_P

!
!  Connecting subroutine for running the MANIC box model
!  from a python script.
!

	use MAP_Model
	use MAP_Python_Initialize
!	use model_routines
	use physical_parameters
	implicit none

! KPP DP - Double precision kind
INTEGER, PARAMETER :: dp_alt = SELECTED_REAL_KIND(14,300)


!!! communication variables
	integer, intent(in) :: N_P_SPEC, traj_length_P, M_P, ND_P
	real(kind=dp_alt), intent(inout) :: C_P(N_P_SPEC)
	real(kind=dp_alt), intent(inout) :: CNP_P(N_P_SPEC)
	real(kind=dp_alt), intent(in) :: CFACTOR_P, DT_P, TSTEPMAX_P, tstart_P, tend_P
	logical, intent(in) :: PDFITE_P, START_P
	real(kind=dp_alt), intent(in), dimension(traj_length_P) :: &
						time_traj_P, temp_traj_P, pres_traj_P, &
						zmbl_traj_P, tide_traj_P, sun_traj_P
	real(kind=dp_alt), intent(in) :: DU_P(M_P,ND_P)
	real(kind=dp_alt), intent(in) :: UEDGE_P(M_P+1,ND_P)
	real(kind=dp_alt), intent(in) :: R_R_P(ND_P)
	real(kind=dp_alt), intent(out) :: tout_P

!!! local variables for controlling integrator system
	integer :: i,j, cou, ntotal
	REAL(kind=dp) :: RSTATE(20)
	real(kind=dp) :: rcntrl(20)
	real(kind=dp) :: t1
	integer :: istatus(20)
	! p equalivalent to C -> molec cm_air^-3
	real(kind=dp), dimension(nsol,m,nd) :: p


!---> initialisation

	C = C_P
	CNP = CNP_P
	CFACTOR = CFACTOR_P
	DT = DT_P
	TSTEPMAX = TSTEPMAX_P
	tstart = tstart_P
	tend = tend_P
	PDFITE = PDFITE_P
	DU = DU_P
	UEDGE = UEDGE_P
	R_R = R_R_P
	
	! initialise trajectory data with negative numbers
	time_traj = -99d0
	temp_traj = -99d0
	pres_traj = -99d0
	zmbl_traj = -99d0
	tide_traj = -99d0
	sun_traj  = -99d0
	
	! load up our real trajectory data
	time_traj(1:traj_length_P) = time_traj_P
	temp_traj(1:traj_length_P) = temp_traj_P
	pres_traj(1:traj_length_P) = pres_traj_P
	zmbl_traj(1:traj_length_P) = zmbl_traj_P
	tide_traj(1:traj_length_P) = tide_traj_P
	sun_traj(1:traj_length_P)  = sun_traj_P

	! set tolerances
	do i = 1,nvar
	  rtol(i) = 1.0d-6
	  atol(i) = 1.0d-21
	end do ! i=1,nc
	rcntrl = 0d0
	ntotal = 0
	rstate = 0d0

	time = tstart

    ! start off the trajectory variables
    temp  = temp_traj(1)
    pres  = pres_traj(1)
    z_mbl = zmbl_traj(1)
    tidal_height = tide_traj(1)
    SUN   = sun_traj(1)
    traj_species = 0


	! Set up the number array
	call Update_NUMBER(C(NVARST:NVAR),C(NFIXST:NSPEC))

	! Initialise kpp model variables (reduced function,
	!         much of this work will now be done in the python 
	!         script calling this program).
	call Python_Initialize(START_P)
	! if this is the first call to MANIC then we want to copy the CNP array back out
	if(START_P)then
		CNP_P = CNP
	end if


	rstate(3) = 1d-5   ! initial timestep, this will be set by solver later
	cou = 0            ! counter for timesteps taken

	istatus = 0


		write(6,*) '  radii1      number1      radii2        number2'
	  	do i=1,m
		    write(6,'(4e12.4)') radius(i,1), Z(i,1), radius(i,2), Z(i,2)
	  	end do
		write(6,*) SPC_NAMES(ind_CH2I2), C(ind_CH2I2)/cfactor !, cnp(15)/cfactor
		write(6,*) SPC_NAMES(ind_IO), C(ind_IO)/cfactor !, cnp(15)/cfactor
		write(6,*) time, tout, tend


!---> time loop
	time_loop : do
		cou = cou + 1

		! determine what tout should be
		tstep = min(tend-time, tstepmax)

		rcntrl(3) = rstate(3)   ! after 1st step use step size suggested by ODE solver

		! call the ODE solver (internal time for solver is just from 0 to tstep,
		!                        time dependent calculations are done outside)
		call integrate(tin = 0d0 ,tout = tstep, & 
			rstatus_u = rstate, rcntrl_u = rcntrl, istatus_u = istatus, &
			icntrl_u = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) )
		Ntotal = Ntotal + ISTATUS(3)

		write(6,*) 'timestep ', cou, ' was ', RSTATE(1)

		! calculate time of next step
		if((RSTATE(1)-time).eq.0d0) write(6,*) 'timestep zero!'
		time = time + RSTATE(1)

	  	write(6,'(i6,e12.4)') cou, RSTATE(1)
		write(6,*) '  radii1      number1      radii2        number2'
	  	do i=1,m
		    write(6,'(4e12.4)') radius(i,1), Z(i,1), radius(i,2), Z(i,2)
	  	end do
		write(6,*) SPC_NAMES(ind_CH2I2), C(ind_CH2I2)/cfactor !, cnp(15)/cfactor
		write(6,*) SPC_NAMES(ind_IO), C(ind_IO)/cfactor !, cnp(15)/cfactor
		write(6,*) time, tout, tend
    	PRINT*,'NSTEPS=',ISTATUS(3),' (',Ntotal,')   (',cou,')'


		! exit at end of time
		if (time >= tend) exit time_loop
	end do time_loop
!---> time loop end

	C_P = C
	tout_P = time


end subroutine manic_python