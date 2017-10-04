! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialization File
! 
! Generated by KPP-2.2 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : simple_Initialize.f90
! Time                 : Tue May  6 12:17:10 2008
! Working directory    : /Users/lowe/work/manchester/chemistry/new-MAP/combined-MAP-code
! Equation file        : simple.kpp
! Output root filename : simple
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE MAP_Initialize

  USE MAP_Parameters, ONLY: dp, NVAR, NFIX
  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Initialize - function to initialize concentrations
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Initialize ( )


  USE MAP_Global

  INTEGER :: i
  REAL(kind=dp) :: x

  CFACTOR = 2.550000e+07_dp

  x = (1.0E-11)*CFACTOR
  DO i = 1, NVAR
    VAR(i) = x
  END DO

  x = (1.0E-11)*CFACTOR
  DO i = 1, NFIX
    FIX(i) = x
  END DO

! constant rate coefficients
  RCONST(3) = 2.2e-10
  RCONST(21) = 4e-12
  RCONST(25) = 2.6e-22
  RCONST(45) = 5e-14
  RCONST(47) = 5e-14
  RCONST(48) = 3e-12
  RCONST(51) = 1e-11
  RCONST(53) = 5.8e-16
  RCONST(55) = 1.4e-15
  RCONST(59) = 2e-12
  RCONST(60) = 1e-13
  RCONST(61) = 4e-13
  RCONST(62) = 4.7e-12
  RCONST(74) = 3.3e-10
  RCONST(77) = 9.5e-15
  RCONST(78) = 1.4e-14
  RCONST(80) = 1e-11
  RCONST(83) = 1.2e-11
  RCONST(84) = 6e-13
  RCONST(86) = 2.2e-12
  RCONST(87) = 3e-13
  RCONST(88) = 5e-11
  RCONST(90) = 8.7e-11
  RCONST(91) = 9e-11
  RCONST(92) = 1e-13
  RCONST(97) = 1.6e-10
  RCONST(100) = 1e-10
  RCONST(102) = 5.7e-11
  RCONST(131) = 5e-14
  RCONST(135) = 4.9e-11
  RCONST(138) = 4.1e-12
  RCONST(139) = 1.6e-12
  RCONST(140) = 1.5e-14
  RCONST(146) = 2.78e-13
  RCONST(155) = 4.5e-10
  RCONST(156) = 2.99e-11
  RCONST(161) = 2e-10
  RCONST(165) = 2.4
  RCONST(167) = 2.1e-10
  RCONST(168) = 1.5e-12
  RCONST(170) = 1.2e-12
  RCONST(184) = 1.5e-11
  RCONST(185) = 1.2e-10
  RCONST(186) = 2.09e-10
  RCONST(188) = 1.1e-15
  RCONST(189) = 3.3e-15
  RCONST(190) = 1.2e-10
  RCONST(191) = 1.2e-11
  RCONST(203) = 0
  RCONST(204) = 0
  RCONST(249) = 0.8
  RCONST(327) = 6e+09
  RCONST(354) = 3.3e+07
  RCONST(355) = 4.2e+06
  RCONST(368) = 2.4e+06
  RCONST(370) = 800000
  RCONST(385) = 3e+06
  RCONST(388) = 1830
  RCONST(486) = 1.8e+06
  RCONST(494) = 160000
  RCONST(498) = 110000
  RCONST(500) = 320
  RCONST(502) = 1e+13
  RCONST(504) = 91000
  RCONST(509) = 1.3e+09
  RCONST(511) = 280000
  RCONST(513) = 3.8e+09
  RCONST(523) = 1.3e+09
  RCONST(525) = 3.5e+08
  RCONST(527) = 280000
  RCONST(529) = 3.8e+09
  RCONST(705) = 0.8
  RCONST(783) = 6e+09
  RCONST(810) = 3.3e+07
  RCONST(811) = 4.2e+06
  RCONST(824) = 2.4e+06
  RCONST(826) = 800000
  RCONST(841) = 3e+06
  RCONST(844) = 1830
  RCONST(942) = 1.8e+06
  RCONST(950) = 160000
  RCONST(954) = 110000
  RCONST(956) = 320
  RCONST(958) = 1e+13
  RCONST(960) = 91000
  RCONST(965) = 1.3e+09
  RCONST(967) = 280000
  RCONST(969) = 3.8e+09
  RCONST(979) = 1.3e+09
  RCONST(981) = 3.5e+08
  RCONST(983) = 280000
  RCONST(985) = 3.8e+09
! END constant rate coefficients

! INLINED initializations



    DO i = 1, NVAR
    	VAR(i) = 0d0
    END DO

    DO i = 1, NFIX
    	FIX(i) = 0d0
    END DO


	! read the model setup and trajectory data
	call model_input
	


	! record the start and end times for the model run
    TSTART = time_traj(1)					! sec
    TEND   = time_traj(traj_length)			! sec
    
    ! start off the trajectory variables
    temp  = temp_traj(1)
    pres  = pres_traj(1)
    z_mbl = zmbl_traj(1)
    tidal_height = tide_traj(1)
    SUN   = sun_traj(1)


	! initialise model with H2O content and the (extra) H+ and OH- ions 
	! from the dissociation of H2O
	call h2o_init


	! record initial chemical concentrations
	cnp   = C

	
	write(6,*) 'temperature = ', temp





	!-------------------------------------
	contains	! subroutines and functions internal to the subroutine initialise
	!-------------------------------------
	
	
	subroutine model_input
	! the main subroutine for initialising the model data
	use MAP_Parameters
	use MAP_Global
	use MAP_Rates
	use Physical_Parameters
	use input
	
	character(len=20) :: word, section
	character(len=20), dimension(20) :: traj_header
	character(len=6) :: bdyfrt	
	real(kind=dp) :: adigit
	integer :: in,jn,nn,idigit
	logical :: eof

	! set the SOURCE species to 1
    FIX(indf_SOURCE) 	= 1d0 !(3.921569e-8)*CFACTOR

	open(9,file='input.dat',status="old")

	! initialise the number of chemical species at zero
	traj_species = 0

	call input_options(error_flag=0,skip_blank_lines=.true.)
	
data_read : do
	call read_line(eof,9)
	if(eof) exit data_read 	! finish reading in data
	call readu(word)
	select case(word)
	case("END") 			! end of data
		exit data_read 		! finish reading in data
	case("GAS","AERO","SET")
		section = word    		! set the section variable
	case("TRAJ") 			! trajectory section heading
		call readu(word)
		select case(word)
		case("HEAD","BODY","SET")	! select which section of trajectory data
			section = word
		case default
			call report("Unrecognised heading name: TRAJ "//trim(word),.true.)
		end select
	case default			! if line isn't section heading then read the data
		select case(section)
		case("GAS")			! gas data
			call reread(-1)
			call reada(word)	! step backwards and read 1st variable again
			call readf(adigit)	! then read the data value
			C(spec2num(word)) = adigit*CFACTOR
		case("AERO")		! aerosol data
			select case(word)
			case("DU")		! lognormal bin width
				call readi(idigit)
				do in=1,m
					call readf(adigit)
					du(in,idigit) = adigit
				end do ! in=1,m
			case("UEDGE")	! lognormal bin edges
				call readi(idigit)
				do in=1,(m+1)
					call readf(adigit)
					uedge(in,idigit) = adigit
				end do ! in=1,(m+1)
			case("R_R")		! reference radii
				call readi(idigit)
				call readf(adigit)
				r_r(idigit) = adigit
			case default	! chemical (and number) concentrations
				call reread(-1)
				call reada(word)	! step backwards and read 1st variable again
				call readi(idigit)
				do in = 1,m
					if(in.lt.10.and.idigit.lt.10) then
	    				write(bdyfrt,'(a2,I1,a2,I1)') '00',in,'00', idigit
	    			elseif(in.lt.100.and.idigit.lt.10) then 
	    				write(bdyfrt,'(a1,I2,a2,I1)') '0',in,'00', idigit
	    			elseif(in.lt.10.and.idigit.lt.100) then
	    				write(bdyfrt,'(a2,I1,a1,I2)') '00',in,'0', idigit
	    			elseif(in.lt.100.and.idigit.lt.100) then
	    				write(bdyfrt,'(a1,I2,a1,I2)') '0',in,'0', idigit
	    			elseif(in.lt.10) then
	    				write(bdyfrt,'(a2,I1,I3)') '00',in, idigit
	    			elseif(idigit.lt.10) then
	    				write(bdyfrt,'(I3,a2,I1)') in,'00', idigit
	    			elseif(in.lt.100) then
	    				write(bdyfrt,'(a1,I2,I3)') '0',in, idigit
	    			elseif(idigit.lt.100) then
	    				write(bdyfrt,'(I3,a1,I2)') in,'0', idigit
	    			else 
	    				write(bdyfrt,'(2I3)') in, idigit
	    			endif
					call readf(adigit)
					select case(word) ! treat number concs differently to chemical concs
					case("NUM")
						C(spec2num('BIN'//trim(bdyfrt)//trim(word))) = adigit
					case default
						C(spec2num('BIN'//trim(bdyfrt)//trim(word))) = adigit
					end select
				end do ! in = 1,(nitems-1)
			end select
		case("HEAD")		! trajectory headers
			traj_header(1) = trim(word)
			do in=2,nitems
				call readu(word)
				select case(word)
				case("TIME","TEMP","PRES","Z_MBL","TIDAL","SUN")
					traj_header(in) = trim(word)
				case default
					idigit = spec2num(word)
					if(idigit.gt.NSPEC.or.idigit.lt.1)then
						write(6,*) "Trajectory header "//trim(word)//" unrecognised -",&
											& " data will not be used"
					else
						if(idigit.le.NVAR)then
							write(6,*) "Warning: "//trim(word)//" is a variable species but",&
											& " is being forced using trajectory data.",&
											& " This *will* cause problems."
						endif
						traj_species = traj_species + 1
						traj_names(traj_species) = trim(word)
						traj_no(traj_species) = idigit
						traj_header(in) = trim(word)
					endif
				end select
			end do ! in=2,nitems
			nn = 1 			! set the line counter for traj data to 1
		case("BODY")		! trajectory data
			call reread(-1)
			time_traj(nn) = -99d0
			temp_traj(nn) = -99d0
			pres_traj(nn) = -99d0
			zmbl_traj(nn) = -99d0
			tide_traj(nn) = -99d0
			sun_traj(nn)  = -99d0
			do in = 1,nitems
				call readf(adigit)
				select case(traj_header(in))
				case("TIME")
					time_traj(nn) = adigit
				case("TEMP")
					temp_traj(nn) = adigit
				case("PRES")
					pres_traj(nn) = adigit
				case("Z_MBL")
					zmbl_traj(nn) = adigit *mtocm
				case("TIDAL")
					tide_traj(nn) = adigit
				case("SUN")
					sun_traj(nn) = adigit
				case default
					do jn = 1,traj_species
						if(traj_header(in).eq.traj_names(jn))then
							species_traj(jn,nn) = adigit*CFACTOR
						endif
					end do ! jn=1,traj_species
				end select
			end do ! in=1,nitems
			nn = nn + 1		! step line counter along one
		case("SET")			! model settings
			select case(word)
			case("DT")
				call readf(adigit)
				dt = adigit 
			case("TSTEPMAX")
				call readf(adigit)
				tstepmax = adigit
			case("PDFITE")
				call readu(word)
				select case(word)
				case("YES")
					pdfite = .true.
				case default
					pdfite = .false.
				end select
			!case("RTOL")
			!case("ATOL")
			end select
		case default		! no section defined
			call report("No section header for data: "//trim(word),.true.)
		end select
	end select
	end do data_read
	
	traj_length = nn-1		! record the number of trajectory data entries


	! update Z from the number variables stored in the C array
	call Update_NUMBER(C(NVARST:NVAR),C(NFIXST:NSPEC))



	close(9)

	end subroutine model_input


	subroutine h2o_init
	! Initialises the model with:
	! a) H2O content of aerosol particles
	! b) H+ and OH- content from dissociation of H2O

	use MAP_Parameters
	use MAP_Global
	use MAP_Rates
	use aerosol_microphysics
	
	real(kind=dp) :: Ctemp

	! ensure that vapour pressures and water content are up-to-date
	call aerosol_bin_sort(C(NVARST:NVAR),C(NFIXST:NSPEC),full=.true.)
	call update_rate_conversion()


	! initialise model with (extra) H+ and OH- ions equal to those 
	! expected from the dissociation of H2O
	    Ctemp = C(ind_BIN001001H2O)*conv_rate(1,1)*k_arr(1.0e-14,-6716e0,temp)
	    C(ind_BIN001001Hplu)  = C(ind_BIN001001Hplu)  + Ctemp**(0.5)
	    C(ind_BIN001001OHmin) = C(ind_BIN001001OHmin) + Ctemp**(0.5)
	    Ctemp = C(ind_BIN001002H2O)*conv_rate(1,2)*k_arr(1.0e-14,-6716e0,temp)
	    C(ind_BIN001002Hplu)  = C(ind_BIN001002Hplu)  + Ctemp**(0.5)
	    C(ind_BIN001002OHmin) = C(ind_BIN001002OHmin) + Ctemp**(0.5)


	end subroutine h2o_init

	
	

ELEMENTAL INTEGER FUNCTION spec2num ( id )

  USE MAP_Monitor, ONLY: SPC_NAMES

  CHARACTER(LEN=*), INTENT(IN) :: id
  INTEGER in

  spec2num = 0
  DO in = 1, SIZE(SPC_NAMES)
    IF (TRIM(SPC_NAMES(in)) == TRIM(id)) THEN
      spec2num = in
      EXIT
    ENDIF
  END DO

END FUNCTION spec2num






! End INLINED initializations

      
END SUBROUTINE Initialize

! End of Initialize function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE MAP_Initialize

