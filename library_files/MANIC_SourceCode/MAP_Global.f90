! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Global Data Module File
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
! File                 : simple_Global.f90
! Time                 : Tue May  6 12:17:10 2008
! Working directory    : /Users/lowe/work/manchester/chemistry/new-MAP/combined-MAP-code
! Equation file        : simple.kpp
! Output root filename : simple
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE MAP_Global

  USE MAP_Parameters, ONLY: dp, NSPEC, NVAR, NFIX, NREACT
  PUBLIC
  SAVE


! Declaration of global variables

! C - Concentration of all species
  REAL(kind=dp) :: C(NSPEC)
! VAR - Concentrations of variable species (global)
  REAL(kind=dp) :: VAR(NVAR)
! FIX - Concentrations of fixed species (global)
  REAL(kind=dp) :: FIX(NFIX)
! VAR, FIX are chunks of array C
      EQUIVALENCE( C(1),VAR(1) )
      EQUIVALENCE( C(234),FIX(1) )  ! NAG doesn't deal well with equivalence - *don't* use FIX!!!
! RCONST - Rate constants (global)
  REAL(kind=dp) :: RCONST(NREACT)
! TIME - Current integration time
  REAL(kind=dp) :: TIME
! SUN - Sunlight intensity between [0,1]
  REAL(kind=dp) :: SUN
! TEMP - Temperature
  REAL(kind=dp) :: TEMP
! RTOLS - (scalar) Relative tolerance
  REAL(kind=dp) :: RTOLS
! TSTART - Integration start time
  REAL(kind=dp) :: TSTART
! TEND - Integration end time
  REAL(kind=dp) :: TEND
! DT - Integration step
  REAL(kind=dp) :: DT
! ATOL - Absolute tolerance
  REAL(kind=dp) :: ATOL(NVAR)
! RTOL - Relative tolerance
  REAL(kind=dp) :: RTOL(NVAR)
! STEPMIN - Lower bound for integration step
  REAL(kind=dp) :: STEPMIN
! STEPMAX - Upper bound for integration step
  REAL(kind=dp) :: STEPMAX
! CFACTOR - Conversion factor for concentration units
  REAL(kind=dp) :: CFACTOR

! INLINED global variable declarations



  !!! MANIC parameters
! nsol - number of soluble species within the system
  integer, parameter :: nsol = 79

  ! m - number of size bins in aerosol model
    integer, parameter :: m = 1
  ! nd - number of divisions for each size bin
    integer, parameter :: nd = 2


  !!! MANIC variable arrays
	real(kind=dp) :: pres, z_mbl, tidal_height

  ! trajectory data arrays
  ! time_traj = time (seconds)
  ! temp_traj = temperature (K)
  ! pres_traj = pressure (mbar)
  ! zmbl_traj = boundary layer height (cm)
  ! tide_traj = tidal height (0-1)
  ! sun_traj  = solar strength (0-1)
	real(kind=dp), dimension(100000) :: time_traj, temp_traj, pres_traj,&
				& zmbl_traj, tide_traj, sun_traj
  ! species_traj = concentrations of (fixed) species (molecules cm_air^-3)
  ! Note: upper limit of 10 species which can be treated this way...
    real(kind=dp), dimension(10,100000) :: species_traj
  ! traj_length  = number of data points in trajectory files
  ! traj_species = number of species stored in species_traj
    integer :: traj_length, traj_species
  ! traj_names = names of the species stored in species_traj
  	character(len=20), dimension(10) :: traj_names
  ! traj_no = positions in C array of species stored in species_traj
    integer, dimension(10) :: traj_no

    
  ! cnp - initial chemical concentrations (molecules cm_air^-3)
    real(kind=dp) :: cnp(NSPEC)

  ! size variables for aerosol
  ! u     - logarithmic particle size
  ! du    - width of logarithmic size bin
  ! ulim  - limits of the logarithmic size scale
    real(kind=dp), dimension(m,nd) :: u
    real(kind=dp), dimension(m,nd) :: du 
    real(kind=dp), dimension(2,nd) :: ulim
  ! uedge - edges of logarithmic size bins
    real(kind=dp), dimension(m+1,nd) :: uedge
  ! r_r - reference radii for the size arrays (m)
    real(kind=dp), dimension(nd) :: r_r
  ! radius - particle radius for that bin (u=ln(radius/r_r) (m)
    real(kind=dp), dimension(m,nd) :: radius

  ! tstep - size of timestep external to ODE intergrator (s)
    real(kind=dp) :: tstep
  ! tout - time of next data output (s)
    real(kind=dp) :: tout

  ! Z     - current particle number in each size bin (# cm^-3)
    real(kind=dp), dimension(m,nd) :: Z

  ! conversion ratio for liquid-phase reaction rates
  ! (from moles l_aq^-1 to molecules cm_air^-3)
    real(kind=dp), dimension(m,nd) :: conv_rate

	! acid - integer array storing data on the last pH state of the aerosol
	!	1 = acid
	!	0 = alkali
	!	-1 = undetermined (chose at first step)
	integer, dimension(m,nd) :: acid

      
  ! act_coeff - activity coefficients for non-ideal 
  !             ion dissociation (dimensionless)
    real(kind=dp), dimension(5,m,nd) :: act_coeff
  ! vap_press - HNO3, NH3 and HCl concentrations from PD-FiTE
    real(kind=dp), dimension(3,m,nd) :: vap_press
  ! pdfite - logical switch choosing if we use PD-FiTE vapour pressures or not
    logical :: pdfite
  
    
  ! v_d - aerosol deposition velocity
  ! v_e - aerosol emission velocity
    real(kind=dp), dimension(m,nd) :: v_d, v_e
  
  ! kelvin - Kelvin effect parameter
    real(kind=dp), dimension(m,nd) :: kelvin
  
  ! water_content - equivalent to total aerosol volume (cm_aq^3 cm_air^-3)
    real(kind=dp), dimension(m,nd) :: water_content
  
  ! wrad - wet radius (m)
    real(kind=dp), dimension(m,nd) :: wrad
  
  ! tmass - total wet mass of one particle (molec)
    real(kind=dp), dimension(m,nd) :: tmass

  ! rh - relative humidity (0-1)
	real(kind=dp) :: rh

  ! tstepmax - maximum timestep allowed
    real(kind=dp) :: tstepmax




! INLINED global variable declarations


END MODULE MAP_Global

