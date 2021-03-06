! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Auxiliary Routines File
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
! File                 : simple_Util.f90
! Time                 : Tue May  6 12:17:10 2008
! Working directory    : /Users/lowe/work/manchester/chemistry/new-MAP/combined-MAP-code
! Equation file        : simple.kpp
! Output root filename : simple
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE MAP_Util

  USE MAP_Parameters
  IMPLICIT NONE

CONTAINS



! User INLINED Utility Functions


	subroutine aero_spc_sort(p)
	! sorts out the aerosol species stored in the C array according to
	! thier bin number
	USE MAP_Global
	! p equalivalent to C -> molec cm_air^-3
	real(kind=dp), dimension(nsol,m,nd) :: p


	   p(1,1,1) = C(ind_BIN001001NO)
	   p(2,1,1) = C(ind_BIN001001NO2)
	   p(3,1,1) = C(ind_BIN001001HONO)
	   p(4,1,1) = C(ind_BIN001001Brmin)
	   p(5,1,1) = C(ind_BIN001001HNO4)
	   p(6,1,1) = C(ind_BIN001001SO4min)
	   p(7,1,1) = C(ind_BIN001001SO2)
	   p(8,1,1) = C(ind_BIN001001DMS)
	   p(9,1,1) = C(ind_BIN001001DMSO)
	  p(10,1,1) = C(ind_BIN001001DMSO2)
	  p(11,1,1) = C(ind_BIN001001O2)
	  p(12,1,1) = C(ind_BIN001001O3)
	  p(13,1,1) = C(ind_BIN001001OH)
	  p(14,1,1) = C(ind_BIN001001HO2)
	  p(15,1,1) = C(ind_BIN001001H2O2)
	  p(16,1,1) = C(ind_BIN001001CO2)
	  p(17,1,1) = C(ind_BIN001001HCOOH)
	  p(18,1,1) = C(ind_BIN001001HCHO)
	  p(19,1,1) = C(ind_BIN001001CH3OO)
	  p(20,1,1) = C(ind_BIN001001CH3OOH)
	  p(21,1,1) = C(ind_BIN001001CH3OH)
	  p(22,1,1) = C(ind_BIN001001DOM)
	  p(23,1,1) = C(ind_BIN001001Cl)
	  p(24,1,1) = C(ind_BIN001001Cl2)
	  p(25,1,1) = C(ind_BIN001001CH3SO3H)
	  p(26,1,1) = C(ind_BIN001001HOCl)
	  p(27,1,1) = C(ind_BIN001001Br)
	  p(28,1,1) = C(ind_BIN001001Br2)
	  p(29,1,1) = C(ind_BIN001001BrCl)
	  p(30,1,1) = C(ind_BIN001001HBr)
	  p(31,1,1) = C(ind_BIN001001HOBr)
	  p(32,1,1) = C(ind_BIN001001IO)
	  p(33,1,1) = C(ind_BIN001001HOI)
	  p(34,1,1) = C(ind_BIN001001I2)
	  p(35,1,1) = C(ind_BIN001001ICl)
	  p(36,1,1) = C(ind_BIN001001IBr)
	  p(37,1,1) = C(ind_BIN001001O2min)
	  p(38,1,1) = C(ind_BIN001001OHmin)
	  p(39,1,1) = C(ind_BIN001001CH3SO2H)
	  p(40,1,1) = C(ind_BIN001001Hplu)
	  p(41,1,1) = C(ind_BIN001001NO2min)
	  p(42,1,1) = C(ind_BIN001001NO3min)
	  p(43,1,1) = C(ind_BIN001001NO4min)
	  p(44,1,1) = C(ind_BIN001001CO3min)
	  p(45,1,1) = C(ind_BIN001001HCO3min)
	  p(46,1,1) = C(ind_BIN001001HCOOmin)
	  p(47,1,1) = C(ind_BIN001001SO3min)
	  p(48,1,1) = C(ind_BIN001001SO32min)
	  p(49,1,1) = C(ind_BIN001001SO42min)
	  p(50,1,1) = C(ind_BIN001001SO5min)
	  p(51,1,1) = C(ind_BIN001001HSO3min)
	  p(52,1,1) = C(ind_BIN001001HSO4min)
	  p(53,1,1) = C(ind_BIN001001HSO5min)
	  p(54,1,1) = C(ind_BIN001001CH2xxmin)
	  p(55,1,1) = C(ind_BIN001001CH3SO2min)
	  p(56,1,1) = C(ind_BIN001001CH3SO3min)
	  p(57,1,1) = C(ind_BIN001001NH4plu)
	  p(58,1,1) = C(ind_BIN001001Br2min)
	  p(59,1,1) = C(ind_BIN001001Br2Clmin)
	  p(60,1,1) = C(ind_BIN001001BrCl2min)
	  p(61,1,1) = C(ind_BIN001001BrOmin)
	  p(62,1,1) = C(ind_BIN001001BrOHmin)
	  p(63,1,1) = C(ind_BIN001001Clmin)
	  p(64,1,1) = C(ind_BIN001001Cl2min)
	  p(65,1,1) = C(ind_BIN001001ClOmin)
	  p(66,1,1) = C(ind_BIN001001ClOHmin)
	  p(67,1,1) = C(ind_BIN001001Imin)
	  p(68,1,1) = C(ind_BIN001001IBr2min)
	  p(69,1,1) = C(ind_BIN001001ICl2min)
	  p(70,1,1) = C(ind_BIN001001IClBrmin)
	  p(71,1,1) = C(ind_BIN001001IO2min)
	  p(72,1,1) = C(ind_BIN001001IO3min)
	  p(73,1,1) = C(ind_BIN001001Naplu)
	  p(74,1,1) = C(ind_BIN001001NO3)
	  p(75,1,1) = C(ind_BIN001001ROOH)
	  p(76,1,1) = C(ind_BIN001001HNO3)
	  p(77,1,1) = C(ind_BIN001001NH3)
	  p(78,1,1) = C(ind_BIN001001HCl)
	  p(79,1,1) = C(ind_BIN001001H2SO4)
	  Z(1,1)    = C(ind_BIN001001NUM)
	   p(1,1,2) = C(ind_BIN001002NO)
	   p(2,1,2) = C(ind_BIN001002NO2)
	   p(3,1,2) = C(ind_BIN001002HONO)
	   p(4,1,2) = C(ind_BIN001002Brmin)
	   p(5,1,2) = C(ind_BIN001002HNO4)
	   p(6,1,2) = C(ind_BIN001002SO4min)
	   p(7,1,2) = C(ind_BIN001002SO2)
	   p(8,1,2) = C(ind_BIN001002DMS)
	   p(9,1,2) = C(ind_BIN001002DMSO)
	  p(10,1,2) = C(ind_BIN001002DMSO2)
	  p(11,1,2) = C(ind_BIN001002O2)
	  p(12,1,2) = C(ind_BIN001002O3)
	  p(13,1,2) = C(ind_BIN001002OH)
	  p(14,1,2) = C(ind_BIN001002HO2)
	  p(15,1,2) = C(ind_BIN001002H2O2)
	  p(16,1,2) = C(ind_BIN001002CO2)
	  p(17,1,2) = C(ind_BIN001002HCOOH)
	  p(18,1,2) = C(ind_BIN001002HCHO)
	  p(19,1,2) = C(ind_BIN001002CH3OO)
	  p(20,1,2) = C(ind_BIN001002CH3OOH)
	  p(21,1,2) = C(ind_BIN001002CH3OH)
	  p(22,1,2) = C(ind_BIN001002DOM)
	  p(23,1,2) = C(ind_BIN001002Cl)
	  p(24,1,2) = C(ind_BIN001002Cl2)
	  p(25,1,2) = C(ind_BIN001002CH3SO3H)
	  p(26,1,2) = C(ind_BIN001002HOCl)
	  p(27,1,2) = C(ind_BIN001002Br)
	  p(28,1,2) = C(ind_BIN001002Br2)
	  p(29,1,2) = C(ind_BIN001002BrCl)
	  p(30,1,2) = C(ind_BIN001002HBr)
	  p(31,1,2) = C(ind_BIN001002HOBr)
	  p(32,1,2) = C(ind_BIN001002IO)
	  p(33,1,2) = C(ind_BIN001002HOI)
	  p(34,1,2) = C(ind_BIN001002I2)
	  p(35,1,2) = C(ind_BIN001002ICl)
	  p(36,1,2) = C(ind_BIN001002IBr)
	  p(37,1,2) = C(ind_BIN001002O2min)
	  p(38,1,2) = C(ind_BIN001002OHmin)
	  p(39,1,2) = C(ind_BIN001002CH3SO2H)
	  p(40,1,2) = C(ind_BIN001002Hplu)
	  p(41,1,2) = C(ind_BIN001002NO2min)
	  p(42,1,2) = C(ind_BIN001002NO3min)
	  p(43,1,2) = C(ind_BIN001002NO4min)
	  p(44,1,2) = C(ind_BIN001002CO3min)
	  p(45,1,2) = C(ind_BIN001002HCO3min)
	  p(46,1,2) = C(ind_BIN001002HCOOmin)
	  p(47,1,2) = C(ind_BIN001002SO3min)
	  p(48,1,2) = C(ind_BIN001002SO32min)
	  p(49,1,2) = C(ind_BIN001002SO42min)
	  p(50,1,2) = C(ind_BIN001002SO5min)
	  p(51,1,2) = C(ind_BIN001002HSO3min)
	  p(52,1,2) = C(ind_BIN001002HSO4min)
	  p(53,1,2) = C(ind_BIN001002HSO5min)
	  p(54,1,2) = C(ind_BIN001002CH2xxmin)
	  p(55,1,2) = C(ind_BIN001002CH3SO2min)
	  p(56,1,2) = C(ind_BIN001002CH3SO3min)
	  p(57,1,2) = C(ind_BIN001002NH4plu)
	  p(58,1,2) = C(ind_BIN001002Br2min)
	  p(59,1,2) = C(ind_BIN001002Br2Clmin)
	  p(60,1,2) = C(ind_BIN001002BrCl2min)
	  p(61,1,2) = C(ind_BIN001002BrOmin)
	  p(62,1,2) = C(ind_BIN001002BrOHmin)
	  p(63,1,2) = C(ind_BIN001002Clmin)
	  p(64,1,2) = C(ind_BIN001002Cl2min)
	  p(65,1,2) = C(ind_BIN001002ClOmin)
	  p(66,1,2) = C(ind_BIN001002ClOHmin)
	  p(67,1,2) = C(ind_BIN001002Imin)
	  p(68,1,2) = C(ind_BIN001002IBr2min)
	  p(69,1,2) = C(ind_BIN001002ICl2min)
	  p(70,1,2) = C(ind_BIN001002IClBrmin)
	  p(71,1,2) = C(ind_BIN001002IO2min)
	  p(72,1,2) = C(ind_BIN001002IO3min)
	  p(73,1,2) = C(ind_BIN001002Naplu)
	  p(74,1,2) = C(ind_BIN001002NO3)
	  p(75,1,2) = C(ind_BIN001002ROOH)
	  p(76,1,2) = C(ind_BIN001002HNO3)
	  p(77,1,2) = C(ind_BIN001002NH3)
	  p(78,1,2) = C(ind_BIN001002HCl)
	  p(79,1,2) = C(ind_BIN001002H2SO4)
	  Z(1,2)    = C(ind_BIN001002NUM)

	end subroutine aero_spc_sort




	subroutine aero_spc_return(p)
	! inserts the aerosol species back into the C array for return
	! to the chemistry calculations
	USE MAP_Global
	! p equalivalent to C -> molec cm_air^-3
	real(kind=dp), dimension(nsol,m,nd) :: p

	  C(ind_BIN001001NO)		  =	   p(1,1,1) 
	  C(ind_BIN001001NO2)		  =	   p(2,1,1) 
	  C(ind_BIN001001HONO)		  =	   p(3,1,1) 
	  C(ind_BIN001001Brmin)	  =	   p(4,1,1)
	  C(ind_BIN001001HNO4)		  =	   p(5,1,1) 
	  C(ind_BIN001001SO4min)	  =	   p(6,1,1)
	  C(ind_BIN001001SO2)		  =	   p(7,1,1) 
	  C(ind_BIN001001DMS)		  =	   p(8,1,1) 
	  C(ind_BIN001001DMSO)		  =	   p(9,1,1) 
	  C(ind_BIN001001DMSO2)	  =    p(10,1,1) 
	  C(ind_BIN001001O2)		  =    p(11,1,1) 
	  C(ind_BIN001001O3)		  =    p(12,1,1) 
	  C(ind_BIN001001OH)		  =    p(13,1,1) 
	  C(ind_BIN001001HO2)		  =    p(14,1,1) 
	  C(ind_BIN001001H2O2)		  =    p(15,1,1) 
	  C(ind_BIN001001CO2)		  =    p(16,1,1) 
	  C(ind_BIN001001HCOOH)	  =    p(17,1,1) 
	  C(ind_BIN001001HCHO)		  =    p(18,1,1) 
	  C(ind_BIN001001CH3OO)	  =    p(19,1,1) 
	  C(ind_BIN001001CH3OOH)	  =    p(20,1,1) 
	  C(ind_BIN001001CH3OH)	  =    p(21,1,1) 
	  C(ind_BIN001001DOM)		  =    p(22,1,1) 
	  C(ind_BIN001001Cl)		  =    p(23,1,1) 
	  C(ind_BIN001001Cl2)		  =    p(24,1,1) 
	  C(ind_BIN001001CH3SO3H)    =    p(25,1,1)
	  C(ind_BIN001001HOCl)		  =    p(26,1,1) 
	  C(ind_BIN001001Br)		  =    p(27,1,1) 
	  C(ind_BIN001001Br2)		  =    p(28,1,1) 
	  C(ind_BIN001001BrCl)		  =    p(29,1,1) 
	  C(ind_BIN001001HBr)		  =    p(30,1,1) 
	  C(ind_BIN001001HOBr)		  =    p(31,1,1) 
	  C(ind_BIN001001IO)		  =    p(32,1,1) 
	  C(ind_BIN001001HOI)		  =    p(33,1,1) 
	  C(ind_BIN001001I2)		  =    p(34,1,1) 
	  C(ind_BIN001001ICl)		  =    p(35,1,1) 
	  C(ind_BIN001001IBr)		  =    p(36,1,1) 
	  C(ind_BIN001001O2min)	  =    p(37,1,1) 
	  C(ind_BIN001001OHmin)	  =    p(38,1,1) 
	  C(ind_BIN001001CH3SO2H)    =    p(39,1,1) 
	  C(ind_BIN001001Hplu)		  =    p(40,1,1) 
	  C(ind_BIN001001NO2min)	  =    p(41,1,1) 
	  C(ind_BIN001001NO3min)	  =    p(42,1,1) 
	  C(ind_BIN001001NO4min)	  =    p(43,1,1) 
	  C(ind_BIN001001CO3min)	  =    p(44,1,1) 
	  C(ind_BIN001001HCO3min)    =    p(45,1,1) 
	  C(ind_BIN001001HCOOmin)    =    p(46,1,1) 
	  C(ind_BIN001001SO3min)	  =    p(47,1,1) 
	  C(ind_BIN001001SO32min)    =    p(48,1,1) 
	  C(ind_BIN001001SO42min)    =    p(49,1,1) 
	  C(ind_BIN001001SO5min)	  =    p(50,1,1) 
	  C(ind_BIN001001HSO3min)    =    p(51,1,1) 
	  C(ind_BIN001001HSO4min)    =    p(52,1,1) 
	  C(ind_BIN001001HSO5min)    =    p(53,1,1) 
	  C(ind_BIN001001CH2xxmin)	  =    p(54,1,1) 
	  C(ind_BIN001001CH3SO2min)  =    p(55,1,1) 
	  C(ind_BIN001001CH3SO3min)  =    p(56,1,1) 
	  C(ind_BIN001001NH4plu)	  =    p(57,1,1) 
	  C(ind_BIN001001Br2min)	  =    p(58,1,1) 
	  C(ind_BIN001001Br2Clmin)   =    p(59,1,1) 
	  C(ind_BIN001001BrCl2min)   =    p(60,1,1) 
	  C(ind_BIN001001BrOmin)	  =    p(61,1,1) 
	  C(ind_BIN001001BrOHmin)    =    p(62,1,1) 
	  C(ind_BIN001001Clmin)	  =    p(63,1,1) 
	  C(ind_BIN001001Cl2min)	  =    p(64,1,1) 
	  C(ind_BIN001001ClOmin)	  =    p(65,1,1) 
	  C(ind_BIN001001ClOHmin)    =    p(66,1,1) 
	  C(ind_BIN001001Imin)		  =    p(67,1,1) 
	  C(ind_BIN001001IBr2min)    =    p(68,1,1) 
	  C(ind_BIN001001ICl2min)    =    p(69,1,1) 
	  C(ind_BIN001001IClBrmin)   =    p(70,1,1) 
	  C(ind_BIN001001IO2min)	  =    p(71,1,1) 
	  C(ind_BIN001001IO3min)	  =    p(72,1,1) 
	  C(ind_BIN001001Naplu)	  =    p(73,1,1) 
	  C(ind_BIN001001NO3)		  =    p(74,1,1)
	  C(ind_BIN001001ROOH)		  =    p(75,1,1)
	  C(ind_BIN001001HNO3)	  	  =    p(76,1,1) 
	  C(ind_BIN001001NH3)		  =    p(77,1,1) 
	  C(ind_BIN001001HCl)		  =    p(78,1,1) 
	  C(ind_BIN001001H2SO4)	  =    p(79,1,1) 
	  C(ind_BIN001001NUM)		      =    Z(1,1)
	  C(ind_BIN001002NO)		  =	   p(1,1,2) 
	  C(ind_BIN001002NO2)		  =	   p(2,1,2) 
	  C(ind_BIN001002HONO)		  =	   p(3,1,2) 
	  C(ind_BIN001002Brmin)	  =	   p(4,1,2)
	  C(ind_BIN001002HNO4)		  =	   p(5,1,2) 
	  C(ind_BIN001002SO4min)	  =	   p(6,1,2)
	  C(ind_BIN001002SO2)		  =	   p(7,1,2) 
	  C(ind_BIN001002DMS)		  =	   p(8,1,2) 
	  C(ind_BIN001002DMSO)		  =	   p(9,1,2) 
	  C(ind_BIN001002DMSO2)	  =    p(10,1,2) 
	  C(ind_BIN001002O2)		  =    p(11,1,2) 
	  C(ind_BIN001002O3)		  =    p(12,1,2) 
	  C(ind_BIN001002OH)		  =    p(13,1,2) 
	  C(ind_BIN001002HO2)		  =    p(14,1,2) 
	  C(ind_BIN001002H2O2)		  =    p(15,1,2) 
	  C(ind_BIN001002CO2)		  =    p(16,1,2) 
	  C(ind_BIN001002HCOOH)	  =    p(17,1,2) 
	  C(ind_BIN001002HCHO)		  =    p(18,1,2) 
	  C(ind_BIN001002CH3OO)	  =    p(19,1,2) 
	  C(ind_BIN001002CH3OOH)	  =    p(20,1,2) 
	  C(ind_BIN001002CH3OH)	  =    p(21,1,2) 
	  C(ind_BIN001002DOM)		  =    p(22,1,2) 
	  C(ind_BIN001002Cl)		  =    p(23,1,2) 
	  C(ind_BIN001002Cl2)		  =    p(24,1,2) 
	  C(ind_BIN001002CH3SO3H)    =    p(25,1,2)
	  C(ind_BIN001002HOCl)		  =    p(26,1,2) 
	  C(ind_BIN001002Br)		  =    p(27,1,2) 
	  C(ind_BIN001002Br2)		  =    p(28,1,2) 
	  C(ind_BIN001002BrCl)		  =    p(29,1,2) 
	  C(ind_BIN001002HBr)		  =    p(30,1,2) 
	  C(ind_BIN001002HOBr)		  =    p(31,1,2) 
	  C(ind_BIN001002IO)		  =    p(32,1,2) 
	  C(ind_BIN001002HOI)		  =    p(33,1,2) 
	  C(ind_BIN001002I2)		  =    p(34,1,2) 
	  C(ind_BIN001002ICl)		  =    p(35,1,2) 
	  C(ind_BIN001002IBr)		  =    p(36,1,2) 
	  C(ind_BIN001002O2min)	  =    p(37,1,2) 
	  C(ind_BIN001002OHmin)	  =    p(38,1,2) 
	  C(ind_BIN001002CH3SO2H)    =    p(39,1,2) 
	  C(ind_BIN001002Hplu)		  =    p(40,1,2) 
	  C(ind_BIN001002NO2min)	  =    p(41,1,2) 
	  C(ind_BIN001002NO3min)	  =    p(42,1,2) 
	  C(ind_BIN001002NO4min)	  =    p(43,1,2) 
	  C(ind_BIN001002CO3min)	  =    p(44,1,2) 
	  C(ind_BIN001002HCO3min)    =    p(45,1,2) 
	  C(ind_BIN001002HCOOmin)    =    p(46,1,2) 
	  C(ind_BIN001002SO3min)	  =    p(47,1,2) 
	  C(ind_BIN001002SO32min)    =    p(48,1,2) 
	  C(ind_BIN001002SO42min)    =    p(49,1,2) 
	  C(ind_BIN001002SO5min)	  =    p(50,1,2) 
	  C(ind_BIN001002HSO3min)    =    p(51,1,2) 
	  C(ind_BIN001002HSO4min)    =    p(52,1,2) 
	  C(ind_BIN001002HSO5min)    =    p(53,1,2) 
	  C(ind_BIN001002CH2xxmin)	  =    p(54,1,2) 
	  C(ind_BIN001002CH3SO2min)  =    p(55,1,2) 
	  C(ind_BIN001002CH3SO3min)  =    p(56,1,2) 
	  C(ind_BIN001002NH4plu)	  =    p(57,1,2) 
	  C(ind_BIN001002Br2min)	  =    p(58,1,2) 
	  C(ind_BIN001002Br2Clmin)   =    p(59,1,2) 
	  C(ind_BIN001002BrCl2min)   =    p(60,1,2) 
	  C(ind_BIN001002BrOmin)	  =    p(61,1,2) 
	  C(ind_BIN001002BrOHmin)    =    p(62,1,2) 
	  C(ind_BIN001002Clmin)	  =    p(63,1,2) 
	  C(ind_BIN001002Cl2min)	  =    p(64,1,2) 
	  C(ind_BIN001002ClOmin)	  =    p(65,1,2) 
	  C(ind_BIN001002ClOHmin)    =    p(66,1,2) 
	  C(ind_BIN001002Imin)		  =    p(67,1,2) 
	  C(ind_BIN001002IBr2min)    =    p(68,1,2) 
	  C(ind_BIN001002ICl2min)    =    p(69,1,2) 
	  C(ind_BIN001002IClBrmin)   =    p(70,1,2) 
	  C(ind_BIN001002IO2min)	  =    p(71,1,2) 
	  C(ind_BIN001002IO3min)	  =    p(72,1,2) 
	  C(ind_BIN001002Naplu)	  =    p(73,1,2) 
	  C(ind_BIN001002NO3)		  =    p(74,1,2)
	  C(ind_BIN001002ROOH)		  =    p(75,1,2)
	  C(ind_BIN001002HNO3)	  	  =    p(76,1,2) 
	  C(ind_BIN001002NH3)		  =    p(77,1,2) 
	  C(ind_BIN001002HCl)		  =    p(78,1,2) 
	  C(ind_BIN001002H2SO4)	  =    p(79,1,2) 
	  C(ind_BIN001002NUM)		      =    Z(1,2)


	end subroutine aero_spc_return




! End INLINED Utility Functions

! Utility Functions from KPP_HOME/util/util
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! UTIL - Utility functions
!   Arguments :
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! ****************************************************************
!                            
! InitSaveData - Opens the data file for writing
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE InitSaveData ()

      USE MAP_Parameters

      open(10, file='simple.dat')

      END SUBROUTINE InitSaveData

! End of InitSaveData function
! ****************************************************************

! ****************************************************************
!                            
! SaveData - Write LOOKAT species in the data file 
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE SaveData ()

      USE MAP_Global
      USE MAP_Monitor

      INTEGER i

      WRITE(10,999) (TIME-TSTART)/3600.D0,  &
                   (C(LOOKAT(i))/CFACTOR, i=1,NLOOKAT)
999   FORMAT(E24.16,100(1X,E24.16))

      END SUBROUTINE SaveData

! End of SaveData function
! ****************************************************************

! ****************************************************************
!                            
! CloseSaveData - Close the data file 
!   Parameters :                                                  
!
! ****************************************************************

      SUBROUTINE CloseSaveData ()

      USE MAP_Parameters

      CLOSE(10)

      END SUBROUTINE CloseSaveData

! End of CloseSaveData function
! ****************************************************************

! ****************************************************************
!                            
! GenerateMatlab - Generates MATLAB file to load the data file 
!   Parameters : 
!                It will have a character string to prefix each 
!                species name with.                                                 
!
! ****************************************************************

      SUBROUTINE GenerateMatlab ( PREFIX )

      USE MAP_Parameters
      USE MAP_Global
      USE MAP_Monitor

      
      CHARACTER(LEN=8) PREFIX 
      INTEGER i

      open(20, file='simple.m')
      write(20,*) 'load simple.dat;'
      write(20,990) PREFIX
990   FORMAT(A1,'c = simple;')
      write(20,*) 'clear simple;'
      write(20,991) PREFIX, PREFIX
991   FORMAT(A1,'t=',A1,'c(:,1);')
      write(20,992) PREFIX
992   FORMAT(A1,'c(:,1)=[];')

      do i=1,NLOOKAT
        write(20,993) PREFIX, SPC_NAMES(LOOKAT(i)), PREFIX, i
993     FORMAT(A1,A6,' = ',A1,'c(:,',I2,');')
      end do
      
      CLOSE(20)

      END SUBROUTINE GenerateMatlab

! End of GenerateMatlab function
! ****************************************************************


! ****************************************************************
!                            
! tag2num - convert equation tags to kpp reaction number
!   Arguments :
!      id        - string with the equation tag
!
! ****************************************************************

ELEMENTAL INTEGER FUNCTION tag2num ( id )

  USE MAP_Monitor, ONLY: EQN_TAGS

  CHARACTER(LEN=*), INTENT(IN) :: id
  INTEGER i

  tag2num = 0 ! mz_rs_20050115
  DO i = 1, SIZE(EQN_TAGS)
    IF (TRIM(EQN_TAGS(i)) == TRIM(id)) THEN
      tag2num = i ! mz_rs_20050115
      EXIT
    ENDIF
  END DO

END FUNCTION tag2num

! End of tag2num function
! ****************************************************************

! End Utility Functions from KPP_HOME/util/util
! End of UTIL function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Shuffle_user2MAP - function to copy concentrations from USER to KPP
!   Arguments :
!      V_USER    - Concentration of variable species in USER's order
!      V         - Concentrations of variable species (local)
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Shuffle_user2MAP ( V_USER, V )

! V_USER - Concentration of variable species in USER's order
  REAL(kind=dp) :: V_USER(NVAR)
! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)

  V(80) = V_USER(1)
  V(170) = V_USER(2)
  V(228) = V_USER(3)
  V(130) = V_USER(4)
  V(29) = V_USER(5)
  V(129) = V_USER(6)
  V(79) = V_USER(7)
  V(166) = V_USER(8)
  V(103) = V_USER(9)
  V(108) = V_USER(10)
  V(227) = V_USER(11)
  V(233) = V_USER(13)
  V(63) = V_USER(14)
  V(152) = V_USER(15)
  V(59) = V_USER(16)
  V(23) = V_USER(17)
  V(25) = V_USER(18)
  V(203) = V_USER(19)
  V(210) = V_USER(20)
  V(106) = V_USER(21)
  V(44) = V_USER(22)
  V(174) = V_USER(23)
  V(60) = V_USER(24)
  V(107) = V_USER(25)
  V(224) = V_USER(26)
  V(111) = V_USER(27)
  V(88) = V_USER(28)
  V(137) = V_USER(29)
  V(26) = V_USER(30)
  V(13) = V_USER(31)
  V(116) = V_USER(32)
  V(120) = V_USER(33)
  V(70) = V_USER(34)
  V(89) = V_USER(35)
  V(104) = V_USER(36)
  V(19) = V_USER(37)
  V(190) = V_USER(39)
  V(180) = V_USER(40)
  V(136) = V_USER(41)
  V(149) = V_USER(42)
  V(124) = V_USER(43)
  V(46) = V_USER(44)
  V(141) = V_USER(45)
  V(115) = V_USER(46)
  V(24) = V_USER(47)
  V(191) = V_USER(48)
  V(179) = V_USER(49)
  V(142) = V_USER(50)
  V(86) = V_USER(51)
  V(40) = V_USER(52)
  V(131) = V_USER(53)
  V(114) = V_USER(54)
  V(144) = V_USER(55)
  V(189) = V_USER(56)
  V(177) = V_USER(57)
  V(123) = V_USER(58)
  V(112) = V_USER(59)
  V(100) = V_USER(60)
  V(57) = V_USER(61)
  V(122) = V_USER(62)
  V(188) = V_USER(63)
  V(92) = V_USER(64)
  V(66) = V_USER(65)
  V(58) = V_USER(66)
  V(1) = V_USER(67)
  V(2) = V_USER(68)
  V(3) = V_USER(69)
  V(4) = V_USER(70)
  V(27) = V_USER(71)
  V(28) = V_USER(72)
  V(148) = V_USER(73)
  V(214) = V_USER(74)
  V(215) = V_USER(75)
  V(5) = V_USER(76)
  V(172) = V_USER(77)
  V(101) = V_USER(78)
  V(199) = V_USER(79)
  V(30) = V_USER(80)
  V(47) = V_USER(81)
  V(135) = V_USER(82)
  V(93) = V_USER(83)
  V(14) = V_USER(84)
  V(194) = V_USER(85)
  V(211) = V_USER(86)
  V(229) = V_USER(87)
  V(109) = V_USER(88)
  V(147) = V_USER(89)
  V(127) = V_USER(90)
  V(158) = V_USER(91)
  V(64) = V_USER(92)
  V(118) = V_USER(93)
  V(15) = V_USER(94)
  V(48) = V_USER(95)
  V(16) = V_USER(96)
  V(17) = V_USER(97)
  V(175) = V_USER(98)
  V(219) = V_USER(99)
  V(49) = V_USER(100)
  V(220) = V_USER(101)
  V(155) = V_USER(102)
  V(171) = V_USER(103)
  V(167) = V_USER(104)
  V(50) = V_USER(105)
  V(201) = V_USER(106)
  V(32) = V_USER(107)
  V(181) = V_USER(108)
  V(94) = V_USER(109)
  V(145) = V_USER(110)
  V(161) = V_USER(111)
  V(217) = V_USER(112)
  V(192) = V_USER(113)
  V(207) = V_USER(115)
  V(139) = V_USER(116)
  V(186) = V_USER(117)
  V(43) = V_USER(118)
  V(119) = V_USER(119)
  V(82) = V_USER(120)
  V(87) = V_USER(121)
  V(76) = V_USER(122)
  V(193) = V_USER(123)
  V(176) = V_USER(124)
  V(185) = V_USER(125)
  V(163) = V_USER(126)
  V(196) = V_USER(127)
  V(143) = V_USER(128)
  V(133) = V_USER(129)
  V(61) = V_USER(130)
  V(41) = V_USER(131)
  V(78) = V_USER(132)
  V(31) = V_USER(133)
  V(6) = V_USER(134)
  V(206) = V_USER(135)
  V(226) = V_USER(136)
  V(67) = V_USER(137)
  V(68) = V_USER(138)
  V(95) = V_USER(139)
  V(102) = V_USER(140)
  V(213) = V_USER(141)
  V(195) = V_USER(142)
  V(96) = V_USER(143)
  V(69) = V_USER(144)
  V(164) = V_USER(145)
  V(33) = V_USER(146)
  V(34) = V_USER(147)
  V(71) = V_USER(148)
  V(125) = V_USER(149)
  V(7) = V_USER(150)
  V(154) = V_USER(151)
  V(84) = V_USER(152)
  V(231) = V_USER(153)
  V(150) = V_USER(154)
  V(55) = V_USER(155)
  V(8) = V_USER(156)
  V(9) = V_USER(157)
  V(216) = V_USER(158)
  V(113) = V_USER(159)
  V(168) = V_USER(160)
  V(35) = V_USER(161)
  V(51) = V_USER(162)
  V(134) = V_USER(163)
  V(90) = V_USER(164)
  V(18) = V_USER(165)
  V(232) = V_USER(166)
  V(209) = V_USER(167)
  V(222) = V_USER(168)
  V(105) = V_USER(169)
  V(146) = V_USER(170)
  V(128) = V_USER(171)
  V(157) = V_USER(172)
  V(65) = V_USER(173)
  V(117) = V_USER(174)
  V(20) = V_USER(175)
  V(52) = V_USER(176)
  V(21) = V_USER(177)
  V(22) = V_USER(178)
  V(169) = V_USER(179)
  V(205) = V_USER(180)
  V(53) = V_USER(181)
  V(198) = V_USER(182)
  V(156) = V_USER(183)
  V(225) = V_USER(184)
  V(187) = V_USER(185)
  V(54) = V_USER(186)
  V(200) = V_USER(187)
  V(37) = V_USER(188)
  V(178) = V_USER(189)
  V(97) = V_USER(190)
  V(160) = V_USER(191)
  V(138) = V_USER(192)
  V(223) = V_USER(193)
  V(208) = V_USER(194)
  V(221) = V_USER(196)
  V(140) = V_USER(197)
  V(197) = V_USER(198)
  V(45) = V_USER(199)
  V(121) = V_USER(200)
  V(83) = V_USER(201)
  V(91) = V_USER(202)
  V(77) = V_USER(203)
  V(184) = V_USER(204)
  V(173) = V_USER(205)
  V(202) = V_USER(206)
  V(162) = V_USER(207)
  V(204) = V_USER(208)
  V(159) = V_USER(209)
  V(132) = V_USER(210)
  V(62) = V_USER(211)
  V(42) = V_USER(212)
  V(81) = V_USER(213)
  V(36) = V_USER(214)
  V(10) = V_USER(215)
  V(212) = V_USER(216)
  V(183) = V_USER(217)
  V(72) = V_USER(218)
  V(73) = V_USER(219)
  V(98) = V_USER(220)
  V(110) = V_USER(221)
  V(218) = V_USER(222)
  V(182) = V_USER(223)
  V(99) = V_USER(224)
  V(74) = V_USER(225)
  V(165) = V_USER(226)
  V(38) = V_USER(227)
  V(39) = V_USER(228)
  V(75) = V_USER(229)
  V(126) = V_USER(230)
  V(11) = V_USER(231)
  V(153) = V_USER(232)
  V(85) = V_USER(233)
      
END SUBROUTINE Shuffle_user2MAP

! End of Shuffle_user2MAP function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Shuffle_MAP2user - function to restore concentrations from KPP to USER
!   Arguments :
!      V         - Concentrations of variable species (local)
!      V_USER    - Concentration of variable species in USER's order
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Shuffle_MAP2user ( V, V_USER )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! V_USER - Concentration of variable species in USER's order
  REAL(kind=dp) :: V_USER(NVAR)

  V_USER(1) = V(80)
  V_USER(2) = V(170)
  V_USER(3) = V(228)
  V_USER(4) = V(130)
  V_USER(5) = V(29)
  V_USER(6) = V(129)
  V_USER(7) = V(79)
  V_USER(8) = V(166)
  V_USER(9) = V(103)
  V_USER(10) = V(108)
  V_USER(11) = V(227)
  V_USER(13) = V(233)
  V_USER(14) = V(63)
  V_USER(15) = V(152)
  V_USER(16) = V(59)
  V_USER(17) = V(23)
  V_USER(18) = V(25)
  V_USER(19) = V(203)
  V_USER(20) = V(210)
  V_USER(21) = V(106)
  V_USER(22) = V(44)
  V_USER(23) = V(174)
  V_USER(24) = V(60)
  V_USER(25) = V(107)
  V_USER(26) = V(224)
  V_USER(27) = V(111)
  V_USER(28) = V(88)
  V_USER(29) = V(137)
  V_USER(30) = V(26)
  V_USER(31) = V(13)
  V_USER(32) = V(116)
  V_USER(33) = V(120)
  V_USER(34) = V(70)
  V_USER(35) = V(89)
  V_USER(36) = V(104)
  V_USER(37) = V(19)
  V_USER(39) = V(190)
  V_USER(40) = V(180)
  V_USER(41) = V(136)
  V_USER(42) = V(149)
  V_USER(43) = V(124)
  V_USER(44) = V(46)
  V_USER(45) = V(141)
  V_USER(46) = V(115)
  V_USER(47) = V(24)
  V_USER(48) = V(191)
  V_USER(49) = V(179)
  V_USER(50) = V(142)
  V_USER(51) = V(86)
  V_USER(52) = V(40)
  V_USER(53) = V(131)
  V_USER(54) = V(114)
  V_USER(55) = V(144)
  V_USER(56) = V(189)
  V_USER(57) = V(177)
  V_USER(58) = V(123)
  V_USER(59) = V(112)
  V_USER(60) = V(100)
  V_USER(61) = V(57)
  V_USER(62) = V(122)
  V_USER(63) = V(188)
  V_USER(64) = V(92)
  V_USER(65) = V(66)
  V_USER(66) = V(58)
  V_USER(67) = V(1)
  V_USER(68) = V(2)
  V_USER(69) = V(3)
  V_USER(70) = V(4)
  V_USER(71) = V(27)
  V_USER(72) = V(28)
  V_USER(73) = V(148)
  V_USER(74) = V(214)
  V_USER(75) = V(215)
  V_USER(76) = V(5)
  V_USER(77) = V(172)
  V_USER(78) = V(101)
  V_USER(79) = V(199)
  V_USER(80) = V(30)
  V_USER(81) = V(47)
  V_USER(82) = V(135)
  V_USER(83) = V(93)
  V_USER(84) = V(14)
  V_USER(85) = V(194)
  V_USER(86) = V(211)
  V_USER(87) = V(229)
  V_USER(88) = V(109)
  V_USER(89) = V(147)
  V_USER(90) = V(127)
  V_USER(91) = V(158)
  V_USER(92) = V(64)
  V_USER(93) = V(118)
  V_USER(94) = V(15)
  V_USER(95) = V(48)
  V_USER(96) = V(16)
  V_USER(97) = V(17)
  V_USER(98) = V(175)
  V_USER(99) = V(219)
  V_USER(100) = V(49)
  V_USER(101) = V(220)
  V_USER(102) = V(155)
  V_USER(103) = V(171)
  V_USER(104) = V(167)
  V_USER(105) = V(50)
  V_USER(106) = V(201)
  V_USER(107) = V(32)
  V_USER(108) = V(181)
  V_USER(109) = V(94)
  V_USER(110) = V(145)
  V_USER(111) = V(161)
  V_USER(112) = V(217)
  V_USER(113) = V(192)
  V_USER(115) = V(207)
  V_USER(116) = V(139)
  V_USER(117) = V(186)
  V_USER(118) = V(43)
  V_USER(119) = V(119)
  V_USER(120) = V(82)
  V_USER(121) = V(87)
  V_USER(122) = V(76)
  V_USER(123) = V(193)
  V_USER(124) = V(176)
  V_USER(125) = V(185)
  V_USER(126) = V(163)
  V_USER(127) = V(196)
  V_USER(128) = V(143)
  V_USER(129) = V(133)
  V_USER(130) = V(61)
  V_USER(131) = V(41)
  V_USER(132) = V(78)
  V_USER(133) = V(31)
  V_USER(134) = V(6)
  V_USER(135) = V(206)
  V_USER(136) = V(226)
  V_USER(137) = V(67)
  V_USER(138) = V(68)
  V_USER(139) = V(95)
  V_USER(140) = V(102)
  V_USER(141) = V(213)
  V_USER(142) = V(195)
  V_USER(143) = V(96)
  V_USER(144) = V(69)
  V_USER(145) = V(164)
  V_USER(146) = V(33)
  V_USER(147) = V(34)
  V_USER(148) = V(71)
  V_USER(149) = V(125)
  V_USER(150) = V(7)
  V_USER(151) = V(154)
  V_USER(152) = V(84)
  V_USER(153) = V(231)
  V_USER(154) = V(150)
  V_USER(155) = V(55)
  V_USER(156) = V(8)
  V_USER(157) = V(9)
  V_USER(158) = V(216)
  V_USER(159) = V(113)
  V_USER(160) = V(168)
  V_USER(161) = V(35)
  V_USER(162) = V(51)
  V_USER(163) = V(134)
  V_USER(164) = V(90)
  V_USER(165) = V(18)
  V_USER(166) = V(232)
  V_USER(167) = V(209)
  V_USER(168) = V(222)
  V_USER(169) = V(105)
  V_USER(170) = V(146)
  V_USER(171) = V(128)
  V_USER(172) = V(157)
  V_USER(173) = V(65)
  V_USER(174) = V(117)
  V_USER(175) = V(20)
  V_USER(176) = V(52)
  V_USER(177) = V(21)
  V_USER(178) = V(22)
  V_USER(179) = V(169)
  V_USER(180) = V(205)
  V_USER(181) = V(53)
  V_USER(182) = V(198)
  V_USER(183) = V(156)
  V_USER(184) = V(225)
  V_USER(185) = V(187)
  V_USER(186) = V(54)
  V_USER(187) = V(200)
  V_USER(188) = V(37)
  V_USER(189) = V(178)
  V_USER(190) = V(97)
  V_USER(191) = V(160)
  V_USER(192) = V(138)
  V_USER(193) = V(223)
  V_USER(194) = V(208)
  V_USER(196) = V(221)
  V_USER(197) = V(140)
  V_USER(198) = V(197)
  V_USER(199) = V(45)
  V_USER(200) = V(121)
  V_USER(201) = V(83)
  V_USER(202) = V(91)
  V_USER(203) = V(77)
  V_USER(204) = V(184)
  V_USER(205) = V(173)
  V_USER(206) = V(202)
  V_USER(207) = V(162)
  V_USER(208) = V(204)
  V_USER(209) = V(159)
  V_USER(210) = V(132)
  V_USER(211) = V(62)
  V_USER(212) = V(42)
  V_USER(213) = V(81)
  V_USER(214) = V(36)
  V_USER(215) = V(10)
  V_USER(216) = V(212)
  V_USER(217) = V(183)
  V_USER(218) = V(72)
  V_USER(219) = V(73)
  V_USER(220) = V(98)
  V_USER(221) = V(110)
  V_USER(222) = V(218)
  V_USER(223) = V(182)
  V_USER(224) = V(99)
  V_USER(225) = V(74)
  V_USER(226) = V(165)
  V_USER(227) = V(38)
  V_USER(228) = V(39)
  V_USER(229) = V(75)
  V_USER(230) = V(126)
  V_USER(231) = V(11)
  V_USER(232) = V(153)
  V_USER(233) = V(85)
      
END SUBROUTINE Shuffle_MAP2user

! End of Shuffle_MAP2user function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! GetMass - compute total mass of selected atoms
!   Arguments :
!      CL        - Concentration of all species (local)
!      Mass      - value of mass balance
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE GetMass ( CL, Mass )

! CL - Concentration of all species (local)
  REAL(kind=dp) :: CL(NSPEC)
! Mass - value of mass balance
  REAL(kind=dp) :: Mass(1)

      
END SUBROUTINE GetMass

! End of GetMass function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE MAP_Util

