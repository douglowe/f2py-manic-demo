#
#  Script for testing out the python interface for MANIC
# 
import sys
sys.path.insert(0, '../library_files/MANIC_Python_Interface')
#sys.path.insert(0, '../library_files/MAP_Library_Python')

import numpy as np
import pandas as pd
import MANIC_Python as MANICP
import MAP_Monitor as MAPM
import MAP_Parameters as MAPP



### set up variables to pass to manic
# define our conversion factor
cfactor_p = 2.55e7

# get number of variables in C chemical array from model parameters
n_p_spec = MAPP.map_parameters.nspec

# initialise the C chemical array (for gas-phase species at the moment!)
#    NOTE: indexes from parameter file are for fortran, so start at 1 not 0!!!!
c_p = np.zeros(n_p_spec,order='F')
c_p[MAPP.map_parameters.ind_no2-1]		= 20.0*cfactor_p
c_p[MAPP.map_parameters.ind_o3-1]		= 20.0e3*cfactor_p
c_p[MAPP.map_parameters.ind_hcho-1]		= 300.0*cfactor_p
c_p[MAPP.map_parameters.ind_pan-1]		= 10.0*cfactor_p
c_p[MAPP.map_parameters.ind_co-1]		= 70.0e3*cfactor_p
c_p[MAPP.map_parameters.ind_hno3-1]		= 5.0*cfactor_p
c_p[MAPP.map_parameters.ind_so2-1]		= 90.0*cfactor_p
c_p[MAPP.map_parameters.ind_ch3sch3-1]	= 60.0*cfactor_p
c_p[MAPP.map_parameters.ind_h2o2-1]		= 600.0*cfactor_p
c_p[MAPP.map_parameters.ind_c2h6-1]		= 500.0*cfactor_p
c_p[MAPP.map_parameters.ind_hcl-1]		= 100.0*cfactor_p
c_p[MAPP.map_parameters.ind_ch3i-1]		= 2.0*cfactor_p
c_p[MAPP.map_parameters.ind_c3h7i-1]	= 1.0*cfactor_p
c_p[MAPP.map_parameters.ind_so2-1]		= 90.0*cfactor_p
c_p[MAPP.map_parameters.ind_no-1]		= 10.0*cfactor_p


# fixed species settings
c_p[MAPP.map_parameters.ind_nh3-1]		= 50.0*cfactor_p
c_p[MAPP.map_parameters.ind_o2-1]		= 2.1e11*cfactor_p
c_p[MAPP.map_parameters.ind_h2o-1]		= 1.25e10*cfactor_p
c_p[MAPP.map_parameters.ind_m-1]		= 1.0e12*cfactor_p
c_p[MAPP.map_parameters.ind_ch4-1]		= 1.8e6*cfactor_p
c_p[MAPP.map_parameters.ind_n2-1]		= 7.9e11*cfactor_p
c_p[MAPP.map_parameters.ind_h2-1]		= 5e5*cfactor_p
c_p[MAPP.map_parameters.ind_co2-1]		= 383e6*cfactor_p
c_p[MAPP.map_parameters.ind_source-1]   = 1e0

# aerosol species initialisation (in molecules cc^-1)
c_p[MAPP.map_parameters.ind_bin001001num-1] 		= 8.99815e+02 
c_p[MAPP.map_parameters.ind_bin001002num-1] 		= 2.49219e+00 
c_p[MAPP.map_parameters.ind_bin001001nh4plu-1] 		= 4.51659e+10 
c_p[MAPP.map_parameters.ind_bin001002nh4plu-1] 		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001001so42min-1]		= 1.00318e+10 
c_p[MAPP.map_parameters.ind_bin001002so42min-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001001hso4min-1]		= 2.30321e+10 
c_p[MAPP.map_parameters.ind_bin001002hso4min-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001001no3min-1]		= 2.07013e+09 
c_p[MAPP.map_parameters.ind_bin001002no3min-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001001naplu-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001002naplu-1]		= 2.38696e+12 
c_p[MAPP.map_parameters.ind_bin001001clmin-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001002clmin-1]		= 2.37798e+12 
c_p[MAPP.map_parameters.ind_bin001001brmin-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001002brmin-1]		= 2.02601e+09 
c_p[MAPP.map_parameters.ind_bin001001imin-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001002imin-1]		= 1.27948e+05 
c_p[MAPP.map_parameters.ind_bin001001io3min-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001002io3min-1]		= 3.39903e+05 
c_p[MAPP.map_parameters.ind_bin001001hco3min-1]		= 0.00000e+00 
c_p[MAPP.map_parameters.ind_bin001002hco3min-1]		= 6.94825e+09 

### duplicate the C array, to create the CNP_P initial array
cnp_p = np.copy(c_p)

# aerosol array size information
m_p 		= 1     # should import this from MAP_Global in future...
nd_p 		= 2
du_p 		= np.array([[2.1e0,3.6e0]],order='F')
uedge_p 	= np.array([[-9e-1,1.2e0],[-1.2e0,2.4e0]],order='F') 
r_r_p		= np.array([[8.8e-8],[1.67e-6]],order='F')



# create the trajectory information
mtocm = 100.0
trajset = pd.read_csv('traj_input.csv',sep='\s+')
time_traj_p = trajset.time
temp_traj_p = trajset.temp
pres_traj_p = trajset.pres
zmbl_traj_p = trajset.z_mbl  * mtocm
tide_traj_p = trajset.tidal
sun_traj_p  = trajset.sun

traj_length_p = len(time_traj_p)


# trajectory information
dt_p = 3600.0
tstepmax_p = 300.0
pdfite_p = 0   # this will become logical within manic?

tstart_p = 0.0
tend_p   = time_traj_p[traj_length_p-1]

# flag to say that we're in the first step (so calculate the water, H+, OH-, content for aerosol)
start_p = 1

print(c_p)

print('HNO3 in = ', c_p[MAPP.map_parameters.ind_hno3-1])

# call the solver, for one step
(c_p_new, cnp_new, tout_p) = MANICP.manic_python( n_p_spec=n_p_spec, c_p=c_p, cnp_p=cnp_p, cfactor_p=cfactor_p, 
										traj_length_p=traj_length_p, time_traj_p=time_traj_p,
										temp_traj_p=temp_traj_p, pres_traj_p=pres_traj_p, 
										zmbl_traj_p=zmbl_traj_p, tide_traj_p=tide_traj_p,
										sun_traj_p=sun_traj_p,
										dt_p=dt_p, tstepmax_p=tstepmax_p, pdfite_p=pdfite_p, 
										m_p=m_p, nd_p=nd_p, du_p=du_p, uedge_p=uedge_p, 
										r_r_p=r_r_p,
										tstart_p=tstart_p, tend_p=tend_p, start_p=start_p)

print(tstart_p, tend_p, tout_p)
print('HNO3 in? = ', c_p[MAPP.map_parameters.ind_hno3-1])
print('HNO3 out = ', c_p_new[MAPP.map_parameters.ind_hno3-1])

print(c_p_new[0:19])



"""
Parameters
----------
c_p : input rank-1 array('f') with bounds (n_p_spec)
cnp_p : input rank-1 array('f') with bounds (n_p_spec)
cfactor_p : input float
traj_length_p : input int
time_traj_p : input rank-1 array('f') with bounds (traj_length_p)
temp_traj_p : input rank-1 array('f') with bounds (traj_length_p)
pres_traj_p : input rank-1 array('f') with bounds (traj_length_p)
zmbl_traj_p : input rank-1 array('f') with bounds (traj_length_p)
tide_traj_p : input rank-1 array('f') with bounds (traj_length_p)
sun_traj_p : input rank-1 array('f') with bounds (traj_length_p)
dt_p : input float
tstepmax_p : input float
pdfite_p : input int
m_p : input int
nd_p : input int
du_p : input rank-2 array('f') with bounds (m_p,nd_p)
uedge_p : input rank-2 array('f') with bounds (m_p + 1,nd_p)
r_r_p : input rank-1 array('f') with bounds (nd_p)
tstart_p : input float
tend_p : input float

Other Parameters
----------------
n_p_spec : input int, optional
    Default: len(c_p)

Returns
-------
c_p : rank-1 array('f') with bounds (n_p_spec)
tout_p : float
"""



