################################################################################3
##############CHECK THIS PARAMETERS EVERY TIME#####################################
#f_T_P ----> The path of the file you are intrested in
#lat_a ----->lattice parameter
# dis_width ---->layers of atoms considered for considering a box round the dislocation
#dis_scale ----->scales the dis_width
# y_prop ------> which parameter you wanna choose : dmin, phop or gmean
# dmin_cut ----> the dmin threshold
#t_yprop_train -----> how often to choose the train data set
#t_yprop_test  -----> how often choose test set

import numpy as np;

import getpass
# getpass.getuser()


################################################################################
# The time step is 1 femtosecond.
# The dump files are written every 10 steps (10 fs).
t_diff = 10; # units: femtoseconds.
# The total time of the simulation.
t_init = 0; t_last = 100000;
# Total number of dump files
num_dump = int((t_last - t_init)/10) + 1;
# The index of the dump-file is the time-step.
# That is, zdump.n.out corresponds to $n^{th}$ femtosecond.

### Directory names


f_T_P			='T100_p100/'

fname0 			= '/home/leila/Leila_sndhard/codes/';
fname1 			= 'dislocation/disl_data/';
fname2 			= 'Lammps_Ni/';
fname3 			= 'Analysis_Data/';

gb_id 			= 'shear/';
#gb_id 			= 'p119/';

fstr1 			= fname0+fname1+fname2+  f_T_P + gb_id;
sim_dump_dir 	= fstr1;
sim_fwc			= 'zdump.';


################################################################################
# Step 2: 	Average atom-positions in dump-files.
#			This is done to reduce  the  highest  frequency  atomic  vibrations
#			The rest of the calculations depend on this parameter
#			(so we create a 't2avg_/' folder)
### Parameters for averaging atom positions
fstr2 			= fname0+fname1+fname2+f_T_P+fname3;
t_avg = 400;
t2_avg = int(t_avg/2.0);
t_stop_avg = t_last - t2_avg;
### Directory names
t2avg_dir 		= fstr2+'t2avg_'+str(t2_avg)+'/';
avg_xcoor_dir	= t2avg_dir+'avg_xcoor/'
avg_fwc			= 'avg_xcoor.';

################################################################################
# Step 3: 	Calculate some properties of Grain Boundary positions.
#			1) Important for calculating mobility
#			2) To determing the last time-step for further analysis
### Parameters

lat_a = 3.62; 
dis_width=1.5*lat_a;
#rollingParam = 10; sc = 0.4;
### Directory names
dis_pos_dir 		= t2avg_dir+'dis_position/';


################################################################################
# Step 4: 	Calculate Parrinello-Behler Fingerprint
### Parameters
__sca__ = lat_a / 4.05   # scales the lattice parameter with aluminum
n_rfunc = 18; n_afunc = 54;
n_vec = n_rfunc + n_afunc;

alpha = np.arange(1, 19, 1) # alpha = [1,2,3...18]
miu1 = np.arange(2.6, 3.9, 0.1);
miu2 = np.arange(4.0, 5.0, 0.2)
miu =  np.dot(__sca__,np.concatenate((miu1,miu2))) # miu = [2.6, 2.7...3.8, 4.0...4.6, 4.8]

cutoff = np.dot(__sca__, 5.0) # d0
sigma = np.dot(__sca__, 0.25) # sigma0 = 0.25

beta = np.arange(1, 55, 1) # b = [1,2,3...53, 54]
zeta = np.arange(0.1, 2.8, 0.1)
zeta = np.concatenate((zeta,zeta))
zeta.sort() # zeta = [0.1, 0.1, 0.2, 0.2, ...2.7, 2.7]
eta = np.dot(1/np.sqrt(__sca__), np.exp(-np.sqrt(beta)/2)) # eta = [0.61, 0.49, 0.42, ...0.026, 0.025]

lamda = np.ones(54)
neg_inx = np.arange(1, 54, 2)
lamda[neg_inx] = -1 # lamda = [1, -1, 1, ...1, -1]


lamda = np.reshape(lamda, (1, 54));
eta = np.reshape(eta, (1, 54));
zeta = np.reshape(zeta, (1, 54));

params = {};
params['miu'] = miu;
params['lamda'] = lamda;
params['eta'] = eta;
params['zeta'] = zeta;
params['sigma'] = sigma;
params['cutoff'] = cutoff;


## GB buffer region
dis_scale = .9;

t_last_pb = 50000;


### Directory names
pb_dir_soap			= t2avg_dir+'pb_pkls/soap_fp/';
pb_dir_str			= t2avg_dir+'pb_pkls/str_fp/'
pbpkls_fwc		= 'pb_fp_';

################################################################################
# Step 5: 	Calculate D2_min
### Parameters
t_r = 600
# # t2_r: Is equal to t_r/2
t2_r = int(t_r/2);
# Number of dump-files to consider for each moving-time window
ndump_mw = int(t2_r/t_diff);
# First and Last time step for p_hop or dmin
t_yprop_init = t2_avg + t2_r;
t_yprop_stop = t_last - t2_avg - t2_r;
### Directory names
t2r_dir1		= t2avg_dir + 't2r_' + str(int(t2_r))+ '/';
dmin_dir		= t2r_dir1 + 'dmin/dmin_pkls/';
dmin_fwc		= 'dmin_';
dmin_dump_dir	= t2r_dir1 + 'dmin/dmin_dumps/';
dmin_dump_fwc	= 'dmin_xcoor';
dmin_svm_dir 	= t2r_dir1 + 'dmin/dmin_svm/';
################################################################################
# Step 5 (v1): 	Calculate P_hop
### Parameters
# t_r: The time duration (length) of the moving time window. (alredy defined fyr D2_min)
### Directory names
phop_dir		= t2r_dir1 + 'phop/phop_pkls/';
phop_fwc		= 'phop_';
phop_dump_dir	= t2r_dir1 + 'phop/phop_dumps/';
phop_dump_fwc	= 'phop_xcoor.';
phop_svm_dir 	= t2r_dir1 + 'phop/phop_svm/';
###############################################################################
# Step 5 (v1): 	Calculate P_hop
### Parameters
# t_r: The time duration (length) of the moving time window. (alredy defined fyr D2_min)
### Directory names
gmean_dir		= t2r_dir1 + 'gmean/gmean_pkls/';
gmean_fwc		= 'gmean_';
gmean_dump_dir	= t2r_dir1 + 'gmean/gmean_dumps/';
gmean_dump_fwc	= 'gmean_xcoor.';
gmean_svm_dir 	= t2r_dir1 + 'gmean/gmean_svm/';
###############################################################################
# For softness-SVM Training
frac_satms = 0.25;

# y_prop = 'phop'; phop_cut = 0.3;
y_prop = 'dmin'; 

dmin_cut = .03;
dmin_cut_crystal=dmin_cut/2;
t_yprop_train = 400
t_yprop_test = 100;
y_soft = 1; y_nsoft = -1;



