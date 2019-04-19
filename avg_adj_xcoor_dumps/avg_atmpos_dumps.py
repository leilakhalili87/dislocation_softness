#!/bin/bash
import sys
import getpass

#dir_path = ('/Users/' + getpass.getuser() + 
#	'/NCSU_Research/Repos/GB_Mobility_Analysis/Python_Analysis_Codes/')
#sys.path.append(dir_path)
import sys
dir_path = ('../');
sys.path.insert(0,dir_path)

from dir_names import *;

import ovito.io as oio;
import ovito.modifiers as ovm;
import ovito.data as ovd;
# from ovito.data import DataCollection, SimulationCell

import numpy as np;
sys.path.insert(0,dir_path)
dir_path = ('../');
import util_funcs as uf;
import pickle as pkl;

##########################################################################################
### A list will contain all the displacements between successive frames
# First time window
# t0 = 1000;

t0 = t2_avg;
# t0 = t_stop;
ndump_avgw = int(t2_avg/t_diff)

pb_dumps = np.arange(t2_avg, t_stop_avg+1, t_diff);
# pb_dumps = np.arange(35600, t_stop_avg+1, t_diff);
pb_dumps = pb_dumps.astype(int);

for tstep in pb_dumps:
# for tstep in pb_dumps[0:10]:
	print(tstep);
	input_dump_wc = sim_dump_dir+sim_fwc;
	filename1 = input_dump_wc+str(tstep)+'.out';
	pipeline = oio.import_file (filename1, sort_particles=True);
	data = pipeline.compute();

	positions = data.particles['Position']; 
	pos_pbd = positions[...];
	num_atoms = np.shape(pos_pbd)[0];

	A_window = np.linspace(tstep, tstep-t2_avg, ndump_avgw+1);
	A_window = A_window.astype(int); A_window = np.delete(A_window, 0);

	B_window = np.linspace(tstep, tstep+t2_avg, ndump_avgw+1);
	B_window = B_window.astype(int); B_window = np.delete(B_window, 0);

	# Calculate displacement vectors
	disp_awin = uf.calc_win_disps(num_atoms, ndump_avgw, A_window, input_dump_wc, pipeline);
	disp_bwin = uf.calc_win_disps(num_atoms, ndump_avgw, B_window, input_dump_wc, pipeline);

	avg_pos = 0*pos_pbd; tot_ndw = 2*ndump_avgw+1;
	for ct1 in range(3):
		avg_disp = ((np.sum(disp_awin[ct1], axis=1) 
		            + np.sum(disp_bwin[ct1], axis=1))/(tot_ndw));
		avg_pos[:,ct1] = (pos_pbd[:,ct1] - avg_disp);

	data.particles.create_property('Position', data = avg_pos);

	##########################
	# dumping the adjusted files
	###########################
	dname = avg_xcoor_dir+avg_fwc + str(tstep) +'.out';
	oio.export_file(data, dname, format = "lammps_dump",
		columns=["Particle Identifier",
				 "Particle Type",
				 "Position.X",
				 "Position.Y",
				 "Position.Z",
				 "c_csym" ], frame=tstep)
	pipeline = [];

sys.path.remove(dir_path)
