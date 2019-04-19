#!/bin/bash
import sys
import getpass

dir_path = ('/home/leila/Leila_sndhard/codes/gb_mobility_test/python_analysis_codes/github/softness_analysis-master');
sys.path.insert(0,dir_path)

from dir_names import *;

import numpy as np;
import ovito.io as oio;
sys.path.insert(0,dir_path)
import util_funcs as uf;
import matplotlib.pyplot as plt;

pb_dumps = np.arange(t2_avg, t_stop_avg+1, 10000); 
pb_dumps = pb_dumps.astype(int);

for pbd in pb_dumps:
# for pbd in pb_dumps[0:1]:
	print(pbd);
	# filename = adj_xcoor_dir+xcadj_fwc + str(pbd) + '.out';
	filename = avg_xcoor_dir+avg_fwc + str(pbd) +'.out';
	node_1 = oio.import_file (filename)
	data = node_1.compute();

	pos_property = data.particle_properties['Position']
	pos = pos_property.array

	## Calculate GB Position and Width
	# op_prop = data.particle_properties['f_gb[1]'];
	op_prop = data.particle_properties['f_gb1'];
	op = op_prop[...];
	x_pos = pos[:,0];

	x_min = np.min(x_pos);
	x_max = np.max(x_pos);

	[x, y] = uf.avg_op_vals(x_pos, op, lat_a, sc);
	y_min = np.min(y); y_max = np.max(y);
	y_wma = uf.weighted_mov_avg(y, rollingParam);
	if rollingParam>0:
		x_wma = x[rollingParam:-rollingParam];
	else:
		x_wma = np.copy(x);


	[gb_pos, coeffs] = uf.calc_gb_pos(x_wma, y_wma);

	gb_min = gb_pos[0]; gb_mid = gb_pos[1]; gb_max = gb_pos[2];
	# print(gb_min);

	## GB buffer region
	gb_width = 2*(gb_mid - gb_min);
	gb_low = gb_mid - gb_width; gb_high = gb_mid + gb_width;

	xpos_low = x_min + 2*cutoff;
	xpos_high = x_max - 2*cutoff;
	xcut_low = np.max([xpos_low, gb_low]);
	xcut_high = np.min([xpos_high, gb_high]);

	y_tanh = uf.tanhyp_func(x, coeffs[0], coeffs[1], coeffs[2], coeffs[3]);

	plt.figure(figsize=(10,10))
	plt.scatter(x, y, s=10, c='b');
	plt.scatter(x_wma, y_wma, s=100, c='r', alpha=0.5);
	plt.plot(x, y_tanh, 'b', linewidth=2);

	# plt.plot([gb_low, gb_low], [y_min, y_max], 'r', linewidth=2, linestyle='dashed');
	# plt.plot([gb_high, gb_high], [y_min, y_max], 'r', linewidth=2, linestyle='dashed');

	# plt.plot([xpos_low, xpos_low], [y_min, y_max], 'b', linewidth=2, linestyle='dashed');
	# plt.plot([xpos_high, xpos_high], [y_min, y_max], 'b', linewidth=2, linestyle='dashed');

	plt.plot([xcut_low, xcut_low], [y_min, y_max], 'b', linewidth=2, linestyle='dashed');
	plt.plot([xcut_high, xcut_high], [y_min, y_max], 'b', linewidth=2, linestyle='dashed');
	plt.plot([gb_mid, gb_mid], [y_min, y_max], 'g', linewidth=2, linestyle='dashed');

	plt.plot([x_min, x_min], [y_min, y_max], 'k', linewidth=3);
	plt.plot([x_max, x_max], [y_min, y_max], 'k', linewidth=3);
	# plt.grid(True);

	fig_name = gb_pos_dir + 'gb_pos_'+str(int(pbd))+'.png';
	plt.savefig(fig_name, dpi=200);
	plt.xlim([x_min, x_max])
	plt.close();

sys.path.remove(dir_path)