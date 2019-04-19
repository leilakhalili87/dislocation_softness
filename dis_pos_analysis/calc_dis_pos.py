######################NOTICE#################3
##########It dumps the dislocation bposition every 100 steps#######################3

#!/bin/bash
import sys
import getpass

dir_path = ('../');
sys.path.insert(0,dir_path)
from dir_names import *;

import ovito.io as oio;
import ovito.modifiers as ovm;
from ovito.modifiers import DislocationAnalysisModifier
from numpy import linalg as LA
import numpy as np;
sys.path.insert(0,dir_path)
import util_funcs as uf;

t_diff=100;

pb_dumps = np.arange(t2_avg, t_stop_avg+1, t_diff); 
pb_dumps = pb_dumps.astype(int);

n_f = np.size(pb_dumps);
# gb_vals = np.zeros((n_f, 3));
gb_vals = {};

#for pbd in pb_dumps[0:2]:
for pbd in pb_dumps:
    print(pbd)
    filename = avg_xcoor_dir+avg_fwc + str(pbd) +'.out';
    node = oio.import_file (filename)
    modifier = DislocationAnalysisModifier()
    modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
    node.modifiers.append(modifier)
    node.compute()
    data= node.output.dislocations
    dis_coor=[]
    dumm=0
    for segment in data.segments:
        dis_coor.append(segment.points)
        dumm=dumm+1;


    x_min=np.mean(dis_coor[0],axis=0)[0:2]
    if dumm>1:
    	x_max=np.mean(dis_coor[1],axis=0)[0:2]
    else:
    	x_max=x_min;
    vec_bound=np.subtract(x_min,x_max)
    length_dis=LA.norm(vec_bound)
    if length_dis<15:
        bound_min_x=np.minimum(x_min[0],x_max[0]);
        bound_max_x=np.maximum(x_min[0],x_max[0]);
    else:
        bound_min_x=np.maximum(x_min[0],x_max[0]);
        bound_max_x=np.minimum(x_min[0],x_max[0]);

    bound_min_y=np.minimum(x_min[1],x_max[1]);
    bound_max_y=np.maximum(x_min[1],x_max[1]);

    gb_pos_x=[bound_min_x,bound_max_x]
    gb_pos_y=[bound_min_y,bound_max_y]
    gb_pos=np.concatenate((gb_pos_x,gb_pos_y), axis=0) #elements are  x_min, x_max,y_min,y_max respectively
    ## Calculate GB Position and Width
    
    gb_vals[str(pbd)] = gb_pos;

    node_1 = [];

import pickle as pkl;
fname = dis_pos_dir+'dis_pos.pkl';
f1 = open(fname, "wb");
pkl.dump(gb_vals, f1, protocol=2);
f1.close()

sys.path.remove(dir_path)
