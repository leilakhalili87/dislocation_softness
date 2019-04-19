import numpy as np
import sys
import ovito.io as oio;
import ovito.modifiers as ovm;
from ovito.modifiers import DislocationAnalysisModifier
from numpy import linalg as LA
import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt;
dir_path = ('../');
sys.path.insert(0,dir_path)
from dir_names import *;


t_diff=500;
pb_dumps = np.arange(t2_avg, t_stop_avg, t_diff);
# input_dump_wc='/home/leila/Leila_sndhard/codes/dislocation/disl_data/Lammps_Ni/T100/shear/';
i=0
num_dump=np.shape(pb_dumps)[0];
x_min=np.zeros(num_dump);
x_max=np.zeros(num_dump);
time=np.zeros(num_dump);
diff_snd=np.zeros(num_dump-1);
diff_first=np.zeros(num_dump-1);
numer_1=0
numer_2=0
for pbd in pb_dumps:
	print(pbd)
	dis=np.zeros(5);
	filename=avg_xcoor_dir+avg_fwc + str(pbd) +'.out';
	node = oio.import_file(filename)
	modifier = DislocationAnalysisModifier()
	modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
	node.modifiers.append(modifier)
	
	data=node.compute()
	cell = data.cell
	box_size = np.diagonal(data.cell.matrix)
	box_length=box_size[0]
    # box_size = numpy.diagonal(node.source.cell.matrix)
    # print("Size X: ", box_size[0])
	dis_coor=[]
	jj=0
	for segment in data.dislocations.segments:
		dis_coor=segment.points
		dis[jj]=np.mean(dis_coor,axis=0)[0]
        #dis[jj]=dis_coor[0][0]
		jj=jj+1

	x_max[i]=np.maximum(dis[0],dis[1]);
	x_min[i]=np.minimum(dis[0],dis[1]);
	if abs(x_max[i]-x_min[i])<box_length/2:
		x_max[i]=np.maximum(dis[0],dis[1]);
		x_min[i]=np.minimum(dis[0],dis[1]);
	else:
		x_min[i]=np.maximum(dis[0],dis[1]);
		x_max[i]=np.minimum(dis[0],dis[1]);
    
	time[i]=pbd*.001;
	i=i+1
snd_partial=np.copy(x_max);
first_partial=np.copy(x_min);
for k in range(num_dump-1):
	diff_snd[k]=snd_partial[k+1]-snd_partial[k]
	diff_first[k]=first_partial[k+1]-first_partial[k]

	if abs(diff_snd[k])>box_length/2:
		snd_partial[k+1:num_dump]=snd_partial[k+1:num_dump]+box_length
		diff_snd[k]=snd_partial[k+1]-snd_partial[k]


	if abs(diff_first[k])>box_length/2:
		first_partial[k+1:num_dump]=first_partial[k+1:num_dump]+box_length
		diff_first[k]=first_partial[k+1]-first_partial[k]

    
# plt.figure(figsize=(10,10))
# plt.scatter(time, diff_snd, s=10, c='b');
# fig_name = input_dump_wc + 'diff'+str(int(pbd))+'.png';
# plt.savefig(fig_name, dpi=200);
# plt.close();
final = {'time': time[0:num_dump], 'diff_first': first_partial, 'diff_snd': snd_partial};
import pickle as pkl;
fname = dis_pos_dir+'diff.pkl';
f1 = open(fname, "wb");
pkl.dump(final,f1);
f1.close()


