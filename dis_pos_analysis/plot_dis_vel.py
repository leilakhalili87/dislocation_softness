import pickle as pkl
import numpy as np
import sys
import matplotlib.pyplot as plt;
dir_path = ('../');
sys.path.insert(0,dir_path)
from dir_names import *;

# input_dump_wc_1='/home/leila/Leila_sndhard/codes/dislocation/disl_data/Lammps_Ni/T100/fig/';
# input_dump_wc='/home/leila/Leila_sndhard/codes/dislocation/disl_data/Lammps_Ni/T100/shear/';
f_max=open(dis_pos_dir+'diff.pkl','rb');

diff=pkl.load(f_max);
time=diff['time'];
diff_first=diff['diff_first'];
diff_snd=diff['diff_snd']

print((diff_first[199]-diff_first[1])/(time[199]-time[1]))
print((diff_snd[199]-diff_snd[1])/(time[199]-time[1]))
plt.figure(figsize=(10,10))
fig_name = sim_dump_dir+'../fig/' + 'diff_first.png';
plt.plot(time, diff_first)
plt.xlabel('time(ps)', fontsize=18)
plt.ylabel('1st_dis_pos(A)', fontsize=16)
plt.savefig(fig_name, dpi=200);

plt.figure(figsize=(10,10))
fig_name = sim_dump_dir+'../fig/' + 'diff_snd.png';
plt.plot(time, diff_snd)
plt.xlabel('time(ps)', fontsize=18)
plt.ylabel('snd_dis_pos(A)', fontsize=16)
plt.savefig(fig_name, dpi=200);


