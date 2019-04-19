from ovito.data import ParticleProperty
from ovito.data import SimulationCell
from scipy.optimize import curve_fit
import numpy as np;
import ovito.modifiers as ovm;

##########################################################
def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

##########################################################
def avg_op_vals(x1, y1, lat_a, sc):
        min_x=np.min(x1) + lat_a*3;
        max_x=np.max(x1) - lat_a*3;

        dx=sc*lat_a;

        x_range=int((max_x-min_x)/dx);

        y=np.zeros(x_range); x=np.zeros(x_range);
        tx1 = min_x;
        for i in range(x_range-1):
                a=np.where((x1>=tx1) & (x1<tx1+dx))
                x[i]=tx1+dx/2
                if np.size(a) ==0:
                    # print("No atoms present in this slice!!")
                    print("No atoms in slice " + str(i));
                    y[i]=float('NaN');
                    print("xmin_win: %f \t xmax_win: %f" %(tx1, tx1+dx));
                else:
                    y[i]=np.average(y1[a])
                tx1=tx1+dx

        i=x_range-1;
        a=np.where((x1>=tx1) & (x1<=tx1+dx))
        x[i]=tx1+dx/2
        if np.size(a) ==0:
            # print("No atoms present in this slice!!")
            y[i] = float('NaN');
            print("xmin_win: %f \t xmax_win: %f" %(tx1, tx1+dx));
        else:
            y[i]=np.average(y1[a])
        tx1=tx1+dx


        nans, y_n= nan_helper(y);
        y[nans]= np.interp(y_n(nans), y_n(~nans), y[~nans])

        # print(y);
        # tind1 = np.where(y == 1e6)[0];
        # nsz = np.size(y);
        # if (np.size(tind1) > 0):
        #     print(tind1);
        #     print(nsz);
        #     for ct1 in tind1:
        #         if (ct1 == 0):
        #             y[ct1] = y[ct1 + 1];
        #         elif (ct1 == nsz-1):
        #             y[ct1] = y[ct1 - 1];
        #         else:
        #             y[ct1] = (y[ct1+1] + y[ct1-1])/2.0;
        #     print(y);

        return [x, y];

#########################################################
def weighted_mov_avg(yvals, k):
        tws = np.hstack((np.arange(1,k+2), np.flip(np.arange(1,k+1),0)))/((2*k+1)*(k+1)-k*(k+1));
        ct1 = 0;
        nsz = np.size(yvals);
        y1 = np.zeros(nsz-2*k, );
        for ct2 in range(k, nsz-k):
                # print(ct2);
                tinds = np.arange(ct2-k, ct2+k+1);
                y1[ct1] = np.sum(tws*yvals[tinds]);
                ct1 = ct1 + 1;
        return y1;

#########################################################
def tanhyp_func(x, a1, a2, a3, a4):
        return (a1*np.tanh(a2*x + a3) + a4);

#########################################################
def calc_invtan_val(y, a1, a2, a3, a4):
        return (np.arctanh((y-a4)/a1) - a3)/a2;


#########################################################
def fit_tanh_func(x_wma, y_wma):
        b1 = min(y_wma); b2 = max(y_wma);
        ### Find x-value where y_wma is equal to (b1+b2)/2
        y2 = (b1+b2)/2.0;
        ind1 = np.where(y_wma < y2);
        tind1 = np.max(ind1); tind2 = np.max(ind1)+1;

        y1 = y_wma[tind1]; y3 = y_wma[tind2];
        x1 = x_wma[tind1]; x3 = x_wma[tind2];
        x2 = x1 + (y2-y1)*(x3-x1)/(y3-y1)

        init_a1 = (b2-b1)/2; init_a4 = (b2+b1)/2;
        init_a2 = 1;
        init_a3 = np.arctanh( ((b1+b2)/2 - init_a4)/init_a1) - x2;

        init_vals = [init_a1, init_a2, init_a3, init_a4]  # for [amp, cen, wid]
        best_vals, covar = curve_fit(tanhyp_func, x_wma, y_wma, p0=init_vals)

        return best_vals;

#########################################################
def calc_gb_pos(x_wma, y_wma):
        best_vals = fit_tanh_func(x_wma, y_wma);
        # print(best_vals);

        y_tanh = tanhyp_func(x_wma, best_vals[0], best_vals[1], best_vals[2], best_vals[3]);
        # plt.plot(x_wma, y_tanh);

        ### Calculate GB position
        a1=best_vals[0]; a2=best_vals[1]; a3=best_vals[2]; a4=best_vals[3];
        b1 = tanhyp_func(-1e08, a1, a2, a3, a4);
        b2 = tanhyp_func(1e08, best_vals[0], best_vals[1], best_vals[2], best_vals[3]);
        b3 = (b1+b2)/2;
        gb_pos = calc_invtan_val(b3, a1, a2, a3, a4);

        trange = 0.99;
        b4 = trange*b1 + (1-trange)*b2; b5 = (1-trange)*b1 + trange*b2;
        gb_min1 = calc_invtan_val(b4, a1, a2, a3, a4);
        gb_max1 = calc_invtan_val(b5, a1, a2, a3, a4);
        # print([b1, b2, b3]); print(gb_pos);

        return [[gb_min1, gb_pos, gb_max1], best_vals];

#########################################################



###############################################################################
def write_phop_dump(filename1, dir_name, phop, t0):
    pipeline = oio.import_file (filename1); ## Particles are already sorted
    data = pipeline.compute();
    data.create_user_particle_property('p_hop', 'float', 1, phop);
    phop_dump = dir_name + "phop_xcoor_%i.out"%(t0);
    oio.export_file(data, phop_dump,format = "lammps_dump",
                columns=["Particle Identifier",
                         "Particle Type",
                         "Position.X",
                         "Position.Y",
                         "Position.Z",
                         "c_csym",
                         "p_hop"], frame=t0)
    pipeline = [];

############################################################################
#Creating the modifier for subtracting center of mass from all x positions
############################################################################
def modify_xcoor(frame, input, output):
    old_pos_property = input.particle_properties['Position']
    old_pos = old_pos_property.array


    p_type_prop = input.particle_properties['Particle Type'];
    p_type = p_type_prop[...];

    ## Using the atoms of type 4
    tind1 = np.where(p_type == 1)[0];
    # Compute new positions by subtracting the same shift vector from all current positions:
    xdisp = np.average(old_pos[tind1,0]);

    # Change simulation cell
    assert(output.cell is input.cell)
    # Make a copy of the simulation cell:
    cell = output.copy_if_needed(output.cell)
    # copy_if_needed() made a deep copy of the simulation cell object.
    # Now the the input and output each point to different objects.
    assert(cell is output.cell)
    assert(cell is not input.cell)

    with cell:
        cell[0,3] = cell[0,3] - xdisp;

    # print(output.cell);

    new_pos=np.array([(old_pos[:,0]-xdisp), old_pos[:,1],old_pos[:,2]])
    new_pos=np.transpose(new_pos)

    # Output the coordinates as 'Position' particle property, effectively replacing old property:
    output.create_particle_property(ParticleProperty.Type.Position, new_pos)

###############################################################################
def calc_win_disps(num_atoms, tsz, A_window, fstr, pipeline):
    ### Displacements for A-window
    dispx_awin = np.zeros((num_atoms, tsz));
    dispy_awin = np.zeros((num_atoms, tsz));
    dispz_awin = np.zeros((num_atoms, tsz));

    for ct1 in range(tsz):
        tstep1 = A_window[ct1];
        filename0 = fstr+str(tstep1)+'.out';

        dmod = ovm.CalculateDisplacementsModifier();
        dmod.reference.load(filename0);
        pipeline.modifiers.append(dmod);

        data = pipeline.compute();
        # disp_ct1 = dfunc.get_sortedIDs_disp(data, pID_ct1);
        displacements = data.particles['Displacement']; 
        disp_ct1 = displacements[...];

        dispx_awin[:, ct1] = disp_ct1[:, 0];
        dispy_awin[:, ct1] = disp_ct1[:, 1];
        dispz_awin[:, ct1] = disp_ct1[:, 2];

        del pipeline.modifiers[0];

    disp_awin = [dispx_awin, dispy_awin, dispz_awin];
    return disp_awin;

##################################################################################
