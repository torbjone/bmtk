from os.path import join
import numpy as np
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt
import h5py

from helpers import load_electrode_params, read_node_types_csv

def load_morphology(morph_file, x, y, x_rot, y_rot, z_rot):
    import LFPy
    cell_params = {"morphology": morph_file}
    cell = LFPy.Cell(**cell_params)

    cell.set_pos(x=x, y=y)
    # cell.set_rotation(y=y_rot, z=z_rot)
    cell.set_rotation(x=x_rot, z=z_rot, y=y_rot)
    return cell


electrode_file = join("..", "biophys_components", "recXelectrodes", "linear_electrode.csv")
morphology_folder = join("..", "biophys_components", "morphologies")
datafolder = "."
network_folder = "network"

pop_names = ['e4Scnn1a', 'e4Rorb', 'e4Nr5a1',
             'i4Pvalb1', 'i4Pvalb2', 'LIFe4Scnn1a',
             'LIFi4Pvalb', "standard_run"]

elec_pos = load_electrode_params(electrode_file)

node_type_ids_dict = read_node_types_csv(join(network_folder, "v1_node_types.csv"))

f_network = h5py.File(join(network_folder, 'v1_nodes.h5'), 'r')
node_type_id = np.array(f_network["nodes"]["v1"]["node_type_id"])

print(node_type_ids_dict)

nodes_0_x = np.array(f_network["nodes"]["v1"]["0"]["x"])
nodes_0_y = np.array(f_network["nodes"]["v1"]["0"]["y"])
nodes_0_z = np.array(f_network["nodes"]["v1"]["0"]["z"])
nodes_0_rot_x = np.array(f_network["nodes"]["v1"]["0"]["rotation_angle_xaxis"])
nodes_0_rot_y = np.array(f_network["nodes"]["v1"]["0"]["rotation_angle_yaxis"])
nodes_0_rot_z = np.array(f_network["nodes"]["v1"]["0"]["rotation_angle_zaxis"])

nodes_1_x = np.array(f_network["nodes"]["v1"]["1"]["x"])
nodes_1_y = np.array(f_network["nodes"]["v1"]["1"]["y"])
nodes_1_z = np.array(f_network["nodes"]["v1"]["1"]["z"])

num_cells = len(node_type_id)
print(num_cells, len(nodes_0_y), len(nodes_1_y))

pop_names = [d_["pop_name"] for d_ in node_type_ids_dict.values()]

num_pops = len(pop_names)
pop_clr = {pop_name: plt.cm.jet(idx / (num_pops - 1))
           for idx, pop_name in enumerate(pop_names)}

plt.close("all")
fig = plt.figure(figsize=[18, 9])
fig.subplots_adjust(right=0.9)

ax1 = fig.add_subplot(121, aspect=1, xlabel="x [$\mu$m]", ylabel="y [$\mu$m]",
                         xlim=[-200, 200])
ax2 = fig.add_subplot(122, aspect=1, xlabel="x [$\mu$m]", ylabel="y [$\mu$m]",
                         xlim=[-400, 2000], ylim=[-1000, 100])
ax1.plot(elec_pos[:, 0], elec_pos[:, 1], 'k.')
ax1.plot(nodes_0_x, nodes_0_y, '.', c='y')
ax1.plot(nodes_1_x, nodes_1_y, 's')

pop_plot_x_shifts = {pop_name: 400 * d for d, pop_name in enumerate(pop_names)}

for cell_idx in range(num_cells):
    pop_id = node_type_id[cell_idx]
    pop_name = node_type_ids_dict[pop_id]["pop_name"]
    model_type = node_type_ids_dict[pop_id]["model_type"]
    print(model_type)
    if model_type == "biophysical":
        morphology_file = node_type_ids_dict[pop_id]["morphology"]
        x = nodes_0_x[cell_idx]
        y = nodes_0_y[cell_idx]
        x_rot = nodes_0_rot_x[cell_idx]
        y_rot = nodes_0_rot_y[cell_idx]
        z_rot = nodes_0_rot_z[cell_idx]
        # print(y_rot, z_rot)
        cell = load_morphology(join(morphology_folder, morphology_file), x, y, x_rot, y_rot, z_rot)

        xshift = pop_plot_x_shifts[pop_name]
        [ax2.plot([cell.xstart[idx] + xshift, cell.xend[idx] + xshift],
                  [cell.ystart[idx], cell.yend[idx]], lw=1, c=pop_clr[pop_name])
                for idx in range(cell.totnsegs)]
        # xshift = 500
        # [ax2.plot([cell.xstart[idx] + xshift, cell.xend[idx] + xshift],
        #           [cell.zstart[idx] - 500, cell.zend[idx] - 500], lw=1, c=pop_clr[pop_name])
        #         for idx in range(cell.totnsegs)]

        fig.savefig("populations_%d.png" % cell_idx)
    else:
        break
