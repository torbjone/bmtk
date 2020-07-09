from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import h5py
from helpers import load_electrode_params, read_node_types_csv

standard_run_data_folder = "output_standard_run"
kernel_data_folder = "."
network_folder = "network"

electrode_file = join("..", "biophys_components", "recXelectrodes", "linear_electrode.csv")
datafolder = "."

# pop_names = ['e4Scnn1a', 'e4Rorb', 'e4Nr5a1', 'i4Pvalb1', 'i4Pvalb2', 'LIFe4Scnn1a', 'LIFi4Pvalb']
pop_names = []

elec_depths = load_electrode_params(electrode_file)[:, 1]
node_type_ids_dict = read_node_types_csv(join(network_folder, "v1_node_types.csv"))
f_network = h5py.File(join(network_folder, 'v1_nodes.h5'), 'r')
node_type_id = np.array(f_network["nodes"]["v1"]["node_type_id"])
f_network.close()
num_cells = len(node_type_id)
# pop_names = [d_["pop_name"] for d_ in node_type_ids_dict.values()]
pop_size_dict = {}

for key, d_ in node_type_ids_dict.items():
    pop_size = len(np.where(node_type_id == key)[0])
    pop_size_dict[d_["pop_name"]] = pop_size

pop_names.append("lgn")
pop_size_dict["lgn"] = 10
num_pops = len(pop_names)

# print(pop_names)
pop_clr = {'lgn': 'g'}#{pop_name: plt.cm.jet(idx / (num_pops - 1))
           #for idx, pop_name in enumerate(pop_names)}

fig = plt.figure(figsize=[18, 9])

# ax_w = 0.55 / (num_pops - 1)
ax_w = 0.55 / (1)
ax_y0_fr = 0.8
ax_y0_ker = 0.42
ax_y0_lfp = 0.1
ax_h_lfp = 0.25
ax_x0_lfp = 0.07
ax_w_lfp = 0.9

ax_h_fr = 0.07
ax_h_ker = 0.3
xlim_ker = [0, 200]
xlim_lfp = [00, 600]
ylim_lfp = [0, -800]

f_lfp = h5py.File(join(datafolder, 'ecp_standard_run.h5'), 'r')
lfp = 1000 * np.array(f_lfp['ecp']["data"]).T#[:, :]
lfp_reconstructed = np.zeros(lfp.shape)
time_lfp = f_lfp['ecp']["time"]
tvec_lfp = np.arange(time_lfp[0], time_lfp[1], time_lfp[2])
num_elecs, num_tsteps = lfp.shape
f_lfp.close()

elec_separation = np.abs(elec_depths[1] - elec_depths[0])
# plot_boost_factor = 16


ax_lfp = fig.add_axes([ax_x0_lfp, ax_y0_lfp, ax_w_lfp, ax_h_lfp], frameon=True,
                      xlabel="time (ms)", xlim=xlim_lfp, ylim=ylim_lfp)

for p_idx, pop_name in enumerate(pop_names):
    fr = np.load(join(kernel_data_folder, "firing_rate_{}.npy".format(pop_name)))

    f = h5py.File(join(kernel_data_folder, 'ecp_{}.h5'.format(pop_name)), 'r')
    channel_id = f['ecp']["channel_id"]
    lfp_kernel = 1000 * np.array(f['ecp']["data"]).T
    time = f['ecp']["time"]
    tvec = np.arange(time[0], time[1], time[2])
    f.close()
    # Kernel must be normalized to show average post-synaptic LFP response
    lfp_kernel /= pop_size_dict[pop_name]
    print(pop_size_dict)
    lfp_kernel /= 10  # Why? Two synapses in each case?

    lfp_kernel -= lfp_kernel[:, -1, None]
    num_elecs, num_tsteps = lfp_kernel.shape
    lfp_kernel[:, :int(num_tsteps/2 - 1)] = 0

    pop_lfp_reconstructed = np.array([np.convolve(lfp_kernel[idx], fr, mode='same')
                                for idx in range(num_elecs)])

    lfp_reconstructed += pop_lfp_reconstructed

    # ax_x0 = 0.07 + p_idx / (num_pops - 1) * 0.8
    ax_x0 = 0.07 + p_idx / (1) * 0.8

    ax_fr = fig.add_axes([ax_x0, ax_y0_fr, ax_w, ax_h_fr], title=pop_name)
    ax_ker = fig.add_axes([ax_x0, ax_y0_ker, ax_w, ax_h_ker], frameon=True, xlabel="time (ms)",
                          xlim=xlim_ker, ylim=ylim_lfp)

    ax_fr.plot(tvec_lfp, fr)
    max_sig = np.max(np.abs(lfp_kernel))
    elec_separation = np.abs(elec_depths[1] - elec_depths[0])
    plot_boost_factor = 4
    normalize_kernel = elec_separation * plot_boost_factor / max_sig
    for elec_idx in range(num_elecs)[::2]:
        ax_ker.plot(tvec, lfp_kernel[elec_idx] * normalize_kernel + elec_depths[elec_idx], c='k')


for elec_idx in range(num_elecs)[::2]:
    max_sig = np.max(np.abs(lfp[:, :]))
    max_sig_re = np.max(np.abs(lfp_reconstructed[:, :]))
    normalize_lfp = elec_separation * plot_boost_factor / max_sig
    normalize_lfp_re = elec_separation * plot_boost_factor / max_sig_re
    ax_lfp.plot(tvec_lfp, lfp_reconstructed[elec_idx] * normalize_lfp + elec_depths[elec_idx], c='r')
    ax_lfp.plot(tvec_lfp, lfp[elec_idx] * normalize_lfp + elec_depths[elec_idx], c='k')

plt.savefig("kernels_and_firing_rates.png")

plt.show()
plt.close("all")
fig = plt.figure(figsize=[18, 9])
ax_lfp = fig.add_subplot(111, frameon=True,
                      xlabel="time (ms)")

elec_idx = 10
max_sig = np.max(np.abs(lfp[elec_idx, :]))
max_sig_re = np.max(np.abs(lfp_reconstructed[elec_idx, :]))

ax_lfp.plot(tvec_lfp, lfp_reconstructed[elec_idx] / max_sig, c='r')
ax_lfp.plot(tvec_lfp, lfp[elec_idx] / max_sig, c='k')

plt.show()