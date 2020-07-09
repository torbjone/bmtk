from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import h5py

from elephant.conversion import BinnedSpikeTrain
from quantities import s, Hz, ms
import neo

from helpers import read_node_types_csv

t_stop = 600
t_start = 0
dt = 0.1
num_lgn_cells = 10
# lgn_delay = 2.0


tvec = np.arange(t_start, t_stop, dt)
fr_bin_size = 1 * dt
tvec_fr = np.arange(t_start, t_stop, fr_bin_size)

data_folder = "output_standard_run"
fr_data_folder = "."
network_folder = "network"
lgn_input_folder = join('..', 'spikes_inputs')

f_spikes = h5py.File(join(data_folder, 'spikes.h5'), 'r')
try:
    spike_node_ids = np.array(f_spikes["spikes"]["v1"]["node_ids"])
    timestamps = np.array(f_spikes["spikes"]["v1"]["timestamps"])
except KeyError:
    spike_node_ids = np.array([])
    timestamps = np.array([])

f_spikes.close()

f_network = h5py.File(join(network_folder, 'v1_nodes.h5'), 'r')
node_type_id = np.array(f_network["nodes"]["v1"]["node_type_id"])
f_network.close()

# Must pick out certain number of LGN neurons
f_lgn = h5py.File(join(lgn_input_folder, 'lgn_spikes.h5'), 'r')
lgn_neuron_ids = np.array(f_lgn['spikes']['lgn']['node_ids'])
lgn_timestamps = np.array(f_lgn['spikes']['lgn']['timestamps'])
use_spike_idxs = np.where(lgn_timestamps < t_stop)

lgn_timestamps = lgn_timestamps[use_spike_idxs]
lgn_neuron_ids = lgn_neuron_ids[use_spike_idxs]
f_lgn.close()


node_type_ids_dict = read_node_types_csv(join(network_folder, "v1_node_types.csv"))

pop_names = [d_["pop_name"] for d_ in node_type_ids_dict.values()]
pop_names.append("lgn")

num_pops = len(pop_names)
print(pop_names)
pop_clr = {pop_name: plt.cm.jet(idx / (num_pops - 1))
           for idx, pop_name in enumerate(pop_names)}

num_cells = len(node_type_id)

fig = plt.figure(figsize=[18, 9])
ax1 = fig.add_subplot(211, xlim=[-10, t_stop], xlabel="time (ms)", ylabel="cell id")
ax2 = fig.add_subplot(212, xlim=[-10, t_stop], xlabel="time (ms)", ylabel="pop. firing rates")

# ax1.scatter(timestamps, spike_node_ids, marker='|')

spike_trains = {pop_name: [] for pop_name in pop_names}
for cell_id in range(num_cells):
    cell_pop_id = node_type_id[cell_id]
    cell_pop_name = node_type_ids_dict[node_type_id[cell_id]]["pop_name"]
    spike_train = np.array(timestamps[np.where(spike_node_ids == cell_id)[0]])
    ax1.scatter(spike_train, np.ones(len(spike_train)) * cell_id, marker='|', color=pop_clr[cell_pop_name])
    spike_trains[cell_pop_name].extend(spike_train)

for cell_id in range(num_lgn_cells):

    spike_train = np.array(lgn_timestamps[np.where(lgn_neuron_ids == cell_id)[0]])
    ax1.scatter(spike_train, np.ones(len(spike_train)) * cell_id, marker='|', color=pop_clr["lgn"])
    spike_trains["lgn"].extend(spike_train)

# lgn_timestamps = lgn_timestamps[np.where(lgn_timestamps < t_stop)]
# st_lgn = neo.SpikeTrain(lgn_timestamps*ms, t_start=t_start*ms, t_stop=t_stop*ms)
# fr_lgn = BinnedSpikeTrain(st_lgn, binsize=fr_bin_size*ms, t_start=t_start*ms,
#                                              t_stop=t_stop*ms).to_array()[0]
#


firing_rates = {}


tot_num_spikes = 0
lines = []
line_names = []
max_fr = 0
for pop_name in pop_names:
    tot_num_spikes += len(spike_trains[pop_name])
    st = neo.SpikeTrain(spike_trains[pop_name]*ms, t_start=t_start*ms, t_stop=t_stop*ms)
    fr = BinnedSpikeTrain(st, binsize=fr_bin_size*ms, t_start=t_start*ms,
                                             t_stop=t_stop*ms).to_array()[0]
    max_fr = np.max([max_fr, np.max(fr)])
    firing_rates[pop_name] = fr


# firing_rates["lgn"] = fr_lgn
for p_idx, pop_name in enumerate(pop_names):
    fr = firing_rates[pop_name]
    l, = ax2.plot(tvec_fr, fr / max_fr + p_idx, c=pop_clr[pop_name])
    lines.append(l)
    line_names.append(pop_name)

    np.save(join(fr_data_folder, "firing_rate_{}.npy".format(pop_name)), fr)


fig.legend(lines, line_names, frameon=False, loc="upper right", ncol=4)
# if not tot_num_spikes == (len(timestamps) + len(lgn_timestamps)):
#     print(tot_num_spikes, len(timestamps), len(lgn_timestamps))
#     raise RuntimeError("Wrong number of spikes found")

plt.savefig("spikes.png")
