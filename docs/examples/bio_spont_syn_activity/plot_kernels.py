from os.path import join
import numpy as np
import matplotlib.pyplot as plt
import h5py
from helpers import load_electrode_params


electrode_file = join("..", "biophys_components", "recXelectrodes", "linear_electrode.csv")
datafolder = "."

pop_names = ["lgn", 'e4Scnn1a', 'e4Rorb', 'e4Nr5a1',
             'i4Pvalb1', 'i4Pvalb2', 'LIFe4Scnn1a', 'LIFi4Pvalb', "standard_run"]

elec_depths = load_electrode_params(electrode_file)[:, 1]

# pop_name = 'e4Rorb'
for pop_name in pop_names:
    # pop_name = 'standard_run'

    f = h5py.File(join(datafolder, 'ecp_{}.h5'.format(pop_name)), 'r')
    print(f["ecp"].keys())

    channel_id = f['ecp']["channel_id"]
    lfp_kernel = 1000 * np.array(f['ecp']["data"]).T#[:, :]
    time = f['ecp']["time"]

    lfp_kernel -= lfp_kernel[:, -1, None]

    print(lfp_kernel.shape)
    num_elecs, num_tsteps = lfp_kernel.shape

    print(time[:])
    tvec = np.arange(time[0], time[1], time[2])

    max_sig = np.max(np.abs(lfp_kernel))
    elec_separation = np.abs(elec_depths[1] - elec_depths[0])
    plot_boost_factor = 4
    normalize_kernel = elec_separation * plot_boost_factor / max_sig

    if pop_name == "standard_run":

        xlim = [0, 2000]
        figname = "stadard_run_lfp"
        plt.close("all")
        fig = plt.figure(figsize=[18, 9])
        fig.subplots_adjust(right=0.9)
        ax = fig.add_subplot(111, frameon=False, xlabel="time (ms)", xlim=xlim, ylim=[-800, 25])
        # ax.axvline(250, ls='--', c='gray')
        ax.plot([xlim[1] + 100] * 2, [0, -elec_separation * plot_boost_factor], lw=3, c='k', clip_on=False)
        ax.text(xlim[1] + 150, -elec_separation * plot_boost_factor / 2, "{:1.2f} $\mu$V".format(max_sig))
        for elec_idx in range(num_elecs)[::2]:
            ax.plot(tvec, lfp_kernel[elec_idx] * normalize_kernel + elec_depths[elec_idx], c='k')
        plt.savefig("standard_run_lfp_{}.png".format(pop_name))

    else:
        xlim = [95, 150]
        plt.close("all")
        fig = plt.figure()
        fig.subplots_adjust(right=0.8)
        ax = fig.add_subplot(111, frameon=False, xlabel="time (ms)", xlim=xlim, ylim=[-800, 50])
        ax.axvline(250, ls='--', c='gray')
        ax.plot([282, 282], [0, -elec_separation * plot_boost_factor], lw=3, c='k', clip_on=False)
        ax.text(283, -elec_separation * plot_boost_factor / 2, "{:1.2f} $\mu$V".format(max_sig))
        for elec_idx in range(num_elecs)[::2]:
            ax.plot(tvec, lfp_kernel[elec_idx] * normalize_kernel + elec_depths[elec_idx], c='k')
        plt.savefig("kernel_lfp_{}.png".format(pop_name))


