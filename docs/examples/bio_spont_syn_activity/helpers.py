import numpy as np

def read_node_types_csv(file_name):
    f = open(file_name)
    header = f.readline().replace("\n", "").split(" ")[1:]
    node_types_dict = {}
    for l in f.readlines():
        l_ = l.split(" ")
        d_ = {header[idx]: l_[idx + 1].replace("\n", "") for idx in range(len(header))}
        node_types_dict[int(l_[0])] = d_
    f.close()
    return node_types_dict

def load_electrode_params(electrode_file):
    elec_depth = np.loadtxt(electrode_file, skiprows=1)[:, 1:]
    return elec_depth

