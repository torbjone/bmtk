import os
import sys


pop_names = ['e4Scnn1a', #'e4Rorb', 'e4Nr5a1',
             'i4Pvalb1',#, 'i4Pvalb2', #'LIFe4Scnn1a', 'LIFi4Pvalb'
             ]

text_to_replace = "INSERT_POP_NAME"
config_file_name = "config.json"

os.system("python3 build_network.py")
os.system("python3 run_bionet.py config_standard_run.json")
os.system("python3 run_bionet.py config_lgn.json")

for pop_name in pop_names:
    f = open(config_file_name, 'r')
    config_text = [l.replace(text_to_replace, pop_name) for l in f.readlines()]
    #print(config_text)
    f_new = open("config_{}.json".format(pop_name), 'w')
    f_new.writelines(config_text)
    f.close()
    f_new.close()
    os.system("python3 run_bionet.py config_{}.json".format(pop_name))



os.system("python3 find_firing_rates.py")
os.system("python3 reconstruct_lfp_from_kernels.py")