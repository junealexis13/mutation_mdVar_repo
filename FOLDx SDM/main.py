import re
import os
import time
import foldx_functions as fx
import subprocess
import datetime
import glob

if __name__ == '__main__':
    #Mutation List that falls within sorting filters
    workspace = '/Users/junealexissantos/Desktop/FOLDX_Modelling'
    model_workspace = os.path.join(workspace, '1.3')
    mutation_workspace = os.path.join(workspace,'1.1')
    cfg_workspace = os.path.join(workspace, '1.2')

    foldx_app = '/Users/junealexissantos/Documents/Executables/foldx/foldx_20221231'
    #7X6A matches all mutations in the list model cant be created
    #modified E505H for wild type mutagenesis
    mut_list = ['T478K','S477N','E484A','Y505H','Q493R','Q498R','N501Y','N440K']

    #foldX prep
    for mutations in mut_list:

        print(f"Creating individual mutation list. - {mutations}"); time.sleep(1)
        #working on mut indi list
        if not os.path.exists(os.path.join(mutation_workspace,mutations)):
            os.makedirs(os.path.join(mutation_workspace,mutations))
        fx.create_mutlist(os.path.join(mutation_workspace, mutations),mutations,'A')

        print(f"Creating individual config list. - {mutations}"); time.sleep(1)
        # working on cfg files
        if not os.path.exists(os.path.join(cfg_workspace,mutations)):
            os.makedirs(os.path.join(cfg_workspace,mutations))
        if not os.path.exists(os.path.join(model_workspace, mutations)):
            os.makedirs(os.path.join(model_workspace, mutations))
        fx.create_cfg(os.path.join(cfg_workspace, mutations),r"/Users/junealexissantos/Desktop/FOLDX_Modelling/omicron_7cwu_RBD.pdb",
        os.path.join(mutation_workspace, mutations, 'individual_list.txt'),out_model=os.path.join(model_workspace, mutations),runs = 10, tags=mutations+'MMod')

    #autoRun
    for files in os.listdir(cfg_workspace):
        if files != '.DS_Store':
            for sfiles in os.listdir(os.path.join(cfg_workspace,files)):
                print('Subprocess Initiated', datetime.datetime.now())
                subprocess.run(["cd",os.path.join(cfg_workspace,files,sfiles)])
                subprocess.run(["chmod","+X", "/Users/junealexissantos/Documents/Executables/foldx/foldx_20221231"])
                # subprocess.run([r"D:Programs\FoldX\foldx.exe",'-f',os.path.join(cfg_workspace,files,sfiles)])
                # print('Done!', f"completed: {sfiles}")

    print('####### INITIALIZING FOLDX RUN - AUTO  ########')
    time.sleep(1.5)
    for sfiles in glob.glob(os.path.join(cfg_workspace,"*","*.cfg")):
        print('Subprocess Initiated', datetime.datetime.now())
        subprocess.run(["cd", "-f", os.path.dirname(sfiles)])
        subprocess.run(["/Users/junealexissantos/Documents/Executables/foldx/foldx_20221231", "-f", sfiles])
        # subprocess.run([r"D:Programs\FoldX\foldx.exe",'-f',os.path.join(cfg_workspace,files,sfiles)])
        print('Done!', f"completed: {sfiles}")