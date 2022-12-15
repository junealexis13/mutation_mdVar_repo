import re
import os
import subprocess 

def create_mutlist(dst: str, mutations: str, chain: str):
    '''Create an individual_list.txt params for FoldX'''

    f=open(os.path.join(dst, 'individual_list.txt'), 'w')
    f.write(f'{mutations[:1] + chain + mutations[1:]};')
    f.close()

def create_cfg(dst: str, pdbfile: str, mutfile: str, out_model:str ,runs=1, tags=''):
    '''Create a configuration file to be used for mutant modelling'''
    f = open(os.path.join(dst, f'config_{tags}.cfg'), 'w')
    f.write(f'command=BuildModel\npdb={os.path.basename(pdbfile)}\npdb-dir={os.path.dirname(pdbfile)}\nmutant-file={mutfile}\noutput-dir={out_model}\nnumberOfRuns={runs}\noutput-file={tags}')
    f.close()

def run_pwrshl(command):
    out = subprocess.run(["powershell", command], capture_output=True)
    return out


# if __name__ == '__main__':
#     d = '/Users/junealexissantos/Desktop/Notes for Publication.docx'