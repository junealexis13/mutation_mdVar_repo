from MDAnalysis.auxiliary.XVG import XVGReader
import matplotlib.pyplot as plt
import MDAnalysis as mda
from os import path
import os
import glob
import shlex
from dataclasses import dataclass
from typing import List
######

@dataclass
class CMD:
    command: str
    arguments: List[str]

class plot_xvg():
    def __init__(self, directory: str, others: list):
        self.directory = path.abspath(directory)
        self.parsed_file = XVGReader(self.directory)
        self.args = others

        self.params = {
                "color" : None,
                "label" : None,
                "plot_size" : 0.1,
                "plot_title" : path.basename(self.directory)
            }

    def update_params(self):
        if self.args:
            for params in list(self.params.keys()):
                if f"--{params}" in self.args or f"-{params}" in self.args:
                    self.params[params]=self.args[self.args.index(f"--{params}")+1]
        elif not self.args :
            print("None args has passed.")

    def plot(self):
        plt.plot([ts.data[0] for ts in self.parsed_file], [ts.data[1] for ts in self.parsed_file], linewidth=float(self.params['plot_size']) ,label= self.params['label'], color = self.params['color'] )
        plt.title(self.params['plot_title'])
        if self.params['label'] is not None:
            plt.legend()

        plt.show()

    def runplot(self):
        self.update_params()
        self.plot()
