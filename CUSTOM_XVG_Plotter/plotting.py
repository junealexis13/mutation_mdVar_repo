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
    def __init__(self, directory, plot_label=None ,title = None, color="blue"):
        self.directory = path.abspath(directory)
        self.title = title
        self.plot_label = plot_label
        self.color = color

        self.parsed_file = XVGReader(self.directory)
    def plot(self):
        if self.plot_label is None:
            self.plot_label = os.path.basename(self.directory)
            plt.legend()
        plt.plot([ts.data[0] for ts in self.parsed_file], [ts.data[1] for ts in self.parsed_file], label= self.plot_label, color = self.color )
        plt.title(self.title)
        plt.show()
