import os, glob
import pandas as pd
import numpy as np

class Frames:
    def __init__(self, spreadsheet: str):
        #load spreadsheet directory
        a = os.basename(spreadsheet).split(".")
        self.dir = spreadsheet

    def load_df(self,):
        #format filter
        if a[1] == 'csv':
            self.df = pd.read_csv(self.dir)
        elif a[1] == 'xls' or a[1] == 'xlsx':
            self.df = pd.read_excel(self.dir, engine='openpyxl')

        return self.df