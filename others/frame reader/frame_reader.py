import os, glob
import pandas as pd
import numpy as np
from io import StringIO
import statistics as stat

class Frames:
    #grabs the MMPBSA Dataset
    STR_SET = ['Complex', 'Receptor', 'Ligand', 'Delta']

    def __init__(self, spreadsheet: str):
        #load spreadsheet directory
        self.file_comp = os.path.basename(spreadsheet).split(".")
        self.dir = spreadsheet
        print(f'FILE LOADED: {os.path.basename(spreadsheet)}')

    def restructure_mmpbsa_output(self):
        #seperate all the aggregated datasets from each other
        with open(self.dir,'r') as rd:
            a = rd.read().split("\n")
            rd.close()
        return a

    def load_df(self, file = False, stringCSV = None):
        #format filter
        if not file:
            if self.file_comp[1] == 'csv':
                self.df = pd.read_csv(self.dir)
            elif self.file_comp[1] == 'xls' or self.file_comp[1] == 'xlsx':
                self.df = pd.read_excel(self.dir, engine='openpyxl')
        elif file:
            csv = StringIO(stringCSV)
            self.df = pd.read_csv(csv, index_col=0)
        return self.df

    def get_array(self,str_code: str):
        #subsequently read the dataframe and find the specific substring of interest and get the cell values below it.
        #this is in coordination with the dataframe structure of MMPBSA output

        #specify if the target is Complex, Receptor, Ligand, or Delta
        
        if str_code.title() in Frames.STR_SET:
            #setting code
            self.str_code = str_code.title() + " Energy Terms"
            #fetching the split data
            mmpbsa_data = self.restructure_mmpbsa_output()
            i = mmpbsa_data.index(self.str_code)
            file = "\n".join(mmpbsa_data[i+1:i+2+334])
            df1 = self.load_df(file=True,stringCSV=file)
            return df1

        else:
            raise Exception("Str code invalid. Choose among the following energy terms.\n 'Complex', 'Receptor', 'Ligand', 'Delta'")

    def get_zero_vdW(self, term = None):
        #Get the frame values where vdW is zero
        df = self.get_array(term)
        df.replace(-0.0, 0, inplace=True)
        return df.loc[(df['VDWAALS'] == 0.00)]

    def fetch_agg(self):
        #get the average values of all data
        colNames = self.df.columns
        agg = [self.df[s].mean() for s in colNames]
        trnspsd =  pd.DataFrame(agg).transpose()
        trnspsd.columns = colNames
        return trnspsd

    def save_summary(self, filename, savedir=os.getcwd()):
        df_agg = self.fetch_agg()
        df_agg.to_csv(os.path.join(savedir,filename+'.csv'), index=True, encoding='utf-8')

if __name__ == "__main__":
    directory = glob.glob(os.path.join('others/frame reader/data',"*","mmpbsa_**PB*"))

    for mmpbsa_data in directory:

            a = Frames(mmpbsa_data)
            d = a.get_zero_vdW('complex')

            #fetch filename 
            fname = os.path.basename(mmpbsa_data).strip("_PB.csv")
            if not os.path.exists("others/frame reader/results"):
                os.makedirs('others/frame reader/results')

            a.save_summary(filename=fname, savedir='others/frame reader/results')


