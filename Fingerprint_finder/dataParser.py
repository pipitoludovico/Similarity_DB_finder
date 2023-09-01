import subprocess
import sys
from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd


class Parser:
    def __init__(self, data):
        self.df2 = None
        if data is None:
            print("Please, input a .csv or excel file!")
        if data.endswith('.csv'):
            self.data = pd.read_csv(data)
            print(self.data)
        if data.endswith('.pdb'):
            try:
                self.obabel_result = subprocess.check_output("obabel -i pdb " + sys.argv[1] + " -o smi", shell=True,
                                                             text=True)
                self.pdbRep = {'SMILES': str(self.obabel_result).split()[0]}
                self.draw = Chem.MolFromSmiles(str(self.obabel_result).split()[0])
                self.data = pd.DataFrame([self.pdbRep])
                Draw.MolToFile(self.draw, f"{str(sys.argv[1]).replace('.pdb', '.png')}")
            except Exception as e:
                if hasattr(e, 'Openbabel has failed: trying with RDKit now. Please check bonds and geometry!'):
                    self.mol = Chem.MolFromPDBFile(data)
                    self.smileFromMol = Chem.MolToSmiles(self.mol)
                    self.pdbRep = {'SMILES': self.smileFromMol}
                    self.data = pd.DataFrame([self.pdbRep])
                    Draw.MolToFile(self.mol, f"{str(sys.argv[1]).replace('.pdb', '.png')}")
                else:
                    print(e)

        if data.endswith('.xlml'):
            self.data = pd.read_excel(data)
        else:
            df_list = []
            for df in pd.read_csv(data, sep=" ", chunksize=100000, low_memory=False):
                df_list.append(df)

    def getBigDf(self):
        return self.df2
