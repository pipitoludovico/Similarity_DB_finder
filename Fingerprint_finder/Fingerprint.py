from rdkit import Chem
from rdkit.Chem import AllChem


class Read_fingerprint:

    def __init__(self):
        self.addH = None
        self.no_H = None
        self.fp = None

    def GetFPmol(self, fingerprint):
        print("Reading PDB file and preparing the moleculer for fingerprint scanning.")
        self.fp = Chem.MolFromPDBFile(str(fingerprint))
        self.no_H = Chem.RemoveHs(self.fp)
        self.addH = Chem.AddHs(self.no_H)
        AllChem.Compute2DCoords(self.addH)
        AllChem.MMFFOptimizeMolecule(self.addH)
        return self.addH
