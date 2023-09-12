from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


class Finder:

    def __init__(self, ref=None, df=None):
        self.finalDF = None
        self.FP2 = FingerprintMols.FingerprintMol(ref, minPath=1, maxPath=7, fpSize=2048,
                                                  bitsPerHash=2, useHs=False, tgtDensity=0.0,
                                                  minSize=128)
        self.DBfingerPrints = [FingerprintMols.FingerprintMol(x, minPath=1, maxPath=7, fpSize=2048,
                                                              bitsPerHash=2, useHs=False, tgtDensity=0.0,
                                                              minSize=128) for x in df['ROMol']]
        self.similarities2 = [DataStructs.FingerprintSimilarity(self.FP2, fp) for fp in self.DBfingerPrints]
        df_copy = df.copy()  # Create a copy of the DataFrame
        df_copy['similarity'] = self.similarities2
        self.finalDF = df_copy.sort_values(['similarity'], ascending=False)

    def getDFwithFP(self):
        return self.finalDF
