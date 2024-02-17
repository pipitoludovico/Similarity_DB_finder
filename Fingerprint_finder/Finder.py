from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols


class Finder:

    def __init__(self, ref=None, df=None):
        self.finalDF = None
        self.reference = FingerprintMols.FingerprintMol(ref, minPath=1, maxPath=8, fpSize=512,
                                                        bitsPerHash=2, useHs=True, tgtDensity=0.05,
                                                        minSize=128)
        self.DBfingerPrints = [FingerprintMols.FingerprintMol(x, minPath=1, maxPath=8, fpSize=512,
                                                              bitsPerHash=2, useHs=True, tgtDensity=0.05,
                                                              minSize=128) for x in df['ROMol']]

        self.similarities = [DataStructs.FingerprintSimilarity(self.reference, fp) for fp in self.DBfingerPrints]
        df_copy = df.copy()  # Create a copy of the DataFrame
        df_copy['similarity'] = self.similarities
        self.finalDF = df_copy.sort_values(['similarity'], ascending=False)

    def getDFwithFP(self):
        return self.finalDF
