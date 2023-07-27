from rdkit.Chem import DataStructs, PandasTools
from rdkit.Chem.Fingerprints import FingerprintMols


class Finder:

    def __init__(self, ref=None, df=None):
        self.finalDF = None
        # self.FP = AllChem.GetMorganFingerprintAsBitVect(ref, radius=3, nBits=2048)
        #
        # self.ECFP4 = [AllChem.GetMorganFingerprintAsBitVect(x, radius=3, nBits=2048) for x in df['ROMol']]
        # self.similarities = DataStructs.BulkAllBitSimilarity(self.FP, self.ECFP4)
        # df['similarity'] = self.similarities
        #
        # df = df.sort_values(['similarity'], ascending=False)

        self.FP2 = FingerprintMols.FingerprintMol(ref, minPath=1, maxPath=7, fpSize=2048,
                                                  bitsPerHash=2, useHs=False, tgtDensity=0.0,
                                                  minSize=128)
        self.DBfingerPrints = [FingerprintMols.FingerprintMol(x, minPath=1, maxPath=7, fpSize=2048,
                                                              bitsPerHash=2, useHs=False, tgtDensity=0.0,
                                                              minSize=128) for x in df['ROMol']]
        self.similarities2 = DataStructs.BulkAllBitSimilarity(self.FP2, self.DBfingerPrints)
        df['similarity2'] = self.similarities2
        self.finalDF = df.sort_values(['similarity2'], ascending=False)

    def getDFwithFP(self):
        return self.finalDF
