from rdkit import Chem
from rdkit.Chem import Descriptors, PandasTools


class FetchDescriptors:
    def __init__(self, dataframe):
        self.dataframe = dataframe
        self.mols = [Chem.MolFromSmiles(i) for i in dataframe['CanonicalSmiles']]
        PandasTools.AddMoleculeColumnToFrame(dataframe, smilesCol='CanonicalSmiles')

    @staticmethod
    def calculate_molecular_descriptors(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        molecular_weight = Descriptors.MolWt(mol)
        formal_charge = Chem.GetFormalCharge(mol)

        return molecular_weight, formal_charge

    def CreateDescriptors(self):
        self.dataframe['molecular_weight'], self.dataframe['formal_charge'] = zip(
            *self.dataframe['CanonicalSmiles'].apply(self.calculate_molecular_descriptors))

    def GetDFwithDescriptors(self):
        return self.dataframe
