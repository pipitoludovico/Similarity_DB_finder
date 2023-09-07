from rdkit import Chem
from rdkit.Chem import Descriptors, PandasTools, Fragments


class FetchDescriptors:
    def __init__(self, dataframe, excludes):
        self.dataframe = dataframe
        self.excludes = excludes
        self.mols = [Chem.MolFromSmiles(i) for i in dataframe['CanonicalSmiles']]
        PandasTools.AddMoleculeColumnToFrame(dataframe, smilesCol='CanonicalSmiles')

    def calculate_molecular_descriptors(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if 'COC' in self.excludes:
            ester_number = Chem.Fragments.fr_ester(mol)
            if ester_number == 0:
                return Descriptors.MolWt(mol), Chem.GetFormalCharge(mol), ester_number
            else:
                return None, None, None
        else:
            return Descriptors.MolWt(mol), Chem.GetFormalCharge(mol), Chem.Fragments.fr_ester(mol)

    def CreateDescriptors(self):
        self.dataframe['molecular_weight'], self.dataframe['formal_charge'], self.dataframe['esters'] = zip(
            *self.dataframe['CanonicalSmiles'].apply(self.calculate_molecular_descriptors))

    def GetDFwithDescriptors(self):
        return self.dataframe
