import traceback
from rdkit import Chem
from rdkit.Chem import Descriptors, PandasTools
import pandas as pd


class FetchDescriptors:
    def __init__(self, dataframe=None, toxicFlag=None):
        self.dataframe = dataframe
        if dataframe is not None:
            PandasTools.AddMoleculeColumnToFrame(dataframe, smilesCol='CanonicalSmiles')
        self.toxicFlag = toxicFlag

    @staticmethod
    def calculate_molecular_descriptors(smiles):
        descriptors = {}
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        else:
            descriptors['molecular_weight'] = Descriptors.MolWt(mol)
            descriptors['formal_charge'] = Chem.GetFormalCharge(mol)
            descriptors['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
            descriptors['NumHDonors'] = Descriptors.NumHDonors(mol)
            descriptors['TPSA'] = Descriptors.TPSA(mol)
            descriptors['NumHeteroatoms'] = Descriptors.NumHeteroatoms(mol)
            descriptors['NumRotatableBonds'] = Descriptors.NumRotatableBonds(mol)
            return descriptors

    def CreateDescriptors(self):
        descriptors_list = self.dataframe['CanonicalSmiles'].apply(self.calculate_molecular_descriptors)
        descriptors_df = pd.DataFrame(descriptors_list.tolist())
        self.dataframe = self.dataframe.reset_index(drop=True)
        self.dataframe = pd.concat([self.dataframe, descriptors_df], axis=1)

    def RemoveToxicAndCreate3D(self, smiles):
        from rdkit import Chem
        from rdkit.Chem import Fragments, AllChem, Descriptors3D
        from rdkit import RDLogger
        RDLogger.DisableLog('rdApp.*')
        public_method_names = [method for method in dir(Descriptors3D) if callable(getattr(Descriptors3D, method)) if
                               not method.startswith('_')]
        descriptors = {}
        toxicFrag = {}
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("mol failed")
            return None
        if self.toxicFlag is True:
            print("Removing toxic fragments from database...")
            try:
                toxicFrag['Esters'] = Chem.Fragments.fr_ester(mol)
                toxicFrag['Ar-NH2'] = Chem.Fragments.fr_Ar_NH(mol)
                toxicFrag["IdroxilAmmine"] = Chem.Fragments.fr_N_O(mol)
                toxicFrag['Thiol'] = Chem.Fragments.fr_SH(mol)
                toxicFrag['Alehydes'] = Chem.Fragments.fr_aldehyde(mol)
                toxicFrag['Carbammate'] = Chem.Fragments.fr_alkyl_carbamate(mol)
                toxicFrag['Allylic Ox'] = Chem.Fragments.fr_allylic_oxid(mol)
                toxicFrag["Azide"] = Chem.Fragments.fr_azide(mol)
                toxicFrag['Azo'] = Chem.Fragments.fr_azo(mol)
                toxicFrag['DiAzo'] = Chem.Fragments.fr_diazo(mol)
                toxicFrag['Alkyl_allide'] = Chem.Fragments.fr_alkyl_halide(mol)
                toxicFrag['Epoxide'] = Chem.Fragments.fr_epoxide(mol)
                toxicFrag['Furan'] = Chem.Fragments.fr_furan(mol)
                toxicFrag['Hydrazine'] = Chem.Fragments.fr_hdrzine(mol)
                toxicFrag['Hydrazone'] = Chem.Fragments.fr_hdrzone(mol)
                toxicFrag['Isocyanates'] = Chem.Fragments.fr_isocyan(mol)
                toxicFrag['IsoThioCyanates'] = Chem.Fragments.fr_isothiocyan(mol)
                toxicFrag['Nitriles'] = Chem.Fragments.fr_nitrile(mol)
                toxicFrag['Nitro'] = Chem.Fragments.fr_nitro(mol)
                toxicFrag['NitroBenzAromatic'] = Chem.Fragments.fr_nitro_arom(mol)
                toxicFrag['NonOrthtoNitroBenz'] = Chem.Fragments.fr_nitro_arom_nonortho(mol)
                toxicFrag['NO2'] = Chem.Fragments.fr_nitroso(mol)
                toxicFrag['Oxime'] = Chem.Fragments.fr_oxime(mol)
                toxicFrag['ThyoCyanate'] = Chem.Fragments.fr_thiocyan(mol)
            except Exception:
                print(traceback.format_exc())
        if any(value >= 1 for value in toxicFrag.values()):
            return None

        try:
            AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
            embed_success = True
        except:
            return None
        if embed_success:
            for method_name in public_method_names:
                method = getattr(Descriptors3D, method_name)
                result = method(mol)
                descriptors[method_name] = result
        return descriptors

    def Get3DDescriptors(self, preFinalDf):
        descriptors = []
        rows_to_remove = []
        for index, row in preFinalDf.iterrows():
            smiles = row['CanonicalSmiles']
            try:
                descriptor = self.RemoveToxicAndCreate3D(smiles)
                if descriptor is not None:
                    descriptors.append(descriptor)
                else:
                    rows_to_remove.append(index)
            except Exception as e:
                print(f"Error processing row {index}: {e}")

        preFinalDf.drop(rows_to_remove, inplace=True)
        descriptors_df = pd.DataFrame(descriptors)
        finalDF = pd.concat([preFinalDf.reset_index(drop=True), descriptors_df.reset_index(drop=True)], axis=1)
        return finalDF

    def GetDFwithDescriptors(self):
        return self.dataframe
