import pandas as pd
from Fingerprint_finder.Parser import ArgParser
from Fingerprint_finder.clean_smiles import SmilesCleaner
from Fingerprint_finder.filter import Filter
from Fingerprint_finder.fetch_descriptors import FetchDescriptors
from Fingerprint_finder.Fingerprint import Read_fingerprint
from Fingerprint_finder.Finder import Finder
from rdkit.Chem import PandasTools
from sys import argv

PARSER = ArgParser()
FingerPrint, Database = PARSER.ParseArgs()

def main():
    df_list = []
    FPreader = Read_fingerprint(FingerPrint)

    for df in (pd.read_csv(f'{Database}', sep=" ", chunksize=10000, low_memory=False)):
        smileCol = [col for col in df.columns if col.lower().startswith('smile')]
        idCol = [col for col in df.columns if "id" in col.lower()]
        initialSmiles = df[smileCol[0]]
        idColName = idCol[0]
        cleaner = SmilesCleaner(initialSmiles)
        canonicalSmiles = cleaner.getCanonicalSmiles()
        df['CanonicalSmiles'] = canonicalSmiles

        filtered = Filter(df, idColName)
        filteredDF = filtered.getFiltered()
        print("Adding molecular descriptors: Weight and Charge.")
        dfWithDesc = FetchDescriptors(filteredDF)
        dfWithDesc.CreateDescriptors()
        df = dfWithDesc.GetDFwithDescriptors()

        finaldf = Finder(FPreader.GetFPmol(), df)
        result = finaldf.getDFwithFP()
        result.drop_duplicates(subset=['CanonicalSmiles'], inplace=True)  # .reset_index(drop=True)
        df_list.append(result)

    df2 = pd.concat(df_list)
    ultimateDF = df2.sort_values(['similarity2'], ascending=False)
    PandasTools.SaveXlsxFromFrame(ultimateDF.head(100),f'{FingerPrint.replace(".pdb", "")}_{Database.replace(".smi", "")}.xlsx', molCol='ROMol')


if __name__ == '__main__':
    main()
