#!/usr/bin/env python3
import pandas as pd
from Fingerprint_finder.Parser import ArgParser
from Fingerprint_finder.clean_smiles import SmilesCleaner
from Fingerprint_finder.Filter import Filter
from Fingerprint_finder.fetch_descriptors import FetchDescriptors
from Fingerprint_finder.Fingerprint import Read_fingerprint
from Fingerprint_finder.Finder import Finder
from rdkit.Chem import PandasTools

PARSER = ArgParser()
FingerPrint, Database, Sort, separator, filterColumn, filterOperator, filterCriterium = PARSER.ParseArgs()


def main():
    df_list = []
    FPreader = Read_fingerprint(FingerPrint)
    print("Separator: ", '"' + separator, '"' + " ,Sorting Similarity Ascending: ", Sort)
    print("Filtering criteria: ", filterColumn, filterOperator, filterCriterium)

    for df in (pd.read_csv(f'{Database}', sep=separator, chunksize=10000, low_memory=False)):
        print(df.head())
        print("...")
        smileCol = [col for col in df.columns if col.lower().startswith('smi')]
        idCol = [col for col in df.columns if "id" in col.lower()]
        initialSmiles = df[smileCol[0]]
        idColName = idCol[0]
        cleaner = SmilesCleaner(initialSmiles)
        canonicalSmiles = cleaner.getCanonicalSmiles()
        df['CanonicalSmiles'] = canonicalSmiles
        try:
            filtered = Filter(df, idColName)
        except:
            print(
                "Wrong separator or database format. Try to define a different separator with -se or check your database format ['id', 'smiles']")
            exit()
        filteredDF = filtered.getFiltered()
        print("Adding molecular descriptors: Weight and Charge.")
        dfWithDesc = FetchDescriptors(filteredDF)
        dfWithDesc.CreateDescriptors()
        df = dfWithDesc.GetDFwithDescriptors()
        finaldf = Finder(FPreader.GetFPmol(), df)
        result = finaldf.getDFwithFP()
        result.drop_duplicates(subset=['CanonicalSmiles'], inplace=True)
        filterQuery = None
        if filterOperator == '==' or filterOperator == "=":
            filterQuery = f"{filterColumn} == {filterCriterium}"
        if filterOperator == '>=':
            filterQuery = f"{filterColumn} >= {filterCriterium}"
        if filterOperator == '<=':
            filterQuery = f"{filterColumn} <= {filterCriterium}"

        if filterOperator == '<':
            filterQuery = f"{filterColumn} < {filterCriterium}"
        if filterOperator == '>':
            filterQuery = f"{filterColumn} > {filterCriterium}"
        filteredResult = result.query(filterQuery)

        df_list.append(filteredResult)

    df2 = pd.concat(df_list)
    df2.drop_duplicates(subset=['CanonicalSmiles'], inplace=True)
    ultimateDF = df2.sort_values(['similarity'], ascending=Sort)
    try:
        PandasTools.SaveXlsxFromFrame(ultimateDF.head(100), f'{FingerPrint.replace(".pdb", "")}_{Database.replace(".smi", "")}.xlsx', molCol='ROMol')
    except:
        print("No data matching your selection. Empty dataframe.")


if __name__ == '__main__':
    main()
