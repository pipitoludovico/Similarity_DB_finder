#!/usr/bin/env python3
from multiprocessing import Pool, Manager
import pandas as pd
from Fingerprint_finder.Parser import ArgParser
from Fingerprint_finder.clean_smiles import SmilesCleaner
from Fingerprint_finder.Filter import Filter
from Fingerprint_finder.fetch_descriptors import FetchDescriptors
from Fingerprint_finder.Fingerprint import Read_fingerprint
from Fingerprint_finder.Finder import Finder
from rdkit.Chem import PandasTools

import time

PARSER = ArgParser()
FingerPrint, Database, Sort, filters, output, sliceStart, sliceEnd, excludes, includes = PARSER.ParseArgs()


def ProcessDF(df, qq, p_FPmol):
    print("...")
    filterQuery = None
    print("Applying filters...")
    if excludes is not None:
        for ex_string in excludes:
            if ex_string == '+':
                ex_string = '\\+'
            df = df[~df['smiles'].str.contains(ex_string)]
            print("\nAfter exclusion:")
            print(df.shape)
    if includes is not None:
        for in_string in includes:
            df = df[df['smiles'].str.contains(in_string)]
            print("\nAfter inclusion:")
            print(df.shape)
    df2 = df.copy()
    del df
    if df2.shape[0] != 0:
        smileCol = [col for col in df2.columns if col.lower().startswith('smi')]
        idCol = [col for col in df2.columns if "id" in col.lower()]
        initialSmiles = df2[smileCol[0]]
        idColName = idCol[0]
        cleaner = SmilesCleaner(initialSmiles)
        canonicalSmiles = cleaner.getCanonicalSmiles()
        df2.drop(columns=['smiles'])
        df2['CanonicalSmiles'] = canonicalSmiles

        try:
            filtered = Filter(df2, idColName)
        except Exception as e:
            print(e)
            print(
                "Wrong separator or database format. Try to define a different separator with -se or check your database format ['id', 'smiles']")
            exit()
        filteredDF = filtered.getFiltered()
        print("Adding molecular descriptors: Weight and Charge.")
        dfWithDesc = FetchDescriptors(filteredDF, excludes)
        print("Creating descriptors...")
        dfWithDesc.CreateDescriptors()
        print("OK TEST")
        df_wd = dfWithDesc.GetDFwithDescriptors()
        finaldf = Finder(p_FPmol, df_wd)
        result = finaldf.getDFwithFP()

        result2 = result.copy()
        if len(filters) != 0:
            for _f in filters:
                print("Filtering:", _f)
                _filter = _f.split()
                if str(_filter[0]).startswith('simi'):
                    column, operator, criterium = _filter[0], _filter[1], float(_filter[2]) / 100
                elif str(_filter[0]).startswith('for'):
                    column, operator, criterium = _filter[0], _filter[1], int(_filter[2])
                else:
                    column, operator, criterium = _filter[0], _filter[1], float(_filter[2])

                if operator == '==' or operator == "=":
                    filterQuery = f"{column} == {criterium}"
                if operator == '>=':
                    filterQuery = f"{column} >= {criterium}"
                if operator == '<=':
                    filterQuery = f"{column} <= {criterium}"
                if operator == '<':
                    filterQuery = f"{column} < {criterium}"
                if operator == '>':
                    filterQuery = f"{column} > {criterium}"
                result2 = result2.copy().query(filterQuery)
            result2.drop_duplicates(subset=['CanonicalSmiles'], inplace=True)
        print("\nAdding filtered dataframe to queue:")
        print(result2.head())
        print('\n')
        qq.put(result2)
    else:
        qq.put(pd.DataFrame())


def main():
    df_list = []
    FPreader = Read_fingerprint()
    FPmol = FPreader.GetFPmol(FingerPrint)
    try:
        manager = Manager()
        q = manager.Queue()
        results = []
        batch = 0
        with Pool() as p:
            for df in pd.read_csv(Database, sep=None, chunksize=100000, engine='python'):
                results.append(p.apply_async(ProcessDF, args=(df, q, FPmol,)))
            for result in results:
                result.wait()
                df_list.append(q.get())
                print("Batch #: ", batch, " completed.")
                batch += 1
            p.close()
            p.terminate()
            p.join()

        df2 = pd.concat(_df for _df in df_list if len(_df) != 0)
        df2.drop_duplicates(subset=['CanonicalSmiles'], inplace=True)

        print("\nTotal number of compounds found: ", df2.shape[0])
        ultimateDF = df2.sort_values(['similarity'], ascending=Sort)
        ultimateDF = ultimateDF.iloc[sliceStart:sliceEnd]
        UltimateIDCol = [uCol for uCol in ultimateDF.columns if "id" in uCol.lower()][0]

        try:
            PandasTools.SaveXlsxFromFrame(ultimateDF.head(output),
                                          f'{FingerPrint.replace(".pdb", "")}_{Database.replace(".smi", "")}.xlsx',
                                          molCol='ROMol')
            dfCSV = pd.DataFrame({'id': ultimateDF[UltimateIDCol], 'smiles': ultimateDF['CanonicalSmiles']})
            dfCSV.head(output).to_csv(f'{FingerPrint.replace(".pdb", "")}_{Database}', index=False)
        except:
            raise Exception("DataFrame must not be empty")
    except pd.errors.ParserError:
        print("Couldn't parse your database. Check if your db has the id and smiles columns")


if __name__ == '__main__':
    start = time.perf_counter()
    main()
    finish = time.perf_counter()
    print("Filtering completed in:", (finish - start))
