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

PARSER = ArgParser()
FingerPrint, Database, Sort, filters, output, sliceStart, sliceEnd, excludes, includes = PARSER.ParseArgs()


def ProcessDF(df, qq, p_FPmol):
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
    except Exception as e:
        print(e)
        print(
            "Wrong separator or database format. Try to define a different separator with -se or check your database format ['id', 'smiles']")
        exit()
    filteredDF = filtered.getFiltered()
    print("Adding molecular descriptors: Weight and Charge.")
    dfWithDesc = FetchDescriptors(filteredDF)
    dfWithDesc.CreateDescriptors()
    df = dfWithDesc.GetDFwithDescriptors()
    finaldf = Finder(p_FPmol, df)
    result = finaldf.getDFwithFP()  # ritorna un dataframe con i descrittori molecolari
    result.drop_duplicates(subset=['CanonicalSmiles'], inplace=True)
    filterQuery = None
    cumulative_result = result.copy()
    if excludes is not None:
        for ex_string in excludes:
            cumulative_result = cumulative_result[~cumulative_result['CanonicalSmiles'].str.contains(ex_string)]
            print("\nAfter exclusion:")
            print(cumulative_result.shape)

    if includes is not None:
        for in_string in includes:
            cumulative_result = cumulative_result[cumulative_result['CanonicalSmiles'].str.contains(in_string)]
            print("\nAfter inclusion:")
            print(cumulative_result.shape)

    for _f in filters:
        _filter = _f.split()
        if str(_filter[0]).startswith('simi'):
            column, operator, criterium = _filter[0], _filter[1], float(_filter[2]) / 100
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
        cumulative_result = cumulative_result.query(filterQuery)
    qq.put(cumulative_result)


def main():
    df_list = []
    FPreader = Read_fingerprint()
    FPmol = FPreader.GetFPmol(FingerPrint)
    try:
        manager = Manager()
        q = manager.Queue()
        results = []
        with Pool() as p:
            for df in pd.read_csv(Database, sep=None, chunksize=1000000, engine='python'):
                results.append(p.apply_async(ProcessDF, args=(df, q, FPmol,)))
                df_list.append(q.get())
        p.join()
        p.close()
        p.terminate()
        df2 = pd.concat(df_list)
        df2 = df2.iloc[sliceStart:sliceEnd]
        df2.drop_duplicates(subset=['CanonicalSmiles'], inplace=True)
        ultimateDF = df2.sort_values(['similarity'], ascending=Sort)
        try:
            print(ultimateDF.shape)
            PandasTools.SaveXlsxFromFrame(ultimateDF.head(output),
                                          f'{FingerPrint.replace(".pdb", "")}_{Database.replace(".smi", "")}.xlsx',
                                          molCol='ROMol')
        except:
            raise Exception("DataFrame must not be empty")
    except pd.errors.ParserError:
        print("Couldn't parse your database. Check if your db has the id and smiles columns")


if __name__ == '__main__':
    main()
