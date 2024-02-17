#!/usr/bin/env python3
import multiprocessing
import os

import numpy as np
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
FingerPrint, Database, Sort, filters, output, sliceStart, sliceEnd, excludes, includes, stereo, toxicFlag = PARSER.ParseArgs()


def ProcessDF(df, qq, p_FPmol, idx):
    log = open(f'reports/chunk_{idx}.log', 'w')
    log.write("...\n")
    filterQuery = None
    log.write("Applying filters...\n")
    if excludes is not None:
        for ex_string in excludes:
            if ex_string == '+':
                ex_string = '\\+'
            df = df[~df['smiles'].str.contains(ex_string)]
            log.write("\nAfter exclusion: ")
            size = str(df.shape)
            log.write(size)
    if includes is not None:
        for in_string in includes:
            df = df[df['smiles'].str.contains(in_string)]
            log.write("\nAfter inclusion: ")
            size = str(df.shape)
            log.write(size)
    df2 = df.copy()
    del df
    if df2.shape[0] != 0:
        smileCol = [col for col in df2.columns if col.lower().startswith('smi')]
        idCol = [col for col in df2.columns if "id" in col.lower()]
        initialSmiles = df2[smileCol[0]]
        idColName = idCol[0]
        cleaner = SmilesCleaner(initialSmiles)
        log.write("\nRemoving stereochemistry and rewriting the SMILES in canonical form...\n")
        canonicalSmiles = cleaner.getCanonicalSmiles(stereo)
        df2.drop(columns=['smiles'])
        df2['CanonicalSmiles'] = canonicalSmiles

        try:
            log.write("Now dropping duplicates from the dataframe...\n")
            filtered = Filter(df2, idColName)
        except Exception as e:
            print(e)
            log.write(
                "Wrong separator or database format. Try to define a different separator with -se or check your database format ['id', 'smiles']\n")
            exit()
        filteredDF = filtered.getFiltered()
        log.write("Adding molecular descriptors...\n")
        dfWithDesc = FetchDescriptors(filteredDF)
        log.write("Creating descriptors...\n")
        dfWithDesc.CreateDescriptors()
        df_wd = dfWithDesc.GetDFwithDescriptors()
        dfHead = str(df_wd.head())
        log.write(dfHead)
        finaldf = Finder(p_FPmol, df_wd)
        result = finaldf.getDFwithFP()
        result2 = result.copy()
        for _f in filters:
            log.write(f"\nFiltering: {_f}")
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
        resHead = str(result2.head())
        log.write(resHead)
        qq.put(result2)
        log.close()
    else:
        log.write("\nAdding the whole data with no filters\n")
        qq.put(pd.DataFrame())
        log.close()


def main():
    os.makedirs('reports', exist_ok=True)
    df_list = []
    FPreader = Read_fingerprint()
    FPmol = FPreader.GetFPmol(FingerPrint)
    halfCPUcount = int((multiprocessing.cpu_count()) / 2)
    try:
        manager = Manager()
        q = manager.Queue()
        results = []
        batch = 0
        with Pool(processes=halfCPUcount) as p:
            for idx, df in enumerate(pd.read_csv(Database, sep=None, chunksize=100000, engine='python')):
                results.append(p.apply_async(ProcessDF, args=(df, q, FPmol, idx,)))
            for result in results:
                result.wait()
                df_list.append(q.get())
                print("Batch #: ", batch, " completed.")
                batch += 1
            p.close()
            p.terminate()
            p.join()
    except pd.errors.ParserError:
        print("Couldn't parse your database. Check if your db has the id and smiles columns")

    try:
        df2 = pd.concat(_df for _df in df_list if not _df.empty)
    except print("No compound is matching your criteria. Please change selection criteria and rerun."):
        exit()
    # sorting on similarity first and dropping left duplicates
    ultimateDF = df2.sort_values(['similarity'], ascending=Sort)
    ultimateDF.drop_duplicates(subset=['CanonicalSmiles'], inplace=True)
    # filtering toxic compounds and getting full 3D descriptors
    finalFilter = FetchDescriptors(toxicFlag=toxicFlag)
    finalDf = finalFilter.Get3DDescriptors(ultimateDF)
    # filtering nans or RDkit's errors:
    finalDf.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
    idCol = [uCol for uCol in finalDf.columns if "id" in uCol.lower()][0]

    try:
        PandasTools.SaveXlsxFromFrame(finalDf.head(output),
                                      f'{FingerPrint.replace(".pdb", "")}_{Database.replace(".smi", "")}.xlsx',
                                      molCol='ROMol')
        dfCSV = pd.DataFrame({'id': finalDf[idCol], 'smiles': finalDf['CanonicalSmiles']})
        dfCSV.head(output).to_csv(f'{FingerPrint.replace(".pdb", "")}_{Database}', index=False)
    except:
        raise Exception("DataFrame must not be empty")


if __name__ == '__main__':
    start = time.perf_counter()
    main()
    finish = time.perf_counter()
    print("Filtering completed in:", (finish - start))
