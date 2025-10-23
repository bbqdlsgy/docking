import os
import subprocess
from tqdm import tqdm

import pandas as pd
import numpy as np

import pubchempy as pcp
from pubchempy import Compound, get_compounds

from rdkit import Chem
from rdkit.Chem import PandasTools

import requests
import urllib
import urllib.request
from urllib.error import HTTPError, URLError

from get_pdb_from_rcsb_01 import get_pdb, gz_extract_pdb

# # pdb download
# pdbs_path = 'pdb_2023_human_X_ray_Prot_resol2/pdb_2023_human_X_ray_Prot_resol2.txt'
# output = 'pdb_2023_human_X_ray_Prot_resol2/rcsb_api_pdb'
# get_pdb (input=pdbs_path, output=output)

# # extract pdb from .gz
# extract_path = 'pdb_2023_human_X_ray_Prot_resol2/rcsb_api_pdb'
# gz_extract_pdb(input=extract_path) # .gz 파일들을 압축풀기

def Extract_MW (csv_path: str,  
                MinMW: int | float, 
                MaxMW: int | float, 
                MW_col: str | None = "Ligand MW",
                PDB_col:str | None = "Entry ID",
                output: str | None = "."):
    """
        Extrac Molecular Weight that user wants from csv 
    Args:
        csv_path (str): csv including pdbid, lig, Molecular Weight of lig
        MinMW (int | float): set minimum Molecular Weight
        MaxMW (int | float): set maximum Molecular Weight
        MW_col (str | None, optional): Molecular Weight column. Defaults to "Ligand MW".
        PDB_col (str | None, optional): PDB id column. Defaults to "Entry ID".
        output (str | None, optional): save path. Defaults to ".".
    """
    pdb_lig_path = csv_path
    pdb_2023_human = pd.read_csv(pdb_lig_path)
    
    # mw = pdb_2023_human[(pdb_2023_human['Ligand MW'] >= 250) & (pdb_2023_human['Ligand MW'] <= 400)]
    mw = pdb_2023_human[(pdb_2023_human[f"{MW_col}"] >= MinMW) & (pdb_2023_human[f"{MW_col}"] <= MaxMW)]
    
    pdb = []
    for pdbid in mw[f"{PDB_col}"].tolist():
         # print(type(pdbid))    
        if type(pdbid) == str:
            pdb.append(pdbid)
            
    mw_min_max = mw[mw[f"{PDB_col}"].isin(pdb)]
    os.makedirs(output, exist_ok=True)
    save_path = output
    mw_min_max.to_csv(f"{save_path}/{MinMW}_{MaxMW}.csv", index=False, index_label=False)
    
    # return pd.read_csv(save_path)

if __name__ == '__main__':
    
    path = '/home/bjw0118/Desktop/Projects/pdb_redocking_summary/pdb_2023_human_X_ray_Prot_resol2/01_pdb_lig_2023_human_X_ray_Prot_resol2.csv'
    output = './pdb_2023_preprocess_test'
    Extract_MW(csv_path=path,
               MinMW=300,
               MaxMW=400,
               MW_col= "Ligand MW",
               PDB_col= "Entry ID",
               output=output)
    