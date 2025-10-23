import sys
sys.path.append('/home/bjw0118/Desktop/aimblelib')

import os
import subprocess
from tqdm import tqdm

import requests
import urllib
import urllib.request
from urllib.error import HTTPError, URLError

import pandas as pd 
from rdkit import Chem
from rdkit.Chem import PandasTools

from aimblelib import common
from aimblelib.pdb import PDB
from aimblelib.ligand import Ligand
from aimblelib.docking import Docking

path = 'pdb_2023_human_X_ray_Prot_resol2/final_output/pubchem_3d_exist/resol_min_pdb_3d_exist.sdf'
sdf = PandasTools.LoadSDF(path)

# pdb, cid, lid of sdf
exist_3d_pdb = sdf['PDB'].tolist()[1:5]
exist_3d_cid = sdf['CID'].tolist()[1:5]
exist_3d_lid = sdf['LID'].tolist()[1:5]
exist_3d_isosmi = sdf['isomeric_smiles'].tolist()[1:5]
    
def DownloadConformer3dSingle (output: str | None = '.'):
    """
        Download single 3d conformer 
    Args:
        output (str | None, optional): Download path. Defaults to '.'.
    """
    
    # save pdb & 3d conformer
    output = output

    # check = []
    for digit, cid, lid, iso_smi in zip(exist_3d_pdb, exist_3d_cid, exist_3d_lid, exist_3d_isosmi):
        os.makedirs(output + f"/{digit}", exist_ok=True)
        
        try:
            print(f'{digit}_{lid}_{cid}', '\n', 
                  f"we are searching 3d conformer with isomeric smiles", {iso_smi})
            url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{iso_smi}/SDF?record_type=3d&response_type=save&response_basename=Conformer3D_COMPOUND_CID_{cid}'
            response = urllib.request.urlopen(url)
            savename = output + f"/{digit}" + f'/3D_COMPOUND_CID_{digit}_{lid}_{cid}.sdf'
            
            # url이 가리키는 주소에 접근해서 해당 자원을 로컬 컴퓨터에 저장하기
            urllib.request.urlretrieve(url, savename)
        except HTTPError as herr:
            print(f'{digit}_{lid}_{cid}' + 'Reason: ', herr.reason, '\n',
                  f"We will search 3d conformer with {cid}")
            url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d&response_type=save&response_basename=Conformer3D_COMPOUND_CID_{cid}'
            response = urllib.request.urlopen(url)
            savename = output + f"/{digit}" + f'/Conformer3D_COMPOUND_{digit}_{lid}_{cid}.sdf'
            
            # url이 가리키는 주소에 접근해서 해당 자원을 로컬 컴퓨터에 저장하기
            urllib.request.urlretrieve(url, savename)

def DownloadConformer3dMultiple (output: str | None = '.'):
    """
        Download multiple 3d conformer 
    Args:
        output (str | None, optional): Download path. Defaults to '.'.
    """
    
    # save pdb & 3d conformer
    # output = 'pdb_2023_human_X_ray_Prot_resol2/final_output/pubchem_3d_exist/final_pdb'
    output = output

    for digit, cid, lid, iso_smi in zip(exist_3d_pdb, exist_3d_cid, exist_3d_lid, exist_3d_isosmi):
        os.makedirs(output + f"/{digit}", exist_ok=True)

        try:
            print(f'{digit}_{lid}_{cid}', '\n', 
                  f"we are searching 3d conformer with cid", '\n', 'cid: ', {cid})
            id_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/conformers/TXT?record_type=3d"
            response = urllib.request.urlopen(id_url)
            saveid = output + f"/{digit}" + f"/{digit}_{lid}_{cid}_3d_id.txt"
            urllib.request.urlretrieve(id_url, saveid)
       
            
            response = requests.get(id_url) # id가 여러개이기 때문에 id를 구분하기 위해 다시 부름.
            for line in response.text.splitlines():
                print(f"{line}_{cid}")
                id_cid_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/conformers/{line}/SDF?response_type=save&response_basename=Conformer3D_COMPOUND_CID_{cid}_{line}"
                response = urllib.request.urlopen(id_cid_url)
                saveid = output + f"/{digit}" + f"/Conformer3D_COMPOUND_{digit}_{lid}_{cid}_{line}.sdf"
                urllib.request.urlretrieve(id_cid_url, saveid)
        except:
            print('except')

def Conformer3dSingleDocking (nummodes: int | None = 10, 
                              energyRange: int | float | None = 5, 
                              exhaustive: int | None = 8, 
                              output: str | None = '.'):
    """
                              
        Docking for Single 3d conformation (~.sdf)
    Args:
        nummodes (int | None, optional): maximum number of binding modes to generate. Defaults to 10.
        energyRange (int | float | None, optional): maximum energy difference between the best binding mode and the worst one displayed (kcal/mol). Defaults to 5.
        exhaustive (int | None, optional): exhaustiveness of the global search (roughly proportional to time). Defaults to 8.
        output (str | None, optional): Download docking result files for multiple 3d conformation. Defaults to '.'.
    """
    
    docking_path = '/home/bjw0118/Desktop/aimblelib'
    save_path = output
    
    
    for digit, cid, sdf_lid, iso_smi in zip(exist_3d_pdb, exist_3d_cid, exist_3d_lid, exist_3d_isosmi):
        
        digit_lower = digit.lower()
        pdb = PDB.constructByID(digit, f"{save_path}/{digit}", makeDir=True)
        liginfo = pdb.getLigandList()
        # print(liginfo)
        for i in range(len(liginfo)):     
            ligid = liginfo[i].monomerName
            chainid = liginfo[i].key.chainID
            SeqNum = liginfo[i].key.sequenceNumber
            
            centerX_, centerY_, centerZ_, sizeX_, sizeY_, sizeZ_ = pdb.getLigandBox(
                                                                                code = ligid,
                                                                                chainID=chainid,
                                                                                sequenceNumber=SeqNum
                                                                                )

            
            if sdf_lid in ligid:
            
                # check 
                # print('pdb: ', digit, '\n',
                # 'sdflig: ', sdf_lid, '\t', 'getlid: ', ligid, '\n',
                # 'chain: ', chainid, '\n',
                # 'seqnum: ', SeqNum, '\n',
                # 'center: ', centerX_, centerY_, centerZ_, '\n',
                # 'size: ', sizeX_, sizeY_, sizeZ_, '\n', '#'*150)
                
                sp = subprocess.Popen(f"python {docking_path}/pdb_2023_docking_cli.py -f {save_path}/{digit}/{digit_lower}.pdb -c {chainid} -b {centerX_} {centerY_} {centerZ_} {sizeX_} {sizeY_} {sizeZ_} \
                                                                -l {save_path}/{digit}/3D_COMPOUND_CID_{digit}_{sdf_lid}_{cid}.sdf --3d -n {nummodes} -e {exhaustive} -r {energyRange} \
                                                                -o {save_path}/{digit}/3D_COMPOUND_CID_{digit}_{sdf_lid}_{cid}"
                                                                , 
                                    shell=True, stdout=subprocess.PIPE,  
                                    stderr=subprocess.PIPE)
                
                s, e = sp.communicate()
                print(s.decode('utf-8'))
                print(e.decode('utf-8'))
                print(sp.returncode)   
                
                
            else:
                continue


def Conformer3dMultipleDocking (nummodes: int | None = 10, 
                                energyRange: int | float | None = 5, 
                                exhaustive: int | None = 8, 
                                output: str | None = '.'):

    """
        Docking for Multiple 3d conformation (~.sdf)
    Args:
        nummodes (int | None, optional): maximum number of binding modes to generate. Defaults to 10.
        energyRange (int | float | None, optional): maximum energy difference between the best binding mode and the worst one displayed (kcal/mol). Defaults to 5.
        exhaustive (int | None, optional): exhaustiveness of the global search (roughly proportional to time). Defaults to 8.
        output (str | None, optional): Download docking result files for multiple 3d conformation. Defaults to '.'.
    """
    
    docking_path = '/home/bjw0118/Desktop/aimblelib'
    save_path = output

    for digit, cid, sdf_lid, iso_smi in zip(exist_3d_pdb, exist_3d_cid, exist_3d_lid, exist_3d_isosmi):
        
        digit_lower = digit.lower()
        pdb = PDB.constructByID(digit, f"{save_path}/{digit}", makeDir=True)
        liginfo = pdb.getLigandList()
        
        for i in range(len(liginfo)):     
            ligid = liginfo[i].monomerName
            chainid = liginfo[i].key.chainID
            SeqNum = liginfo[i].key.sequenceNumber
            
            centerX_, centerY_, centerZ_, sizeX_, sizeY_, sizeZ_ = pdb.getLigandBox(code = ligid,
                                                                                chainID=chainid,
                                                                                sequenceNumber=SeqNum)

            
            if sdf_lid in ligid:
            
                # check 
                # print('pdb: ', digit, '\n',
                # 'sdflig: ', sdf_lid, '\t', 'getlid: ', ligid, '\n',
                # 'chain: ', chainid, '\n',
                # 'seqnum: ', SeqNum, '\n',
                # 'center: ', centerX_, centerY_, centerZ_, '\n',
                # 'size: ', sizeX_, sizeY_, sizeZ_, '\n', '#'*150)
            
                id_path = os.listdir(f"{save_path}/{digit}")
                for file in id_path:
                    sdf = file.endswith('.sdf')
                    if sdf == True:
                        fn, format_ = file.split('.')
                        id = fn.split('_')[-1]
                        sp = subprocess.Popen(f"python {docking_path}/pdb_2023_docking_cli.py -f {save_path}/{digit}/{digit_lower}.pdb -c {chainid} -b {centerX_} {centerY_} {centerZ_} {sizeX_} {sizeY_} {sizeZ_} \
                                                                        -l {save_path}/{digit}/Conformer3D_COMPOUND_{digit}_{sdf_lid}_{cid}_{id}.sdf --3d -n {nummodes} -e {exhaustive} -r {energyRange} \
                                                                        -o {save_path}/{digit}/Conformer3D_COMPOUND_{digit}_{sdf_lid}_{cid}_{id}_result/", 
                                                                                                                                        
                                            shell=True, stdout=subprocess.PIPE,  
                                            stderr=subprocess.PIPE)
                        
                        s, e = sp.communicate()
                        print(s.decode('utf-8'))
                        print(e.decode('utf-8'))
                        print(sp.returncode)   
                    else: 
                        pass
                
                
            else:
                continue


if __name__ == '__main__':
    
    ## single
    conformer3d_SingleDown_path = 'pdb_2023_human_X_ray_Prot_resol2/final_output/pubchem_3d_exist/pdb_2023_single_final_test4'
    DownloadConformer3dSingle(output=conformer3d_SingleDown_path)
    
    # conformer3d_Single_docking_Down_path = 'pdb_2023_human_X_ray_Prot_resol2/final_output/pubchem_3d_exist/pdb_2023_single_final_test4'
    # Conformer3dSingleDocking(output=conformer3d_Single_docking_Down_path, nummodes=10, exhaustive=24)
    
    
    ## multiple
    # conformer3d_MultiDown_path = 'pdb_2023_human_X_ray_Prot_resol2/final_output/pubchem_3d_exist/pdb_2023_multiple_final_test1'
    # DownloadConformer3dMultiple (output = conformer3d_MultiDown_path)
    
    # conformer3d_multiple_docking_Down_path = 'pdb_2023_human_X_ray_Prot_resol2/final_output/pubchem_3d_exist/pdb_2023_multiple_final_test1'
    # Conformer3dMultipleDocking (output = conformer3d_multiple_docking_Down_path, nummodes=10, exhaustive=24)
    
    
    
