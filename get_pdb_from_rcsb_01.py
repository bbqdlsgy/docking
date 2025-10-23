import os
import subprocess


def get_pdb (input: str, output: str | None = '.'):
    """
        Download pdb.gz files from RCSB using batch_download.sh \n
        url="https://www.rcsb.org/scripts/batch_download.sh"
    Args:
        input (str): text file with PDB ID(s) [4-digit code (e.g. 8te6 or 8TE6)] \n 
                        (url:https://www.rcsb.org/stats/all-released-structures)
        output (str | None, optional): . Defaults to '.'.
    """
    os.makedirs(output, exist_ok=True)
    
    # shell=True: 앞에 사용했던 명령어 그대로 사용하여 실행.
    # stdout, stderr: 잘 수행되는지 체크하기 위함.
    sp = subprocess.Popen(f"bash batch_download.sh -f {input} -p -o {output}", 
                          shell=True, stdout=subprocess.PIPE,  
                          stderr=subprocess.PIPE)
    
    s, e = sp.communicate()
    print(s.decode('utf-8'))
    print(e.decode('utf-8'))
    print(sp.returncode)

def gz_extract_pdb (input: str):
    """
        Extract pdb.gz files in {input} directory
    Args:
        input (str): Extract pdb files
    """
    sp = subprocess.Popen(f"gunzip {input}/*.gz", 
                          shell=True, stdout=subprocess.PIPE, 
                          stderr=subprocess.PIPE)
    
    s, e = sp.communicate()
    print(s.decode('utf-8'))
    print(e.decode('utf-8'))
    print(sp.returncode)
    
# if __name__ == "__main__":
#     test_path = 'pdb_2023_human_X_ray_Prot_resol2/pdb_2023_human_X_ray_Prot_resol2.txt'
#     # get_pdb(input=test_path, output='test')
#     gz_extract_pdb(input='test')