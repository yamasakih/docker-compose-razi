import argparse
from pathlib import Path
from os.path import splitext

from joblib import Parallel, delayed
import numpy as np
from rdkit import Chem

from molvs import standardize_smiles


def desalt(line):
    smiles, name = line.strip().split('\t')
    try:
        smiles = standardize_smiles(smiles)
    except ValueError as e:
        return None
    if not Chem.MolFromSmiles(smiles):
        return None
    if len(smiles.split('.')) >= 2:        
        smiles_set = smiles.split('.')
        num_atoms = np.array([Chem.MolFromSmiles(s).GetNumAtoms() for s in smiles_set])
        return smiles_set[num_atoms.argmax()], name
    return (smiles, name)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input-file',
                                action='store',
                                default=None,
                                type=str,
                                required=True,
                                help='input SMILES file')

    parser.add_argument('-o', '--output-file',
                                action='store',
                                default=None,
                                type=str,
                                required=False,
                                help='output file name (defalut: "input SMILES file name"_desalt.smi')

    parser.add_argument('-n', '--n-jobs',
                                action='store',
                                default=1,
                                type=int,
                                required=False,
                                help='Number of CPU to desalt compounds. (default: 1)')
    
    args = parser.parse_args()    
    
    input_file = args.input_file
    n_jobs = args.n_jobs
    
    p = Path(args.input_file)
    if not p.exists():
        raise ValueError(f"{args.input_file} doesn't exist.")
    if not p.is_file():
        raise ValueError(f"{args.input_file} is not file.")
        
    if args.output_file is None:
        root, ext = splitext(args.input_file)
        output_file = f'{root}_desalt.smi'
    else:
        output_file = args.output_file
        
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    lines = Parallel(n_jobs=n_jobs, verbose=1)( [delayed(desalt)(line) for line in lines] )
    
    with open(output_file, 'w') as f:
        for line in lines:
            if line is not None:
                smiles, name = line
                f.write(f'{smiles}\t{name}\n')
