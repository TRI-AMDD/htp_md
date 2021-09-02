import os
from rdkit import Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors
import pandas as pd
import numpy as np
from math import *
import pickle
from sklearn.ensemble import RandomForestRegressor

def process_smiles(smiles, form_ring, has_H):
    """Form a ring molecule for monomer."""
    mol = Chem.MolFromSmiles(smiles)
    if has_H:
        mol = Chem.AddHs(mol)
    if form_ring:
        rxn = AllChem.ReactionFromSmarts('([Cu][*:1].[*:2][Au])>>[*:1]-[*:2]')
        results = rxn.RunReactants([mol])
        assert len(results) == 1 and len(results[0]) == 1, smiles
        mol = results[0][0]
    Chem.SanitizeMol(mol)
    return mol

def random_forests_prediction(smiles, prop, form_ring=1, has_H=0):
    # The function to predict transport properties based on pre-trained random forest models. 
    # Input: smiles: the list of smile structures of the monomer; 
    # prop: the property to predict, including "conductivity", "li_diff", "tfsi_diff", "poly_diff" and "transference"
    # form_ring: whether form a ring structure based on the monomer; has_H: whether the ring structure contain H or not
    calc = Calculator(descriptors, ignore_3D=True) #initialize descriptors calculator
    mols = []
    for smile in smiles:
        smile = process_smiles(smile, form_ring, has_H)#pre-process simile structures
        mols.append(smile)
    df = calc.pandas(mols)
    df = df.apply(pd.to_numeric, errors='coerce')#force all molecules have same dimension of features
    df = df.select_dtypes(include=['int','float32','float64'])
    df = df.fillna(0)
    rf = pickle.load(open('pre-trained-models-random-forests/rf_%s.sav'%(prop), 'rb'))
    output = rf.predict(df)
    return output


output = random_forests_prediction(["CCN(CCCC(C)O[Cu])CCOC(=O)[Au]","CC(C)C(O[Cu])C(=O)N(C)CCCOC(=O)[Au]"], "conductivity")
print (output)