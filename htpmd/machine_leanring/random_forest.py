import json
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

def random_forests_prediction(smiles, prop, form_ring=True, has_H=True):
    # The function to predict transport properties based on pre-trained random forest models. 
    # Input: smiles: the smile structure of the monomer; 
    # prop: the property to predict, including "conductivity", "li_diff", "tfsi_diff", "poly_diff" and "transference"
    # form_ring: whether form a ring structure based on the monomer; has_H: whether the ring structure contain H or not
    smiles = process_smiles(smiles, form_ring, has_H) #pre-process smiles structure
    calc = Calculator(descriptors, ignore_3D=True) #initialize descriptors calculator
#    features=calc(Chem.MolFromSmiles(smiles))
    benzene = Chem.MolFromSmiles("c1ccccc1")
    features = calc(benzene)
    print (features)
    rf = pickle.load(open('pre-trained-models-random-forests/rf_%s.sav'%(prop), 'rb'))
    output = rf.predict(features)
    return output


random_forests_prediction("O=C([Au])OCCOCCNCCOCCO[Cu]", "conductivity")
