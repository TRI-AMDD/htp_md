import csv
import functools
import os.path as osp
import torch
import numpy as np
from torch.utils.data import Dataset
from torch_geometric.data import Data, DataLoader
from rdkit import Chem
from rdkit.Chem import AllChem


x_map = {
    'atomic_num': list(range(0, 119)),
    'degree': list(range(0, 11)),
    'formal_charge': list(range(-5, 7)),
    'num_hs': list(range(0, 9)),
    'hybridization': [
        'UNSPECIFIED',
        'S',
        'SP',
        'SP2',
        'SP3',
        'SP3D',
        'SP3D2',
        'OTHER',
    ],
    'is_aromatic': [False, True],
    'is_in_ring': [False, True],
}


e_map = {
    'bond_type': [
        'misc',
        'SINGLE',
        'DOUBLE',
        'TRIPLE',
        'AROMATIC',
    ],
    'stereo': [
        'STEREONONE',
        'STEREOZ',
        'STEREOE',
        'STEREOCIS',
        'STEREOTRANS',
        'STEREOANY',
    ],
    'is_conjugated': [False, True],
}


def onehot(feature_list, cur_feature):
    assert cur_feature in feature_list
    vector = [0] * len(feature_list)
    index = feature_list.index(cur_feature)
    vector[index] = 1
    return vector


def process_smiles(smiles, form_ring, has_H):
    """Form a ring molecule for monomer."""
    mol = Chem.MolFromSmiles(smiles)
    if has_H:
        mol = Chem.AddHs(mol)
    if form_ring:
        rxn = AllChem.ReactionFromSmarts('([Cu][*:1].[*:2][Au])>>[*:1]-[*:2]')
        results = rxn.RunReactants([mol])
        if not (len(results) == 1 and len(results[0]) == 1):
            # Double bond
            rxn = AllChem.ReactionFromSmarts(
                '([Cu]=[*:1].[*:2]=[Au])>>[*:1]=[*:2]')
            results = rxn.RunReactants([mol])
        assert len(results) == 1 and len(results[0]) == 1, smiles
        mol = results[0][0]
    Chem.SanitizeMol(mol)
    return mol


class PolymerDataset(Dataset):
    """Polymer conducitivty dataset.

    args:
        root_dir: str, directory that stores smiles csv
        label: str, train, val, or test
        split: int, the part that is used for validation
    """

    def __init__(self, smiles_list, log10=True,
                 form_ring=True, has_H=True):
        self.raw_data = []
        for smiles in smiles_list:
            self.raw_data.append(['test', smiles, 1.])
        self.log10 = log10
        self.form_ring = form_ring
        self.has_H = has_H

    def __len__(self):
        return len(self.raw_data)

    @functools.lru_cache(maxsize=None)
    def __getitem__(self, idx):
        poly_id, smiles, target = self.raw_data[idx]

        mol = process_smiles(smiles, form_ring=self.form_ring,
                             has_H=self.has_H)
        target = float(target)
        if self.log10:
            target = np.log10(target)
        target = torch.tensor(target).float()

        xs = []
        for atom in mol.GetAtoms():
            x = []
            x += onehot(x_map['atomic_num'], atom.GetAtomicNum())
            x += onehot(x_map['degree'], atom.GetTotalDegree())
            x += onehot(x_map['formal_charge'], atom.GetFormalCharge())
            x += onehot(x_map['num_hs'], atom.GetTotalNumHs())
            x += onehot(x_map['hybridization'], str(atom.GetHybridization()))
            x += onehot(x_map['is_aromatic'], atom.GetIsAromatic())
            x += onehot(x_map['is_in_ring'], atom.IsInRing())
            xs.append(x)

        x = torch.tensor(xs).to(torch.float)

        edge_indices, edge_attrs = [], []
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            e = []
            e += onehot(e_map['bond_type'], str(bond.GetBondType()))
            e += onehot(e_map['stereo'], str(bond.GetStereo()))
            e += onehot(e_map['is_conjugated'], bond.GetIsConjugated())

            edge_indices += [[i, j], [j, i]]
            edge_attrs += [e, e]

        edge_index = torch.tensor(edge_indices)
        edge_index = edge_index.t().to(torch.long).view(2, -1)
        edge_attr = torch.tensor(edge_attrs).to(torch.float)

        # Sort indices.
        if edge_index.numel() > 0:
            perm = (edge_index[0] * x.size(0) + edge_index[1]).argsort()
            edge_index, edge_attr = edge_index[:, perm], edge_attr[perm]

        data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, y=target,
                    smiles=smiles, poly_id=poly_id)
        return data
