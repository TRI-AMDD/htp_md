import argparse
import sys
import os
import os.path as osp
import csv

import numpy as np
import torch
import torch.nn.functional as F
from torch_geometric.data import DataLoader
from htpmd.ml_models.polymernet.data import PolymerDataset
from htpmd.ml_models.polymernet.model import PolymerNet
from torch_geometric.utils import softmax


PROP_CONFIGS = {
    'conductivity': {'mean': -4.262819, 'std': 0.222358, 'log10': True},
    'li_diffusivity': {'mean': -7.81389, 'std': 0.205920, 'log10': True},
    'poly_diffusivity': {'mean': -7.841585, 'std': 0.256285, 'log10': True},
    'tfsi_diffusivity': {'mean': -7.60879, 'std': 0.217374, 'log10': True},
    'molarity': {'mean': 1.411856, 'std': 0.0581892, 'log10': False},
    'transference_number': {'mean': 0.0623139, 'std': 0.281334, 'log10': False},
}


def predict(smiles_list, property):
    # setting for the pre-trained models
    has_H, form_ring = False, True
    log10 = PROP_CONFIGS[property]['log10']
    mean, std = PROP_CONFIGS[property]['mean'], PROP_CONFIGS[property]['std']
    fea_len, n_layers, n_h = 16, 4, 2

    pred_dataset = PolymerDataset(
        smiles_list, log10=log10, form_ring=form_ring, has_H=has_H)
    pred_loader = DataLoader(
        pred_dataset, batch_size=128, shuffle=False)

    data_example = pred_dataset[0]

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = PolymerNet(
        data_example.num_features, data_example.num_edge_features,
        fea_len, n_layers, n_h).to(device)

    # Predict on pred dataset
    cur_path = os.path.dirname(__file__)
    model.load_state_dict(torch.load(
        osp.join(cur_path, f'pre_trained_gnns/{property}.pth'),
        map_location='cpu'))
    model.eval()
    poly_ids = []
    preds = []
    targets = []
    smiles = []
    for data in pred_loader:
        data = data.to(device)
        pred = model(data)
        # De-normalize prediction
        pred = pred * std + mean
        preds.append(pred.cpu().detach().numpy())
        targets.append(data.y.cpu().detach().numpy())
        poly_ids += data.poly_id
        smiles += data.smiles
    preds = np.concatenate(preds)
    targets = np.concatenate(targets)

    if log10:
        preds = 10**preds
    return preds


if __name__ == '__main__':
    preds = predict(['CC(CCNCC(C)OC(=O)[Au])CCO[Cu]'], 'poly_diffusivity')
    print(preds)
