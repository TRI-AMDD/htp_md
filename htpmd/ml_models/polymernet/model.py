import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU, GRU

from torch_geometric.nn import CGConv, GlobalAttention, NNConv, Set2Set


class PolymerNet(torch.nn.Module):
    """CGConv + Global attention pooling."""

    def __init__(self, node_in_len, edge_in_len, fea_len, n_layers, n_h):
        super(PolymerNet, self).__init__()
        self.node_embed = Linear(node_in_len, fea_len)
        self.edge_embed = Linear(edge_in_len, fea_len)
        self.cgconvs = nn.ModuleList([
            CGConv(fea_len, fea_len, aggr='mean', batch_norm=True)
            for _ in range(n_layers)])

        self.pool = GlobalAttention(
            gate_nn=Sequential(Linear(fea_len, fea_len), Linear(fea_len, 1)),
            nn=Sequential(Linear(fea_len, fea_len), Linear(fea_len, fea_len)))
        self.hs = nn.ModuleList(
            [Linear(fea_len, fea_len) for _ in range(n_h - 1)])
        self.out = Linear(fea_len, 1)

    def forward(self, data):
        out = F.leaky_relu(self.node_embed(data.x))
        edge_attr = F.leaky_relu(self.edge_embed(data.edge_attr))

        for cgconv in self.cgconvs:
            out = cgconv(out, data.edge_index, edge_attr)

        out = self.pool(out, data.batch)

        for hidden in self.hs:
            out = F.leaky_relu(hidden(out))
        out = self.out(out)
        return torch.squeeze(out, dim=-1)
