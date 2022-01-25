import pytest
from htpmd.ml_models import random_forest


@pytest.mark.parametrize(
    'prop',
    [
        'conductivity',
        'li_diffusivity',
        'poly_diffusivity',
        'tfsi_diffusivity',
        'transference_number',
    ])
def test_gnn(prop):
    preds = random_forest.random_forests_prediction(
        ['CC(CCNCC(C)OC(=O)[Au])CCO[Cu]'], prop)
    assert preds is not None
