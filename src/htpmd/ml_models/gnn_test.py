import pytest
from htpmd.ml_models import gnn


@pytest.mark.parametrize(
    'prop',
    [
        'conductivity',
        'li_diffusivity',
        'poly_diffusivity',
        'tfsi_diffusivity',
        'molarity',
        'transference_number',
    ])
def test_gnn(prop):
    preds = gnn.predict(['CC(CCNCC(C)OC(=O)[Au])CCO[Cu]'], prop)
    assert preds is not None
