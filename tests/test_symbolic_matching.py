import sympy as sp
from collections import OrderedDict

from mod2nml import nml_helpers as nml
import mod2nml.symbolic as sym


def test_match_exp():
    beta = 0.125 * sp.exp(-(nml.v + 65) / 80)
    r = nml.Exponential.match(beta)
    assert r.rate == 0.125
    assert r.midpoint == -65
    assert r.scale == -80


def test_match_explinear():
    alpha = -0.01 * (nml.v + 60) / (sp.exp(-(nml.v + 60) / 10) - 1)
    r = nml.ExpLinear.match(alpha)
    assert r.rate == 0.1
    assert r.midpoint == -60
    assert r.scale == 10


def test_match_sigmoid():
    minf = 1 / (1 + sp.exp((-nml.v - 48) / 10))
    r = nml.Sigmoidal.match(minf)
    assert r.rate == 1
    assert r.midpoint == -48
    assert r.scale == 10


def test_match_standard_rates():
    for ex in [
        ".125*exp(-(v+65)/80)",
        "1 / ( 1 + exp( ( - v - 48 ) / 10 ) )",
        "-.01*(v+55)/(exp(-(v+55)/10)-1)",
    ]:
        expr = sp.parse_expr(ex, local_dict={"v": nml.v})
        matched = nml.match_standard_rates(expr).as_symbolic()
        assert sp.simplify(expr / matched) == 1


# mAlpha = (0.182 * (v- -32))/(1-(exp(-(v- -32)/6)))
# mBeta  = (0.124 * (-v -32))/(1-(exp(-(-v -32)/6)))
# hAlpha = (-0.015 * (v- -60))/(1-(exp((v- -60)/6)))
# hBeta  = (-0.015 * (-v -60))/(1-(exp((-v -60)/6)))
# minf  = 1 / ( 1 + exp( ( - v - 48 ) / 10 ) )


def test_match_simple_odes():
    n = sp.Symbol("n", real=True)  # state var

    dnti = sp.parse_expr("(ninf-n)/ntau", local_dict={"n": n})
    m = nml.match_alpha_beta_tau_inf(dnti, n)
    assert m["inf"].name == "ninf"
    assert m["tau"].name == "ntau"

    dnab = sp.parse_expr("(1-n)*alphan-n*betan", local_dict={"n": n})
    m = nml.match_alpha_beta_tau_inf(dnab, n)
    assert m["alpha"].name == "alphan"
    assert m["beta"].name == "betan"


def test_match_multiexpr_odes():
    n = sp.Symbol("n", real=True)  # state var
    ctxt = dict([("n", n), ("v", nml.v)])
    exprs = [
    ("alphan", "-.01*(v+55)/(exp(-(v+55)/10)-1)"),
    ("betan", ".125*exp(-(v+65)/80)"),
    ("sum", "alphan + betan"),
    ("taun", "1/sum"),
    ("ninf", "alphan/sum"),
    ("dn", "(ninf-n)/taun")]

    replaced, replacements = nml.replace_standards_in_sequence(exprs, ctxt)
    m = nml.match_alpha_beta_tau_inf(replaced["dn"], n)
    assert m["alpha"] == replaced["alphan"]
    assert m["beta"] == replaced["betan"]
    assert type(replacements["betan"][0]) == nml.Exponential
    assert type(replacements["alphan"][0]) == nml.ExpLinear

def test_q10():
    tadj = sym.match_q10('3^((celsius - 6.3)/10)', {})
    assert tadj['q10'] == 3
    assert tadj['temp'] == sp.S('celsius')
    assert tadj['exp_temp'] == 6.3

