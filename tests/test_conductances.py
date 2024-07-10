from collections import OrderedDict
import sympy as sp

from mod2nml import nmodl_nml as m2n


def test_find_current():
    mod = """
    NEURON {
        USEION na READ ena WRITE ina
        RANGE gna
    }
    BREAKPOINT {
        LOCAL drf
        drf = (v - ena)
        ina = gna*drf
    }
    """
    ina = m2n.find_currents(m2n.parse_mod(mod))['ina'][1]
    assert str(ina['ina']) == 'gna*(-ena + v)'

def test_match_gates():
    mod = """
    NEURON {
        USEION na READ ena WRITE ina
        RANGE gna
    }
    STATE { m h }
    BREAKPOINT {
        ina = gna*m*m*m*h*(v - ena)
    }
    """
    ast = m2n.parse_mod(mod)
    ina = m2n.find_currents(ast)['ina'][1]
    gna = m2n.conductance('ina', ina)
    assert str(gna) == 'gna*h*m**3'
    states = m2n.get_state_vars(ast)
    gates = m2n.match_cond_states(gna, states)
    gate_pow = {str(g): str(p) for g, p in gates.items()}
    assert gate_pow['m'] == 'm**3'
    assert gate_pow['h'] == 'h'

def test_gate_dynamics():
    mod = """
    NEURON {
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        RANGE gkbar, gk
        RANGE gna
        GLOBAL alpha, beta
    }
    PARAMETER { gkbar = .036 }
    STATE { n }
    BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkbar*n*n*n*n
	    ik = gk*(v - ek)
    }
    DERIVATIVE states {
        alpha = .01*-(v+55)/(exp(-(v+55)/10)-1)
        beta = .125*exp(-(v+65)/80)
        n' = alpha*(1-n) - beta*n
    }
    """
    ast = m2n.parse_mod(mod)
    states = m2n.get_state_vars(ast)
    dxs = m2n.get_gate_dynamics(ast)
    simp_dxs = m2n.simplify_derivatives(dxs,states)
    print(simp_dxs)


def test_gate_dynamics_tauinf():
    mod = """
    NEURON {
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        RANGE gkbar, gk
        RANGE gna
        GLOBAL inf, tau
    }
    PARAMETER { gkbar = .036 }
    STATE { n }
    BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkbar*n*n*n*n
	    ik = gk*(v - ek)
    }
    DERIVATIVE states {
        tau = 10
        inf = 1/(1+exp(-(v+65)/80))
        n' = (inf-n)/tau
    }
    """
    ast = m2n.parse_mod(mod)
    states = m2n.get_state_vars(ast)
    dxs = m2n.get_gate_dynamics(ast)
    simp_dxs = m2n.simplify_derivatives(dxs,states)
    print(simp_dxs)

def test_rates_in_function():
    with open("sample_mods/k_hh.mod") as f:
        mod = f.read()
    ast = m2n.parse_mod(mod)
    m2n.ast_inline_fold(ast)

    states = m2n.get_state_vars(ast)
    dxs = m2n.get_gate_dynamics(ast)
    simp_dxs = m2n.simplify_derivatives(dxs,states)
    fwd = simp_dxs['n'].forward[0]
    rev = simp_dxs['n'].reverse[0]

    # from modfile
    #    alpha = .01*alpha_f(-(v+55),10)
    #    beta = .125*exp(-(v+65)/80)

    assert isinstance(fwd, m2n.nml.ExpLinear)
    assert isinstance(rev, m2n.nml.Exponential)

    assert fwd.rate == 0.1
    assert fwd.midpoint == -55
    assert fwd.scale == 10

    assert rev.rate == 0.125
    assert rev.midpoint == -65
    assert rev.scale == -80


def test_reused_vars_funccall():
    with open("./sample_mods/na_hh.mod") as f:
        mod = f.read()
    ast = m2n.parse_mod(mod)
    m2n.ast_inline_fold(ast)
    m2n.analyse_currents(ast)


def test_multiple_conds():
    with open("./sample_mods/hh_clean.mod") as f:
        mod = f.read()
    ast = m2n.parse_mod(mod)
    m2n.ast_inline_fold(ast)
    for c in m2n.analyse_currents(ast):
        print(c)

def test_q10():
    with open("./sample_mods/k_hh_q10.mod") as f:
        mod = f.read()
    ast = m2n.parse_mod(mod)
    m2n.ast_inline_fold(ast)
    for c in m2n.analyse_currents(ast):
        print(c)

