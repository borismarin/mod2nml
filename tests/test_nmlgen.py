from collections import OrderedDict
import sympy as sp

from mod2nml import nmodl_nml as m2n

def test_compile_snippet():
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

def test_compile_file():
    target = open('./sample_mods/k_hh.nml').read()
    nmldoc = str(m2n.compile_mod("sample_mods/k_hh.mod"))
    compiled = nmldoc.replace('NeuroMLDocument', 'neuroml') #TODO: how ugly...

    print(compiled)
    print(target)
    from .xmlcomp import xml_compare
    assert xml_compare(target, compiled)
