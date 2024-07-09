from collections import OrderedDict
import sympy as sp

v, rate, midpoint, scale = sp.symbols("v,rate,midpoint,scale", real=True)


def extract_exp_arg(expr, inv=False):
    ret = ()
    exps = expr.atoms(sp.exp)
    if exps:
        argexp = exps.pop().args[0]  # TODO: think about cases with >1 exp
        vv = (v - midpoint) / scale
        if inv:
            vv = -vv
        s = sp.solve(argexp - vv, midpoint, scale)
        V = sp.Symbol("V", real=True)
        cann = expr.subs(vv.subs(s), -V if inv else V)
        ret = cann, V, s
    return ret


def sympy_context_for_var(var, blk_vars, blk_assignments):
    # TODO: do that only for assignments that are relevant to var
    # (visit NAME nodes in its def)

    # locals / globals
    ctxt = OrderedDict([(v, sp.symbols(v, real=True))
                        for v in blk_vars
                        if v not in [v for v, k in blk_assignments]])

    for v, eq in blk_assignments:
        if v in ctxt:
            print('## Variable', v, 'rewritten') # TODO: create test for that
        ctxt[v] = sp.sympify(eq, ctxt)
        if v == var:
            break
    return ctxt


def match_alpha_beta_tau_inf(expr, statevar):
    n = statevar
    a = sp.Wild("a", exclude=[n])
    b = sp.Wild("b", exclude=[n])

    ddn = expr.diff(n)  # d(dn/dt)/dn = -1/tau = -1/(alpha+beta) for 1st order kinetics
    eqs = {}
    if m := (
        sp.solve(expr, n)[0].match(a / (a + b))
    ):  # if that matches, we have inf = alpha/(alpha + beta)
        eqs["alpha"] = m[a]
        eqs["beta"] = m[b]
    elif s := sp.solve(expr, n):
        eqs["tau"] = -1 / ddn
        eqs["inf"] = s[0]

    return eqs

def match_q10(expr):
    q10 = sp.Wild("q10")
    temp = sp.Wild("temp")
    exptemp = sp.Wild("exptemp")
    return sp.S(expr).match(q10 ** ((temp - exptemp)/10))

