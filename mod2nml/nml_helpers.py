import sympy as sp

from .symbolic import *

import neuroml


class StdExpr:
    def __init__(self, r, m, s):
        self.rate = r
        self.midpoint = m
        self.scale = s
    def __repr__(self):
        return f'{type(self).__name__}[rate:{self.rate}, '\
                f'midpoint:{self.midpoint}, scale:{self.scale}]'

class Exponential(StdExpr):
    @classmethod
    def match(cls, expr):
        ret = None
        if c := extract_exp_arg(expr):
            cann, V, s = c
            (r,) = sp.solve(cann - rate * sp.exp(V), rate)
            if len(r.free_symbols) < 1:  # no spurious solutions found
                ret = cls(r, s[midpoint], s[scale])
        return ret
        # V = Wild('V', real=True)
        # r = Wild('r', real=True)
        # vv = (v-midpoint)/scale
        # ret = {}
        # if m := expr.match(r*exp(V)):
        #    sol = solve(m[V]-vv, midpoint, scale)
        #    ret['rate'] = m[r]
        #    ret['midpoint'] = sol[midpoint]
        #    ret['scale'] = sol[scale]
        # return ret

    def as_symbolic(self):
        x = (v - self.midpoint) / self.scale
        return self.rate * sp.exp(x)


class Sigmoidal(StdExpr):
    @classmethod
    def match(cls, expr):
        ret = None
        if c := extract_exp_arg(expr, inv=True):
            cann, V, s = c
            (r,) = sp.solve(cann - rate / (1 + sp.exp(-V)), rate)
            if len(r.free_symbols) < 1:  # no spurious solutions found
                ret = cls(r, s[midpoint], s[scale])
        return ret
        # V = Wild('V', real=True)
        # r = Wild('r', real=True)
        # vv = -(v-midpoint)/scale
        # ret = {}
        # if m := expr.match(r/(1+exp(V))):
        #    sol = solve(m[V]-vv, midpoint, scale)
        #    ret['rate'] = m[r]
        #    ret['midpoint'] = sol[midpoint]
        #    ret['scale'] = sol[scale]
        # return ret

    def as_symbolic(self):
        x = (v - self.midpoint) / self.scale
        return self.rate / (1 + sp.exp(-x))


class ExpLinear(StdExpr):
    @classmethod
    def match(cls, expr):
        ret = None
        if c := extract_exp_arg(expr, inv=True):
            cann, V, s = c
            vv = (v - s[midpoint]) / s[scale]
            try:
                (r,) = sp.solve(cann + vv.subs(s) * rate / (sp.exp(-V) - 1), rate)
            except ValueError: # TODO: better handling
                return None
            if len(r.free_symbols) < 1:  # no spurious solutions found
                ret = cls(r, s[midpoint], s[scale])
        return ret
        # V = Wild('V', real=True)
        # lin = Wild('lin', real=True)
        # vv = (v-midpoint)/scale
        # ret = {}
        # if m := expr.match(lin/(exp(V)-1)):
        #    sol = solve(m[V]-vv, midpoint, scale)
        #    m[lin] = -m[lin]
        # elif m := expr.match(lin/(1-exp(V))):
        #    sol = solve(m[V]-vv, midpoint, scale)
        # if m:
        #    sol[scale] = -sol[scale]
        #    ret['rate'] = solve(m[lin]/rate-vv.subs(sol),rate)[0]
        #    ret['midpoint'] = sol[midpoint]
        #    ret['scale'] = sol[scale]
        # return ret

    def as_symbolic(self):
        x = (v - self.midpoint) / self.scale
        return self.rate * x / (1 - sp.exp(-x))




def match_dynamics(exprs, statevar, replacements={}):
    n = statevar
    a = sp.Wild("a", exclude=[n])
    b = sp.Wild("b", exclude=[n])
    expr = exprs[statevar.name+'\\prime']
    ddn = expr.diff(n)  # d(dn/dt)/dn = -1/tau = -1/(alpha+beta) for 1st order kinetics

    ret = None
    if m := (
        sp.solve(expr, n)[0].match(a / (a + b))
    ):  # if that matches, we have inf = alpha/(alpha + beta)
        alpha = m[a]
        beta = m[b]
        for r,v in replacements.items():
            if exprs[r] == alpha:
                alpha = v
            if exprs[r] == beta:
                beta = v
        ret = GateHHrates(alpha, beta)
    elif s := sp.solve(expr, n):
        tau = -1/ddn
        inf = s[0]
        for r,v in replacements.items():
            if exprs[r] == tau:
                tau = v
            if exprs[r] == inf:
                inf = v
        ret = GateHHtauInf(tau, inf)

    return ret

class GateHHrates:
    def __init__(self, alpha, beta):
        self.forward = alpha
        self.reverse = beta
    def __repr__(self):
        return f'{type(self).__name__}[forward:{self.forward}, '\
                f'reverse:{self.reverse}]'

class GateHHtauInf:
    def __init__(self, tau, inf):
        self.time_course = tau
        self.steady_state = inf
    def __repr__(self):
        return f'{type(self).__name__}[time_course:{self.time_course}, '\
                f'steady_state:{self.steady_state}]'


def match_standard_rates(expr, subs={}):
    m = None
    for r in [cls for cls in StdExpr.__subclasses__()]:
        if m := r.match(expr):
            break
    return m

def pprint(x):
    ret = x
    if x.is_Float:
        if x == int(x):
            ret = sp.Integer(x)
        else:
            ret = f'{float(x):g}'
    return ret

std2nml_rates = {Exponential: 'HHExpRate', ExpLinear: 'HHExpLinearRate', Sigmoidal:'HHSigmoidRate'}
def std_rates(expression):
    try:
        expr, symbol = expression
    except TypeError:
        expr = expression
    match expr:
        case Exponential() | ExpLinear() | Sigmoidal():
            return neuroml.HHRate(type=std2nml_rates[type(expr)],
                                  rate=f'{pprint(expr.rate)}per_ms',
                                  midpoint=f'{pprint(expr.midpoint)}mV',
                                  scale=f'{pprint(expr.scale)}mV')
        case _:
            return neuroml.HHRate(type='Rate does not match standard form!')

def replace_reused_symbols(seq, ctxt):
    found = []
    replaced = []
    for name, expr in seq:
        if name in found:
            print('found reused name', name, 'with value', expr)
            replaced.append(name)
        found.append(name)
        syex = sp.S(expr, ctxt)
        ctxt[name] = syex
    return ctxt, replaced


def replace_standards_in_sequence(seq, ctxt):
    def absorb_scalar(name, expr): # TODO: factor out
        stds_to_reps = {v[1]:v[0] for v in replacements.values()}
        stds = [r[1] for r in replacements.values()]
        a = sp.Wild('a', properties=[lambda k: k.is_Number])
        s = sp.Wild('s', properties=[lambda k: k.is_Atom])
        if expr.has(*stds):
            if m := expr.match(a*s):
                std = stds_to_reps[m[s]]
                std.rate *= m[a]
                sym = next(syms)
                replacements[name] = (std, sym)
                ctxt[name] = sym
        return
    replacements = {}
    syms = sp.numbered_symbols("std",real=True)
    for name, expr in seq:
        syex = sp.S(expr, ctxt)
        if m := match_standard_rates(syex):
            syex = next(syms)
            replacements[name] = (m, syex)
        ctxt[name] = syex
        absorb_scalar(name, syex) #TODO: other optimizations possible here
    return ctxt, replacements


def match_hh_rates(conductance, gates_exponents, dynamics, q10s, nml_chan):
    for gating_var in (m for m in gates_exponents if m != 'gbar'):
        varname = gating_var.name
        ge = gates_exponents[gating_var]
        n_particles = ge.exp if isinstance(ge, sp.Pow) else 1

        if q10s[varname]:
            exptemp = f"{q10s[varname][1]['exp_temp']:.3g}degC"
            q10 = q10s[varname][1]['q10']
            q10el = neuroml.Q10Settings(experimental_temp=exptemp, q10_factor=q10)

        match dyn:= dynamics[varname]:
            case GateHHrates():
                gate = neuroml.GateHHRates(id=varname, instances=n_particles)
                gate.forward_rate = std_rates(dyn.forward)
                gate.reverse_rate = std_rates(dyn.reverse)
                if q10s[varname]:
                    gate.q10_settings = q10el
                nml_chan.gate_hh_rates.append(gate)
            case GateHHtauInf():
                gate = neuroml.GateHHTauInf(id=varname, instances=n_particles)
                gate.steady_state = std_rates(dyn.steady_state)
                gate.time_course = std_rates(dyn.time_course)
                if q10s[varname]:
                    gate.q10_settings = q10el
                nml_chan.gate_hh_tau_infs.append(gate)
            case _:
                gate = neuroml.Gate(id='Could not match gate dynamics to known forms!!')
                nml_chan.gates.append(gate)
