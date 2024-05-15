TITLE hh_rates.mod   squid sodium, potassium, and leak channels -- rate constants and no q10
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	(S) = (siemens)
}

? interface
NEURON {
        SUFFIX hhtauinf
        REPRESENTS NCIT:C17145   : sodium channel
        REPRESENTS NCIT:C17008   : potassium channel
        USEION na READ ena WRITE ina REPRESENTS CHEBI:29101
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, gna, gk
        GLOBAL alpha_m, alpha_n, beta_m, beta_n, alpha_h, beta_h
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gnabar = .12 (S/cm2)	<0,1e9>
        gkbar = .036 (S/cm2)	<0,1e9>
        gl = .0003 (S/cm2)	<0,1e9>
        el = -54.3 (mV)
}
 
STATE {
        m h n
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ena (mV)
        ek (mV)

	gna (S/cm2)
	gk (S/cm2)
        ina (mA/cm2)
        ik (mA/cm2)
        il (mA/cm2)
        :minf hinf ninf
        alpha_m beta_m alpha_h beta_h alpha_n beta_n
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m*m*m*h
	ina = gna*(v - ena)
        gk = gkbar*n*n*n*n
	ik = gk*(v - ek)      
        il = gl*(v - el)
}
 
 
INITIAL {
	rates(v)
	m = alpha_m/(alpha_m+beta_m)
	h = alpha_h/(alpha_h+beta_h)
	n = alpha_n/(alpha_n+beta_n)
}

? states
DERIVATIVE states {  
        rates(v)
        m' = (1-m) * alpha_m - m * beta_m 
        h' = (1-h) * alpha_h - h * beta_h
        n' = (1-n) * alpha_n - n * beta_n
}
 


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
        LOCAL  alpha, beta, sum, q10

UNITSOFF
                :"m" sodium activation system
        alpha = .1 * vtrap(-(v+40),10)
        beta =  4 * exp(-(v+65)/18)
                :"h" sodium inactivation system
        alpha = .07 * exp(-(v+65)/20)
        beta = 1 / (exp(-(v+35)/10) + 1)
                :"n" potassium activation system
        alpha = .01*vtrap(-(v+55),10)
        beta = .125*exp(-(v+65)/80)
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON
