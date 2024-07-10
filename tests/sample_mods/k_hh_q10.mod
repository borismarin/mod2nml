TITLE k_hh.mod   
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	    (S) = (siemens)
}

? interface
NEURON {
        SUFFIX k_hh
        USEION k READ ek WRITE ik REPRESENTS CHEBI:29103
        RANGE gkbar, gk
        GLOBAL ninf, ntau
	    THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
        gkbar = .036 (S/cm2)	<0,1e9>
}
 
STATE {
        n
}
 
ASSIGNED {
        v (mV)
        celsius (degC)
        ek (mV)

	    gk (S/cm2)
        ik (mA/cm2)
        ninf
	    ntau (ms)
}
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
        gk = gkbar*n*n*n*n
	    ik = gk*(v - ek)      
}
 
 
INITIAL {
	rates(v)
	n = ninf
}

? states
DERIVATIVE states {  
        rates(v)
        n' = (ninf-n)/ntau
}
 


? rates
PROCEDURE rates(v(mV)) {  
        LOCAL  alpha, beta, sum, q10

UNITSOFF
        q10 = 3^((celsius - 6.3)/10)

        :"n" potassium activation system
        alpha = .01*alpha_f(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
	    sum = alpha + beta
        ntau = 1/(q10 * sum)
        ninf = alpha/sum
}

FUNCTION alpha_f(x,y) {  
        alpha_f =  x/(exp(x/y) - 1)
}
 
 
UNITSON
