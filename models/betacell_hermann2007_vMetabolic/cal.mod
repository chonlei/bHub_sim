TITLE CaL
: ICaL current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX cal
	USEION ca READ eca,cai,ica,cao WRITE ica
	RANGE  g, gbar, Vm, km, mtau, Vh, kh, htau, ccal, ncal
	GLOBAL minf, hinf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	eca		(mV)            : must be explicitly def. in hoc
	Vm 	(mV)
	km		(mV)
	mtau	(ms) 
	Vh 	(mV)
	kh 	(mV)
	htau 	(ms) 
	ccal 	(mM) 
	ncal
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(mM) = (milli/liter)
}

ASSIGNED {
	v 		(mV)
	cai		(mM)
	ica 		(mA/cm2)
	g		(S/cm2) 		
	minf
	hinf
	cao		(mM)
}
 

STATE { m  h }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = m*h * (1 - hill(cai, ccal, ncal))
	ica = gbar* g * (v - eca)
}

INITIAL {
	trates(v)
	m = minf
	h = hinf
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau
}

FUNCTION hill(x, xh, n) {
	hill = (x^n)/(x^n + xh^n)
}


PROCEDURE trates(v) {
	minf = 1/(1 + exp((v - Vm)/km))
	hinf = 1/(1 + exp((v - Vh)/kh))
}
