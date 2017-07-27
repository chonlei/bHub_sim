TITLE KCA
: IKCa current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX kca
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE  gbar, Vm, km, mtau, nkca, a, b, i
	GLOBAL minf, cinf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	ek		(mV)            : must be explicitly def. in hoc
	Vm 	(mV)
	km 	(mV)
	mtau 	(ms)
	nkca
	a 	(mV)
	b 	(mV)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
:	(mM) = (millimolars)
}

ASSIGNED {
	v 		(mV)
	cai		(mM)
	ik 		(mA/cm2)
	i 		(mA/cm2)
	g		(S/cm2) 	
	minf
	cinf		(mM)
}

STATE { m  c }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar * m * hill(cai, c, nkca)
	i = g * (v - ek)
	ik = i
}

INITIAL {
	trates(v)
	m = minf
	c = cinf
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
        c' = (cinf-c)/mtau
}

FUNCTION hill(x, xh, n) {
	hill = (x^n)/(x^n + xh^n)
}

PROCEDURE trates(v) {
	minf = 1/(1 + exp((v - Vm)/km))
	cinf = 1e-3 (mM) * exp((a - v)/b)	: 1e-3 uM -> mM
}
