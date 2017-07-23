TITLE KBK
: IKBK current from Pedesen (2010)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX kbk
	USEION ca READ ica
	USEION k READ ek WRITE ik
	RANGE  gbar, Vm, km, mtau, bbk
	GLOBAL minf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	ek		(mV)            : must be explicitly def. in hoc
	Vm = 0		(mV)
	km = -10	(mV)
	mtau = 2	(ms)
	bbk = 2e-2	(mA/cm2)	: TO BE CHECKED!!
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

ASSIGNED {
	v 		(mV)
	ik 		(mA/cm2)
	ica 		(mA/cm2)
	g		(S/cm2) 	
	minf
}

STATE { m }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*m
	ik = g * (v - ek) * ( bbk - ica )
}

INITIAL {
	trates(v)
	m = minf
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
}

PROCEDURE trates(v) {
	minf = 1/(1 + exp((v - Vm)/km))
}
