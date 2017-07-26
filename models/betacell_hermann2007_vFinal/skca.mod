TITLE SKCA
: IsKCa current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX skca
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE  gbar, Vm, km, mtau
	GLOBAL minf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	ek		(mV)            : must be explicitly def. in hoc
	Vm (mM)
	km (mM)
	mtau  (ms)
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
	ik 		(mA/cm2)
	g		(S/cm2) 	
	minf
}

STATE { m }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar * m
	ik = g * (v - ek)
}

INITIAL {
	trates(cai)
	m = minf
}

DERIVATIVE states {   
        trates(cai)      
        m' = (minf-m)/mtau
}

PROCEDURE trates(cai) {
	minf = 1/(1 + exp((cai - Vm)/km))
}
