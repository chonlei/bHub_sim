TITLE CaT
: ICaT current from Pedesen (2010)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX cat
	USEION ca READ eca WRITE ica
	RANGE  gbar, Vm, km, Vh, kh, htau
	GLOBAL minf, hinf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	eca		(mV)            : must be explicitly def. in hoc
	Vm = -40	(mV)
	km = -4		(mV)
	Vh = -64	(mV)
	kh = 8		(mV)
	htau = 7	(ms) 	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

ASSIGNED {
	v 		(mV)
	ica 		(mA/cm2)
	g		(S/cm2) 		
	minf
	hinf
}
 

STATE { h }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*minf*h
	ica = g * (v - eca)
}

INITIAL {
	trates(v)
	h = hinf
}

DERIVATIVE states {   
        trates(v)
	h' = (hinf-h)/htau
}

PROCEDURE trates(v) {
	minf = 1/(1 + exp((v - Vm)/km))
	hinf = 1/(1 + exp((v - Vh)/kh))
}
