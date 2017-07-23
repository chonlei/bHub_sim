TITLE CaL
: ICaL current from Pedesen (2010)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX cal
	USEION ca READ eca WRITE ica
	RANGE  gbar, Vm, km, Vh, kh, htau, Phi
	GLOBAL minf, hinf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	eca		(mV)            : must be explicitly def. in hoc
	Vm = -25	(mV)
	km = -6		(mV)
	Vh = -42	(mV)
	kh = 6		(mV)
	htau = 20	(ms) 	
	Phi = 57	(mV)
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
	: TBU cap hinf by [0,1]
	: hinf = 1/(1 + exp((v - Vh)/kh))??
	hinf = 1 + (minf*(v - eca)/Phi)
	if ( hinf > 1 ) {
		hinf = 1
	}
	if ( hinf < 0 ) {
		hinf = 0
	}
}
