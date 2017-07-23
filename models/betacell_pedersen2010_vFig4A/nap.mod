TITLE Na
: INa current from Pedesen (2010)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX nap
	USEION na READ ena WRITE ina
	RANGE  gbar, Vm, km, Vh, kh, htau
	GLOBAL minf, hinf
}

PARAMETER {
:CHECK UNIT OF gbar (done)
	gbar		(S/cm2)	    : must be explicitly def. in hoc
	ena		(mV)            : must be explicitly def. in hoc
	Vm = -18	(mV)
	km = -5		(mV)
	Vh = -42	(mV)
	kh = 6		(mV)
	htau = 2	(ms) 	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

ASSIGNED {
	v 		(mV)
	ina 		(mA/cm2)
	g		(S/cm2) 		
	minf
	hinf
}
 

STATE { h }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*minf*h
	ina = g * (v - ena)
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
