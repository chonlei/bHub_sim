TITLE HERG
: IHERG current from Pedesen (2010)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX herg
	USEION k READ ek WRITE ik
	RANGE  gbar, Vm, km, mtau, Vh, kh, htau
	GLOBAL minf, hinf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	ek		(mV)            : must be explicitly def. in hoc
	Vm = -30	(mV)
	km = -10	(mV)
	mtau = 100	(ms)
	Vh = -42	(mV)
	kh = 17.5	(mV)
	htau = 50	(ms)
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
	g		(S/cm2) 	
	minf
	hinf
}

STATE { m  h }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*m*h
	ik = g * (v - ek)
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

PROCEDURE trates(v) {
	minf = 1/(1 + exp((v - Vm)/km))
	hinf = 1/(1 + exp((v - Vh)/kh))
}
