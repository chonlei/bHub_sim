TITLE KATP
: IKATP current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX katp
	USEION k READ ek WRITE ik
	RANGE  m,gbar, gamma, gammakatp, nkatp, mtau, minf,gammatoset
}

PARAMETER {
	: must be explicitly def. in hoc
	gbar		(S/cm2)	
	ek		(mV) 
	gammakatp (mM)
	nkatp  (mM)
	mtau (ms)
	gammatoset (mM)
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
	gamma	(mM)
	minf
}

STATE {m}

BREAKPOINT {
	: control glucose level here.
	if (t>3e3) {gamma = gammatoset}
	:if (t>5e4) {gamma = 1}
	
	SOLVE states METHOD cnexp
	ik = gbar * (1 - m) * (v - ek)
}

INITIAL {
	gamma = 1 :gammakatp
	trates(gamma)
	m = minf
}

DERIVATIVE states {   
        trates(gamma)
	m' = (minf-m)/mtau
}

PROCEDURE trates(gamma) {
	minf = 1/(1 + exp((gamma - gammakatp)/nkatp))
}
