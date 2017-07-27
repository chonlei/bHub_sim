TITLE PMCA
: IPMCA current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX pmca
	USEION ca READ cai WRITE ica
	RANGE  h, gbar, cpmca, npmca
}

PARAMETER {
	: must be explicitly def. in hoc
	gbar		(mA/cm2)	: should be ibar (in paper)
	cai		(mM)
	cpmca		(mM)
	npmca 
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
	ica 		(mA/cm2)
	h
}


BREAKPOINT {
	h = hill(cai, cpmca, npmca)
        ica = gbar * h
}

FUNCTION hill(x, xh, n) {
	hill = (x^n)/(x^n + xh^n)
}
