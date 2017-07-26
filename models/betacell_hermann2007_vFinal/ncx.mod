TITLE NCX
: INCX current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX ncx
	USEION na WRITE ina
	USEION ca READ cai WRITE ica
	RANGE  h, gbar, cncx, nncx, i, alpha, zca
}

PARAMETER {
	: must be explicitly def. in hoc
	gbar		(mA/cm2)	: should be ibar (in paper)
	cai		(mM)
	cncx		(mM)
	nncx 
	alpha 
	zca 
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
	ina 		(mA/cm2)
	ica 		(mA/cm2)
	i 		(mA/cm2)
	h
}


BREAKPOINT {
	h = hill(cai, cncx, nncx)
        i = gbar * h
	ina = alpha * i
	ica = (-1.0/3.0) * zca * alpha * i
}

FUNCTION hill(x, xh, n) {
	hill = (x^n)/(x^n + xh^n)
}
