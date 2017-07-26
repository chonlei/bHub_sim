TITLE NaK
: INaK current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX nak
	USEION na READ nai WRITE ina
	USEION k READ ki WRITE ik
	RANGE  gbar:, knak, nnak, nanak, ntnak, i, alpha
}

PARAMETER {
	: must be explicitly def. in hoc
	gbar		(mA/cm2)	: should be ibar (in paper)
	nai		(mM)
	ki		(mM)
	knak 	(mM)
	nnak
	nanak 	 (mM)
	ntnak 
	alpha 
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
	ik 		(mA/cm2)
	i 		(mA/cm2)
}


BREAKPOINT {
        i = gbar * (1 - hill(ki, knak, nnak)) * hill(nai, nanak, ntnak)
	ina = 2 * alpha * i
	ik = (-4.0/3.0) * alpha * i
}

FUNCTION hill(x, xh, n) {
	hill = (x^n)/(x^n + xh^n)
}
