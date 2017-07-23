TITLE KATP
: IKATP current from Pedesen (2010)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX katp
	USEION k READ ek WRITE ik
	RANGE  gbar
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	ek		(mV)            : must be explicitly def. in hoc
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
}


BREAKPOINT {
	ik = gbar * (v - ek)
}

