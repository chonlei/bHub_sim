TITLE CaPQ
: ICaPQ current from Pedesen (2010)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX capq
	USEION ca READ eca WRITE ica
	RANGE  gbar, Vm, km
	GLOBAL minf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	eca		(mV)            : must be explicitly def. in hoc
	Vm = -10	(mV)
	km = -10	(mV)
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
}


BREAKPOINT {
        trates(v)
        g = gbar*minf
	ica = g * (v - eca)
}

INITIAL {
	trates(v)
}

PROCEDURE trates(v) {
	minf = 1/(1 + exp((v - Vm)/km))
}
