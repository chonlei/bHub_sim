TITLE Kv
: IKv current from Pedesen (2010)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX kv
	USEION k READ ek WRITE ik
	RANGE  gbar, Vm, km, mtau, mtau0
	GLOBAL minf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	ek		(mV)            : must be explicitly def. in hoc
	Vm = 18		(mV)
	km = -10	(mV)
	mtau0 = 2	(ms)		: Herrington et al. 2005 reported 10ms
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
	mtau		(ms)		
	minf
}

STATE { m }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*m
	ik = g * (v - ek)
}

INITIAL {
	trates(v)
	m = minf
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau
}

PROCEDURE trates(v) {
	minf = 1/(1 + exp((v - Vm)/km))
	UNITSOFF
	if ( v > -26.6 ) {
		mtau = mtau0 + 10 * exp((-20 - v)/6)
	} else { 
		mtau = mtau0 + 30
	}
	UNITSON
}
