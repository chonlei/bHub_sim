TITLE K-DR
: IK current from Jo et al (2005) and Sherman (1996)
: L.Briant Jan 2017

NEURON {
    THREADSAFE
	SUFFIX kdr
	USEION k READ ek WRITE ik
	RANGE gbar, Vm, km, mtau0
	GLOBAL minf
}

PARAMETER {
	gbar	(S/cm2)		: must be explicitly def. in hoc
	ek		(mV)		: must be explicitly def. in hoc
	v 		(mV)
	Vm=-17	(mV)
	km=5.6	(mV)
	mtau0 = 11	(ms)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

ASSIGNED {
	ik 		(mA/cm2)
	minf
	mtau (ms)
}

STATE { m}

BREAKPOINT {
        SOLVE states METHOD cnexp
	ik = gbar*m*(v - ek)
} 

INITIAL {
	trates(v)
	m=minf  
}

DERIVATIVE states {   
        trates(v)      
        m' = (minf-m)/mtau0
}

PROCEDURE trates(v) {  
	minf = 1/(1 + exp((Vm-v)/km))
}
