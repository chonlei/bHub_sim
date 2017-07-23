TITLE Ca
: ICa current from Jo et al (2005) and Sherman (1996)
: L.Briant Jan. 2017

NEURON {
    THREADSAFE
	SUFFIX cal
	USEION ca READ eca WRITE ica
	RANGE  gbar, Vm, km, mtau0
	GLOBAL minf
}

PARAMETER {
	gbar	(S/cm2)			: must be explicitly def. in hoc
	eca		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	Vm = -20	(mV)
	km = 12		(mV)
	mtau0 = 0.01	(ms) 	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

ASSIGNED {
	ica 		(mA/cm2)
	gca			(S/cm2) 		
	minf
}
 

STATE {m}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gca = gbar*m
		ica = gca * (v - eca)
}

INITIAL {
	trates(v)
	m=minf
}

DERIVATIVE states {   
        trates(v)
        m' = (minf-m)/mtau0
}

PROCEDURE trates(vm) {
	minf = 1/(1 + exp((Vm - v)/km))
}
