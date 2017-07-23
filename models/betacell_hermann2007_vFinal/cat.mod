TITLE CaT
: ICaT current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX cat
	USEION ca READ eca,ica WRITE ica
	RANGE  g
	GLOBAL minf, hinf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	eca		(mV)            : must be explicitly def. in hoc
	Vm 	(mV)
	km	(mV)
	mtau	(ms) 
	Vh	(mV)
	kh 	(mV)
	htau	(ms) 
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
	hinf
}
 

STATE { m  h }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = m*h
	ica = gbar * g * (v - eca)
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
