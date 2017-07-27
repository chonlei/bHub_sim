TITLE Kv
: IKv current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX kv
	USEION k READ ek WRITE ik
	RANGE  gbar, Vm, km, mtau, mtau0, Vh, kh, htau
	GLOBAL minf, hinf
}

PARAMETER {
	gbar		(S/cm2)		: must be explicitly def. in hoc
	ek		(mV)            : must be explicitly def. in hoc
	Vm	(mV)
	km 	(mV)
	mtau0 	(ms)		: Herrington et al. 2005 reported 10ms
	Vh 		(mV)
	kh 	(mV)
	htau 	(ms)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
:	(mM) = (millimolar)
}

ASSIGNED {
	v 		(mV)
	ik 		(mA/cm2)
	g		(S/cm2) 
	mtau		(ms)		
	minf
	hinf
}

STATE { m  h }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*m*h
	ik = g * (v - ek)
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
	mtau = mtau0/(exp((v - ek)/(65 (mV))) + exp((ek - v)/(20 (mV))))
	minf = 1/(1 + exp((v - Vm)/km))
	hinf = 1/(1 + exp((v - Vh)/kh))
}
