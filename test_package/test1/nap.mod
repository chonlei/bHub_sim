TITLE Na
: INa current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX nap
	USEION na READ ena WRITE ina
	:RANGE  gbar, Vm, km, mtau, mtau0, Vh, kh, htau
	:GLOBAL minf, hinf
}

PARAMETER {
	gbar		(S/cm2)	    : must be explicitly def. in hoc
	ena		(mV)            : must be explicitly def. in hoc
	Vm (mV)
	km 	(mV)
	mtau0  (ms)
	Vh 	(mV)
	kh 	(mV)
	htau	(ms) 
	V0	(mV)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

ASSIGNED {
	v 		(mV)
	ina 		(mA/cm2)
	g		(S/cm2) 		
	mtau
	minf
	hinf
}
 

STATE { m  h }

BREAKPOINT {
        SOLVE states METHOD cnexp
        g = gbar*m*h
	ina = g * (v - ena)
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
	mtau = mtau0/(exp((v - V0)/(40 (mV))) + exp((V0 - v)/(50 (mV))))
	minf = 1/(1 + exp((v - Vm)/km))
	hinf = 1/(1 + exp((v - Vh)/kh))
}
