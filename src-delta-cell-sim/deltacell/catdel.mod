COMMENT

catdel.mod

Author: Linford Briant

Based on alpha-cell CaT as in Watts and Sherman (2014).
	
ENDCOMMENT

NEURON {
    THREADSAFE
	SUFFIX catdel
	USEION ca READ eca WRITE ica
	RANGE  gbar, Vm, km, Vh, kh, th1, th2
	GLOBAL minf, hinf
}

PARAMETER {
	gbar	(S/cm2)			: must be explicitly def. in hoc
	eca		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	Vm = -49	(mV)
	km = 4		(mV)	
	Vh = -52	(mV)
	kh = -5		(mV)
	th1 = 20	(ms) 
	th2	= 5		(ms)
	tm1 = 15	(ms)
	tm2 = 0
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
	hinf
	tauh
	taum
}
 

STATE {m h}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gca = gbar*m*m*m*h
		ica = gca*(v - eca)
}

INITIAL {
	trates(v)
	m=minf
	h=hinf
}

DERIVATIVE states {   
        trates(v)
        m' = (minf-m)/taum
		h' = (hinf-h)/tauh
}

FUNCTION alptaum(v(mV)) {
  alptaum = (tm1/(exp(-(v+50)/12)+exp((v+50)/12)))+tm2
}

FUNCTION alptauh(v(mV)) {
  alptauh = (th1/(exp(-(v+50)/15)+exp((v+50)/15)))+th2
}

PROCEDURE trates(vm) {
	minf = 1/(1+exp(-(v-Vm)/km))
	hinf = 1/(1+exp(-(v-Vh)/kh))
	taum = alptaum(v)
	tauh = alptauh(v)
}
