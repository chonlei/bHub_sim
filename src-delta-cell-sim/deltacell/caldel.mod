COMMENT

caldel.mod

Author: Linford Briant

Based on alpha-cell CaL as in Watts and Sherman (2014).
	
ENDCOMMENT

NEURON {
    THREADSAFE
	SUFFIX caldel
	USEION ca READ eca WRITE ica
	RANGE  gbar, Vm, km, Vh, kh, th1, th2
	GLOBAL minf, hinf, taum, tauh
}

PARAMETER {
	gbar	(S/cm2)			: must be explicitly def. in hoc
	eca		(mV)            : must be explicitly def. in hoc
	v 		(mV)
	Vm = -30	(mV)
	km = 10		(mV)	
	Vh = -33	(mV)
	kh = -5		(mV)
	th1 = 10	(ms) 
	th2	= 0		(ms)
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
        gca = gbar*m*m*h
		ica = gca * (v - eca)
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
  alptaum = (1/(exp(-(v+23)/20)+exp((v+23)/20)))+0.05
}

FUNCTION alptauh(v(mV)) {
  alptauh = (th1/(exp(-(v+0)/20)+exp((v+0)/20)))+th2
}

PROCEDURE trates(vm) {
	minf = 1/(1+exp(-(v-Vm)/km))
	hinf = 1/(1+exp(-(v-Vh)/kh))
	taum = alptaum(v)
	tauh = alptauh(v)
}
