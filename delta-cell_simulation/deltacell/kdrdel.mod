COMMENT

kdrdel.mod

Author: Linford Briant

Based on alpha-cell KDR as in Watts and Sherman (2014).
	
ENDCOMMENT

NEURON {
	SUFFIX kdrdel
	USEION k READ ek WRITE ik
	RANGE m, gk, gbar, sm, vm, vtaum, staum, tm1, tm2
	GLOBAL minf,taum
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

PARAMETER {
	gbar = 0.03   	(S/cm2)
	v 				(mV)
	vm=-25			(mV)
	sm=23			(mV)
	vtaum=10		(mV)
	staum=25		(ms)
	tm1=1.5		(ms)
	tm2=15			(ms)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gk		(S/um2)
	ek		(mV)
	minf
	taum 	(ms)
}
 

STATE { m }

INITIAL { 
	rates(v)
	m = minf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
	gk = gbar*m*m*m*m
	ik = gk*(v-ek)
} 

FUNCTION alpm(v(mV)) {
  alpm = exp(-(v-vm)/sm)
}

FUNCTION alptaum(v(mV)) {
  alptaum = (tm1/(exp(-(v+vtaum)/staum)+exp((v+vtaum)/staum)))+tm2
}

DERIVATIVE states { 
        rates(v)
        m' = (minf - m)/taum
}

PROCEDURE rates(v (mV)) { :callable from hoc
		taum = alptaum(v)
        minf = 1/(1 + alpm(v))
}