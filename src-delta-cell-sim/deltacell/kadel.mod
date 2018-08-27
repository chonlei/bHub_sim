COMMENT

kadel.mod

Author: Linford Briant

Based on alpha-cell KA as in Watts and Sherman (2014).
	
ENDCOMMENT

NEURON {
	SUFFIX kadel
	USEION k READ ek WRITE ik
	RANGE gka,gbar,km,Vm,kh,Vh,taum,th1,th2
	GLOBAL minf,hinf,tauh
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
	Vm=-45			(mV)
	km=10			(mV)
	Vh=-68			(mV)
	kh=-10			(mV)
	taum=0.1		(ms)
	th1=60			(ms)
	th2=5			(ms)
} 

ASSIGNED {
	ik 		(mA/cm2)
	gk		(S/um2)
	ek		(mV)
	minf	(1)
	tauh	(ms)
	hinf	(1)
	gka
}
 

STATE {
	m
	h
}

INITIAL {
        rates(v)
        m=minf
        h=hinf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gka = gbar*m*h
	ik = gka*(v-ek)
}

FUNCTION alpm(v(mV)) {
  alpm = exp(-(v-Vm)/km)
}

FUNCTION alph(v(mV)) {
  alph = exp(-(v-Vh)/kh)
}

FUNCTION alptauh(v(mV)) {
  alptauh = (th1/(exp(-(v-5)/20)+exp((v-5)/20)))+th2
}

DERIVATIVE states { 
        rates(v)
        m' = (minf - m)/taum
        h' = (hinf - h)/tauh
}

PROCEDURE rates(v (mV)) { :callable from hoc
		tauh = alptauh(v)
        minf = 1/(1 + alpm(v))
		hinf = 1/(1 + alph(v))
}
