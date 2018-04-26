TITLE delta-cell na+ current E121114-1

COMMENT	
	Author: Linford Briant
ENDCOMMENT

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v 					(mV)
	ena 				(mV)
	gbar=0.012 			(S/cm2)
	V2m=-28.8749674		(mV)
	sm=-5.45289598		(mV)
	V2h=-45.388663		(mV)
	sh=4.99005762		(mV)
	f1=1				(1)
	a1taum=0.0001072	(s)
	b1taum=12			(mV)
	c1taum=40			(mV)
	a2taum=0.000152		(s)
	b2taum=-20			(mV)
	c2taum=22.69		(mV)
	a1tauh=0.0001636	(s)
	b1tauh=-15			(mV)
	c1tauh=8.6856		(mV)
	a2tauh=0.0001857	(s)
	b2tauh=5			(mV)
	c2tauh=35.35		(mV)
}

NEURON {
	SUFFIX nadel
	USEION na READ ena WRITE ina
        RANGE gbar,gna,V2m,sm,V2h,sh,f1
        GLOBAL minf,hinf,tauh,taum
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

ASSIGNED {
	ina (mA/cm2)
        minf
        hinf      
        tauh
        taum
        gna
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gbar*m*m*m*m*m*(f1*h)
	ina = gna*(v-ena)
}

FUNCTION alpm(v(mV)) {
  alpm = exp((v-V2m)/sm)
}

FUNCTION alph(v(mV)) {
  alph = exp((v-V2h)/sh)
}

FUNCTION alptauh(v(mV)) {
  alptauh=1e3*(a1tauh*exp(-((v-b1tauh)/c1tauh)^2)+a2tauh*exp(-((v-b2tauh)/c2tauh)^2))
}

FUNCTION alptaum(v(mV)) {
  alptaum=1e3*(a1taum*exp(-((v-b1taum)/c1taum)^2)+a2taum*exp(-((v-b2taum)/c2taum)^2))
}

DERIVATIVE states { 
        rates(v)
        m' = (minf - m)/taum
        h' = (hinf - h)/tauh
}

PROCEDURE rates(v (mV)) { :callable from hoc
		taum = alptaum(v)
		tauh = alptauh(v)
        minf = 1/(1 + alpm(v))
		hinf = 1/(1 + alph(v))
}
