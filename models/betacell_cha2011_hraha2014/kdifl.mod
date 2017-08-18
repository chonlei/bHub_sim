TITLE DKDT
: d[K]/dt from Hermann (2007)
: C.L.Lei Mar. 2017

NEURON {
	SUFFIX kdifl
	USEION k READ ik WRITE ki
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(mM) = (milli/liter)
	F = (faraday) (coulomb)		::: REMEMBER to *1e-3
	R = (k-mole) (joule/degC)             ::: CHECK
}

PARAMETER {
	voli
}

ASSIGNED {
	ik 	(mA/cm2)
}

STATE {
	ki	 (mM)
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
	~ ki << ( -1/voli*(ik)/(F*1e-3) )
}
