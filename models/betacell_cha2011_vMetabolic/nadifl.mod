TITLE DNADT
: d[Na]/dt from Hermann (2007)
: C.L.Lei Mar. 2017

NEURON {
	SUFFIX nadifl
	USEION na READ ina WRITE nai
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(mM) = (milli/liter)
	F = (faraday) (coulomb)		::: REMEMBER to *1e-3
	R = (k-mole) (joule/degC) 
}

PARAMETER {
	voli
}

ASSIGNED {
	ina 	(mA/cm2)
}

STATE {
	nai	 (mM)
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
	~ nai << ( -1/voli*(ina)/(F*1e-3) )
}
