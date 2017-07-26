TITLE DNADT
: d[Na]/dt from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
	SUFFIX nadifl
	USEION na READ ina,ena,nao WRITE nai
	RANGE D
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	D = .6 (um2/ms)
	ena (mV)
	nao (mM)
}

ASSIGNED {
	ina (milliamp/cm2)
	diam (um)
}

STATE {
	nai (mM)
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
	~ nai << (-ina/(FARADAY)*3/diam*(1e4))
}
