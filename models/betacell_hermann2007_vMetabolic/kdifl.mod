TITLE DKDT
: d[K]/dt from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
	SUFFIX kdifl
	USEION k READ ik WRITE ki
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
}

ASSIGNED {
	ik (milliamp/cm2)
	diam (um)
}

STATE {
	ki (mM)
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
	~ ki << (-ik/(FARADAY)*3/diam*(1e4))
}
