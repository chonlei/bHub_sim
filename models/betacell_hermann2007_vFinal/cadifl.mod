TITLE DCADT
: d[Ca]/dt from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
	SUFFIX cadifl
	USEION ca READ ica,cai WRITE cai
	RANGE D, zca, xc, c0, kc
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	D = .6 (um2/ms)
	zca
	c0 (mM)
	kc (mM)
}

ASSIGNED {
	ica (milliamp/cm2)
	diam (um)
	xc
	eca (mV)
	temp (mV)
}

STATE {
	cai (mM)
}

BREAKPOINT {
	xc = c0*kc/(cai+kc)^2
	
	SOLVE conc METHOD sparse
}

KINETIC conc {
	~ cai << (-ica/(FARADAY)*3/diam*(1e4)/(1+xc)/zca)
}
