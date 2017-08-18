TITLE DCADT
: d[Ca]/dt from Cha (2011)
: C.L.Lei Mar. 2017

NEURON {
	SUFFIX cadifl
	USEION ca READ ica WRITE cai
	USEION Jcaer READ Jcaeri VALENCE 0 :: NOT REALLY ION - should use POINTER
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
	fi
	voli
}

ASSIGNED {
	ica 	(mA/cm2)
	Jcaeri	(mA/cm2)
}

STATE {
	cai (mM)
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
	~ cai << ( fi*(-ica/(2.*(F*1e-3))+Jcaeri)/voli )
}
