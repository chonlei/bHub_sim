TITLE ECAUPDATE
: eca update from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX ecaup
	USEION ca READ cai,cao WRITE eca
}

PARAMETER {
	T
	zca
	CorrVCa		(mV)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(mM) = (milli/liter)
	FARADAY = (faraday) (coulomb)
}

ASSIGNED {
	cai		(mM)
	cao 		(mM)
	eca		(mV)
}
 
INITIAL {
	eca = geteca(cai,cao,T)
}

BREAKPOINT {
	eca = geteca(cai,cao,T)
}


FUNCTION geteca(cai,cao,T) {
        UNITSOFF
	geteca = 8.31*T/zca/(FARADAY)*log(cao/cai) * (1e3) - CorrVCa
	UNITSON
}

