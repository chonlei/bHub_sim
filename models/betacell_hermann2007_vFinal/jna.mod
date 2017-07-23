TITLE JNA
: JNa current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX jna
	USEION na READ ina, ena,nai,nao WRITE ina
	RANGE  JNA
}

PARAMETER {
	ena	(mV)
	
	V0	(mV)
	V_NaV	(mV)
	kappa_NaV	(mV)
	W_NaV	(mV)
	lamb_NaV	(mV)
	g_NaV_bar	(S/cm2)
	V_Na_bar	(mV)
	
	I_NaK_bar	(mA/cm2)
	K0
	K_NaK
	n_NaK
	N0
	K_NaK_til
	n_NaK_til
	
	I_NCX_bar	(mA/cm2)
	C0
	C_NCX
	n_NCX
	
	rho_NaV
	rho_NaK
	rho_NCX
	alpha_NCX
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(mM) = (milli/liter)
}

ASSIGNED {
	v 		(mV)
	ina 		(mA/cm2)
	JNA		(mA/cm2)	
	Iss_NaV	(mA/cm2)
	Iss_NaK	(mA/cm2)
	Iss_NCX	(mA/cm2)
nai (mM)
nao (mM)
}


BREAKPOINT {
	ina = -JNA
}

INITIAL {
	Iss_NaV = sigmoid_act(V0,V_NaV,kappa_NaV)*sigmoid_inact(V0,W_NaV,lamb_NaV)*g_NaV_bar*(V0 - ena)
	Iss_NaK = I_NaK_bar*(1. - hill(K0,K_NaK,n_NaK))*hill(N0,K_NaK_til,n_NaK_til)
	Iss_NCX = I_NCX_bar*hill(C0,C_NCX,n_NCX)
	
	JNA = rho_NaV*Iss_NaV + 3.*rho_NaK*Iss_NaK + rho_NCX*Iss_NCX*alpha_NCX
:printf("nxna: %f\n",rho_NCX*Iss_NCX*alpha_NCX)
}


FUNCTION hill(x, xh, n) {
	hill = (x^n)/(x^n + xh^n)
}

FUNCTION sigmoid_act(v1,vab,kab) {
	sigmoid_act = 1/(1 + exp((vab - v1)/kab))
}

FUNCTION sigmoid_inact(v1,vab,kab) {
	sigmoid_inact = 1/(1 + exp((v1 - vab)/kab))
}
