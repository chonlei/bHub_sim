TITLE JK
: JK current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX jk
	USEION k READ ik,ek WRITE ik
	RANGE  JK
}

PARAMETER {
	ek	(mV)
	
	I_NaK_bar	(mA/cm2)
	K0
	K_NaK
	n_NaK
	N0
	K_NaK_til
	n_NaK_til
	
	ass = 45.	(mV)
	bss = 30.	(mV)
	gamma0	(mM)
	gamma_KATP	(mM)
	kappa_KATP	(mM)
	g_KATP_bar	(S/cm2)
	
	V0 	(mV)
	V_K_bar	(mV)
	V_KV	(mV)
	kappa_KV	(mV)
	W_KV	(mV)
	lamb_KV	(mV)
	g_KV_bar	(S/cm2)
	
	C0	(mM)
	C_sKCa	(mM)
	kappa_sKCa	(mM)
	g_sKCa_bar	(S/cm2)
	
	g_KCa_bar	(S/cm2)
	n_KCa
	V_KCa	(mM)
	kappa_KCa	(mM)
	
	rho_KATP
	rho_KV
	rho_sKCa
	rho_KCa
	rho_NaK
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
	ik 		(mA/cm2)
	JK		(mA/cm2)
	Iss_NaK	(mA/cm2)	
	Iss_KATP	(mA/cm2)
	Iss_KV	(mA/cm2)
	Iss_sKCa	(mA/cm2)
	Iss_KCa	(mA/cm2)
	C_KCa_V0 (mM)
}


BREAKPOINT {
	ik = -JK
}

INITIAL {
	Iss_NaK = I_NaK_bar*(1. - hill(K0,K_NaK,n_NaK))*hill(N0,K_NaK_til,n_NaK_til)
	Iss_KATP = (1. - sigmoid_act(gamma0,gamma_KATP,kappa_KATP))*g_KATP_bar*(V0 - ek)
	Iss_KV = sigmoid_act(V0,V_KV,kappa_KV)*sigmoid_inact(V0,W_KV,lamb_KV)*g_KV_bar*(V0 - ek)
	Iss_sKCa = sigmoid_act(C0,C_sKCa,kappa_sKCa)*g_sKCa_bar*(V0 - ek)
	C_KCa_V0 = 1e-3 (mM) * exp((ass - V0)/bss)
	Iss_KCa = g_KCa_bar*(V0 - ek)*hill(C0,C_KCa_V0,n_KCa)*sigmoid_act(V0,V_KCa,kappa_KCa)
	
	JK = rho_KATP*Iss_KATP + rho_KV*Iss_KV + rho_sKCa*Iss_sKCa + rho_KCa*Iss_KCa - 2.*rho_NaK*Iss_NaK
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
