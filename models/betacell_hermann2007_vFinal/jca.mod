TITLE JCA
: JCa current from Hermann (2007)
: C.L.Lei Feb. 2017

NEURON {
    THREADSAFE
	SUFFIX jca
	USEION ca READ ica,eca WRITE ica
	RANGE rho_CaL,rho_CaT,rho_NCX,rho_PMCA
}

PARAMETER {
	eca	(mV)
	
	I_NCX_bar (mA/cm2)
	C0	(mM)
	C_NCX	(mM)
	n_NCX
	V0		(mV)
	
	V_CaL		(mV)
	kappa_CaL		(mV)
	W_CaL	(mV)
	lamb_CaL	(mV)
	g_CaL_bar	(S/cm2)
	V_Ca_bar		(mV)
	C_CaL	(mM)
	n_CaL
	
	V_CaT	(mV)
	kappa_CaT	(mV)
	g_CaT_bar	(S/cm2)
	W_CaT	(mV)
	lamb_CaT	(mV)
	
	I_PMCA_bar (mA/cm2)
	C_PMCA	(mM)
	n_PMCA
	
	rho_CaL
	rho_CaT
	rho_NCX
	rho_PMCA
	z_Ca
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
	ica 		(mA/cm2)
	Iss_NCX	(mA/cm2)
	Iss_CaL	(mA/cm2)
	Iss_CaT	(mA/cm2)
	Iss_PMCA	(mA/cm2)
	JCA		(mA/cm2)	
}


BREAKPOINT {
	ica = -JCA
}

INITIAL {
	Iss_NCX = I_NCX_bar*hill(C0,C_NCX,n_NCX)
	Iss_CaL = sigmoid_act(V0,V_CaL,kappa_CaL)*sigmoid_inact(V0,W_CaL,lamb_CaL)*g_CaL_bar*(V0 - eca)*(1. - hill(C0,C_CaL,n_CaL))
	Iss_CaT = sigmoid_act(V0,V_CaT,kappa_CaT)*g_CaT_bar*(V0 - eca)*sigmoid_inact(V0,W_CaT,lamb_CaT)
	Iss_PMCA = I_PMCA_bar*hill(C0,C_PMCA,n_PMCA)
	
	JCA = rho_CaL*Iss_CaL + rho_CaT*Iss_CaT - z_Ca*rho_NCX*Iss_NCX + rho_PMCA*Iss_PMCA
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
