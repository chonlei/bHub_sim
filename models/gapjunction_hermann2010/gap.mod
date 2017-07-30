NEURON {
   POINT_PROCESS Gap
   POINTER vgap, naj, kj, caj
   RANGE g, tau
:   NONSPECIFIC_CURRENT i
   USEION na READ nai WRITE ina
   USEION k READ ki WRITE ik
   USEION ca READ cai WRITE ica
}

PARAMETER { 
   g = 1.0 (S/cm2) : g = 1/6*rho_gap*g_bar_gap
   tau = 500.0 (ms)
   celsius = 36.7842 (degC)
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
:	(mM) = (millimolars)
	F = (faraday) (coulomb)
	R = (k-mole) (joule/degC)             ::: CHECK
}

ASSIGNED {
   v (mV)
   vgap (mV)
   nai		(mM)
   ki		(mM)
   cai		(mM)
   ina 		(mA/cm2)
   ik 		(mA/cm2)
   ica 		(mA/cm2)
   vinfna	(mV)
   vinfk	(mV)
   vinfca	(mV)
}

STATE { vgapna  vgapk  vgapca }

BREAKPOINT { 
   SOLVE states METHOD cnexp
   ina = g * vgapna
   ik = g * vgapk
   ica = g * vgapca
   :i = (v - vgap)*g 
}

INITIAL {
   trates(v)
   vgapna = vinfna
   vgapk = vinfk
   vgapca = vinfca
}

DERIVATIVE states {   
   trates(v)      
   vgapna' = (vinfna-vgapna)/tau
   vgapk' = (vinfk-vgapk)/tau
   vgapca' = (vinfca-vgapca)/tau
}

PROCEDURE trates(v) {
	Tem=(celsius+273.15)
	RTF=R*Tem/(F*1e-3)
	RTF2=RTF/2
	vinfna = v - vgap - RTF*log(naj/nai)
	vinfk = v - vgap - RTF*log(kj/ki)
	vinfca = v - vgap - RTF2*log(caj/cai)
}
