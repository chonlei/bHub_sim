NEURON {
   POINT_PROCESS Gap
   POINTER vgap, naj, kj, caj
   RANGE g, i
:   NONSPECIFIC_CURRENT i
   USEION na READ nai WRITE ina
   USEION k READ ki WRITE ik
   USEION ca READ cai WRITE ica
}

PARAMETER { g = 1.0 (microsiemens) }

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
:	(mM) = (millimolars)
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



:--------------------------
BREAKPOINT { i = (v - vgap)*g }
