: Graded synaptic transmission

COMMENT
Model of SST-GIRK channel activation

This synaptic g does not saturate.
Therefore it does not have a gmax.  Instead the relevant parameter is k, 
the "transfer function slope", which is in units of (uS/mM3).

Modified from the model described by De Schutter et al. J. Neurophysiol 69:1225-35, 1993.

ENDCOMMENT

NEURON {
	POINT_PROCESS sst
	POINTER cp	: presynaptic [Ca] in (mM)
	RANGE e, k, g, i
	NONSPECIFIC_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (microsiemens)
		(molar) = (1/liter)
		(mM) = (millimolar)
}

PARAMETER {
	e = -75  (mV) : reversal potential
	k = 1e3  (uS/mM3)
}

ASSIGNED {
	v  (mV)
	cp (mM)
	g  (uS)
	i  (nA)
}

BREAKPOINT {
	g = k*cp^3
	i = g*(v - e)
}
