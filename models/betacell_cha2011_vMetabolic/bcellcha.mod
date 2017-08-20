TITLE BIGFILE
: ALL current from Cha (2011)
: C.L.Lei Mar. 2017


NEURON {
    THREADSAFE
	SUFFIX bcellcha
	USEION na READ nai,ena WRITE ina
	USEION k READ ki,ek WRITE ik
	USEION ca READ cai,eca WRITE ica
	USEION Jcaer WRITE Jcaeri VALENCE 0 :: NOT REALLY ION - should use POINTER
        RANGE gammatoset, gammaapplytime, GKto, gKATP, PCaER, Pleak, KRe, Kfa, Pop, totalATP
}

PARAMETER {
	:******************Initial Conditions******************
	Voltage0
	Nai0
	Ki0
	Cai0
	CaER0
	ATP0
	MgADP0
	Re0
	d_CaL0
	U_CaL0
	fus0
	r_KDr0
	q_KDr0
	m_Kto0
	h_Kto0
	E1_tota0
	I10
	I20
	
	
	:******************MAIN PARAMS******************
        gammatoset
        gammaapplytime
	
	:******************ThermodynamicConstants******************
	Tem
	RTF
	RTF2
	
	:******************Cellular Parameters******************
	Naout=140.0
	Kout=5.4
	Caout=2.6
	Cm=6.158
	voli=764.0
	volER=280.0
	fi=0.01
	fer=0.025 
	
	:***********ICaV Parameters (ICaL)*************
	RCaLNa=0.0000185
	RCaLK=0.000367
	PCaL=48.9
	
	:***********IKslow Parameters*************
	PKslow=0.2
	nKslow=2.2
	KdKslow=0.00074
	
	:***********IKDr Parameters*************
	pKDr=2.1
	
	:***********ICapump Parameters*************
	P_PMCA=1.56
	K_PMCA=0.00014
	
	:***********ISOC Parameters*************
	PSOC=0.00764
	KCaer=0.003
	RNa_K_SOC=0.8
	pIbNSC=0.00396
	
	:***********ITRPM Parameters*************
	pTRPM=0.0234
	KTRPM=0.00076
	RNa_K_TRPM=0.8
	
	:***********IKATP Parameters*************
	kdd=0.01 
	ktt=0.05
	ktd=0.026
	gKATP:=2.31
	
	:***********IKto Parameters*************
	GKto:=2.13
	
	:***********Na/Ca exchange Parameters*************
	KdNao=87.5
	KdCao=1.38
	KdNai=20.75
	KdCai=0.0184
	k3=1.0
	k4=1.0
	AmpINaCa=204.0
	
	:***********Na/Kpump Parameters*************
	Pii=1.9
	Proton=0.0001
	Kd_MgATP=0.6
	Kd_Nao0=26.8
	Kd_Nai0=5.0
	Kd_Ko0=0.8
	Kd_Ki0=18.8
	delta_Nao=0.44
	delta_Nai=-0.14
	delta_Ko=0.23
	delta_Ki=-0.14
	k1_plus=1.253
	k2_plus=0.139
	k3_plus=6.96
	k4_plus=0.52 
	k1_minus=0.139
	k2_minus=0.0139
	k3_minus=13900.0
	k4_minus=0.348
	PNaK=350.0
	
	:**************Metabolism****************
	Nt=10.0
	totalATP:=4.0
	:ERcalciumdynamics(Jserca,Jout)
	PCaER:=0.096
	KCarp=0.0005
	
	Pleak:=0.46
	
	:GlycolysisAndOxidativephospholylation(ATP,MgADP,Re)
	KmATP=0.5
	hgl=2.5
	Kg=13.
	
	Pop:=0.0005
	Kop=0.02
	
	KRe:=0.000126
	Kfa:=0.0000063
	Stoichi=2.5
	Rvol=2.5
	kATPCa=0.187
	kATP=0.000062
	kADPf=0.0002
	kADPb=0.00002
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
	(mM) = (milli/liter)
	F = (faraday) (coulomb)
	R = (k-mole) (joule/degC)             ::: CHECK
}

ASSIGNED {
	v 		(mV)
	ina		(mA)
	ik		(mA)
	ica		(mA)
	nai		(mM)
	ki		(mM)
	cai		(mM)
	Jcaeri
	celsius		(degC)
        
        :*************Hub Param****************
	Glucose
	
	:******************Cellular Parameters******************
	Caer
	
	:**************Metabolism***************
	Jserca
	Jout
	fGlu
	JOP
	
	:*********ConstantField(v,NaCF,KCF,CaCF,ENa,EK,ECa)***********
	Denom1
	NaCF
	KCF
	Denom2
	CaCF
	
	:**********reversal potentials************
	ek
	ena
	eca
	IbNSC1
	IbNSC2
	IbNSC0
	dalpha
	dbeta
	VpOpen
	
	::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	:************calcium-dependent gate***********
	SingleiCaL
	Ualpha
	Ubeta
	
	:************ultra-slowgate***********
	usalpha
	usbeta
	
	:**********ICaL**********
	RundownATP
	pO
	ICaL1
	ICaL2
	ICaL3
	ICaL0
	
	:::calculateIKDr()
	alphar
	betar
	alphaq
	betaq
	IKDr2
	IKDr0
	
	:::calculateIKto()
	alpham
	betam
	alphah
	betah
	IKto2
	IKto0
	
	:::calculateITRPM()
	PoTRPM
	ITRPM1
	ITRPM2
	ITRPM0
	
	:::calculateIKATP()
	pOatp
	IKATP2
	IKATP0
	
	:::calculateINaK()
	fVm
	Kd_Nao
	Kd_Nai
	Kd_Ko
	Kd_Ki
	Nai_
	Naout_
	Ki_
	Kout_
	MgATP_
	a1_plus_
	a2_plus_
	a3_plus_
	a4_plus_
	a1_minus_
	a2_minus_
	a3_minus_
	a4_minus_
	denomi
	numer
	iglc
	vcyc
	INaK0
	INaK1
	INaK2
	
	:::calculateINaCa(IsFixSlow)
	pE1Na
	pE1Ca
	pE2Na
	pE2Ca
	k1
	k2
	fCa
	alpha1
	beta1
	alpha2
	beta2
	kf
	kb
	E2_tot
	INaCa0
	INaCa1
	INaCa3
	
	:::calculateIPMCA()
	IPMCA0
	IPMCA1
	IPMCA3
	
	:::calculateIKslow()
	PoKslow
	IKslow2
	IKslow0
	
	:::calculateISOC()
	PoSOC
	ISOC1
	ISOC2
	ISOC3
	ISOC0
	
	:::calculatedydt()
	Itot
	INatot
	IKtot
	ICatot

	ADPb
	
}
 

STATE { 
	ATP
	MgADP
	Re
	d_CaL
	U_CaL
	fus
	r_KDr
	q_KDr
	m_Kto
	h_Kto
	E1_tota
	I1
	I2 
}

BREAKPOINT {
        if (t > gammaapplytime) {Glucose = gammatoset}
        SOLVE states METHOD cnexp
	ina = INatot
	ik = IKtot
	ica = ICatot
}


INITIAL {
	ATP	=	ATP0
	MgADP	=	MgADP0
	Re	=	Re0
	d_CaL	=	d_CaL0
	U_CaL	=	U_CaL0
	fus	=	fus0
	r_KDr	=	r_KDr0
	q_KDr	=	q_KDr0
	m_Kto	=	m_Kto0
	h_Kto	=	h_Kto0
	E1_tota	=	E1_tota0
	I1	=	I10
	I2	=	I20
        trates(v)
}

DERIVATIVE states {   
        trates(v)
	
	ATP' = JOP-((INaK0+IPMCA0)/(F*1e-3)+Jserca/2.)/voli-(kATP+kATPCa*cai)*ATP
	MgADP' = -0.55*(JOP-((INaK0+IPMCA0)/(F*1e-3)+Jserca/2.)/voli-(kATP+kATPCa*cai)*ATP)+0.55*kADPb*ADPb-kADPf*MgADP
	
	Re' = (KRe*fGlu+Kfa)*(Nt-Re)-JOP*Rvol/Stoichi
	d_CaL' = dalpha*(1-d_CaL)-dbeta*d_CaL
	U_CaL' = Ualpha*(1-U_CaL)-Ubeta*U_CaL
	fus' = usalpha*(1-fus)-usbeta*fus
	r_KDr' = alphar*(1-r_KDr)-betar*r_KDr
	q_KDr' = alphaq*(1-q_KDr)-betaq*q_KDr
	m_Kto' = alpham*(1-m_Kto)-betam*m_Kto
	h_Kto' = alphah*(1-h_Kto)-betah*h_Kto
	E1_tota' = E2_tot*kf+I1*beta1+I2*beta2-E1_tota*(kb+alpha1+alpha2)
	I1' = E1_tota*alpha1-I1*beta1
	I2' = E1_tota*alpha2-I2*beta2
}


PROCEDURE trates(v) {
	:******************ThermodynamicConstants******************
	Tem=(celsius+273.15)
	RTF=R*Tem/(F*1e-3)
	RTF2=RTF/2

	:******************Cellular Parameters******************
	Caer=CaER0+(fer*voli/2./volER)*(Cm/(F*1e-3)/voli*(v-Voltage0)-(nai-Nai0)-(ki-Ki0)-2./fi*(cai-Cai0))
	
	:**************Metabolism****************
	Jserca=PCaER*cai^2/(cai^2+KCarp^2)
	Jout=Pleak*(Caer-cai)
	fGlu=ATP/(KmATP+ATP)*Glucose^hgl/(Kg^hgl+Glucose^hgl)
	JOP=Pop*Re*MgADP^2/(MgADP^2+Kop^2)
	
	:::Export to calculate dCa/dt
	Jcaeri = -Jserca+Jout
	
	:*********ConstantField(v,NaCF,KCF,CaCF,ENa,EK,ECa)***** ******
	Denom1=(1.-exp(-v/RTF))
	NaCF=(v/RTF)/Denom1*(nai-Naout*exp(-v/RTF))
	KCF=(v/RTF)/Denom1*(ki-Kout*exp(-v/RTF)) 
	Denom2=(1.-exp(-v/RTF2))
	CaCF=(v/RTF2)/Denom2*(cai-Caout*exp(-v/RTF2))
	
	:**********reversal potentials************
	:::calculateIbNSC()
	IbNSC1=pIbNSC*NaCF
	IbNSC2=0.01*KCF
	IbNSC0=IbNSC1+IbNSC2
	:::calculateICaL()
	dalpha=1./(0.88*exp(-(v-3.)/50.)+0.09*exp(-(v-3.)/600.))
	dbeta=1./(5.48*exp((v-3.)/12.)+1.245*exp((v-3.)/30.))
	VpOpen=d_CaL^2
	
	::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	:************calcium-dependent gate***********
	SingleiCaL=0.0676*CaCF
	Ualpha=0.0042*2.
	Ubeta=0.1159*(-1.15*SingleiCaL*VpOpen+cai)*2.
	
	:************ultra-slowgate***********
	usalpha=1./(75000.*exp(v/34.))
	usbeta=1./(5000.*exp(-v/19.)+500.*exp(-v/100.))
	
	:**********ICaL**********
	RundownATP=1./(1.+(1.4/ATP)^3)
	pO=(VpOpen*U_CaL*(0.4+0.6*fus))*RundownATP 
	ICaL1=RCaLNa*PCaL*pO*NaCF
	ICaL2=RCaLK*PCaL*pO*KCF
	ICaL3=PCaL*pO*CaCF
	ICaL0=ICaL1+ICaL2+ICaL3
	
	:::calculateIKDr()
	alphar=1.1/(25.*exp(-(v-3.)/8.)+1.*exp(-(v-3.)/100.))
	betar=1.1/(25.*exp(v/100.))
	alphaq=1./800.
	betaq=1./(1000.*exp(-v/8.)+100.*exp(-v/100.))
	IKDr2=pKDr*r_KDr*(0.6*q_KDr+0.4)*KCF
	IKDr0=IKDr2
	
	:::calculateIKto()
	alpham=0.4/(5.46*exp(-v/20.))
	betam=0.4/(2.48*exp(v/60.))
	alphah=1.7/(969.*exp(v/500.))
	betah=1.7/(13.2*exp(-v/9.)+6.93*exp(-v/1000.))
	IKto2=GKto*m_Kto*h_Kto*(v-ek)
	IKto0=IKto2
	
	:::calculateITRPM()
	PoTRPM=1./(1.+(KTRPM/cai)^1.7)
	ITRPM1=pTRPM*RNa_K_TRPM*NaCF*PoTRPM
	ITRPM2=pTRPM*KCF*PoTRPM 
	ITRPM0=ITRPM1+ITRPM2
	
	:::calculateIKATP()
	pOatp=(0.08*(1.+2.*MgADP/kdd)+0.89*(MgADP/kdd)^2)/(1.+MgADP/kdd)^2/(1.+0.45*MgADP/ktd+ATP/ktt)
	IKATP2=gKATP*pOatp*(v-ek)
	IKATP0=IKATP2
	
	:::calculateINaK()
	fVm=(F*1e-3)*v/(R*Tem)
	Kd_Nao=Kd_Nao0*exp(delta_Nao*fVm)
	Kd_Nai=Kd_Nai0*exp(delta_Nai*fVm)
	Kd_Ko=Kd_Ko0*exp(delta_Ko*fVm)
	Kd_Ki=Kd_Ki0*exp(delta_Ki*fVm)
	Nai_=nai/Kd_Nai
	Naout_=Naout/Kd_Nao
	Ki_=ki/Kd_Ki
	Kout_=Kout/Kd_Ko
	MgATP_=ATP/Kd_MgATP
	a1_plus_=(k1_plus*Nai_^3.0)/((1.+Nai_)^3+(1.+Ki_)^2-1.)
	a2_plus_=k2_plus
	a3_plus_=k3_plus*Kout_^2/((1.+Naout_)^3+(1.+Kout_)^2-1.)
	a4_plus_=k4_plus*MgATP_/(1+MgATP_)
	a1_minus_=k1_minus*MgADP
	a2_minus_=k2_minus*Naout_^3/((1.+Naout_)^3+(1.+Kout_)^2-1.)
	a3_minus_=k3_minus*Pii*Proton/(1.+MgATP_)
	a4_minus_=k4_minus*Ki_^2/((1.+Nai_)^3+(1.+Ki_)^2-1.)
	denomi=(a1_minus_+a1_plus_)*a2_minus_*a3_minus_+a1_plus_*a2_plus_*(a3_plus_+a3_minus_)+a2_plus_*a3_plus_*(a4_plus_+a4_minus_)+(a2_plus_+a2_minus_)*a3_minus_*a4_minus_+(a1_minus_+a1_plus_)*a3_plus_*a4_plus_+a1_minus_*(a3_plus_+a3_minus_)*a4_minus_+a1_plus_*(a2_plus_+a2_minus_)*a4_plus_+a1_minus_*a2_minus_*(a4_plus_+a4_minus_)
	numer=a1_plus_*a2_plus_*a3_plus_*a4_plus_-a1_minus_*a2_minus_*a3_minus_*a4_minus_
	iglc=(0.4+0.6*exp(-Glucose/5.84))
	vcyc=(numer/denomi)*iglc
	INaK0=PNaK*vcyc
	INaK1=3*INaK0
	INaK2=-2*INaK0
	
	:::calculateINaCa(IsFixSlow)
	pE1Na=1./(1+(KdNai/nai)^3*(1+cai/KdCai))
	pE1Ca=1./(1+(KdCai/cai)*(1+(nai/KdNai)^3))
	pE2Na=1./(1+(KdNao/Naout)^3*(1+Caout/KdCao))
	pE2Ca=1./(1+(KdCao/Caout)*(1+(Naout/KdNao)^3))
	k1=exp(0.32*v/RTF)
	k2=exp((0.32-1)*v/RTF)
	fCa=cai/(cai+0.004)
	alpha1=pE1Na*(fCa*0.002+(1-fCa)*0.0015)
	beta1=fCa*0.0012+(1-fCa)*0.0000005
	alpha2=fCa*0.00003+(1-fCa)*0.01
	beta2=fCa*0.09+(1-fCa)*0.0001
	kf=k2*pE2Na+k4*pE2Ca
	kb=k1*pE1Na+k3*pE1Ca 
	E2_tot=1-E1_tota-I1-I2
	INaCa0=AmpINaCa*(k1*pE1Na*E1_tota-k2*pE2Na*E2_tot)
	INaCa1=3*INaCa0
	INaCa3=-2*INaCa0
	
	:::calculateIPMCA()
	IPMCA0=P_PMCA*cai^2/(cai^2+K_PMCA^2)
	IPMCA1=-IPMCA0
	IPMCA3=2*IPMCA0
	
	:::calculateIKslow()
	PoKslow=1./(1.+(KdKslow/cai)^nKslow)
	IKslow2=PKslow*PoKslow*KCF
	IKslow0=IKslow2
	
	:::calculateISOC()
	PoSOC=1./(1.+exp((Caer-KCaer)/0.003))
	ISOC1=PSOC*RNa_K_SOC*PoSOC*NaCF
	ISOC2=PSOC*PoSOC*KCF
	ISOC3=PSOC*PoSOC*CaCF*20
	ISOC0=ISOC1+ISOC2+ISOC3
	
	:::calculatedydt()
	Itot=IbNSC0+IKDr0+IKto0+IKATP0+ITRPM0+ICaL0+INaK0+INaCa0+IPMCA0+IKslow0+ISOC0
	INatot=IbNSC1+ITRPM1+ICaL1+INaK1+INaCa1+IPMCA1+ISOC1 
	IKtot=IbNSC2+IKDr2+IKto2+IKATP2+ITRPM2+ICaL2+INaK2+IKslow2+ISOC2
	ICatot=ICaL3+INaCa3+IPMCA3+ISOC3

	ADPb=totalATP-ATP-MgADP/0.55
}

