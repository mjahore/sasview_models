static double form_volume(double R, double rg1, double poly_sig, double rg2, double v1, double v2) {
	double vol;
	double Ng;

	vol = 1.00; // Complex structures can be normalized w/ scale.
	return vol;
}

static double radius_effective(int mode, double R, double poly_sig, double rg1, double rg2) {
    switch(mode) {
	// Core radius
	case 1:
		return R;
		break;

	// Outer radius
	case 2:
		return (R + rg1 + rg2);
		break;
    }
}


static double Iq(double q, double I0, double m, double sld_c, double sld1, double sld2, double sld_solvent, double R, double rc, double poly_sig, double rg1, double rg2, double nu1, double nu2, double v1, double v2) {

	// Q regions:
        double Q1, Q2, Q3;
	double P1, P2, P3;

	// Number of grafted chains.
	double Ng = 4.00 * M_PI * pow(0.1*R, 2.0) * poly_sig;

	// Parameters for polymer form factors/amplitudes:
	double onu1, o2nu1, onu2, o2nu2, Usub1, Usub2;
	double vc = M_4PI_3 * pow(R, 3.0);

	// Misc terms:
	double Fs, Fp1, Fp2, Pp1, Pp2, E1, E2;
	double inten;

	// There are 9 terms in the intensity.
	double term1, term2, term3, term4, term5, term6, term7, term8, term9;

	// Calculate Q's:
        Q1 = 1.0/R * sqrt(5.0*m/2.0);
        Q2 = 1.0/R * sqrt(3.0*m/4.0);
        Q3 = 1.0/rc * sqrt(3.0*m/4.0);
        P1 =      exp(-pow(Q1,2.0)*pow(R,  2.0)/5) * pow(Q1, m);
        P2 =      exp(-pow(Q2,2.0)*pow(R,  2.0)/6) * pow(Q2, 0.25*m);
	P3 =      exp(-pow(Q3,2.0)*pow(rc, 2.0)/6) * pow(Q3, 0.25*m);

	// Exponents/Pre-factors for incomplete gamma function.
	onu1  = 1.0/nu1;
	onu2  = 1.0/nu2;
	o2nu1 = 1.0/2.0/nu1;
	o2nu2 = 1.0/2.0/nu2;
	Usub1 = (q * rg1) * (q * rg1) * (2.0*nu1 + 1.0) * (2.0*nu1 + 2.0) / 6.0;
	Usub2 = (q * rg2) * (q * rg2) * (2.0*nu2 + 1.0) * (2.0*nu2 + 2.0) / 6.0;

	// Form factor amplitude for core:
	if (q < Q1) {
		Fs = (sld_c - sld_solvent) * vc * exp(-pow(q, 2.0)*pow(R, 2.0)/10.0);
	} else {
		Fs = (sld_c - sld_solvent) * vc * P1 * pow(q, -0.5*m);
	}
	
	// Core propagators:
	if (q < Q2) {
		E1 = exp(-pow(q, 2.0)*pow(R, 2.0)/6.0);
	} else {
		E1 = P2 * pow(q, -0.25*m);
	}

	if (q < Q3) {
		E2 = exp(-pow(q, 2.0)*pow(rc, 2.0)/6.0);
	} else {
		E2 = P3 * pow(q, -0.25*m);
	}

	// Form factor amplitudes for polymers:
	Fp1 = v1 * (sld1 - sld_solvent)*(o2nu1 * pow(Usub1, -o2nu1) * sas_gamma(o2nu1) * sas_gammainc(o2nu1, Usub1));
	Fp2 = v2 * (sld2 - sld_solvent)*(o2nu2 * pow(Usub2, -o2nu2) * sas_gamma(o2nu2) * sas_gammainc(o2nu2, Usub2));
	
	// Form factors for polymers:
	Pp1 = v1*v1*pow((sld1 - sld_solvent), 2.0) * (onu1 * pow(Usub1, -o2nu1)*sas_gamma(o2nu1)*sas_gammainc(o2nu1, Usub1) - onu1 * pow(Usub1, -onu1)*sas_gamma(onu1)*sas_gammainc(onu1, Usub1));
	Pp2 = v2*v2*pow((sld2 - sld_solvent), 2.0) * (onu2 * pow(Usub2, -o2nu2)*sas_gamma(o2nu2)*sas_gammainc(o2nu2, Usub2) - onu2 * pow(Usub2, -onu2)*sas_gamma(onu2)*sas_gammainc(onu2, Usub2));

	// Term 1: Nanoparticle Core
	term1 = Fs*Fs;

	// Term 2: Polymer Block Self Term
	term2 = Ng * (Pp1 + Pp2);

	// Term 3: Block 1/Nanoparticle Crossterm
	term3 = 2.0 * Ng * Fs * E1 * Fp1;

	// Term 4: Block 2/Nanoparticle Crossterm
	term4 = 2.0 * Ng * Fs * E2 * Fp2;

	// Term 5: Block 1/Block 1 Crossterm
	term5 = Ng * (Ng - 1.0) * Fp1 * E1 * E1 * Fp1;

	// Term 6: Block 2/Block 2 Crossterm
	term6 = Ng * (Ng - 1.0) * Fp2 * E2 * E2 * Fp2;

	// Term 7: Block 2/Block 1 Crossterm
	term7 = Ng * Ng * Fp1 * E1 * E2 * Fp2;

	// Final intensity:
	inten = 1.0e-4 * (term1 + term2 + term3 + term4 + term5 + term6 + term7);
	return inten;

 

}

