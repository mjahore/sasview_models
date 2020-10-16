static double form_volume(double radius, double i_shell, double poly_sig, double v, double I0) {
	return 1.00;
}

static double radius_effective(int mode, double radius, double i_shell, double poly_sig, double v, double I0) {
    switch(mode) {
	// Core radius
	case 1:
		return (radius + i_shell);
		break;

	// Outer radius
	case 2:
		return (radius + i_shell);
		break;
    }
}

static double Iq(double q, double volf, double sld_c, double sld_s, double sld1, double sld2, double sld_solvent, double radius, double i_shell, double poly_sig, double rc, double Cinfty, double M0, double M1, double M2, double nu1, double nu2, double v, double I0, double rg3, double nu3) {

	// Number of grafted chains.
	double Ng   = 4.00 * M_PI * pow(0.1*(radius+i_shell), 2.0) * poly_sig;

	// Bond angle (radians)
	double theta0 = 68.0 * M_PI/180.0;

	// Parameters for polymer form factors/amplitudes:
	double onu1, o2nu1, onu2, o2nu2, Usub1, Usub2, Usub3;
	double onu3, o2nu3;

	// Misc terms:
	double r_coreshell, Fs, Fp1, Fp2, Pp1, Pp2, E1, E2;
	double vcore, vcoreshell, vtotal, inten, v1, v2, v3;

	// There are 9 terms in the intensity.
	double term1, term2, term3, term4, term5, term6, term7, term8;

	// Calculate Kuhn Length
	double b = Cinfty * 1.54 / cos(theta0/2.0);

	// Calculate degree of polymerizations, given Cinfty:
	double N1 = (M1/M0) * pow(cos(theta0/2.0), 2) / Cinfty;
	double N2 = (M2/M0) * pow(cos(theta0/2.0), 2) / Cinfty;

	// Exponents/Pre-factors for incomplete gamma function.
	onu1  = 1.0/nu1;
	onu2  = 1.0/nu2;
	onu3  = 1.0/nu3;
	o2nu1 = 1.0/2.0/nu1;
	o2nu2 = 1.0/2.0/nu2;
	o2nu3 = 1.0/2.0/nu3;
	Usub1 = pow(q*b, 2.0) * pow(N1, 2.0*nu1) / 6.0;
	Usub2 = pow(q*b, 2.0) * pow(N2, 2.0*nu2) / 6.0;
	Usub3 = (q * rg3) * (q * rg3) * (2.0*nu3 + 1.0) * (2.0*nu3 + 2.0) / 6.0;

	// Form factor amplitude for core:
	v1 = N1 * v;
	v2 = N2 * v;

	// Volumes, etc.
	r_coreshell = radius + i_shell;
	vcore       = M_4PI_3 * cube(radius);
	vcoreshell  = M_4PI_3 * cube(r_coreshell);
	vtotal      = vcoreshell + Ng * (v1 + v2);
	Fs = (sld_c - sld_s) * vcore * sas_3j1x_x(q*radius) + (sld_s - sld_solvent) * vcoreshell * sas_3j1x_x(q*r_coreshell);
	
	// Phase factors:
	//
	// E1: Core phase factor:
	E1 = sas_sinx_x(q*r_coreshell);

	// E2: Diblock chain propagator:
	double rg1 = b * pow(N1, nu1)/sqrt((2.0*nu1+1.0)*(2.0*nu1+2.0));
	E2 = exp(-pow(q*rc, 2.0));

	// Form factor amplitudes for polymers:
	Fp1 = v1 * (sld1 - sld_solvent)*(o2nu1 * pow(Usub1, -o2nu1) * sas_gamma(o2nu1) * sas_gammainc(o2nu1, Usub1));
	Fp2 = v2 * (sld2 - sld_solvent)*(o2nu2 * pow(Usub2, -o2nu2) * sas_gamma(o2nu2) * sas_gammainc(o2nu2, Usub2));
	
	// Form factors for polymers:
	Pp1 = v1*v1*pow((sld1 - sld_solvent), 2.0) * (onu1 * pow(Usub1, -o2nu1)*sas_gamma(o2nu1)*sas_gammainc(o2nu1, Usub1) - onu1 * pow(Usub1, -onu1)*sas_gamma(onu1)*sas_gammainc(onu1, Usub1));
	Pp2 = v2*v2*pow((sld2 - sld_solvent), 2.0) * (onu2 * pow(Usub2, -o2nu2)*sas_gamma(o2nu2)*sas_gammainc(o2nu2, Usub2) - onu2 * pow(Usub2, -onu2)*sas_gamma(onu2)*sas_gammainc(onu2, Usub2));

	// Term 1: Nanoparticle Core Term
	term1 = Fs*Fs;

	// Term 2: Diblock Copolymer Term
	term2 =  Ng * (Pp1 + Pp2 + 2.0 * Fp1 * Fp2);

	// Term 3: Polymer Block A <--> Core Cross Term
	term3 = 2.0 * Ng * (Fs * E1 * Fp1);

	// Term 4: Polymer Block B <--> Core Cross Term
	term4 = 2.0 * Ng * (Fs * E1 * E2 * Fp2);

	// Term 5: Polymer Block A <--> Polymer Block A Cross Term
	term5 = Ng * (Ng - 1) * Fp1 * E1 * E1 * Fp1;

	// Term 6: Polymer Block A <--> Polymer Block B Cross Term
	term6 = Ng * Ng * Fp1 * E1 * E1 * E2 * Fp2;

	// Term 7: Polymer Block B <--> Polymer Block B Cross Term
	term7 = Ng * (Ng - 1.0) * Fp2 * E2 * E1 * E1 * E2 * Fp2;

	// Term 8: Free chains (if any)
	term8 = (onu3 * pow(Usub3, -o2nu3)*sas_gamma(o2nu3)*sas_gammainc(o2nu3, Usub3) - onu3 * pow(Usub3, -onu3)*sas_gamma(onu3)*sas_gammainc(onu3, Usub3));
	
	// Final intensity:
	inten = 1.0e-4 * volf * (term1 + term2 + term3 + term4 + term5 + term6 + term7)/vtotal + I0*1.0e-4*term8;
	return inten;

 

}
