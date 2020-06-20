static double form_volume(double radius, double i_shell, double rg1, double poly_sig, double rg2, double v1, double v2, double I0) {
	double vol;
	double Ng;

	Ng = 4.00 * M_PI * pow((radius+i_shell)*0.1, 2.0) * poly_sig;
        vol = M_4PI_3 * cube(radius + i_shell) + Ng * (v1 + v2) + I0 * v2;
	vol = 1.00; // Complex structures can be normalized w/ scale.
	return vol;
}

static double radius_effective(int mode, double radius, double i_shell, double poly_sig, double rg1, double rg2, double v1, double v2, double I0) {
    switch(mode) {
	// Core radius
	case 1:
		return (radius + i_shell);
		break;

	// Outer radius
	case 2:
		return (radius + i_shell + rg1 + rg2);
		break;
    }
}

static double Iq(double q, double volf, double sld_c, double sld_s, double sld1, double sld2, double sld_solvent, double radius, double i_shell, double rc, double poly_sig, double rg1, double rg2, double nu1, double nu2, double v1, double v2, double I0, double rg3) {

	// Number of grafted chains.
	double Ng = 4.00 * M_PI * pow(0.1*(radius+i_shell), 2.0) * poly_sig;

	// Parameters for polymer form factors/amplitudes:
	double onu1, o2nu1, onu2, o2nu2, Usub1, Usub2, Usub3;

	// Misc terms:
	double r_coreshell, Fs, Fp1, Fp2, Pp1, Pp2, E1, E2;
	double vcore, vcoreshell, vtotal, inten, v3;

	// There are 9 terms in the intensity.
	double term1, term2, term3, term4, term5, term6, term7, term8, term9;

	// Exponents/Pre-factors for incomplete gamma function.
	onu1  = 1.0/nu1;
	onu2  = 1.0/nu2;
	o2nu1 = 1.0/2.0/nu1;
	o2nu2 = 1.0/2.0/nu2;
	Usub1 = (q * rg1) * (q * rg1) * (2.0*nu1 + 1.0) * (2.0*nu1 + 2.0) / 6.0;
	Usub2 = (q * rg2) * (q * rg2) * (2.0*nu2 + 1.0) * (2.0*nu2 + 2.0) / 6.0;
	Usub3 = (q * rg3) * (q * rg3) * (2.0*nu2 + 1.0) * (2.0*nu2 + 2.0) / 6.0;

	// Form factor amplitude for core:
	r_coreshell = radius + i_shell;
	vcore       = M_4PI_3 * cube(radius);
	vcoreshell  = M_4PI_3 * cube(r_coreshell);
	vtotal      = vcoreshell + Ng * (v1 + v2) + I0*v2;
	Fs = (sld_c - sld_s) * vcore * sas_3j1x_x(q*radius) + (sld_s - sld_solvent) * vcoreshell * sas_3j1x_x(q*r_coreshell);
	
	// Phase factors:
	E1 = sas_sinx_x(q*(r_coreshell));
	E2 = sas_sinx_x(q*(rc));

	// Form factor amplitudes for polymers:
	Fp1 = v1 * (sld1 - sld_solvent)*(o2nu1 * pow(Usub1, -o2nu1) * sas_gamma(o2nu1) * sas_gammainc(o2nu1, Usub1));
	Fp2 = v2 * (sld2 - sld_solvent)*(o2nu2 * pow(Usub2, -o2nu2) * sas_gamma(o2nu2) * sas_gammainc(o2nu2, Usub2));
	
	// Form factors for polymers:
	Pp1 = v1*v1*pow((sld1 - sld_solvent), 2.0) * (onu1 * pow(Usub1, -o2nu1)*sas_gamma(o2nu1)*sas_gammainc(o2nu1, Usub1) - onu1 * pow(Usub1, -onu1)*sas_gamma(onu1)*sas_gammainc(onu1, Usub1));
	Pp2 = v2*v2*pow((sld2 - sld_solvent), 2.0) * (onu2 * pow(Usub2, -o2nu2)*sas_gamma(o2nu2)*sas_gammainc(o2nu2, Usub2) - onu2 * pow(Usub2, -onu2)*sas_gamma(onu2)*sas_gammainc(onu2, Usub2));

	// Term 1: Nanoparticle Core
	term1 = Fs*Fs;

	// Term 2: Polymer Self Term
	// term2 = 2.0 * Ng * Fp1 * Fp2;
	term2 = 0.0;

	// Term 3: Polymer Block Self Term
	term3 = Ng * (Pp1 + Pp2);

	// Term 4: Block 1/Nanoparticle Crossterm
	term4 = 2.0 * Ng * Fs * E1 * Fp1;

	// Term 5: Block 2/Nanoparticle Crossterm
	term5 = 2.0 * Ng * Fs * E2 * Fp2;

	// Term 6: Block 1/Block 1 Crossterm
	term6 = Ng * (Ng - 1.0) * Fp1 * E1 * E1 * Fp1;

	// Term 7: Block 2/Block 2 Crossterm
	term7 = Ng * (Ng - 1.0) * Fp2 * E2 * E2 * Fp2;

	// Term 8: Block 2/Block 1 Crossterm
	 term8 = Ng * (Ng - 1.0) * Fp1 * E1 * E2 * Fp2;

	// Term 9: Free chains (if any)
	term9 = pow(sld2-sld_solvent, 2.0) * v2 * (onu2 * pow(Usub3, -o2nu2)*sas_gamma(o2nu2)*sas_gammainc(o2nu2, Usub3) - onu2 * pow(Usub3, -onu2)*sas_gamma(onu2)*sas_gammainc(onu2, Usub3));
	

	// Final intensity:
	inten = 1.0e-4 * volf * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8)/vtotal + I0*1.0e-4*term9;
	return inten;

 

}



/////
// 2D Version
/////
static double Iq_xy(double qx, double qy, double sld_c, double sld_s, double sld1, double sld2, double sld_solvent, double radius, double i_shell, double rc, double poly_sig, double rg1, double rg2, double nu1, double nu2, double v1, double v2, double I0, double rg3) {

	double q = sqrt(qx*qx + qy*qy);

	// Number of grafted chains.
	double Ng = 4.00 * 3.14159 * (radius * 0.1) * (radius * 0.1) * poly_sig;

	// Parameters for polymer form factors/amplitudes:
	double onu1, o2nu1, onu2, o2nu2, Usub1, Usub2, Usub3;

	// Misc terms:
	double r_coreshell, Fs, Fp1, Fp2, Pp1, Pp2, E1, E2;
	double vcore, vcoreshell, vtotal, inten;

	// There are 9 terms in the intensity.
	double term1, term2, term3, term4, term5, term6, term7, term8, term9;

	// Exponents/Pre-factors for incomplete gamma function.
	onu1  = 1.0/nu1;
	onu2  = 1.0/nu2;
	o2nu1 = 1.0/2.0/nu1;
	o2nu2 = 1.0/2.0/nu2;
	Usub1 = (q * rg1) * (q * rg1) * (2.0*nu1 + 1.0) * (2.0*nu1 + 2.0) / 6.0;
	Usub2 = (q * rg2) * (q * rg2) * (2.0*nu2 + 1.0) * (2.0*nu2 + 2.0) / 6.0;
	Usub3 = (q * rg3) * (q * rg3) * (2.0*nu2 + 1.0) * (2.0*nu2 + 2.0) / 6.0;

	// Form factor amplitude for core:
	r_coreshell = radius + i_shell;
	vcore       = M_4PI_3 * cube(radius);
	vcoreshell  = M_4PI_3 * cube(r_coreshell);
	vtotal      = vcoreshell + Ng * (v1 + v2);
	Fs = (sld_c - sld_s) * vcore * sas_3j1x_x(q*radius) + (sld_s - sld_solvent) * vcoreshell * sas_3j1x_x(q*r_coreshell);
	
	// Phase factors:
	E1 = sas_sinx_x(q*(r_coreshell+rg1));
	E2 = sas_sinx_x(q*(rc+rg2));

	// Form factor amplitudes for polymers:
	Fp1 = v1 * (sld1 - sld_solvent)*(o2nu1 * pow(Usub1, -o2nu1) * sas_gamma(o2nu1) * sas_gammainc(o2nu1, Usub1));
	Fp2 = v2 * (sld2 - sld_solvent)*(o2nu2 * pow(Usub2, -o2nu2) * sas_gamma(o2nu2) * sas_gammainc(o2nu2, Usub2));
	
	// Form factors for polymers:
	Pp1 = v1*v1*pow((sld1 - sld_solvent), 2.0) * (onu1 * pow(Usub1, -o2nu1)*sas_gamma(o2nu1)*sas_gammainc(o2nu1, Usub1) - onu1 * pow(Usub1, -onu1)*sas_gamma(onu1)*sas_gammainc(onu1, Usub1));
	Pp2 = v2*v2*pow((sld2 - sld_solvent), 2.0) * (onu2 * pow(Usub2, -o2nu2)*sas_gamma(o2nu2)*sas_gammainc(o2nu2, Usub2) - onu2 * pow(Usub2, -onu2)*sas_gamma(onu2)*sas_gammainc(onu2, Usub2));

	// Term 1: Nanoparticle Core
	term1 = Fs*Fs;

	// Term 2: Polymer Self Term
	// term2 = 2.0 * Ng * Fp1 * Fp2;
	term2 = 0.00;

	// Term 3: Polymer Block Self Term
	term3 = Ng * (Pp1 + Pp2);

	// Term 4: Block 1/Nanoparticle Crossterm
	term4 = 2.0 * Ng * Fs * E1 * Fp1;

	// Term 5: Block 2/Nanoparticle Crossterm
	term5 = 2.0 * Ng * Fs * E2 * Fp2;

	// Term 6: Block 1/Block 1 Crossterm
	term6 = Ng * (Ng - 1.0) * Fp1 * E1 * E1 * Fp1;

	// Term 7: Block 2/Block 2 Crossterm
	term7 = Ng * (Ng - 1.0) * Fp2 * E2 * E2 * Fp2;

	// Term 8: Block 2/Block 1 Crossterm
	 term8 = Ng * (Ng - 1.0) * Fp1 * E1 * E2 * Fp2;

	// Term 9: Free chains (if any)
	//term9 = pow(sld2-sld_solvent, 2.0) * I0 * I0 * (onu2 * pow(Usub3, -o2nu2)*sas_gamma(o2nu2)*sas_gammainc(o2nu2, Usub3) - onu2 * pow(Usub3, -onu2)*sas_gamma(onu2)*sas_gammainc(onu2, Usub3));
	term9 = 0.00;

	// Final intensity:
	inten = 1.0e-4 * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9)/vtotal;
	return inten;

 

}

