#include <complex.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimBlackHoleRingdown.h>
#include <lal/LALSimIMR.h>

#include <LALSimIMRPhenomHM.c>

int main(void) {

    // COMPLEX16FrequencySeries *htilde = NULL;
    COMPLEX16FrequencySeries *hptilde = NULL;
    COMPLEX16FrequencySeries *hctilde = NULL;

    // REAL8 m1_in = 80.0 * LAL_MSUN_SI;
    // REAL8 m2_in = 10.0 * LAL_MSUN_SI;
    REAL8 m1_in = 0.5 * LAL_MSUN_SI;
    REAL8 m2_in = 0.5 * LAL_MSUN_SI;
    REAL8 chi1z_in = 0.0;
    REAL8 chi2z_in = 0.0;
    REAL8 deltaF = 1.0;
    REAL8 f_min = 1.0;
    REAL8 f_max = 1000.0;
    REAL8 fRef_in = 0.0;
    REAL8 phi0 = 0.0;
    REAL8 inclination = 1.57;
    REAL8 distance = 1e6 * LAL_PC_SI;

    // int ret = XLALIMRPhenomHMMultiModeStrain(&hptilde, &hctilde, m1_in,m2_in,chi1z_in,chi2z_in,deltaF,f_min,f_max,fRef_in,phi0,inclination,distance);
    // LALDict *extraParams=NULL;
    // int ret = XLALSimIMRPhenomDGenerateFD(&htilde,phi0,fRef_in,deltaF,m1_in,m2_in,chi1z_in,chi2z_in,f_min,f_max,distance,extraParams );

    printf("testing main function:\n");
    int ret = XLALIMRPhenomHMMultiModeStrain(&hptilde, &hctilde, m1_in,m2_in,chi1z_in,chi2z_in,deltaF,f_min,f_max,fRef_in,phi0,inclination,distance);


    if (ret != XLAL_SUCCESS){
        printf("ERROR");
    }

    printf("hptilde->length = %i\n", hptilde->data->length);

    printf("hptilde->data->data[2000] = %.16g\n", creal((hptilde->data->data)[2000]));
    // printf("htilde->data->data[2000] = %.9g\n", creal((htilde->data->data)[2000]));

  //   // for(int i=0;i<hptilde->data->length;i++){
  //   //     printf("hptilde->data->data[%i] = %f\n", i, hptilde->data->data[i]);
  //   // }
  //
  //
  //   REAL8 a = 0.;
  //   REAL8 b = 0.;
  //   REAL8 fi = 0.;
  //   REAL8 fr = 0.;
  //   REAL8 f1 = 0.;
  //   REAL8 f2lm = 0.;
  //
  //   // ret = XLALIMRPhenomHMFreqDomainMapParams(&a, &b, &fi, &fr, &f1, &f2lm, 0.07, 2, 1, 0.25, 0, 0, 1);
  //   ret = XLALIMRPhenomHMFreqDomainMapParams(&a, &b, &fi, &fr, &f1, &f2lm, 0.07, 2, 1, 0.25, 0, 0, 0);
  //   if (ret != XLAL_SUCCESS){
  //       printf("ERROR");
  //   }
  //   printf("a = %f\n", a);
  //   printf("b = %f\n", b);
  //   printf("fi = %f\n", fi);
  //   printf("fr = %f\n", fr);
  //   printf("f1 = %f\n", f1);
  //   printf("f2lm = %f\n", f2lm);
  //
  //   printf("&a = %p\n", &a);
  //
  //   REAL8 Mf22 = XLALSimIMRPhenomHMFreqDomainMap(0.07, 2, 1, 0.25, 0, 0, 0);
  //   printf("Mf22 = %g\n", Mf22);
  //
  //
  //   // Test HMPhasePreComp struct is working
  //   INT4 ell = 4;
  //   INT4 mm = 3;
  //   REAL8 eta = 0.16;
  //   REAL8 chi1z = 0.0;
  //   REAL8 chi2z = 0.0;
  //   // Convention m1 >= m2
  //   // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
  //   REAL8 Seta = sqrt(1.0 - 4.0*eta);
  //   REAL8 m1 = 0.5 * (1.0 + Seta);
  //   REAL8 m2 = 0.5 * (1.0 - Seta);
  //
  //   const REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z);
  //   if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");
  //
  //
  //   HMPhasePreComp z;
  //   ret = XLALSimIMRPhenomHMPhasePreComp(&z, ell, mm, eta, chi1z, chi2z, finspin);
  //   if (ret != XLAL_SUCCESS){
  //       XLALPrintError("XLAL Error - XLALSimIMRPhenomHMPhasePreComp failed\n");
  //       XLAL_ERROR(XLAL_EDOM);
  //   }
  //
  // printf("\n\n\n\n");
  //
  // printf("z.ai = %f\n", z.ai);
  // printf("z.bi = %f\n", z.bi);
  // printf("z.a2lm = %f\n", z.a2lm);
  // printf("z.b2lm = %f\n", z.b2lm);
  // printf("z.ar = %f\n", z.ar);
  // printf("z.br = %f\n", z.br);
  // printf("z.fi = %f\n", z.fi);
  // printf("z.f2lm = %f\n", z.f2lm);
  // printf("z.fr = %f\n", z.fr);
  // printf("z.PhDBconst = %f\n", z.PhDBconst);
  // printf("z.PhDCconst = %f\n", z.PhDCconst);
  // printf("z.PhDDconst = %f\n", z.PhDDconst);
  // printf("z.PhDBAterm = %f\n", z.PhDBAterm);
  //
  //
  //   // Test XLALSimIMRPhenomHMPhasePreComp
  //
  //
  //   // Testing amplitude Opt function
  //   ell = 2;
  //   mm = 2;
  //
  //
  //   /*compute phenomD amp and phase coefficients*/
  //   IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1z, chi2z, finspin);
  //   if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
  //   AmpInsPrefactors amp_prefactors;
  //   int errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
  //   XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");
  //
  //   REAL8 Mf_wf = 0.01;
  //   REAL8 ampOpt = XLALSimIMRPhenomHMAmplitudeOpt( Mf_wf, eta, chi1z, chi2z, ell, mm,pAmp, &amp_prefactors );
  //
  //   REAL8 amp = XLALSimIMRPhenomHMAmplitude( Mf_wf, eta, chi1z, chi2z, ell, mm );
  //
  //
  //   printf("ampOpt = %.24f\n", ampOpt);
  //   printf("amp = %.24f\n", amp);
  //
  //
  //   // testing IMRPhenomHMSingleModehlmOpt
  //   // REAL8 Mfref = Mf_wf;
  //   // COMPLEX16 hlmopt = IMRPhenomHMSingleModehlmOpt(eta, chi1z, chi2z, ell, mm, Mf_wf, Mfref, phi0, &z, pAmp, &amp_prefactors);
  //   // COMPLEX16 hlmopt = IMRPhenomHMSingleModehlmOpt(eta, chi1z, chi2z, ell, mm, Mf_wf, Mfref, phi0, &z, pAmp);
  //   // printf("hlmopt  =  %.9f + i %.9f\n", creal(hlmopt), cimag(hlmopt));
  //
  //   // COMPLEX16 hlm = IMRPhenomHMSingleModehlm(eta, chi1z, chi2z, ell, mm, Mf_wf, Mfref, phi0, &z);
  //
  //   // printf("hlm  =  %.9f + i %.9f\n", creal(hlm), cimag(hlm));
  //


    printf("Hello World!\n");
    return 0;
}
