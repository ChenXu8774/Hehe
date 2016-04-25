double  JpsiMass =  3.096916;
double  PsiMass  =  3.68609;
double  MinMs    =  0.5;	//1.5	//2.0 
double  MaxMs    =  9.5;//7.5	//5.0 
double  NBins    =  180;	//120	//60. 
double  MsRange  =  9.0;	//6.0	//3.0     // MaxMs - MinMs
double  SqrtTwoPI=  2.50662827463100050241;

enum PAR {P0, P1, P2, P3, Njpsi, Mjpsi, Sjpsi, Npsi, Mpsi, Spsi};

gStyle->SetOptFit();

void fit()
{

	//use TChain to read in pDST
	TChain T("T");
	//T.Add("/direct/phenix+spin/phnxsp01/yuhw/taxi/Run13pp510Muon/3268/data/38677*.root");
	T.Add("/gpfs/mnt/gpfs02/phenix/spin3/spin2/chenxu/spinAnalyzer/crosscheck/Run15_dimuons.root");


	//extract mass spectrom using TTree::Draw()
	T.Draw("mass>>h","mass<5&&mass>2&&charge==0&&Tr0_pz*Tr1_pz>0");


	//set up fitting function
	TF1 *fdimu;
	fdimu = new TF1("fdimuExpoNorth", GausPoissonFunc, 2, 5, 8);
	fdimu->SetParNames("P0","P1","P2","P3","Njpsi","Mjpsi","Sjpsi",
			"Npsi","Mpsi","Spsi");

	fdimu->SetParameter(P0, 1000.);
	fdimu->SetParameter(P1, 1);
	fdimu->SetParameter(P2, 0.);
	fdimu->SetParameter(P3, 0.);

	fdimu->SetParameter(Njpsi, 10000.);
	fdimu->SetParameter(Mjpsi, 3.095);
	fdimu->SetParLimits(Mjpsi, 3.000,3.200);
	fdimu->SetParameter(Sjpsi, 0.1436);

	fdimu->SetParameter(Npsi, 89.);
	fdimu->SetParLimits(Npsi, 0., 9999999.);


	//fit !!!
	h->Fit(fdimu,"RES");
}

double GausPoissonFunc(double *m, double *par)
{
	double x = m[0];

	double p0 = par[P0];
	double p1 = par[P1];
	double p2 = par[P2];
	double p3 = par[P3];

	double N1 = par[Njpsi];
	double M1 = par[Mjpsi];
	double S1 = par[Sjpsi];

	double N2 = par[Npsi];
	double M2 = M1*(PsiMass/JpsiMass);
	double S2 = S1*(PsiMass/JpsiMass);
	//double S2 = S1*(M2/M1);

	double f =  p0*TMath::Poisson(x,p1) +
		(MsRange/NBins)*(N1/S1/SqrtTwoPI)*exp(-0.5*(x-M1)*(x-M1)/S1/S1)+
		(MsRange/NBins)*(N2/S2/SqrtTwoPI)*exp(-0.5*(x-M2)*(x-M2)/S2/S2);

	return f;
}
