#include <iostream>
#include <complex>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "fir1.h"

using namespace std;
#define PI	M_PI	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)

FILE *fp;
double W;
double rho;

double a0;
double a1;

vector<double> x1;
vector<double> x2;
vector<double> x3;


double W_ref;

double amplitude;

double    firlp0f1[]={9.3022e-05, 1.4470e-03, 4.8311e-03, 1.1993e-02, 2.4043e-02, 4.0934e-02,
   6.1231e-02,   8.2284e-02,   1.0076e-01,   1.1342e-01,   1.1793e-01,   1.1342e-01,
   1.0076e-01,   8.2284e-02,   6.1231e-02,   4.0934e-02,   2.4043e-02,   1.1993e-02,
   4.8311e-03,   1.4470e-03,   9.3022e-05};

double    firlp0f05[] ={ 0.0050906,   0.0071284,   0.0125997,   0.0216072,   0.0336952,   0.0478761,
						 0.0627470,   0.0766788,   0.0880505,   0.0954887,   0.0980756,   0.0954887,
						 0.0880505,   0.0766788,   0.0627470,   0.0478761,   0.0336952,   0.0216072,
						 0.0125997,   0.0071284,   0.0050906};

double firlpf03[]=  
 {-5.1311e-04,  -5.2000e-04,  -5.3406e-04,  -5.5439e-04,  -5.7943e-04,  -6.0696e-04,  -6.3414e-04,
  -6.5747e-04,  -6.7288e-04,  -6.7576e-04,  -6.6103e-04,  -6.2323e-04,  -5.5662e-04,  -4.5527e-04,
  -3.1317e-04,  -1.2437e-04,   1.1688e-04,   4.1606e-04,   7.7817e-04,   1.2077e-03,   1.7083e-03,
   2.2832e-03,   2.9342e-03,   3.6625e-03,   4.4681e-03,   5.3498e-03,   6.3055e-03,   7.3316e-03,
   8.4234e-03,   9.5751e-03,   1.0780e-02,   1.2029e-02,   1.3314e-02,   1.4625e-02,   1.5951e-02,
   1.7281e-02,   1.8603e-02,   1.9904e-02,   2.1172e-02,   2.2394e-02,   2.3560e-02,   2.4655e-02,
   2.5670e-02,   2.6593e-02,   2.7414e-02,   2.8125e-02,   2.8718e-02,   2.9186e-02,   2.9524e-02,
   2.9728e-02,   2.9796e-02,   2.9728e-02,   2.9524e-02,   2.9186e-02,   2.8718e-02,   2.8125e-02,
   2.7414e-02,   2.6593e-02,   2.5670e-02,   2.4655e-02,   2.3560e-02,   2.2394e-02,   2.1172e-02,
   1.9904e-02,   1.8603e-02,   1.7281e-02,   1.5951e-02,   1.4625e-02,   1.3314e-02,   1.2029e-02,
   1.0780e-02,   9.5751e-03,   8.4234e-03,   7.3316e-03,   6.3055e-03,   5.3498e-03,   4.4681e-03,
   3.6625e-03,   2.9342e-03,   2.2832e-03,   1.7083e-03,   1.2077e-03,   7.7817e-04,   4.1606e-04,
   1.1688e-04,  -1.2437e-04,  -3.1317e-04,  -4.5527e-04,  -5.5662e-04,  -6.2323e-04,  -6.6103e-04,
  -6.7576e-04,  -6.7288e-04,  -6.5747e-04,  -6.3414e-04,  -6.0696e-04,  -5.7943e-04,  -5.5439e-04,
  -5.3406e-04,  -5.2000e-04,  -5.1311e-04};
  
void IQ_gen (double x)          // 2.4 IQ-Signal Generation
{
    static double x_1 = 0.0;
    double I = x;
    double Q = a0 * x + a1 * x_1;
    x_1 = x;

    x1.clear();
    x1.push_back(I);
    x1.push_back(Q);
}
void REF_IQ_gen (double x)      // 2.4 IQ-Signal Generation
{
    static double x_1 = 0.0;
    double I = x;
    double Q = a0 * x + a1 * x_1;
    x_1 = x;

    x2.clear();
    x2.push_back(I);
    x2.push_back(Q);
}

void REF_gen (double x_ref)     // unused
{
    double I = sin(x_ref);
    double Q = -cos(x_ref);

    x2.clear();
    x2.push_back(I);
    x2.push_back(Q);
}

void vector_op()                //2.3 Vector Operation
{

    double betrag_x1     = sqrt( x1[0]*x1[0] + x1[1]*x1[1] );
    amplitude=betrag_x1;
    double betrag_x2     = sqrt( x2[0]*x2[0] + x2[1]*x2[1] );

    double skalarprodukt =  x1[0] * x2[0]  +  x1[1] * x2[1];
    double determinante  =  x1[1] * x2[0]  -  x1[0] * x2[1];
    double res           = /*acos*/ ( skalarprodukt /( betrag_x1 * betrag_x2) );
    double res2          = /*asin*/ ( determinante /( betrag_x1 * betrag_x2) );

    x3.clear();
    x3.push_back(res);
    x3.push_back(res2);

 }


double s_1 = 0.0;
double s_2 = 0.0;
double s_3 = 0.0;
double s_4 = 0.0;
double s_5 = 0.0;


double Fmeas (double s_0, double fsample)   // 2.2 Frequency Measurement
{
	double K = 2.0;

	double LS1; //= s_1 * s_2  -  s_0  * s_3;
	double LS2; // = s_0 * s_3  -  s_1  * s_4;
	
	LS1=s_1*s_4-s_2*s_3;
	LS2=s_0*s_5-s_1*s_4;
	
	double f = (fsample/(TWOPI * K)) *acos(LS2/(2*LS1));

	fprintf(fp,"%10f %10.8f %10.8f %10f %10f %10f %10f %10f\n",f, LS1, LS2, s_0,s_1,s_2,s_3,s_4);
	
	s_5 = s_4;
	s_4 = s_3;
	s_3 = s_2;
	s_2 = s_1;
	s_1 = s_0;

	return f;


}

#define NUM_OF_AVG 50
#define NUM_OF_MEM (5+NUM_OF_AVG)

double sample[NUM_OF_MEM];
double Fmeas_uniavg (double s_0, double fsample)   // 2.2 Frequency Measurement
{
  double K = 2.0;

  double LS1=0;
  double LS2=0;
  sample[0]=s_0;
  int i;
  for(i=0;i<NUM_OF_AVG;i++){
    LS1=LS1+(sample[0+i]*sample[3+i]-sample[1+i]*sample[2+i]);
    LS2=LS2+(sample[0+i]*sample[5+i]-sample[1+i]*sample[4+i]);
  }
  
  double f = (fsample/(TWOPI * K)) *acos(LS2/(2*LS1));

  //fprintf(fp,"%10f %10.8f %10.8f %10.8f %10f %10f %10f %10f %10f %10f\n",f, (LS2/(2*LS1)), LS1, LS2, sample[0],sample[1],sample[2],sample[3],sample[4],amplitude);
  
  for(i=(NUM_OF_MEM-1);i>0;i--){
    sample[i]=sample[i-1];
  }
  
  return f;


}




/***********************************************************************************************************************************\
              |--> DUT (filter,etc) --sin(nW + phi)-> PseudoIQ ---(sin(nW + phi),-cos(nW + phi))--->|-----------|
  sin(nW)---->|                                                                                     | Vector_op |----> phasor--> phi
              |-------------------------------------> PseudoIQ ---(sin(nW),-cos(nW))--------------->|-----------|
                                                     (Referenz)
\***********************************************************************************************************************************/



void analysePhasor()
{

	



    fp = fopen( "sniff.txt", "w+");

    double fsample = 1000.0;
    double fsignal = 50.5;
    double fnominal = 50.5;
    double fcut = 50.0;
    double f=0;

    W = TWOPI * fsignal / fsample; // Î©
    W_ref =  TWOPI * fnominal / fsample;
    
    rho = W - PI/2.0;

    a0 = tan(rho);
    a1 = 1.0 / cos(rho);


	//Init FIR
	Fir1 fir(firlpf03,101);
	fir.reset();


	//compare with averaging appraoch
	
    double amp=(rand()%26000)*0.01;
    
	printf("Sample  Measurement    Phase  Freq\n");
	
    for (int n = 0; n <= 2000; n++)
    {
        // Signal
        double x = amp*(sin(  n * W )+ 1000*(rand()%1000)*0.000000001);
        double xf =fir.filter(x);
        double xref = sin(  n * W_ref);
        
		f = Fmeas(xf,fsample);
		
       // Filter
    
        // Pseudo IQ
        IQ_gen(x);

        // Referenz
        REF_IQ_gen(xref);

        // Vector Operation
        vector_op();
        double w_rad = acos(x3[0]);
        double w_deg = w_rad / PI * 180.0;

		double f_unifilter=Fmeas_uniavg(x,fsample);
        // Print Result(s)
        printf("%5i %10f   %10f  %10f   %10f vs%10f A:%f\n",n, x,xf, w_deg,f_unifilter,f,amplitude);

    }
    fclose(fp);
}

int main()
{
    //analyseFFT();

    analysePhasor();

    return 0;
}
