#include <iostream>
#include <complex>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

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

#define NUM_OF_AVG 20
#define NUM_OF_MEM (5+NUM_OF_AVG)

double sample[NUM_OF_MEM];
float Fmeas_uniavg (float s_0, float fsample)   // 2.2 Frequency Measurement
{
  double K = 2.0;

  float LS1=0;
  float LS2=0;
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

int map_position(int rbp, int shift, int i){
	int map=rbp-(5-shift+i);
	if(map<0)map=map+NUM_OF_MEM;
//	printf("m(%i)",map);
	return map;
}


int fm(int v){
	if(v<0){
		v=v+NUM_OF_MEM;
	}
	return v;
}

int ringbuf_position=0;
float grid_sample[NUM_OF_MEM];
float Fmeas_ringbuffer(float s_0, float fsample){
	if(ringbuf_position>(NUM_OF_MEM-1)){
		ringbuf_position=0;
	}
	printf("f:%i\n",ringbuf_position);
	grid_sample[ringbuf_position]=s_0;


	double K = 2.0;

	float LS1=0;
	float LS2=0;
	double f;
	int i;
	
	for(i=0;i<NUM_OF_AVG;i++){
		LS1=LS1+(grid_sample[fm(-0+ringbuf_position-i)]*grid_sample[fm(-3+ringbuf_position-i)]-grid_sample[fm(-1+ringbuf_position-i)]*grid_sample[fm(-2+ringbuf_position-i)]);
		LS2=LS2+(grid_sample[fm(-0+ringbuf_position-i)]*grid_sample[fm(-5+ringbuf_position-i)]-grid_sample[fm(-1+ringbuf_position-i)]*grid_sample[fm(-4+ringbuf_position-i)]);
		printf("LS1:%3i*%3i+%3i*%3i\n",fm(ringbuf_position-i),fm(ringbuf_position-3-i),fm(ringbuf_position-1-i),fm(ringbuf_position-2-i));
		printf("LS2:%3i*%3i+%3i*%3i\n",fm(ringbuf_position-i),fm(ringbuf_position-5-i),fm(ringbuf_position-1-i),fm(ringbuf_position-4-i));
	}
	f = (fsample/(TWOPI * K)) *acos(LS2/(2*LS1));
	//fprintf(fp,"%10f %10.8f %10.8f %10.8f %10f %10f %10f %10f %10f %10f\n",f, (LS2/(2*LS1)), LS1, LS2, sample[0],sample[1],sample[2],sample[3],sample[4],amplitude);


	ringbuf_position=ringbuf_position+1;

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
    double fsignal = 50.515;
    double fnominal = 50.515;
    double fcut = 50.0;
    double f=0;

    W = TWOPI * fsignal / fsample; // Î©
    W_ref =  TWOPI * fnominal / fsample;
    
    rho = W - PI/2.0;

    a0 = tan(rho);
    a1 = 1.0 / cos(rho);


	//compare with averaging appraoch
	
    double amp=(rand()%26000)*0.01;
    
	printf("Sample  Measurement    Phase  Freq\n");
	
    for (int n = 0; n <= 2000; n++)
    {
        // Signal
        double x = amp*(sin(  n * W ));//+ 1000*(rand()%1000)*0.000000001);
        double xref = sin(  n * W_ref);
        
		f = Fmeas_ringbuffer(x,fsample);
		
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
        printf("%5i %10.2f   %10f   o:%-10f r:%-10f A:%f\n",n%NUM_OF_MEM, x, w_deg,f_unifilter,f,amplitude);

    }
    fclose(fp);
}

int main()
{
    //analyseFFT();

    analysePhasor();

    return 0;
}
