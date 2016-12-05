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

#define FM(v) ((v<0)? (v+NUM_OF_MEM): (v))

int ringbuf_position=0;
float grid_sample[NUM_OF_MEM];
float Fmeas_ringbuffer(float s_0, float fsample){
	ringbuf_position=ringbuf_position+1;
	if(ringbuf_position>(NUM_OF_MEM-1)){
		ringbuf_position=0;
	}
	//printf("f:%i\n",ringbuf_position);
	grid_sample[ringbuf_position]=s_0;
	//printf("Bufferwrite:%f\n",grid_sample[ringbuf_position]);

	double K = 2.0;

	float LS1=0;
	float LS2=0;
	double f;
	int i;
	int lowerend=-NUM_OF_AVG+ringbuf_position;
	for(i=ringbuf_position;i>(5);i--){
		LS1=LS1+(grid_sample[(-(5-0)+i)]*grid_sample[(-(5-3)+i)]-grid_sample[(-(5-1)+i)]*grid_sample[(-(5-2)+i)]);
		LS2=LS2+(grid_sample[(-(5-0)+i)]*grid_sample[(-(5-5)+i)]-grid_sample[(-(5-1)+i)]*grid_sample[(-(5-4)+i)]);
		//printf("   LS1:%3i*%3i+%3i*%3i\n",(i-5),(i-(5-3)),(+i-(5-1)),(-(5-2)+i));
		//printf("   LS2:%3i*%3i+%3i*%3i\n",(i-5),(i-(5-5)),(-(5-1)+i),(-(5-4)+i));
	}
	int fixcorner=lowerend;
	if(fixcorner<-1) fixcorner=-1;
	for(i;i>(fixcorner);i--){
		LS1=LS1+(grid_sample[FM(-(5-0)+i)]*grid_sample[FM(-(5-3)+i)]-grid_sample[FM(-(5-1)+i)]*grid_sample[FM(-(5-2)+i)]);
		LS2=LS2+(grid_sample[FM(-(5-0)+i)]*grid_sample[FM(-(5-5)+i)]-grid_sample[FM(-(5-1)+i)]*grid_sample[FM(-(5-4)+i)]);
//printf("fc:LS1:%3i*%3i+%3i*%3i\n",FM(i-5),FM(i-(5-3)),FM(+i-(5-1)),FM(-(5-2)+i));
		//printf("fc:LS2:%3i*%3i+%3i*%3i\n",FM(i-5),FM(i-(5-5)),FM(-(5-1)+i),FM(-(5-4)+i));
	}
	for(i;i>(lowerend);i--){
		LS1=LS1+(grid_sample[(NUM_OF_MEM-(5-0)+i)]*grid_sample[(NUM_OF_MEM-(5-3)+i)]-grid_sample[(NUM_OF_MEM-(5-1)+i)]*grid_sample[(NUM_OF_MEM-(5-2)+i)]);
		LS2=LS2+(grid_sample[(NUM_OF_MEM-(5-0)+i)]*grid_sample[(NUM_OF_MEM-(5-5)+i)]-grid_sample[(NUM_OF_MEM-(5-1)+i)]*grid_sample[(NUM_OF_MEM-(5-4)+i)]);
		//printf("c1:LS1:%3i*%3i+%3i*%3i\n",fm(i-5),fm(i-(5-3)),fm(+i-(5-1)),fm(-(5-2)+i));
		//printf("   LS1:%3i*%3i+%3i*%3i\n",NUM_OF_MEM+(i-5),NUM_OF_MEM+(i-(5-3)),NUM_OF_MEM+(+i-(5-1)),NUM_OF_MEM+(-(5-2)+i));
		//printf("c2:LS2:%3i*%3i+%3i*%3i\n",fm(i-5),fm(i-(5-5)),fm(-(5-1)+i),fm(-(5-4)+i));
		//printf("   LS2:%3i*%3i+%3i*%3i\n",NUM_OF_MEM+(i-5),NUM_OF_MEM+(i-(5-5)),NUM_OF_MEM+(-(5-1)+i),NUM_OF_MEM+(-(5-4)+i));	
	}
	f = (fsample/(TWOPI * K)) *acos(LS2/(2*LS1));
	//fprintf(fp,"%10f %10.8f %10.8f %10.8f %10f %10f %10f %10f %10f %10f\n",f, (LS2/(2*LS1)), LS1, LS2, sample[0],sample[1],sample[2],sample[3],sample[4],amplitude);


	return f;
  
}

float grid_rb_i;
float grid_rb_q;

void IQ_ringbuffer(float f,float fsample){
	float W = (TWOPI / fsample) * f; // Î©
	
	float rho = W - PI/2.0;

	a0 = tan(rho);
	a1 = 1.0 / cos(rho);

	float avg=NUM_OF_MEM-1;
	float I,Isq=0;
	float Q,Qsq=0;
	int lowerend=-NUM_OF_MEM+1+ringbuf_position;
	int i;
	float v0,v1;
	for(i=ringbuf_position;i>1;i--){
		v0=grid_sample[(i)];
		v1=grid_sample[(i-1)];
		I = v0;
		Isq = Isq + I*I;
		Q = a0 * (v0) + a1 * v1;
		Qsq = Qsq+Q*Q;
		//printf("fir:%3i(%f)*%3i(%f)\n",(i),v0,(i-1),v1);
	}
	int fixcorner=lowerend;
	if(fixcorner<-1) fixcorner=-1;
	for(i;i>(fixcorner);i--){
		v0=grid_sample[FM(i)];
		v1=grid_sample[FM(i-1)];
		I = v0;
		Isq = Isq + I*I;
		Q = a0 * (v0) + a1 * v1;
		Qsq = Qsq+Q*Q;
		//printf("fix:%3i(%f)*%3i(%f)\n",FM(i),v0,FM(i-1),v1);	
		//printf("I:%f|Q:%f\n",I,Q);
	}
	for(i;i>(lowerend);i--){
		v0=grid_sample[FM(i)];
		v1=grid_sample[FM(i-1)];
		I = v0;
		Isq = Isq + I*I;
		Q = a0 * (v0) + a1 * v1;
		Qsq = Qsq+Q*Q;
		//printf("sec:%3i(%f)*%3i(%f)\n",(i+NUM_OF_MEM),v0,(i-1+NUM_OF_MEM),v1);
	}
	grid_rb_i=Isq/avg;
	grid_rb_q=Qsq/avg;
}

float IQ_ringbuffer_getamplitude(void){
	float amplitude;
	amplitude=sqrt(grid_rb_i+grid_rb_q);
	return amplitude;
}

float IQ_get_phasor_oldphase;

float IQ_ringbuffer_getphasor(float f,float fsample){
	//x1: Signal x2:Reference
	//0=Q - 1=I
	//skalarprodukt =  x1[0] * x2[0]  +  x1[1] * x2[1];
	//I=AMP*cos(phi)
	//Q=AMP*sin(phi)
	//We fix phase at zero, thus I=1,Q=0
	float skalarprodukt=grid_sample[ringbuf_position];
	float normed= skalarprodukt /(IQ_ringbuffer_getamplitude());
	float phase_rad = acosf (normed);
	float v0,v1;

	v0=grid_sample[FM(ringbuf_position-0)];
	v1=grid_sample[FM(ringbuf_position-1)];
	int oldbiggerthennew=(v0<v1)?1:-1;
	printf("ObN:%i - normed: %f\n",oldbiggerthennew, normed);
	

	if(oldbiggerthennew<0){
		phase_rad=-phase_rad+2*M_PI;
	}

	if(fabs(normed)>(1-f/(2*fsample))) {
		float fix=IQ_get_phasor_oldphase;
		float estimatedphase=IQ_get_phasor_oldphase+(f/fsample)*2*PI;
		if(estimatedphase>TWOPI)estimatedphase=estimatedphase-TWOPI;
		printf("Measurement might be invalid. Phase estimated: %f\n",180*(estimatedphase/PI));
		if(fabs(estimatedphase-phase_rad)>(f/(4*fsample))*TWOPI){
			printf("use estimated value.\n");
			phase_rad=estimatedphase;
		}
	}
	if(fabs(normed)>1){
		float fix=IQ_get_phasor_oldphase;
		float estimatedphase=IQ_get_phasor_oldphase+(f/fsample)*2*PI;
		phase_rad=estimatedphase;
		printf("over unity fix");
	}
	
	IQ_get_phasor_oldphase=phase_rad;
	return phase_rad;
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

    double fsample = 1000.00;
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
	
    double amp=200;
    printf("Amplitude: %lf\n",amp);
	printf("Sample  Measurement    Phase  Freq\n");
	
	float iq_amp;
    for (int n = 0; n <= 2000; n++)
    {
        // Signal
        double x = amp*(cos(  n * W )+ 1000*(rand()%1000)*0.000000001);
        double xref = sin(  n * W_ref);
        
		f = Fmeas_ringbuffer(x,fsample);
		
		IQ_ringbuffer(f,fsample);
		iq_amp=IQ_ringbuffer_getamplitude();
       // Filter
    
        // Pseudo IQ
        IQ_gen(x);

        // Referenz
        REF_IQ_gen(xref);

        // Vector Operation
        vector_op();
        float w_rad = IQ_ringbuffer_getphasor(f,fsample);
        float w_deg = w_rad / PI * 180.0;

		double f_unifilter=Fmeas_uniavg(x,fsample);
        // Print Result(s)
        
        
        float truephase=((n*W) / PI) * 180.0;
        while(truephase>360) {truephase=truephase-360;}
        if(n>NUM_OF_MEM){
			printf("%5i v:%10.2lf r:%-10f A:%f ph:%f tph:%f°\n",n%NUM_OF_MEM,x, f,iq_amp,w_deg,truephase);
		}
    }
    fclose(fp);
}

int main()
{
    //analyseFFT();

    analysePhasor();

    return 0;
}
