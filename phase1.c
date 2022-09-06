
#include<stdio.h>
# include <stdlib.h>
# include <math.h>
#include <string.h>

float objective_function(float x);
float df(float x);
float ddf(float x);
float Bounding_Phase(float a, float b);
float Newton_Rapson(float a, float b);

int main()
{	
	float a, b; /*Ranges of x*/
	printf("Enter\n a = lower limit of x\n b = upper limit on x\n");
	scanf("%f %f",&a,&b);
	
	Bounding_Phase(a,b);

}



float objective_function(float x)
{
	
	//return (-(2*x-5)*(2*x-5)*(2*x-5)*(2*x-5) +(x*x-1)*(x*x-1)*(x*x-1));
    //return (-8-x*x*x+2*x+2*exp(x));
    //return (-4*x*sin(x));
    //return (2*(x-3)*(x-3)+exp(0.5*x*x));
    //return (x*x-10*exp(0.1*x));
    //return (-20*sin(x)+15*x*x);

}

float df(float x){
	
	
	//return (-8*(2*x-5)*(2*x-5)*(2*x-5)+6*x*(x*x-1)*(x*x-1));
    //return (-3*x*x+2+2*exp(x));
    //return (-4*sin(x)-4*x*cos(x));
    //return (4*(x-3)+x*exp(0.5*x*x));
    //return (2*x-exp(0.1*x));
    //return (-20*cos(x)+30*x);
}


float ddf(float x){
	
	
	//return (-48*(2*x-5)*(2*x-5)+ 6*(x*x-1)+24*x*x*(x*x-1));
    //return (-6*x+2*exp(x));
   // return (-8*cos(x)+4*x*sin(x));
    //return (4+exp(0.5*x*x)+x*x*exp(0.5*x*x));
    //return (2-(0.1)*exp(0.1*x));
    //return (20*sin(x)+30);
	
}




float Bounding_Phase(float a, float b)
{
	int i=1/*Number of iterations*/, feval,j=0;
	float delta;
	float x0, f_1, f0, f1,x[500],f[500];
	FILE *out;
	
	printf("Enter the initial guess=\t");
	scanf("%f",&x0);
	printf("Enter the increment ?=\t");
	scanf("%f",&delta);
	/*Step 1*/
	feval = 0; /*function evaluation*/
	f_1 = objective_function(x0-delta); /*Calculate objective_function*/
	f0 = objective_function(x0);
	f1 = objective_function(x0+delta);
	if(f_1>f0 && f0>f1) delta=delta;
	else if(f_1<f0 && f0<f1) delta=-delta;
	x[0]=x0;
	f[0]=f0;

	 
	feval = feval + 3;
	out = fopen("bounding_phase.out","w");/*Output file*/
	
	fprintf(out,"#It\tx\tf(x)\n");
	do{  j=j+1;
		fprintf(out,"%d\t%0.2f\t%0.2f\n",i,x[j-1],f[j-1]);
		x[j]=x[j-1]+pow(2,j-1)*delta;
		f[j]=objective_function(x[j]);
		
		if(f[j]>f[j-1])
		break;
		feval++;
		i = i + 1;
		
	
	}while( f[j]<f[j-1] || f[j]==f[j-1]);/*Step 3*/
	printf("The minimum point lies between (%f,%f)",x[j-2],x[j]);
	printf("\n#Total number of function evaluations: %d",feval);
	/*Store in the file*/
	fprintf(out,"\n#The minimum point lies between (%f,%f)",x[j-2],x[j]);
	fprintf(out,"\n#Total number of function evaluations: %d",feval);
	fclose(out);
	
	Newton_Rapson(x[j-2],x[j]);
	
	
	return 0;
}


float Newton_Rapson(float a, float b){
	
	
	//Step1
	float x[100];
	printf("\nEnter initial guess for newton rapson method between (%f,%f) : ",a,b);
	scanf("%f",&x[0]);
	float e;
	printf("Enter allowable error  : ");
	scanf("%f",&e);
	float df1,df2,f1;
	df1=df(x[0]);
	f1=df1;
	//Step2
	df2=ddf(x[0]);
	
	//Step3
	
	if(f1<0){
		f1=-1*f1;
	}
	int k=0;
	FILE *out;
	out = fopen("newton_rapson.out","w"); /*Output file*/
	fprintf(out,"#It\tx\tf(x)\tdf(x)\tddf(x)\n");	
	while(f1>e){
		
		k=k+1;
		x[k]=x[k-1]-(df1)/(df2);
		df1=df(x[k]);
		df2=ddf(x[k]);
		f1=df1;
		if(f1<0){
			f1=-1*f1;
		}
		
		fprintf(out,"%d\t%f\t%f\t%f\t%f\n",k,x[k],objective_function(x[k]),df1,df2);
		
	}
	
	
	printf("\n%f",x[k]);
	
	
	
	return 0;
}
