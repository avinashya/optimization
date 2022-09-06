#include<stdio.h>
#include<conio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define epsilon pow(10, -3)
double Newton_Rapson(double z[100],double S[100],double a,int n);
double Bounding_Phase(double a[100], double S[100], int n);
double objective_function(double a[100],double S[100],double x, int n);
double df(double z[100],double S[100],double a,int n);
double ddf(double z[100],double S[100],double a,int n);
double diff_1(double c[100],int n);
double diff_2(double c[100],int n);
double func(double x[100],int n);
 
double R=0.1;
double pi=3.14159;

double* powell(double x[100],int n){
 	 double *p ,xp[200],alpha,d[200],dabs=0,t=0,dummy[200][1];
    double S[200][200];
    int i,j,Ns=0;
    
 	for(i=1;i<=n;i++)
    {
    for(j=1;j<=n;j++)
    {
    	
    if(i==j)
    S[j][i]=1;
    
    else
    S[j][i]=0;
    
    }
    }
    double dum[n+1], store[n+1][2];
    
    for(i=1;i<=n;i++)
    xp[i]=x[i];
    
	do{
		Ns++;
	   dabs=0;
	   
       for(i=1;i<=n;i++)
       {  
        dum[i]=S[i][1];
       }
       
        alpha= Bounding_Phase(xp,dum,n);
        
        for(i=1;i<=n;i++)
        {
        xp[i]= xp[i]+ alpha*S[i][1];
        store[i][1]=xp[i];
  
    }
        for(j=2;j<=n;j++)
        {
    	
         for(i=1;i<=n;i++)
		 {
		 dum[i]=S[i][j];
         store[i][2]=xp[i];
            }
         alpha= Bounding_Phase(xp,dum,n);
    
         for(i=1;i<=n;i++)
		 {
        xp[i]= store[i][2]+ alpha*S[i][j];
    
        }
        }
        for(i=1;i<=n;i++)
        store[i][2]=xp[i];
        
        for(i=1;i<=n;i++)
            dum[i]=S[i][1];
        alpha= Bounding_Phase(xp,dum,n);
        
        for(i=1;i<=n;i++){
        xp[i]= store[i][2]+ alpha*S[i][1];

        }
    for(i=1;i<=n;i++)
    d[i]= (xp[i]-store[i][1]);
   
    for(i=1;i<=n;i++)
    t = t + pow( (xp[i]-store[i][1]), 2);
   
    dabs=sqrt(t);
   
    for(j=n;j>=2;j--)   
    {
    for(i=1;i<=n;i++)
    S[i][j]=S[i][j-1];
    }
   
    for(i=1;i<=n;i++)
    S[i][1]=(d[i])/dabs;
   
    }while(dabs>epsilon && Ns<150);
   p=malloc(sizeof(double)*(n+10));
   
   for(i=1;i<=n;i++)
   p[i]=xp[i];
    
   return p;
 }
 
double objective_function(double a[100]/*the starting point*/,double S[100],double x, int n)
{ 
int i;
double c[n+1];
for(i=1;i<=n;i++)
{
    c[i]=a[i]+x*S[i];
}

	return func(c,n);
}

double Bounding_Phase(double a[100] , double S[100], int n)/*the starting point*/
{
     
	int i=1 , feval,j=0;/*Number of iterations*/
	double delta=0.1;
	double x0=0.2, f_1, f0, f1,x[50000],f[50000];
	 
	/*Step 1*/
	feval = 0; /*function evaluation*/
	f_1 = objective_function(a,S,x0-delta,n); /*Calculate objective_function*/
	f0 = objective_function(a,S,x0,n);
	f1 = objective_function(a,S,x0+delta,n);
	
	if(f_1>f0 && f0>f1) 
	delta=delta;
	else if(f_1<f0 && f0<f1) 
	delta=-delta;
	x[0]=x0;
	f[0]=f0;
	
	feval = feval + 3;
	
	do{  j=j+1;
		
		x[j]=x[j-1]+pow(2,j-1)*delta;
		f[j]=objective_function(a,S,x[j],n);
		
		if(f[j]>f[j-1])
		break;
		feval++;
		i = i + 1;
		  
	}while( (f[j]<f[j-1] || f[j]==f[j-1]) && i <1000);/*Step 3*/

	
	return Newton_Rapson(a,S,(x[j-2]+x[j])/2,n);
}

double Newton_Rapson(double z[100],double S[100],double a,int n){
	
	
	//Step1
	double x[10000];
	
	double e=0.001;
	x[0]=a;
	
	double df1,df2;
	df1=df(z,S,x[0],n);

	//Step2
	df2=ddf(z,S,x[0],n);
	
	//Step3

	int k=0;

	while(abs(df1)>e && k<1000){
		
		k=k+1;
		x[k]=x[k-1]-(df1)/(df2);
		df1=df(z,S,x[k],n);
		df2=ddf(z,S,x[k],n);

	}

	return x[k];
}

double df(double z[100],double S[100],double a,int n){// differentiation as a fun of alpha
	
	double c[100];
	int i;
	
    for(i=1;i<=n;i++)
    
    c[i]=z[i]+a*S[i];

    return (diff_1(c,n));
}

double ddf(double z[100],double S[100],double a,int n){// double diff of fun of alpha
    
    double c[100];
	int i;

    for(i=1;i<=n;i++)
    
    c[i]=z[i]+a*S[i];

    return (diff_2(c,n));
	
}

double constraint(double x[100],int n){// constraint func
 	int i,j;
 	double y1,y2,y,g1,g2,g3,g4,g5,g6,g;
 	
 /*	y1=pow(x[1]-5,2)+pow(x[2]-5,2)-100;
 	y2=-1*pow(x[1]-6,2)-1*pow(x[2]-5,2)+82.81;*/
 	
 	   y1=-1*x[1]*x[1]+x[2]-1;
 	 y2=-1+x[1]-pow(x[2]-4,2);  
 	
 	
 	if(y1<0){
 		y1=y1*y1;
	 }
 	
 	else{
	 y1=0;
	 }
 	
 	
 	if(y2<0){
 		y2=y2*y2;
	 }
 	
 	else{
 		y2=0;
	 }
 	
 	y=y1+y2;
	 
	 return	y;
	 
 } 
 
 	 /*g1=1-0.0025*(x[4]+x[6]);
 	 g2=1-0.0025*(-1*x[4]+x[5]+x[7]);
 	 g3=1-0.01*(x[8]-x[6]);
 	 g4=-100*x[1]+x[1]*x[6]-833.33252*x[4]+83333.333;
 	 g5=-x[2]*x[4]+x[2]*x[7]+1250*x[4]-1250*x[5];
 	 g6=-x[3]*x[5]+x[3]*x[8]+2500*x[5]-1250000;
 	 
 	 	if(g1<0){
 		g1=g1*g1;
	 }
 	
 	else{
	 g1=0;
	 }
 	
 	
 	if(g2<0){
 		g2=g2*g2;
	 }
 	
 	else{
 		g2=0;
	 }
 	 
 	 	if(g3<0){
 		g3=g3*g3;
	 }
 	
 	else{
	 g3=0;
	 }
 	
 	
 	if(g4<0){
 		g4=g4*g4;
	 }
 	
 	else{
 		g4=0;
	 }
	 
	 if(g5<0){
 		g5=g5*g5;
	 }
 	
 	else{
	 g5=0;
	 }
 	
 	if(g6<0){
 		g6=g6*g6;
	 }
 	
 	else{
 		g6=0;
	 }
	 	g=g1+g2+g3+g4+g5+g6;
	 
	 return	g;
}*/
	
  

double func(double x[100],int n){
    int i,j;
    double value=0;
    
   // value=pow(x[1]-10,3)+pow(x[2]-20,3);
    
    value=-1*(pow(sin(2*pi*x[1]),3)*sin(2*pi*x[2]))/(pow(x[1],3)*(x[1]+x[2]));
     //value= x[1]+x[2]+x[3];
    
    return (value+R*constraint(x,n));
}

double diff_1(double x[100],int n)
{
 
	int i,j;
	double y1,y2;
	
	double h=0.001;
	double x1[100],x2[100];
	for(i=1;i<=n;i++){
		x1[i]=x[i]+h;
	}
	y1=func(x1,n);
	
	for(i=1;i<=n;i++){
		x2[i]=x[i]-h;
	}
	y2=func(x2,n);
	
	return (y1-y2)/(2*h);
}

double diff_2(double x[100],int n)
{
	int i,j;
	double y1,y2,y3;
	double h=0.001;
	double x1[100],x2[100];
	
	y3=func(x,n);
	
	for(i=1;i<=n;i++){
		x1[i]=x[i]+h;
	}
	y1=func(x1,n);
	
	for(i=1;i<=n;i++){
		x2[i]=x[i]-h;
	}
	y2=func(x2,n);
	
	return (y1+y2-2*y3)/(h*h);
	
}

int main()
{
    double x[200], xp[200],alpha,d[200],dabs=0,t=0,dummy[200][1];
    double S[200][200];
    int i,n,j,Ns=0;
    double fr1,fr2,*p;
    FILE *avi;
    avi=fopen("output_convergence.txt","w");
    
    printf("Enter the number of variables\n");
    scanf("%d",&n);
    printf("Enter the starting point\n");
    
    p=malloc(sizeof(double)*(n+10));

    for(i=1;i<=n;i++)
    scanf("%lf",&x[i]);
    
    p= powell(x,n);
    
    for(i=1;i<=n;i++){
   xp[i]=p[i];
   printf("%lf\t",xp[i]);
   fprintf(avi,"%lf\t",xp[i]);

   }
    printf("%lf\t",func(xp,n));
   	fprintf(avi,"%lf\t",func(xp,n));
   	
   printf("\n");
   fprintf(avi,"\n");
   fr1=func(xp,n);
   R=R*10;
   fr2=1000;
  	int count=0;
   while(fabs(fr1-fr2)>=0.1  && count++ <50){
   	p= powell(xp,n);
    for(i=1;i<=n;i++){
   		xp[i]=p[i];
      	printf("%lf\t",xp[i]);
      	fprintf(avi,"%lf\t",xp[i]);
   
   	}
   	   printf("%lf\t",func(xp,n));
   	fprintf(avi,"%lf\t",func(xp,n));
   	
   printf("\n");
   fprintf(avi,"\n");
   fr2=func(xp,n);
   R=R*10;
   fr1=fr2;
   fr2=func(xp,n);
   }
   
    for(i=1;i<=n;i++)
    {
    	printf("%lf\t",xp[i]);
    	fprintf(avi,"%lf\t",xp[i]);
     
	}
		printf("%lf\t",func(xp,n));
    	fprintf(avi,"%lf\t",func(xp,n));

   return 0;
}
