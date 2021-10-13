/*
Implementation of the bayesian algorithm for inference DNA replication origin statistics from NGS data with an array of missing values
 
datasets
chr10_{456,567}
chr10_wt_{456,567}
chr10_rat1_{456,567}

Created on: 28 Feb 2018
Author: Alina Bazarova

*/


#include <boost/program_options.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/any.hpp>

#include <iostream>
#include <random> //for random number generator
#include <fstream>
#include <ctime>
#define PI 3.1415

#define EPS 0.000001 // small constant for various purposes

#define EPSO 0.04 //a constant for normalising the proposal in the reversible jump
//#include "iostream.h"
using namespace std;


int compare_doubles (const void *X, const void *Y) //for qsort for doubles
{
       double x = *((double *)X);
       double y = *((double *)Y);

       if (x > y)
       {
               return 1;
       }
       else
       {
               if (x < y)
               {
                       return -1;
               }
               else
               {
                       return 0;
               }
       }
}

int compare_ints (const void *X, const void *Y)//for qsort for ints
{
       int x = *((int *)X);
       int y = *((int *)Y);

       if (x > y)
       {
               return 1;
       }
       else
       {
               if (x < y)
               {
                       return -1;
               }
               else
               {
                       return 0;
               }
       }
}



int ind (double a){//indicator function strict

	if (a>0) return 1;
	else return 0;

}

int indw (double a){//indicator function non-strict

	if (a>=0) return 1;
	else return 0;

}



double sum(double a[],int k,int n){//sum of the elements of the array of doubles

	int i;
	double suma=0.0;
	if (k<=n){
	for (i=k-1;i<n;++i) {

		suma=suma+a[i];
	}
	}
	return suma;
}

double sumint(int a[],int k,int n){//sum of the elements of the array of ints

	int i;
	int suma=0;
	if (k<=n){
	for (i=k-1;i<n;++i) suma+=a[i];
	}
	return suma;
}

double mean(double a[],int k, int n){//mean value of the elements of the array of doubles
	return sum(a,k,n)/(n-k+1);
}




double dotprod(double a[],double b[],int k,int n){//scalar product of two vectors written in arrays of doubles

	int i;

	double dotpr=0.0;

for	(i=k-1;i<n;++i) dotpr+=(a[i]*b[i]);

return dotpr;
}

double dotprodp(double a[],int b[],int n){//scalar product of two vectors one written in array of doubles, another written in array of ints
int i;
double dotpr;
dotpr=0;
for (i=0;i<n;++i){
dotpr=dotpr+a[i]*b[i];

}
return dotpr;
}

double meanp(double a[], int n,int fl[],int fr[]){//for computing mean over active only origins; fl[] and fr[] indicate whether left or right origin is active


double k;
int p;
int f[n];

for (int i=0;i<n;++i) f[i]=fl[i]*fr[i];

p=sumint(f,1,n);

k=dotprodp(a,f,n);

	return k/p;

}

int dotprodpp(int f1[],int f2[],int f3[],int n,int j,int k){
//scalar product which takes into account only active origins; f1[],f2[],f3[] indicate whether origins are active or not


int i;
double dotpr;
dotpr=0;

if (j==2){
	if (k==1){
for (i=0;i<n;++i){
dotpr=dotpr+f2[i]*f1[i];

}
	}
	else{
		for (i=0;i<n;++i){
		dotpr=dotpr+f2[i]*f3[i];

		}


	}
}
else{

	for (i=0;i<n;++i){
	dotpr=dotpr+f1[i]*f3[i];

	}


}
return dotpr;
}

double var(double a[],int k,int n){
	//variance of an array of doubles

return (dotprod(a,a,k,n)-(n-k+1)*mean(a,k,n)*mean(a,k,n))/(n-k);
}


double varp(double a[], int n, int fl[], int fr[]){
	//for computing variance over active only origins; fl[] and fr[] indicate whether left or right origin is active


double v;
double vec[8000];
int i,p,f[8000];


for (int j=0;j<n;++j) f[j]=fl[j]*fr[j];

p=sumint(f,1,n);



for (i=0;i<n;++i){
	vec[i]=a[i]*f[i];
}


if (p>1) v=(dotprod(vec,vec,1,n)-p*meanp(a,n,fl,fr)*meanp(a,n,fl,fr))/(p-1);

else {
	v=dotprod(vec,vec,1,n)-p*meanp(a,n,fl,fr)*meanp(a,n,fl,fr);
}


	return v;
}


double cov(double a[], double b[],int k,int n){
	//covariance of two arrays of double

return dotprod(a,b,k,n)-(n-k+1)*mean(a,k,n)*mean(b,k,n);
}

double xst(double t21diff, int N){ //collision point if both origin fires, double

	return (double) (N+t21diff)/2.0;
}



int fire(double t123sum, double t21diff, double t23diff, int diff1,int diff2,int i){
	//determines whether the origin i is obscured (0) or not (1) given the current time parameters and differences between first two (diff1) and last two (diff2) origins

	int r;

	if (i==1){

		if (((xst(t21diff,diff1)<1)&&(xst(-t23diff,diff2)>=1))) r=0;

		else {

			if ((xst(-t23diff,diff2)<1)&&(xst(t21diff-t23diff,diff1+diff2)<1)) r=0;

			else r=1;
		}


	}
	else{
		if (i==2){

			if (((xst(t21diff,diff1)>=diff1)||(xst(-t23diff,diff2)<1))) r=0;

			else r=1;

		}
		else{

			if ((xst(-t23diff,diff2)>=diff2)&&(xst(t21diff,diff1)<diff1)) r=0;

			else{

				if ((xst(t21diff,diff1)>=diff1)&&(xst(t21diff-t23diff,diff1+diff2)>=diff1+diff2)) r=0;

				else r=1;
			}

		}
	}
	return r;
}




void xstar(double t123sum,double t21diff, double t23diff,int N1, int N2, int &b1, int &b2, int fi1, int fi2, int fi3){

//collision points for a given profile with time parameters determined by the parameters 1-3, distances 4-5, licensing indicators 8-10


	double a1,a2;

if (fi1*fi2*fi3==1){


a1=xst(t21diff,N1);
a2=xst(-t23diff,N2);
//a3=xst(t3,t1,N);

if (fire(t123sum,t21diff,t23diff,N1,N2,1)==0){

	a1=0;

	if (t23diff>N2-1){

		a2=0;

		}
	else{

		if (t23diff<-N2) a2=N2;
	}
}

if (fire(t123sum,t21diff,t23diff,N1,N2,2)==0){

	if (N1-N2<t21diff-t23diff){

		a2=-N1+xst(t21diff-t23diff,N1+N2);

		a1=N1;

		if (a2>N2)	a2=N2;
	}
	else{

		a2=0;

		a1=xst(t21diff-t23diff,N1+N2);

		if (a1<0) a1=0;
	}

}
if (fire(t123sum,t21diff,t23diff,N1,N2,3)==0){

a2=N2;

if (t21diff>N1) a1=N1;

else{

	if (t21diff<-N1+1) a1=0;
}


}
}
else{
	if (fi1+fi2+fi3==2){

	if (fi1==0){

		a1=0;

		a2=xst(-t23diff,N2);

		if (a2<0) a2=0;

		if (a2>N2) a2=N2;
	}
	else{

		if (fi2==0){

			if (N1-N2<t21diff-t23diff){

				a2=-N1+xst(t21diff-t23diff,N1+N2);

				a1=N1;

				if (a2>N2) a2=N2;

			}
			else{

				a2=0;

				a1=xst(t21diff-t23diff,N1+N2);

				if (a1<0) a1=0;

			}
		}
		else{

			a2=N2;

			a1=xst(t21diff,N1);

			if (a1<0) a1=0;

			if (a1>N1) a1=N1;

		}

	}
}
	else{

		if (fi1==1){

			a1=N1;
			a2=N2;

		}
		else{

			if (fi2==1){

				a1=0;
				a2=N2;

			}
			else{

				a1=0;
				a2=0;
			}
		}

	}

}
b1=floor(a1);
b2=floor(a2);

}




void obsc(double t123sum[],double t21diff[], double t23diff[], double t31diff[], int diff1, int diff2, int M, int i,double temp[]){

	//replaces firing times for a given time difference by the obscuring time in case obscuring took place;not used at the moment

	int j;



	for (j=0;j<M;++j){
		if (i==1){

			if (t21diff[j]>=diff1) {

				if (t23diff[j]>=diff2-1){

					if (-t23diff[j]+t21diff[j]<0) temp[j]=t123sum[j]-t23diff[j]+diff1-1;
					else temp[j]=t123sum[j]-t21diff[j]+diff1;

				}
				else{
					temp[j]=t123sum[j]-t21diff[j]+diff1;
				}
			}
			else{
				if (t23diff[j]>=diff2-1){

					temp[j]=t123sum[j]-t23diff[j]+diff1-1;
				}
				else{

					if (t21diff[j]<=-diff1+1) temp[j]=t123sum[j]+t21diff[j]-diff1+1;
					else{

						if (t23diff[j]<=-diff1) temp[j]=t123sum[j]+t23diff[j]-diff1;
						else temp[j]=t123sum[j];
					}
				}


			}



		}
		else{

			if (i==2){

if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,1)==0) {


	if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,2)==0) temp[j]=-diff1;

	else temp[j]=-diff1+1;
}

else{

	if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,2)==0) {


			if (t31diff[j]>=diff1-diff2) temp[j]=diff1;


			else temp[j]=diff2-1+t31diff[j];


	}

	else temp[j]=t21diff[j];

}


			}
			else{

				if (i==3){

				if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,3)==0) {

									temp[j]=-diff2;

				}

				else{

					if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,2)==0){

					if (t31diff[j]<diff1-diff2) temp[j]=diff2-1;

					else {

						temp[j]=diff1-t31diff[j];
					}

					}
					else{

						temp[j]=t23diff[j];

					}
				}

				}
				else{

					if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,2)==0){

						if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,3)==0) temp[j]=diff1+diff2;

						else{

							if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,1)==0) temp[j]=-diff1-diff2+1;

							else temp[j]=t31diff[j];

						}


					}
					else{

						if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,1)==0) {

							if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,3)==0) temp[j]=1;

							else temp[j]=-diff1+1-t23diff[j];
						}

						else{

						if (fire(t123sum[j],t21diff[j],t23diff[j],diff1,diff2,3)==0) temp[j]=t21diff[j]+diff2;

						else temp[j]=t31diff[j];

						}


					}




				}
			}

}
}
}




void multrnd1 (double q0, double q1, double q2,int &fi0, int &fi1,int &fi2,int &fi3){

	//multinomial random variable;not used at the moment


	std::random_device gen;

	std::uniform_real_distribution <double>unifrnd(0.0,1.0);




		double u=unifrnd(gen);


		if (u<q0){

		fi0=1;
		fi1=0;
		fi2=0;

		}
		else{

			if (u<q0+q1){

				fi0=0;
				fi1=1;
				fi2=0;

			}
			else{

				if (u<q0+q1+q2){

					fi0=0;
					fi1=0;
					fi2=1;
				}
				else{

					fi0=0;
					fi1=0;
					fi2=0;

				}
			}
		}

fi3=1-fi0-fi1-fi2;



}




int sign(int a){//sgn function

	int r=0;

	if (a>0) r=1;
	else if (a<0) r=-1;

	return r;
}

double betarnd(int a, int b){//beta random variable; not used at the moment

	std::random_device gen;

	std::gamma_distribution<double> gammarndx(a,1.0);

	std::gamma_distribution<double> gammarndy(b,1.0);

	double x=gammarndx(gen);

	double y=gammarndy(gen);

	return x/(x+y);

}

void dirrnd(int sumfi0, int sumfi1,int sumfi2,double &q0,double &q1, double &q2,int M){

	//Dirichlet random variable; not used at the moment

	std::random_device gen;

	std::gamma_distribution<double> gammarnd0(sumfi0+1,1.0);

	std::gamma_distribution<double> gammarnd1(sumfi1+1,1.0);

	std::gamma_distribution<double> gammarnd2(sumfi2+1,1.0);

	std::gamma_distribution<double> gammarnd3(M-sumfi0-sumfi1-sumfi2+1,1.0);

	double y0=gammarnd0(gen);

	double y1=gammarnd1(gen);

	double y2=gammarnd2(gen);

	double y3=gammarnd3(gen);

q0=y0/(y0+y1+y2+y3);

q1=y1/(y0+y1+y2+y3);

q2=y2/(y0+y1+y2+y3);

}



void multrnd (int M, double q0, double q1,double q2,int fi0[], int fi1[],int fi2[],int fi3[]){

	//another function which is not used

	std::random_device gen;

	std::uniform_real_distribution <double>unifrnd(0.0,1.0);


	for (int i=0;i<M;++i){

		double u=unifrnd(gen);


		if (u<q0){

		fi0[i]=1;
		fi1[i]=0;
		fi2[i]=0;

		}
		else{

			if (u<q0+q1){

				fi0[i]=0;
				fi1[i]=1;
				fi2[i]=0;

			}
			else{

				if (u<q0+q1+q2){

					fi0[i]=0;
					fi1[i]=0;
					fi2[i]=1;
				}
				else{

					fi0[i]=0;
					fi1[i]=0;
					fi2[i]=0;

				}
			}
		}
		fi3[i]=1-fi0[i]-fi1[i]-fi2[i];
	}

	cout<<sumint(fi0,1,M)<<" "<<sumint(fi1,1,M)<<" "<<sumint(fi2,1,M)<<" "<<sumint(fi3,1,M)<<"\n";



}

void sim(double t1m, double t2m, double t3m, double t1s, double t2s, double t3s, int it, int N1,int N2,double xj[],double b,double tau, double q1, double q2, double q3){

	//parameters 1-3 are mean, 4-6 sd firing times, it - amount of iterations to average out
	//N1,N2 distances between origins
	//b,tau,qi - parameters of the model

int i,j,k;

//the following lines for especially large values of it

//double *t1 = (double*)malloc(1000000*sizeof(double));
//double *t2 = (double*)malloc(1000000*sizeof(double));
//double *t3 = (double*)malloc(1000000*sizeof(double));
//int *xst1 = (int*)malloc(1000000*sizeof(int));
//int *xst2 = (int*)malloc(1000000*sizeof(int));


double t1[100000], t2[100000], t3[100000];
int xst1[100000],xst2[100000];

int fi1,fi2,fi3;
std::random_device gen;
std::uniform_real_distribution <double>unifrnd(0.0,1.0);



ofstream f51,f52,f53,f54,f55,f56;
f51.open("filedata1t11000.txt");//files to write in values of firing times
f52.open("filedata1t21000.txt");
f53.open("filedata1t31000.txt");

f54.open("filedata1f11000.txt");//files to write in values of licensing indicators
f55.open("filedata1f21000.txt");
f56.open("filedata1f31000.txt");

cout<<"hi!\n";

for	(i=0;i<it;++i){
	std::normal_distribution<double>normrnd1(t1m,t1s);
	std::normal_distribution<double>normrnd2(t2m,t2s);
	std::normal_distribution<double>normrnd3(t3m,t3s);


	fi1=0;
	fi2=0;
	fi3=0;

	t1[i]=normrnd1(gen);
	t2[i]=normrnd2(gen);
	t3[i]=normrnd3(gen);

	xstar(t1[i]+t2[i]+t3[i],t2[i]-t1[i],t2[i]-t3[i],N1,N2,xst1[i],xst2[i],fi1,fi2,fi3);

	while (((((((xst1[i]==N1)&&(xst2[i]==N2))||((xst1[i]==0)&&(xst2[i]==0)))||((fi2==0)&&(fi1+fi2+fi3==1)))||(fi1+fi2+fi3==0))||((t3[i]-t1[i]>N1+N2)&&(fi2==0)))){

		double u=unifrnd(gen);

		if (u<q1) fi1=1;
		else fi1=0;

		u=unifrnd(gen);

		if (u<q2) fi2=1;
		else fi2=0;

		u=unifrnd(gen);

		if (u<q3) fi3=1;
		else fi3=0;

		xstar(t1[i]+t2[i]+t3[i],t2[i]-t1[i],t2[i]-t3[i],N1,N2,xst1[i],xst2[i],fi1,fi2,fi3);


	}


	cout<<i<<" hi!\n";
	f54<<fi1<<"\n";

	f55<<fi2<<"\n";

	f56<<fi3<<"\n";



	f51<<t1[i]<<"\n";

	f52<<t2[i]<<"\n";

	f53<<t3[i]<<"\n";

}

//if using large it

//free(t1);
//free(t2);
//free(t3);

f51.close();
f52.close();
f53.close();

f54.close();
f55.close();
f56.close();



int ind;
double a;
	for (k=0;k<N1+N2;++k){
		std::normal_distribution<double> normrnd(0,tau);
	printf("%d\n",k); //prints the position number for which the value is being computed

for (j=0;j<it;++j){
		if (((k>=xst1[j])&&(k<N1))||(k>=xst2[j]+N1)) {

			ind=1;
		}
		else ind=0;

	xj[k]=xj[k]+(((1-b)*ind+b)/it);

	if (xj[k]<-1000000) exit(1);
	}

	a=normrnd(gen);

	xj[k]*=exp(a);

}

//if using large it

//	free(xst1);
//	free(xst2);
}

void simr(int it, int N1,int N2,double xjr[],double b,double tau){

	//simulates data on the reverse strand given profile on the forward strand

int i,j,k;



double xj[40000];

double t1[100000], t2[100000], t3[100000];
int xst1[100000],xst2[100000];

int fi1,fi2,fi3;
std::random_device gen;





ifstream f51,f52,f53,f54,f55,f56;
f51.open("filedata1t11000.txt");//files to read off the firing times
f52.open("filedata1t21000.txt");
f53.open("filedata1t31000.txt");

f54.open("filedata1f11000.txt");//files to read off the licensing indicators
f55.open("filedata1f21000.txt");
f56.open("filedata1f31000.txt");


for	(i=0;i<it;++i){


	f54>>fi1;
	f55>>fi2;
	f56>>fi3;

	f51>>t1[i];
	f52>>t2[i];
	f53>>t3[i];

	xstar(t1[i]+t2[i]+t3[i],t2[i]-t1[i],t2[i]-t3[i],N1,N2,xst1[i],xst2[i],fi1,fi2,fi3);


}



//free(t1);
//free(t2);
//free(t3);

f51.close();
f52.close();
f53.close();

f54.close();
f55.close();
f56.close();


int ind;
double a;

	for (k=0;k<N1+N2;++k){
		std::normal_distribution<double> normrnd(0,tau);


for (j=0;j<it;++j){
		if (((k>=xst1[j])&&(k<N1))||(k>=xst2[j]+N1)) {

			ind=1;

		}
		else ind=0;



	xj[k]=xj[k]+((ind+b)/(double) it);


	}




xjr[k]=1-xj[k]+2*b;

	a=normrnd(gen);
	xjr[k]*=exp(a);


}


}



void obsc_samp(double t123sum[], double t21diff[],double t23diff[],int n, int obscpr[], int obscpr1[], int obscpr3[],int fi1[],int fi2[],int fi3[],int M,int N1,int N2){
	//computes obscuring number for each origin


	for (int i=0;i<n;++i) obscpr[i]=0;

	for (int j=0;j<M;++j) {

		if (fi1[j]+fi2[j]+fi3[j]==3){

if (fire(t123sum[j],t21diff[j],t23diff[j],N1,N2,1)==0){


			if (fire(t123sum[j],t21diff[j],t23diff[j],N1,N2,3)==0) obscpr[4]+=1;

			else obscpr[1]+=1;

		}
else{

	if (fire(t123sum[j],t21diff[j],t23diff[j],N1,N2,3)==0) obscpr[3]+=1;

	else {

		if (fire(t123sum[j],t21diff[j],t23diff[j],N1,N2,2)==0) obscpr[2]+=1;

		else obscpr[0]+=1;

	}

}
}
		else{
			if ((fi1[j]==0)&&(fi3[j]==1)){

				if ((fire(t123sum[j],t21diff[j],t23diff[j],N1,N2,3)==0)&&(t21diff[j]<N1)) obscpr1[1]+=1;

				else{

					if (t23diff[j]<-N2) obscpr1[1]+=1;

					else obscpr1[0]+=1;

			}
			}
			else{


				if ((fi1[j]==1)&&(fi3[j]==0)){

					if ((fire(t123sum[j],t21diff[j],t23diff[j],N1,N2,1)==0)&&(t23diff[j]<N2-2)) obscpr3[1]+=1;

					else{

						if (t21diff[j]<-N1+2) obscpr3[1]+=1;

						else obscpr3[0]+=1;
			}
			}
		}

	}

}
}



int obsc_cur(double t123,double t21, double t23, int f1, int f2, int f3, int N1, int N2){

	//determines whether obscuring takes place for a given profile;

	//	0 - no, i- ith origin obscured, 4 - both end origins are obscured; if all origins are licensed

	// 0 no obscuring; 1 obscuring in case on of the end origins is not licensed

	//	obscuring status of the current profile

	int res;

	if (f1*f2*f3==1){

		if (fire(t123,t21,t23,N1,N2,1)==0){

			if (fire(t123,t21,t23,N1,N2,3)==0) res=4;

			else res=1;
		}
		else{

			if (fire(t123,t21,t23,N1,N2,3)==0) res=3;

			else {

				if (fire(t123,t21,t23,N1,N2,2)==0) res=2;

				else res=0;
			}


		}


	}
	else{

		if (f3==0){

if ((fire(t123,t21,t23,N1,N2,1)==0)&&(t23<N2-2)) res=1;

else{

	if (t21<-N1+2) res=1;

	else res=0;

}



		}
		else{

			if (f1==0){

				if ((fire(t123,t21,t23,N1,N2,3)==0)&&(t21<N1)) res=1;

				else{

					if (t23<-N2) res=1;

					else res=0;
			}
		}

	}



}

	return res;
}


void obsc_count(int fi1[], int fi2[], int fi3[],int &obscM, int &obsc1M, int &obsc3M, int M){

//computes how many times all 3 origins are active, how many times only left (1)/right one (3) is inactive


	obscM=0;
	obsc1M=0;
	obsc3M=0;

	for (int i=0; i<M ;++i) {

		if (fi1[i]*fi2[i]*fi3[i]==1) obscM+=1;

		else {

			if ((fi1[i]==0)&&(fi3[i]==1)) obsc1M+=1;

			else {

				if ((fi3[i]==0)&&(fi1[i]==1)) obsc3M+=1;
			}
		}


}

	return;
}



void MCMC(double xj[],double xjr[],int M, int N1,int N2,int it,int minNA,int maxNA){
	//MCMC function

	  clock_t time;

	  time=clock();

double t123sum[4992],t21diff[4992],t23diff[4992],t31diff[4992];// firing time parameters

int fi1[4992],fi2[4992],fi3[4992];//licensing indicators

ofstream f1;
ofstream f2;
ofstream f3;
ofstream f4;
ofstream f5;
ofstream f6;
ofstream f7,f8;

ifstream f21,f123,f23;
ifstream ff1,ff2,ff3;


f1.open("fileb.txt");//output for parameter b
f2.open("filetau.txt");//output for parameter tau
f4.open("filemu.txt");//output for parameter mu
f5.open("filesig.txt");//output for parameter sigma
f6.open("filet.txt");//output for firing times parameters (t1+t2+t3, t2-t1, t2-t3)
f7.open("filef.txt");//output for licensing indicator parameters (f1, f2, f3)
f8.open("fileq.txt");//output for lisencing probability parameters (q1, q2, q3)



std::random_device gen;

std::uniform_real_distribution <double>unifrnd(0.0,1.0); //uniform[0,1]


double diff21m,diff23m, sum123m,diff21s,sum123s, diff23s, diff31m=0,diff31s,s21=1.0,s23=1.0,s31=1.0,s=1.0;

//distributions for priors

std::normal_distribution <double>normrndpmu(0,10000);

std::normal_distribution <double>normrndpb(0,sqrt(10)*1000);

std::gamma_distribution <double>gammarndpsig(3.000001,60000);

std::gamma_distribution <double>gammarndptau(0.01,0.01);

//initialisation of parameters tau, b, q, sig, mu

double tau=gammarndptau(gen),b;

while ((b<=0)&&(b>=1)) b=normrndpb(gen);



//double tau=0.000001,b=0.000005; //if the parameters are fixed

//double q1=1.0,q2=1.0,q3=1.0; //if the parameters are fixed
double q1y,q2y,q3y;

double q1=unifrnd(gen),q2=unifrnd(gen),q3=unifrnd(gen);

//double sig1=pow(10000,-2),sig2=pow(10000,-2),sig3=pow(10000,-2),sig21=pow(200*sqrt(2),-2),sig32=pow(2000*sqrt(2),-2); //if the parameters are fixed

double sig1=gammarndpsig(gen),sig2=gammarndpsig(gen),sig3=gammarndpsig(gen);

double sig1alpha=3.00001,sig1beta=60000,sig2alpha=3.00001,sig2beta=60000,sig3alpha=3.000001,sig3beta=60000.01,taub=0.0000001;
double taualpha=0.01,taubeta=0.01;

//mu1=0,mu2=0,mu3=0; //if the parameters are fixed

double taumu2=0.00000001,taumu1=0.00000001,taumu3=0.000000001,mu12=-23000.0/3.0,mu21=23000.0;

double mu1=normrndpmu(gen),mu2=normrndpmu(gen),mu3=normrndpmu(gen);

int res1=0,res2=0,res3=0,res1a=0,res2a=0,res3a=0; //counters for tuning the sd in random walk

int i,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,i20,i21,i22,i23,i24;//various counters

int xst1[4992],xst2[4992]; //collision points, M=4992

//initialisation of firing times and licensing indicators

for (i=0;i<M;++i){

	std::uniform_real_distribution <double>unifrnd0(-3*N2,3*N2); //uniform[-3N2,3N2]
	std::uniform_real_distribution <double>unifrnd01(-N1-N2,N1+N2); //uniform[-N1-N2,N1+N2]

	t123sum[i]=0;
	t23diff[i]=unifrnd0(gen);
	t31diff[i]=unifrnd01(gen);
	t21diff[i]=t31diff[i]+t23diff[i];

	while (((t23diff[i]<=-N2)&&(t21diff[i]<=-N1+1))||((t21diff[i]<-3*N2)||(t21diff[i]>3*N2))){

		t23diff[i]=unifrnd0(gen);
		t31diff[i]=unifrnd01(gen);
		t21diff[i]=t31diff[i]+t23diff[i];
	}


	//different normalising factors depending on the identifiability constraint in the current region

	double normq10=q1*q2*q3+(1-q1)*q2*q3+q1*q2*(1-q3)+(1-q1)*q2*(1-q3);

	double normq20=q1*q2*q3+(1-q1)*q2*q3+q1*q2*(1-q3)+q1*(1-q2)*q3+(1-q1)*q2*(1-q3);

	double normq30=q1*q2*q3+(1-q2)*q1*q3+(1-q1)*q2*(1-q3);

	double normq40=q1*q2*q3+(1-q2)*q1*q3+(1-q1)*q2*(1-q3)+(1-q1)*q2*q3;

	double normq50=q1*q2*q3+(1-q2)*q1*q3+(1-q1)*q2*(1-q3)+q1*q2*(1-q3);

	double uf0=unifrnd(gen);

	//lisencing indicators

	if ((fabs(t31diff[i])>N1+N2)){// |t3-t1|>N1+N2

		if (uf0<q1*q2*q3/normq10){

			fi1[i]=1;
			fi2[i]=1;
			fi3[i]=1;
		}
		else{

			if (uf0<((1-q1)*q2*q3+q1*q2*q3)/normq10){

				fi1[i]=0;
				fi2[i]=1;
				fi3[i]=1;


			}
			else{

				if (uf0<((1-q1)*q2*q3+q1*q2*q3+q1*q2*(1-q3))/normq10){

				fi1[i]=1;
				fi2[i]=1;
				fi3[i]=0;
				}
				else{

					fi1[i]=0;
					fi2[i]=1;
					fi3[i]=0;
				}

			}
		}
	}
	else{


		if (t21diff[i]>N1){

			if (t23diff[i]>N2-2){// t2-t1>N1 and t2-t3>N2-2

				if (uf0<q1*q2*q3/normq30){

					fi1[i]=1;
					fi2[i]=1;
					fi3[i]=1;

				}
				else{

					if (uf0<((1-q2)*q1*q3+q1*q2*q3)/normq30){

						fi1[i]=1;
						fi2[i]=0;
						fi3[i]=1;
					}
					else{

						fi1[i]=0;
						fi2[i]=1;
						fi3[i]=0;

					}


				}



			}
			else{//t2-t1>N1 and t2-t3<=N2-2

				if (uf0<q1*q2*q3/normq40){

					fi1[i]=1;
					fi2[i]=1;
					fi3[i]=1;

				}
				else{

					if (uf0<((1-q2)*q1*q3+q1*q2*q3)/normq40){

						fi1[i]=1;
						fi2[i]=0;
						fi3[i]=1;
					}
					else{

						if (uf0<((1-q2)*q1*q3+q1*q2*q3+(1-q1)*q2*(1-q3))/normq40){

						fi1[i]=0;
						fi2[i]=1;
						fi3[i]=0;
						}
						else{


							fi1[i]=0;
							fi2[i]=1;
							fi3[i]=1;


						}
					}


				}


			}

		}
		else{



			if (t23diff[i]>N2-2){//t2-t1<=N1 and t2-t3>N2-2

				if (uf0<q1*q2*q3/normq50){

					fi1[i]=1;
					fi2[i]=1;
					fi3[i]=1;

				}
				else{

					if (uf0<((1-q2)*q1*q3+q1*q2*q3)/normq50){

						fi1[i]=1;
						fi2[i]=0;
						fi3[i]=1;
					}
					else{

						if (uf0<((1-q2)*q1*q3+q1*q2*q3+(1-q1)*q2*(1-q3))/normq50){

						fi1[i]=0;
						fi2[i]=1;
						fi3[i]=0;
						}
						else{


							fi1[i]=1;
							fi2[i]=1;
							fi3[i]=0;


						}
					}


				}



			}

			else{//t2-t1<=N1 and t2-t3<=N2-2


				if (uf0<q1*q2*q3/normq20){

					fi1[i]=1;
					fi2[i]=1;
					fi3[i]=1;
				}
				else{

					if (uf0<((1-q1)*q2*q3+q1*q2*q3)/normq20){

						fi1[i]=0;
						fi2[i]=1;
						fi3[i]=1;


					}
					else{

						if (uf0<((1-q3)*q2*q1+q1*q2*q3+(1-q1)*q2*q3)/normq20){

						fi1[i]=1;
						fi2[i]=1;
						fi3[i]=0;

						}

						else{

							if (uf0<((1-q3)*q2*q1+q1*q2*q3+(1-q1)*q2*q3+(1-q2)*q1*q3)/normq20){

							fi1[i]=1;
							fi2[i]=0;
							fi3[i]=1;

							}
							else{

								fi1[i]=0;
								fi2[i]=1;
								fi3[i]=0;
							}
						}
					}
				}


			}
		}




	}

	//collision points for ith profile

	xstar(t123sum[i],t21diff[i],t23diff[i],N1,N2,xst1[i],xst2[i],fi1[i],fi2[i],fi3[i]);

	//check they are feasible

	if (((xst1[i]<0)||(xst2[i]<0))||((xst1[i]>N1)||(xst2[i]>N2))){


		cout<<t123sum[i]<<" t123 "<<t21diff[i]<<" t21 "<<t23diff[i]<<" t23 "<<xst1[i]<<" xst1 "<<xst2[i]<<" xst2 "<<fi1[i]<<" f1 "<<fi2[i]<<" f2 "<<fi3[i]<<" f3 in\n";


		exit(1);
	}

}




int L=M/16; //initial number of the profiles

double temp2[4992],temp3[4992],temp4[4992]; //not used

//amount of active origins out of L

int sumfi3=sumint(fi3,1,M/16);

int sumfi1=sumint(fi1,1,M/16);

int sumfi2=sumint(fi2,1,M/16);

//not used

obsc(t123sum,t21diff,t23diff,t31diff,N1,N2,M/16,2,temp2);

obsc(t123sum,t21diff,t23diff,t31diff,N1,N2,M/16,3,temp3);

obsc(t123sum,t21diff,t23diff,t31diff,N1,N2,M/16,4,temp4);

//mean and sd of the t1+t2+t3 move

sum123m=mean(t123sum,1,M/16);

sum123s=0;

//initial value for the t1+t2+t3 move

double sum1230=4000;

//computes how many times one of the three types of profiles taking place (please see obsc_count() function)

int obscM,obsc1M,obsc3M;// how many times all three origins are active

obsc_count(fi1,fi2,fi3,obscM,obsc1M,obsc3M,L);

//obscuring probabilities of the three types of profiles taking place (please see obsc_samp() function)

int obscpr[5],obscpr1[2],obscpr3[2];



obsc_samp(t123sum,t21diff,t23diff,5,obscpr,obscpr1,obscpr3,fi1,fi2,fi3,L,N1,N2);


cout<<obscpr[0]<<" "<<obscpr[1]<<" "<<obscpr[2]<<" "<<obscpr[3]<<" "<<obscpr[4]<<" obscpr \n";//prints obscuring probabilities for the case when all origins are active

double logxj[4000],logxjr[4000]; //for logged data, size corresponds to the length of the data array

double F1[4000]; //for profile

double logvec2=0; //for squared profiles term


//computation of the initial value of the profile in two parts

for (i6=0;i6<N1;++i6){



F1[i6]=0;

if ((i6<minNA)||(i6>maxNA)){//taking into account the unsequencable region

	logxj[i6]=log(xj[i6]);

	logxjr[i6]=log(xjr[i6]);
}
else{

	logxj[i6]=0;

	logxjr[i6]=0;
}


	for (i4=0;i4<L;++i4){

		F1[i6]=F1[i6]+(1.0/(double) L)*indw(i6-xst1[i4]);

	}

	if ((i6<minNA)||(i6>maxNA)){//taking into account the unsequencable region

	logvec2=logvec2+logxj[i6]*logxj[i6]-2*logxj[i6]*log((1-b)*F1[i6]+0.5*b)+log((1-b)*F1[i6]+0.5*b)*log((1-b)*F1[i6]+0.5*b)+logxjr[i6]*logxjr[i6]-2*logxjr[i6]*log((1-b)*(1-F1[i6])+0.5*b)+log((1-b)*(1-F1[i6])+0.5*b)*log((1-b)*(1-F1[i6])+0.5*b);


	}

}




	for (i6=N1;i6<N1+N2;++i6){



		F1[i6]=0;


		if ((i6<minNA)||(i6>maxNA)){

		logxj[i6]=log(xj[i6]);

		logxjr[i6]=log(xjr[i6]);
		}
		else{

			logxj[i6]=0;

			logxjr[i6]=0;
		}


		for (i4=0;i4<L;++i4){


			F1[i6]=F1[i6]+(1.0/(double) L)*indw(i6-xst2[i4]-N1);


		}



		if ((i6<minNA)||(i6>maxNA)){

		logvec2=logvec2+logxj[i6]*logxj[i6]-2*logxj[i6]*log((1-b)*F1[i6]+0.5*b)+log((1-b)*F1[i6]+0.5*b)*log((1-b)*F1[i6]+0.5*b)+logxjr[i6]*logxjr[i6]-2*logxjr[i6]*log((1-b)*(1-F1[i6])+0.5*b)+log((1-b)*(1-F1[i6])+0.5*b)*log((1-b)*(1-F1[i6])+0.5*b);


		}


	}







	diff21m=mean(temp2,1,M/16); //not used

	diff23m=mean(temp3,1,M/16);

	diff31m=mean(temp4,1,M/16);

	//mean and sd which takes into account licensed origins only


	double	diff21m1=meanp(t21diff,M/16,fi2,fi1);

		diff21s=sqrt(varp(t21diff,M/16,fi2,fi1));

double		diff23m1=meanp(t23diff,M/16,fi2,fi3);

		diff23s=sqrt(varp(t23diff,M/16,fi2,fi3));

double		diff31m1=meanp(t31diff,M/16,fi3,fi1);

		diff31s=sqrt(varp(t31diff,M/16,fi3,fi1));

		//computes the amount of time each of the firing time differences is licensed

int dp21=dotprodpp(fi1,fi2,fi3,M/16,2,1);

int dp23=dotprodpp(fi1,fi2,fi3,M/16,2,3);

int dp31=dotprodpp(fi1,fi2,fi3,M/16,3,1);

//set counters to zero

int resb=0,resba=0,restau=0,restaua=0;

//initial values for sds in random walk moves

double bs=1.0;

double taus=10.0;

//burn-in for the initial amount of realisations (L=M/16)

int burnin0=100000;

int ii=0;

int burnina=0;//auxiliary value

for (i1=0;i1<it;++i1){

printf("%d\n",i1);

if (ii<4){

	//doubles the amount of realisations after burnina+burnin0/(pow(2,ii)) iterations

if (i1==burnina+burnin0/(pow(2,ii))){

	burnina+=burnin0/(pow(2,ii));

	L=2*L;

	//other parameters also have to be recomputed after doubling L

	obscM=2*obscM;

	obsc1M=2*obsc1M;

	obsc3M=2*obsc3M;

	for (int iob=0;iob<5;++iob) obscpr[iob]=2*obscpr[iob];

	for (int iob1=0;iob1<2;++iob1){

		obscpr1[iob1]=2*obscpr1[iob1];

		obscpr3[iob1]=2*obscpr3[iob1];

	}

	dp21=2*dp21;
	dp23=2*dp23;
	dp31=2*dp31;

	sumfi1=2*sumfi1;
	sumfi2=2*sumfi2;
	sumfi3=2*sumfi3;



	diff21s=sqrt( ((double) dp21-2)/ ((double)dp21-1))*diff21s;

	diff23s=sqrt(((double) dp23-2)/( (double) dp23-1))*diff23s;

	diff31s=sqrt(((double) dp31-2)/( (double) dp31-1))*diff31s;

	sum123s=sqrt(((double) L-2)/((double) L-1))*sum123s;


	//assigning values to firing times, licensing indicators, collision points


	for (int ib=L/2;ib<L;++ib) {

		t123sum[ib]=t123sum[ib-L/2];

		t21diff[ib]=t21diff[ib-L/2];
		t23diff[ib]=t23diff[ib-L/2];
		t31diff[ib]=t31diff[ib-L/2];

		fi1[ib]=fi1[ib-L/2];
		fi2[ib]=fi2[ib-L/2];
		fi3[ib]=fi3[ib-L/2];

		temp4[ib]=temp4[ib-L/2]; //not used
		temp2[ib]=temp2[ib-L/2];
		temp3[ib]=temp3[ib-L/2];

		xst1[ib]=xst1[ib-L/2];
		xst2[ib]=xst2[ib-L/2];
	}

	ii+=1;
}


}

//tuning the sd for random walk


if ((i1+1)%25==0){

	printf("rat1 %lf\n",(double) res1a/(double) res1);// acceptance rate for t1+t2+t3 move
	printf("rat2 %lf\n",(double) res2a/(double) res2);//acceptance rate for tie difference move
	printf("ratb %lf\n",(double) resba/(double) resb);//acceptance rate for b
	printf("rattau %lf\n",(double) restaua/(double) restau);//acceptance rate for tau

	if (i1<97375){


if ((double) resba/(double) resb<0.2) bs=0.8*bs;

else{

	if ((double) resba/(double) resb>0.4) bs=1.2*bs;
}


if ((double) restaua/(double) restau<0.2) taus=0.9*taus;

else{

	if ((double) restaua/(double) restau>0.4) taus=1.1*taus;
}
	}

	cout<<s21*diff21s<<" "<<s23*diff23s<<" "<<s31*diff31s<<"\n";

	cout<<bs<<"\n";



 	res1a=0;
 	res2a=0;
 	res3a=0;
 	res1=0;
 	res2=0;
 	res3=0;
 	resb=0;
 	resba=0;
 	restau=0;
 	restaua=0;

 }



//sampling b

std::normal_distribution<double>normrndb(b,bs);

double by=normrndb(gen);

while((by<0)||(by>1)){

	by=normrndb(gen);

}

resb+=1;


//MH step computations for b (acceptance probability)

double pow_alpha_b=0;

double logvec_b=0;

double check=0.0;

for (i12=0;i12<N1+N2;++i12){

	if ((i12<minNA)||(i12>maxNA)){

	double incr=(log((1-by)*F1[i12]+0.5*by)-log((1-b)*F1[i12]+0.5*b))*(2*logxj[i12]-log((1-b)*F1[i12]+0.5*b)-log((1-by)*F1[i12]+0.5*by))+(log((1-by)*(1-F1[i12])+0.5*by)-log((1-b)*(1-F1[i12])+0.5*b))*(2*logxjr[i12]-log((1-b)*(1-F1[i12])+0.5*b)-log((1-by)*(1-F1[i12])+0.5*by));

	logvec_b+=incr;

	pow_alpha_b+=(incr+(1.0/tau)*(log((1-by)*F1[i12]+0.5*by)-log((1-b)*F1[i12]+0.5*b)+log((1-by)*(1-F1[i12])+0.5*by)-log((1-b)*(1-F1[i12])+0.5*b)));

	}

}



pow_alpha_b*=(tau/2.0);

pow_alpha_b-=(0.5*taub*(by-b)*(by+b));

//acceptance-rejection for b

if (pow_alpha_b>0){

	b=by;

	logvec2=logvec2-logvec_b;

	resba+=1;

}
else{

	double ub=unifrnd(gen);

	if (ub<exp(pow_alpha_b)){

		b=by;

		logvec2=logvec2-logvec_b;

		resba+=1;
	}

}

if (logvec2<0) {//check that square term is not negative

	cout<<" b logvec2\n";

	exit(1);
}

f1<<b<<"\n";//write b to a file

cout<<b<<" b\n";//print b


//sampling tau

std::normal_distribution<double>normrndtau(tau,taus);

double tauy=normrndtau(gen);

while((tauy<=0)){

	tauy=normrndtau(gen);

}

restau+=1;

//acceptance probability


double pow_alpha_tau=0;

pow_alpha_tau+=((N1+N2-(maxNA-minNA+1)+taualpha-1)*log(tauy/tau)+0.5*logvec2*(tau-tauy)+0.25*(N1+N2-(maxNA-minNA+1))*((1.0/tau)-(1.0/tauy))+taubeta*(tau-tauy));

//acceptance-rejection for tau


if (pow_alpha_tau>0){

	tau=tauy;

	restaua+=1;

}
else{

	double utau=unifrnd(gen);

	if (utau<exp(pow_alpha_tau)){

		tau=tauy;

		restaua+=1;

	}
}

 cout<<tau<<"\n"; //print tau

 f2<<tau<<"\n"; //write tau to a file

 //sampling mu1,mu2,mu3, precisions sig1,sig2,sig3

std::normal_distribution<double> normrndp1((1.0/3.0)*(mean(t123sum,1,L)-2*mean(t21diff,1,L)+mean(t23diff,1,L))*sig1*L/(sig1*L+taumu1),1.0/sqrt(sig1*L+taumu1));
std::normal_distribution<double> normrndp2((1.0/3.0)*(mean(t123sum,1,L)+mean(t21diff,1,L)+mean(t23diff,1,L))*sig2*L/(sig2*L+taumu2),1.0/sqrt(sig2*L+taumu2));
std::normal_distribution<double> normrndp3((1.0/3.0)*(mean(t123sum,1,L)-2*mean(t23diff,1,L)+mean(t21diff,1,L))*sig3*L/(sig3*L+taumu3),1.0/sqrt(sig3*L+taumu3));





mu1=normrndp1(gen);
mu2=normrndp2(gen);
mu3=normrndp3(gen);

//auxiliary values for precisions


double sig11=0.5*((1.0/9.0)*(dotprod(t123sum,t123sum,1,L)+4.0*dotprod(t21diff,t21diff,1,L)+dotprod(t23diff,t23diff,1,L)-4.0*dotprod(t123sum,t21diff,1,L)-4.0*dotprod(t21diff,t23diff,1,L)+2.0*dotprod(t123sum,t23diff,1,L))-(2.0/3.0)*mu1*(sum(t123sum,1,L)-2*sum(t21diff,1,L)+sum(t23diff,1,L))+L*mu1*mu1+2*sig1beta);
double sig21=0.5*((1.0/9.0)*(dotprod(t123sum,t123sum,1,L)+dotprod(t21diff,t21diff,1,L)+dotprod(t23diff,t23diff,1,L)+2.0*dotprod(t123sum,t21diff,1,L)+2.0*dotprod(t21diff,t23diff,1,L)+2.0*dotprod(t123sum,t23diff,1,L))-(2.0/3.0)*mu2*(sum(t123sum,1,L)+sum(t21diff,1,L)+sum(t23diff,1,L))+L*mu2*mu2+2*sig2beta);
double sig31=0.5*((1.0/9.0)*(dotprod(t123sum,t123sum,1,L)+dotprod(t21diff,t21diff,1,L)+4.0*dotprod(t23diff,t23diff,1,L)+2.0*dotprod(t123sum,t21diff,1,L)-4.0*dotprod(t21diff,t23diff,1,L)-4.0*dotprod(t123sum,t23diff,1,L))-(2.0/3.0)*mu3*(sum(t123sum,1,L)-2*sum(t23diff,1,L)+sum(t21diff,1,L))+L*mu3*mu3+2*sig3beta);

std::gamma_distribution<double> gammarndp1(L/2.0+sig1alpha,1.0/sig11);
std::gamma_distribution<double> gammarndp2(L/2.0+sig2alpha,1.0/sig21);
std::gamma_distribution<double> gammarndp3(L/2.0+sig3alpha,1.0/sig31);


sig1=gammarndp1(gen);
sig2=gammarndp2(gen);
sig3=gammarndp3(gen);

//print sd

cout<<sqrt(1.0/sig3)<<" sig3 \n";
cout<<sqrt(1.0/sig2)<<" sig2 \n";
cout<<sqrt(1.0/sig1)<<" sig1 \n";

//every 10th write iteration values of the mu and sd into a file

if ((i1+1)%10==0) f4<<mu1<<" "<<mu2<<" "<<mu3<<"\n";

if ((i1+1)%10==0) f5<<sqrt(1.0/sig1)<<" "<<sqrt(1.0/sig2)<<" "<<sqrt(1.0/sig3)<<"\n";


//how many times each of the origins is active

cout<<sumfi1<<" sumfi1 "<<sumfi2<<" sumfi2 "<<sumfi3<<" sumfi3\n";

//probability of being obscured given certain amount of origins is licensed

cout<<obscpr[0]/(double) obscM<<" "<<obscpr[1]/(double) obscM<<" "<<obscpr[2]/(double) obscM<<" "<<obscpr[3]/(double) obscM<<" "<<obscpr[4]/(double) obscM<<"\n";

cout<<obscpr1[0]/(double) obsc1M<<" "<<obscpr1[1]/(double) obsc1M<<" "<<obscpr3[0]/(double) obsc3M<<" "<<obscpr3[1]/(double) obsc3M<<"\n";


//sampling q

q1y=betarnd(sumfi1+1,L-sumfi1+1);

cout<<q1<<" q1 "<<q1y<<" q1y\n"; //q1 and proposal for it

q1=q1y;



q2y=betarnd(sumfi2+1,L-sumfi2+1);

cout<<q2<<" q2 "<<q2y<<" q2y\n"; //q2 and proposal for it

q2=q2y;



q3y=betarnd(sumfi3+1,L-sumfi3+1);

cout<<q3<<" q3 "<<q3y<<" q3y\n"; //q3 and proposal for it

q3=q3y;




cout<<q1<<" q1 "<<q2<<" q2 "<<q3<<" q3\n"; //new q

f8<<q1<<" "<<q2<<" "<<q3<<"\n"; //write q to a file


//updating time parameters

for (i5=0;i5<L;++i5){




	double y123,y21,y23;



res1+=1;
res2+=1;
res3+=1;

//updating t1+t2+t3


if (i1>0) {
	std::normal_distribution<double> normrnd1(sum123m,sum123s);

	y123=normrnd1(gen);
}

else {

	std::normal_distribution<double> normrnd1(sum123m,sum1230);

	y123=normrnd1(gen);
}




double sum123my,sum123sy,diff21my,diff21sy,diff23my,diff23sy,diff31my,diff31sy,tempy21,tempy23,tempy31;


//mean and variance of the proposal

sum123my=sum123m+(1.0/L)*(y123-t123sum[i5]);


sum123sy=sqrt(sum123s*sum123s+(1.0/(L-1))*(y123*y123-t123sum[i5]*t123sum[i5]+L*(sum123m*sum123m-sum123my*sum123my)));






//ratio of conditional proposals


double pow1=0;

if (i1>0)
 	 pow1=(-pow(t123sum[i5]-sum123my,2))/(2.0*pow(sum123sy,2))+pow(y123-sum123m,2.0)/(2.0*pow(sum123s,2))+log(sum123s)-log(sum123sy);

else         pow1=(-pow(t123sum[i5]-sum123my,2))/(2.0*pow(sum1230,2))+pow(y123-sum123m,2.0)/(2.0*pow(sum1230,2));




//collision points for new value of parameter

int xsty1,xsty2;


xstar(y123,t21diff[i5],t23diff[i5],N1,N2,xsty1,xsty2,fi1[i5],fi2[i5],fi3[i5]);

//check collision point values are feasible


if (((xsty1<0)||(xsty2<0))||((xsty1>N1)||(xsty2>N2))){


	cout<<y123<<" y123 "<<t21diff[i5]<<" t21 "<<t23diff[i5]<<" t23 "<<xsty1<<" xsty1 "<<xsty2<<" xsty2 "<<fi1[i5]<<" f1 "<<fi2[i5]<<" f2 "<<fi3[i5]<<" f3 1\n";

	exit(1);

}

//MH step computation

double sumfdiff1=0;

double sumfdiff2=0;


int i35;

int minf1=min(xsty1,xst1[i5]);

int maxf1=max(xsty1,xst1[i5]);

int minf2=min(xsty2,xst2[i5]);

int maxf2=max(xsty2,xst2[i5]);




int sgn1=sign(xst1[i5]-xsty1),sgn2=sign(xst2[i5]-xsty2);






double pow_alpha;

//acceptance probability

pow_alpha=(-tau/2.0)*(sumfdiff1+sumfdiff2);


pow_alpha+=(pow1-(sig1/6.0)*(y123-t123sum[i5])*((1.0/3.0)*(y123-2*t21diff[i5]+t23diff[i5])+(1.0/3.0)*(t123sum[i5]-2*t21diff[i5]+t23diff[i5])-2*mu1)-(sig2/6.0)*(y123-t123sum[i5])*((1.0/3.0)*(y123+t21diff[i5]+t23diff[i5])+(1.0/3.0)*(t123sum[i5]+t21diff[i5]+t23diff[i5])-2*mu2)-(sig3/6.0)*(y123-t123sum[i5])*((1.0/3.0)*(y123-2*t23diff[i5]+t21diff[i5])+(1.0/3.0)*(t123sum[i5]-2*t23diff[i5]+t21diff[i5])-2*mu3));


int sk1=0;


//criteria for flat profile elimination elimination

if ((((xsty1==0)&&(xsty2==0))||((xsty1==N1)&&(xsty2==N2)))) sk1=1;


//acceptance-rejection step

if ((pow_alpha>=0)&&(sk1==0)){

	double diff=(y123-t123sum[i5])/((double) L);

	t123sum[i5]=y123;

	//shift to preserve symmetry

	for (int i51=0;i51<L;++i51){

		t123sum[i51]=t123sum[i51]-diff;

	}
	 res1a+=1;




	 logvec2+=(sumfdiff1+sumfdiff2);

	 xst1[i5]=xsty1;
	 xst2[i5]=xsty2;

	 sum123m=sum123my-diff;
	 sum123s=sum123sy;



}

else{//same when one has to sample from uniform[0,1]

	double alpha=exp(pow_alpha);
	 double u=unifrnd(gen);



	 if ((u<=alpha)&&(sk1==0)){


			double diff=(y123-t123sum[i5])/((double) L);

			t123sum[i5]=y123;

			for (int i51=0;i51<L;++i51){

				t123sum[i51]=t123sum[i51]-diff;

			}
			 res1a+=1;




			 logvec2+=(sumfdiff1+sumfdiff2);



			 xst1[i5]=xsty1;
			 xst2[i5]=xsty2;

			 sum123m=sum123my-diff;
			 sum123s=sum123sy;

	 }

}


//updating time differences

double uo=unifrnd(gen);

//flag to determine the region we are moving in

int flago,flago1,flago3;

//proposals to determine probabilities of moving into a certain region (notation analogous to obscpr* values)

int obscpry[5],obscpr1y[2],obscpr3y[2];

//normalising factor for the region where all the origins are licensed

double norma=fmax(obscpr[2]/((double) obscM),EPSO)+fmax(obscpr[1]/((double) obscM),EPSO)+fmax(obscpr[3]/((double) obscM),EPSO)+fmax(obscpr[0]/((double) obscM),EPSO)+fmax(obscpr[4]/((double) obscM),EPSO);

//probabilities of moving within a region when all the origins are licensed

double a010=fmax(obscpr[2]/((double) obscM),EPSO)/norma,a100=fmax(obscpr[1]/((double) obscM),EPSO)/norma,a001=fmax(obscpr[3]/((double) obscM),EPSO)/norma,a000=fmax(obscpr[0]/((double) obscM),EPSO)/norma,a101=1-a000-a100-a001-a010;

//normalising factor for the region where right (3)/left (1) origin is not licensed

double norma3=fmax(obscpr3[1]/((double) obsc3M),EPSO)+fmax(obscpr3[0]/((double) obsc3M),EPSO);

double norma1=fmax(obscpr1[1]/((double) obsc1M),EPSO)+fmax(obscpr1[0]/((double) obsc1M),EPSO);

//probabilities of moving within a region where right (3)/left (1) origin is not licensed

double a110_1=fmax(obscpr3[1]/((double) obsc3M),EPSO)/norma3,a110_0=1-a110_1,a011_1=fmax(obscpr1[1]/((double) obsc1M),EPSO)/norma1,a011_0=1-a011_1;

//choosing which region to move into according to the probabilities computed above and recording the "proposal" probability obscpry*

if (fi1[i5]+fi2[i5]+fi3[i5]==3){//all origins are licensed

for (int io=0;io<5;++io) obscpry[io]=obscpr[io];

int obscc=obsc_cur(t123sum[i5],t21diff[i5],t23diff[i5],1,1,1,N1,N2);


if (uo<a010) {

	flago=1;
	flago1=0;
	flago3=0;

	if (obscc!=2) {


		obscpry[obscc]-=1;

		obscpry[2]+=1;
	}

	//check that proposal count does not go below zero

	if (obscpry[obscc]<0) {

		cout<<obscpr[obscc]<<" "<<obscpry[obscc]<<" "<<obscc<<" "<<a010<<" 2\n";
		exit(1);

	}


}
else {


	flago=0;


	if (uo<a010+a000){

		flago1=0;
		flago3=0;

		if (obscc!=0) {

			obscpry[obscc]-=1;

			obscpry[0]+=1;
		}

		//check that proposal count does not go below zero

		if (obscpry[obscc]<0) {

			cout<<obscpr[obscc]<<" "<<obscpry[obscc]<<" "<<obscc<<" "<<a000<<" 0\n";

			exit(1);

		}

	}
	else{
		if (uo<a010+a000+a001){

		flago1=0;
		flago3=1;

		if (obscc!=3) {

			obscpry[obscc]-=1;

			obscpry[3]+=1;
		}

		//check that proposal count does not go below zero

		if (obscpry[obscc]<0) {

			cout<<obscpr[obscc]<<" "<<obscpry[obscc]<<" "<<obscc<<" "<<a001<<" 3\n";

			exit(1);

		}

		}
		else{

			if (uo<a010+a000+a001+a100){

				flago1=1;
				flago3=0;

				if (obscc!=1) {

					obscpry[obscc]-=1;

					obscpry[1]+=1;
				}

				//check that proposal count does not go below zero

				if (obscpry[obscc]<0) {

					cout<<obscpr[obscc]<<" "<<obscpry[obscc]<<" "<<obscc<<" "<<a100<<" 1\n";
					exit(1);

				}

			}
			else{

				if (obscc!=4) {


					obscpry[obscc]-=1;

					obscpry[4]+=1;
				}

				//check that proposal count does not go below zero

				if (obscpry[obscc]<0) {

					cout<<obscpr[obscc]<<" "<<obscpry[obscc]<<" "<<obscc<<" "<<a101<<" 4\n";
					exit(1);

				}

				flago1=1;
				flago3=1;

			}
		}
	}


}


}
else{

	if ((fi1[i5]==0)&&(fi3[i5]==1)){//only left origin is not licensed

for (int io1=0;io1<2;++io1) {

        obscpr1y[io1]=obscpr1[io1];

}

		int obscc=obsc_cur(t123sum[i5],t21diff[i5],t23diff[i5],0,1,1,N1,N2);


		if (uo<a011_0){

			if (obscc!=0) {

				obscpr1y[1]=obscpr1[1]-1;

				obscpr1y[0]=obscpr1[0]+1;

			}


			flago3=0;
			//check that proposal count does not go below zero

			if (obscpr1y[obscc]<0) {

				cout<<fi1[i5]<<" "<<fi2[i5]<<" "<<fi3[i5]<<"\n";

				cout<<t21diff[i5]<<" "<<t23diff[i5]<<" "<<N1<<" "<<N2<<"\n";

				cout<<fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,3)<<" obsc 3\n";

				cout<<obscpr1[obscc]<<" "<<obscpr1y[obscc]<<" "<<obscc<<" "<<a011_0<<" 011 0\n";
				exit(1);

			}
		}
		else {

			if (obscc!=1) {

				obscpr1y[1]=obscpr1[1]+1;

				obscpr1y[0]=obscpr1[0]-1;

			}


			flago3=1;

			//check that proposal count does not go below zero

			if (obscpr1y[obscc]<0) {

				cout<<obscpr1[obscc]<<" "<<obscpr1y[obscc]<<" "<<obscc<<" "<<a011_1<<" 011 1\n";
				exit(1);

			}
		}


	}
	else{

		if ((fi3[i5]==0)&&(fi1[i5]==1)){//only right origin is not licensed

for (int io1=0;io1<2;++io1) {

        obscpr3y[io1]=obscpr3[io1];

}

			int obscc=obsc_cur(t123sum[i5],t21diff[i5],t23diff[i5],1,1,0,N1,N2);

			if (uo<a110_0){

				if (obscc!=0) {

					obscpr3y[1]=obscpr3[1]-1;

					obscpr3y[0]=obscpr3[0]+1;

				}

				flago1=0;

				//check that proposal count does not go below zero

				if (obscpr3y[obscc]<0) {

					cout<<obscpr3[obscc]<<" "<<obscpr3y[obscc]<<" "<<obscc<<" "<<a110_0<<" 110 0\n";
					exit(1);

				}
			}
			else {

				if (obscc!=1) {

					obscpr3y[1]=obscpr3[1]+1;

					obscpr3y[0]=obscpr3[0]-1;

				}

				flago1=1;

				//check that proposal count does not go below zero

				if (obscpr3y[obscc]<0) {

					cout<<obscpr3[obscc]<<" "<<obscpr3y[obscc]<<" "<<obscc<<" "<<a110_1<<" 110 1\n";
					exit(1);

				}
			}

		}
	}
}


//proposal normalising factor for the region where all the origins are licensed

double normay=fmax(obscpry[2]/((double) obscM),EPSO)+fmax(obscpry[1]/((double) obscM),EPSO)+fmax(obscpry[3]/((double) obscM),EPSO)+fmax(obscpry[0]/((double) obscM),EPSO)+fmax(obscpry[4]/((double) obscM),EPSO);

//proposal probabilities of moving within a region when all the origins are licensed

double a010y=fmax(obscpry[2]/((double) obscM),EPSO)/normay,a100y=fmax(obscpry[1]/((double) obscM),EPSO)/normay,a001y=fmax(obscpry[3]/((double) obscM),EPSO)/normay,a000y=fmax(obscpry[0]/((double) obscM),EPSO)/normay,a101y=1-a000y-a100y-a001y-a010y;

//proposal normalising factor for the region where right (3)/left (1) origin is not licensed

double normay3=fmax(obscpr3y[1]/((double) obsc3M),EPSO)+fmax(obscpr3y[0]/((double) obsc3M),EPSO);

double normay1=fmax(obscpr1y[1]/((double) obsc1M),EPSO)+fmax(obscpr1y[0]/((double) obsc1M),EPSO);

//proposal probabilities for the region where right (3)/left (1) origin is not licensed

double a110y_1=fmax(obscpr3y[1]/((double) obsc3M),EPSO)/normay3,a110y_0=1-a110y_1,a011y_1=fmax(obscpr1y[1]/((double) obsc1M),EPSO)/normay1,a011y_0=1-a011y_1;



//parameters and generators for sampling from truncated normal distribution

double y31,z,up,rho;


//generators for normal distribution

std::normal_distribution<double> normrnd31(diff31m1,s31*diff31s);
std::normal_distribution<double> normrnd23(diff23m1,s23*diff23s);
std::normal_distribution<double> normrnd21(diff21m1,s21*diff21s);

//auxiliary values for truncated normal sampler

double mumin21=(N1-diff21m1)/(s21*diff21s), mumin23=(N2-2-diff23m1)/(s23*diff23s);


double mumin21m=(-N1+2-diff21m1)/(s21*diff21s), mumin23m=(-N2-diff23m1)/(s23*diff23s);

//sampling from truncated normal/normal depending on the region

if (fi1[i5]+fi2[i5]+fi3[i5]==3){//all origins are licensed

if (flago==1){//middle origin obscured

	y31=normrnd31(gen); //proposal for t3-t1

	if (y31<N1-N2){//middle origin obscured by the right one

		if (diff23m1-N2+2>EPS){

			y23=0;


			while(y23<N2-2){


				y23=normrnd23(gen); //proposal for t2-t3

			}


		}
		else{

			double  al232=(mumin23+sqrt(mumin23*mumin23+4))/2.0;

			std::exponential_distribution<double> exp232(al232);

			up=1;





			rho=0;

			while (rho<up){

			z=exp232(gen)+mumin23;

			rho=exp(-pow(al232-z,2.0)/2.0);

			up=unifrnd(gen);


			y23=z;
			}

			y23=s23*diff23s*y23+diff23m1; //proposal for t2-t3


		}



	y21=y23+y31; //proposal for t2-t1

	}
	else{//middle origin is obscured by the left one

		if (diff21m1-N1>EPS){

		y21=0;

		while(y21<N1){

			y21=normrnd21(gen); //proposal for t2-t1
		}

		}

		else{

			double al212=(mumin21+sqrt(mumin21*mumin21+4))/2.0;

			std::exponential_distribution<double> exp212(al212);

			up=1;

			rho=0;

			while (rho<up){

			z=exp212(gen)+mumin21;

			rho=exp(-pow(al212-z,2.0)/2.0);

			up=unifrnd(gen);

			y21=z;
			}

			y21=s21*diff21s*y21+diff21m1; //proposal for t2-t1

		}

		y23=y21-y31; //proposal for t2-t3
	}

}

else{



	if (flago3==0) {//right origin is not obscured

		if (mumin23*mumin23m>=0){

			std::uniform_real_distribution<double> unifrnd23(mumin23m,mumin23);

		up=1;

		rho=0;

		while (rho<up){

			z=unifrnd23(gen);

			if (mumin23<=0)  rho=exp(0.5*((mumin23-z)*(mumin23+z)));

			if (mumin23m>=0)  rho=exp(0.5*((mumin23m-z)*(mumin23m+z)));

			up=unifrnd(gen);

			y23=z;

	}

	y23=s23*diff23s*y23+diff23m1; //proposal for t2-t3





	}
	else{

		y23=N2+2;


		while ((y23>=N2-2)||(y23<=-N2)){


			y23=normrnd23(gen); //proposal for t2-t3

		}

	}


}
	else{//right origin is obscured

		if (diff23m1<=-N2){

			y23=0;

			while (y23>-N2){

				y23=normrnd23(gen); //proposal for t2-t3
			}

		}
		else{

			double al23m=(-mumin23m+sqrt(mumin23m*mumin23m+4))/2.0;

			std::exponential_distribution<double> exp23m(al23m);

			up=1;

			rho=0;

			while (rho<up){

			z=exp23m(gen)-mumin23m;

			rho=exp(-pow(al23m-z,2.0)/2.0);


			up=unifrnd(gen);

			y23=-z;

			}

			y23=s23*diff23s*y23+diff23m1; //proposal for t2-t3


	}

	}



	if (flago1==0) {//left origin is not obscured

		if (mumin21*mumin21m>=0){


			std::uniform_real_distribution<double> unifrnd21(mumin21m,mumin21);

			up=1;

			rho=0;

			while (rho<up){

				z=unifrnd21(gen);

				if (mumin21<=0)  rho=exp(0.5*((mumin21-z)*(mumin21+z)));

				if (mumin21m>=0)  rho=exp(0.5*((mumin21m-z)*(mumin21m+z)));

				up=unifrnd(gen);

				y21=z;

		}

		y21=s21*diff21s*y21+diff21m1; //proposal for t2-t1





	}
	else{

		y21=N1+2;


		while ((y21>=N1)||(y21<=-N1+2)){


			y21=normrnd21(gen); //proposal for t2-t1

		}

	}


}
	else{//left origin is obscured

		if (diff21m1<=-N1+1){

			y21=0;

			while (y21>-N1+2){

				y21=normrnd21(gen); //proposal for t2-t1
			}

		}
		else{

			double al21m=(-mumin21m+sqrt(mumin21m*mumin21m+4))/2.0;

			std::exponential_distribution<double> exp21m(al21m);

			up=1;

			rho=0;

			while (rho<up){

			z=exp21m(gen)-mumin21m;

			rho=exp(-pow(al21m-z,2.0)/2.0);


			up=unifrnd(gen);

			y21=-z;

			}

			y21=s21*diff21s*y21+diff21m1; //proposal for t2-t1


			}

	}



	y31=y21-y23; //proposal for t3-t1

}
}
else{

	if (fi2[i5]==0){//middle origin not licensed

		y31=normrnd31(gen);//proposal for t3-t1

		if (dp21>dp23){//t2-t1 is realised more frequently than t2-t3

			y21=normrnd21(gen); //proposal for t2-t1

			y23=y21-y31; //proposal for t2-t3
		}
		else{//t2-t3 is realised more frequently than t2-t3

			y23=normrnd23(gen); //proposal for t2-t3

			y21=y23+y31; //proposal for t2-t1
		}


	}
	else{

		if (fi1[i5]==0){//left origin is not licensed

			y21=normrnd21(gen); //proposal for t2-t1
		}

		else{//left origin is licensed

			if (flago1==0){//left origin is not obscured

			if (mumin21*mumin21m>=0){

				std::uniform_real_distribution<double> unifrnd21(mumin21m,mumin21);

			up=1;

			rho=0;

			while (rho<up){

				z=unifrnd21(gen);

				if (mumin21<=0)  rho=exp(0.5*((mumin21-z)*(mumin21+z)));

				if (mumin21m>=0)  rho=exp(0.5*((mumin21m-z)*(mumin21m+z)));

				up=unifrnd(gen);

				y21=z;

		}

		y21=s21*diff21s*y21+diff21m1; //proposal for t2-t1



		}
		else{

			y21=N1+2;


			while ((y21>=N1)||(y21<=-N1+2)){


				y21=normrnd21(gen);//proposal for t2-t1

			}

		}
		}
			else{//left origin is obscured

				if (diff21m1<=-N1+2){

					y21=0;

					while (y21>-N1+2){

						y21=normrnd21(gen); //proposal for t2-t1
					}

				}
				else{

					double al21m=(-mumin21m+sqrt(mumin21m*mumin21m+4))/2.0;

					std::exponential_distribution<double> exp21m(al21m);

					up=1;

					rho=0;

					while (rho<up){

					z=exp21m(gen)-mumin21m;

					rho=exp(-pow(al21m-z,2.0)/2.0);


					up=unifrnd(gen);

					y21=-z;

					}

					y21=s21*diff21s*y21+diff21m1; //proposal for t2-t1


					}


			}


		}

		if (fi3[i5]==0){//right origin is not licensed

			y23=normrnd23(gen); //proposal for t2-t3

		}
		else{//right origin is licensed

			if (flago3==0){//right origin is not obscured

			if (mumin23*mumin23m>=0){


					std::uniform_real_distribution<double> unifrnd23(mumin23m,mumin23);

				up=1;

				rho=0;

				while (rho<up){

					z=unifrnd23(gen);

					if (mumin23<=0)  rho=exp(0.5*((mumin23-z)*(mumin23+z)));

					if (mumin23m>=0)  rho=exp(0.5*((mumin23m-z)*(mumin23m+z)));

					up=unifrnd(gen);

					y23=z;

			}

			y23=s23*diff23s*y23+diff23m1; //proposal for t2-t3





		}
		else{

			y23=N2+2;


			while ((y23>=N2-2)||(y23<=-N2)){


				y23=normrnd23(gen); //proposal for t2-t3

			}

		}

			}
			else{//right origin is obscured

				if (diff23m1<=-N2){

					y23=0;

					while (y23>-N2){

						y23=normrnd23(gen);//proposal for t2-t3
					}

				}
				else{

					double al23m=(-mumin23m+sqrt(mumin23m*mumin23m+4))/2.0;

					std::exponential_distribution<double> exp23m(al23m);

					up=1;

					rho=0;

					while (rho<up){

					z=exp23m(gen)-mumin23m;

					rho=exp(-pow(al23m-z,2.0)/2.0);


					up=unifrnd(gen);

					y23=-z;

					}

					y23=s23*diff23s*y23+diff23m1; //proposal for t2-t3


			}


			}


	}
		y31=y21-y23; //proposal for t3-t1

}
}



int sk2=0;

//proposal collision points

xstar(t123sum[i5],y21,y23,N1,N2,xsty1,xsty2,fi1[i5],fi2[i5],fi3[i5]);

//criteria to eliminate flat profiles

if ((((xsty1==0)&&(xsty2==0))||((xsty1==N1)&&(xsty2==N2)))) sk2=1;

//check that the collision times values are feasible

if (((xsty1<0)||(xsty2<0))||((xsty1>N1)||(xsty2>N2))){


	cout<<t123sum[i5]<<" t123 "<<y21<<" y21 "<<y23<<" y23 "<<xsty1<<" xsty1 "<<xsty2<<" xsty2 "<<fi1[i5]<<" f1 "<<fi2[i5]<<" f2 "<<fi3[i5]<<" f3 2\n";

	cout<<diff21s<<" diff21s "<<diff23s<<" diff23s "<<diff31s<<" diff31s\n";

	cout<<dp21<<" dp21 "<<dp23<<" dp23 "<<dp31<<" dp31\n";

	exit(1);
}



double pow2=0;

//not used

if (fire(t123sum[i5],y21,y23,N1,N2,2)==0){

	if (fire(t123sum[i5],y21,y23,N1,N2,1)==0){

		tempy21=-N1;
		tempy23=N2-1;
	}
	else{

		if (y31<N1-N2) {

			tempy23=N2-1;
			tempy21=y31+N2;

		}
		else{

			tempy21=N1;

			if (fire(t123sum[i5],y21,y23,N1,N2,3)==1) tempy23=N1-y31;

			else tempy23=-N2;

		}

	}
}
else{

	if (fire(t123sum[i5],y21,y23,N1,N2,1)==0){

		if	(fire(t123sum[i5],y21,y23,N1,N2,3)==0){

			tempy21=-N1+1;

			tempy23=-N2;

		}
		else{

			tempy21=-N1+1;

			tempy23=y23;
		}

	}
	else{

		if	(fire(t123sum[i5],y21,y23,N1,N2,3)==0){

		tempy21=y21;

		tempy23=-N2;

		}
		else{

			tempy21=y21;

			tempy23=y23;
		}
	}
}


tempy31=tempy21-tempy23;

if ((fabs(tempy21)>N1)||(fabs(tempy23>N2))){

	cout<<tempy21<<" tempy21 "<<tempy23<<" tempy23 "<<tempy31<<" tempy31 "<<flago1<<" flago1 "<<flago<<" flago "<<flago3<<" flago3 "<<fi1[i5]<<" fi1 "<<fi2[i5]<<" fi2 "<<fi3[i5]<<" fi3\n";

	cout<<y21<<" y21 "<<y23<<" y23 "<<y31<<" y31\n";

	exit(1);

}


//parameters of the proposal distribution

//mean value t2-t1

double diff21my1,diff23my1,diff31my1;

if (dp21>0) diff21my1=diff21m1+(1.0/dp21)*(y21-t21diff[i5])*fi2[i5]*fi1[i5];

else diff21my1=diff21m1;

//sd t2-t1

if (dp21>1) diff21sy=sqrt(diff21s*diff21s+(1.0/(dp21-1.0))*fi2[i5]*fi1[i5]*(y21*y21-t21diff[i5]*t21diff[i5])+((double) dp21/(dp21-1.0))*(diff21m1*diff21m1-diff21my1*diff21my1));

else {

	diff21sy=sqrt(fi2[i5]*fi1[i5]*(y21*y21-t21diff[i5]*t21diff[i5])+dp21*(diff21m1*diff21m1-diff21my1*diff21my1));
}

//mean value t2-t3

if (dp23>0) diff23my1=diff23m1+(1.0/dp23)*fi2[i5]*fi3[i5]*(y23-t23diff[i5]);

else diff23my1=diff23m1;

//sd t2-t3

if (dp23>1) diff23sy=sqrt(diff23s*diff23s+(1.0/(dp23-1.0))*fi2[i5]*fi3[i5]*(y23*y23-t23diff[i5]*t23diff[i5])+((double) dp23/(dp23-1.0))*(diff23m1*diff23m1-diff23my1*diff23my1));

else diff23sy=sqrt(fi2[i5]*fi3[i5]*(y23*y23-t23diff[i5]*t23diff[i5])+dp23*(diff23m1*diff23m1-diff23my1*diff23my1));

//mean value t3-t1

if (dp31>0) diff31my1=diff31m1+(1.0/dp31)*fi3[i5]*fi1[i5]*(y31-t31diff[i5]);

else diff31my1=diff31m1;

//sd t3-t1

if (dp31>1) diff31sy=sqrt(diff31s*diff31s+(1.0/(dp31-1.0))*fi3[i5]*fi1[i5]*(y31*y31-t31diff[i5]*t31diff[i5])+((double) dp31/(dp31-1.0))*(diff31m1*diff31m1-diff31my1*diff31my1));

else diff31sy=sqrt(fi3[i5]*fi1[i5]*(y31*y31-t31diff[i5]*t31diff[i5])+dp31*(diff31m1*diff31m1-diff31my1*diff31my1));



//not used

if (dp21>0) diff21my=diff21m+(1.0/dp21)*(tempy21-temp2[i5])*fi2[i5]*fi1[i5];

else diff21my=diff21m;


if (dp23>0) diff23my=diff23m+(1.0/dp23)*fi2[i5]*fi3[i5]*(tempy23-temp3[i5]);

else diff23my=diff23m;

if (dp31>0) diff31my=diff31m+(1.0/dp31)*fi3[i5]*fi1[i5]*(tempy31-temp4[i5]);

else diff31my=diff31m;


//for quantiles of standard normal distribution of the proposal

double mumin21y=(N1-diff21my1)/(s21*diff21sy);
double mumin23y=(N2-2-diff23my1)/(s23*diff23sy);
double mumin21my=(-N1+2-diff21my1)/(s21*diff21sy);
double mumin23my=(-N2-diff23my1)/(s23*diff23sy);

//conditional proposal density

if (fi1[i5]*fi2[i5]*fi3[i5]==1){//all origins are licensed

	if (flago==1){//middle origin of the proposal is obscured

		pow2=pow(y31-diff31m1,2.0)/(2.0*pow(s31*diff31s,2.0))+log(s31*diff31s);

		pow2-=log(a010);

		if (y31<N1-N2){//middle origin o the proposal is obscured by the right one

			pow2+=(pow(y23-diff23m1,2.0)/(2.0*pow(s23*diff23s,2.0))+log(s23*diff23s)+log(0.5-0.5*erf(mumin23/sqrt(2.0))));

		}
		else{//middle origin of the proposal is obscured by the left one

			pow2+=(pow(y21-diff21m1,2.0)/(2.0*pow(s21*diff21s,2.0))+log(s21*diff21s)+log(0.5-0.5*erf(mumin21/sqrt(2.0))));

		}
	}
	else{

		if (flago1==1){//left origin of the proposal is obscured

			pow2=pow(y21-diff21m1,2.0)/(2.0*pow(s21*diff21s,2.0))+log(s21*diff21s)+log(0.5+0.5*erf(mumin21m/sqrt(2.0)));


			if (flago3==0) pow2-=log(a100); //right origin of the proposal is not obscured

			else pow2-=log(a101);//right origin of the proposal is obscured


		}
		else{//left origin of the proposal is not obscured


			pow2=pow(y21-diff21m1,2.0)/(2.0*pow(s21*diff21s,2.0))+log(s21*diff21s)+log(0.5*erf(mumin21/sqrt(2.0))-0.5*erf(mumin21m/sqrt(2.0)));

			if (flago3==0) pow2-=log(a000);//right origin of the proposal is not obscured

			else pow2-=log(a001);//right origin of the proposal is obscured

		}

			if (flago3==1){//right origin of the proposal is obscured

				pow2+=(pow(y23-diff23m1,2.0)/(2.0*pow(s23*diff23s,2.0))+log(s23*diff23s)+log(0.5+0.5*erf(mumin23m/sqrt(2.0))));
			}

				else{//right origin of the proposal is not obscured

					pow2+=(pow(y23-diff23m1,2.0)/(2.0*pow(s23*diff23s,2.0))+log(s23*diff23s)+log(0.5*erf(mumin23/sqrt(2.0))-0.5*erf(mumin23m/sqrt(2.0))));


				}



			}



	if (fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,2)==0){//current middle origin is obscured


		pow2-=(pow(t31diff[i5]-diff31my1,2.0)/(2.0*pow(s31*diff31sy,2.0))+log(s31*diff31sy));


		pow2+=log(a010y);

		if (t31diff[i5]<N1-N2){//current middle origin is obscured by the left one

			pow2-=(pow(t23diff[i5]-diff23my1,2.0)/(2.0*pow(s23*diff23sy,2.0))+log(s23*diff23sy)+log(0.5-0.5*erf(mumin23y/sqrt(2.0))));

		}
		else{//current middle origin is obscured by the right one

			pow2-=(pow(t21diff[i5]-diff21my1,2.0)/(2.0*pow(s21*diff21sy,2.0))+log(s21*diff21sy)+log(0.5-0.5*erf(mumin21y/sqrt(2.0))));

		}



	}
	else{

		if (fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,1)==0){//current left origin is obscured

			if (fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,3)==0) pow2+=log(a101y); //current right origin is obscured

			else pow2+=log(a100y);


			pow2-=(pow(t21diff[i5]-diff21my1,2.0)/(2.0*pow(s21*diff21sy,2.0))+log(s21*diff21sy)+log(0.5+0.5*erf(mumin21my/sqrt(2.0))));

		}
			else{//current right origin is not obscured

				if (fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,3)==0) pow2+=log(a001y);

				else pow2+=log(a000y);


				pow2-=(pow(t21diff[i5]-diff21my1,2.0)/(2.0*pow(s21*diff21sy,2.0))+log(s21*diff21sy)+log(0.5*erf(mumin21y/sqrt(2.0))-0.5*erf(mumin21my/sqrt(2.0))));


			}


			if (fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,3)==0){//current right origin is not obscured

				pow2-=(pow(t23diff[i5]-diff23my1,2.0)/(2.0*pow(s23*diff23sy,2.0))+log(s23*diff23sy)+log(0.5+0.5*erf(mumin23my/sqrt(2.0))));
			}

			else{//current right origin is not obscured

				pow2-=(pow(t23diff[i5]-diff23my1,2.0)/(2.0*pow(s23*diff23sy,2.0))+log(s23*diff23sy)+log(0.5*erf(mumin23y/sqrt(2.0))-0.5*erf(mumin23my/sqrt(2.0))));


			}



}
}
else{

if (fi2[i5]==0){//middle origin is not licensed

	pow2=pow(y31-diff31m1,2.0)/(2.0*pow(s31*diff31s,2.0))+log(s31*diff31s);


if (dp21>dp23) pow2+=(pow(y21-diff21m1,2.0)/(2.0*pow(s21*diff21s,2.0))+log(s21*diff21s)); //t2-t1  is realised more frequently than t2-t3

else pow2+=(pow(y23-diff23m1,2.0)/(2.0*pow(s23*diff23s,2.0))+log(s23*diff23s)); //t2-t3  is realised more frequently than t2-t3


	pow2-=(pow(t31diff[i5]-diff31my1,2.0)/(2.0*pow(s31*diff31sy,2.0))+log(s31*diff31sy));



	if (dp21>dp23)	pow2-=(pow(t21diff[i5]-diff21my1,2.0)/(2.0*pow(s21*diff21sy,2.0))+log(s21*diff21sy)); //t2-t1  is realised more frequently than t2-t3

	else pow2-=(pow(t23diff[i5]-diff23my1,2.0)/(2.0*pow(s23*diff23sy,2.0))+log(s23*diff23sy)); //t2-t3  is realised more frequently than t2-t3
}
else{

	if (fi1[i5]==0){//left origin is not licensed

		pow2=pow(y21-diff21m1,2.0)/(2.0*pow(s21*diff21s,2.0))+log(s21*diff21s);

		pow2-=(pow(t21diff[i5]-diff21my1,2.0)/(2.0*pow(s21*diff21sy,2.0))+log(s21*diff21sy));
	}
	else{//left origin is licensed

		//left origin of the prooposal is not obscured

		if (flago1==0) pow2=pow(y21-diff21m1,2.0)/(2.0*pow(s21*diff21s,2.0))+log(s21*diff21s)+log(0.5*erf(mumin21/sqrt(2.0))-0.5*erf(mumin21m/sqrt(2.0)))-log(a110_0);

		//left of the proposal origin is obscured

		else pow2=pow(y21-diff21m1,2.0)/(2.0*pow(s21*diff21s,2.0))+log(s21*diff21s)+log(0.5+0.5*erf(mumin21m/sqrt(2.0)))-log(a110_1);

		//current left origin is not obscured

		if (((fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,1)==1)&&(t21diff[i5]>-N1+2))||((fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,1)==0)&&(t23diff[i5]>N2-2))) pow2-=(pow(t21diff[i5]-diff21my1,2.0)/(2.0*pow(s21*diff21sy,2.0))+log(s21*diff21sy)+log(0.5*erf(mumin21y/sqrt(2.0))-0.5*erf(mumin21my/sqrt(2.0)))-log(a110y_0));

		//current left origin is obscured

		else pow2-=(pow(t21diff[i5]-diff21my1,2.0)/(2.0*pow(s21*diff21sy,2.0))+log(s21*diff21sy)+log(0.5+0.5*erf(mumin21my/sqrt(2.0)))-log(a110y_1));

	}

	if (fi3[i5]==1){//right origin is licensed


		//right origin of the proposal is obscured

		if (flago3==0) pow2+=(pow(y23-diff23m1,2.0)/(2.0*pow(s23*diff23s,2.0))+log(s23*diff23s)+log(0.5*erf(mumin23/sqrt(2.0))-0.5*erf(mumin23m/sqrt(2.0)))-log(a011_0));

		//right origin of the proposal is not obscured

		else pow2+=(pow(y23-diff23m1,2.0)/(2.0*pow(s23*diff23s,2.0))+log(s23*diff23s)+log(0.5+0.5*erf(mumin23m/sqrt(2.0)))-log(a011_1));

		//current right origin is not obscured

		if (((fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,3)==1)&&(t23diff[i5]>-N2))||((fire(t123sum[i5],t21diff[i5],t23diff[i5],N1,N2,3)==0)&&(t21diff[i5]>N1))) pow2-=(pow(t23diff[i5]-diff23my1,2.0)/(2.0*pow(s23*diff23sy,2.0))+log(s23*diff23sy)+log(0.5*erf(mumin23y/sqrt(2.0))-0.5*erf(mumin23my/sqrt(2.0)))-log(a011y_0));

		//current right origin is obscured

		else pow2-=(pow(t23diff[i5]-diff23my1,2.0)/(2.0*pow(s23*diff23sy,2.0))+log(s23*diff23sy)+log(0.5+0.5*erf(mumin23my/sqrt(2.0)))-log(a011y_1));
	}
	else{//right origin is not licensed

		pow2+=(pow(y23-diff23m1,2.0)/(2.0*pow(s23*diff23s,2.0))+log(s23*diff23s));


		pow2-=(pow(t23diff[i5]-diff23my1,2.0)/(2.0*pow(s23*diff23sy,2.0))+log(s23*diff23sy));


	}


}




}

//MH step computation


int i36;

minf1=min(xsty1,xst1[i5]);

maxf1=max(xsty1,xst1[i5]);

minf2=min(xsty2,xst2[i5]);

maxf2=max(xsty2,xst2[i5]);








sgn1=sign(xst1[i5]-xsty1);

sgn2=sign(xst2[i5]-xsty2);

double Fmin1=0;

double Fmin2=0;




double prodfdiff1=0.0;

double prodfdiff2=0.0;

double prodfdiff1r=0.0;

double prodfdiff2r=0.0;

double prodfdiff1_log=0.0;

double prodfdiff1r_log=0.0;

double prodfdiff2_log=0.0;

double prodfdiff2r_log=0.0;

//finding collision points which will be affected by the proposal

int xsta1[4992];

xsta1[0]=minf1;

int count1=0;

for (i23=0;i23<L;++i23) {




	if ((xst1[i23]>minf1)&&(xst1[i23]<maxf1)) {







		count1+=1;

		xsta1[count1]=xst1[i23];



	}
}

//identifying parts of the profile affected by these collision points

Fmin1=F1[minf1];

Fmin2=F1[minf2+N1];


xsta1[count1+1]=maxf1;


//sorting collision points affected by the proposal


qsort((void*)xsta1,count1+2,sizeof(int),compare_ints);


//computing the likelihood ratio based on the collision points from above


//term which depends on the collision points only

//we take into account unsequencable region


for (i2=0;i2<count1+1;++i2){

	if (xsta1[i2]<minNA){

		int NA=min(xsta1[i2+1],minNA-1);

double incr1=log(((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))*((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))+(1-b)*(1.0/L)*((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))*sgn1)*(-NA+xsta1[i2])*log(1+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2));

double incr1r=log(((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))*((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))-(1-b)*(1.0/L)*((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))*sgn1)*(-NA+xsta1[i2])*log(1-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i2));


prodfdiff1_log+=incr1;

prodfdiff1r_log+=incr1r;

prodfdiff1+=(incr1-(1.0/tau)*(-NA+xsta1[i2])*log(1+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2)));

prodfdiff1r+=(incr1r-(1.0/tau)*(-NA+xsta1[i2])*log(1-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i2)));
	}

	if (xsta1[i2+1]>maxNA){

		int NA=max(xsta1[i2],maxNA+1);

double incr1=log(((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))*((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))+(1-b)*(1.0/L)*((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))*sgn1)*(-xsta1[i2+1]+NA)*log(1+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2));

double incr1r=log(((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))*((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))-(1-b)*(1.0/L)*((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))*sgn1)*(-xsta1[i2+1]+NA)*log(1-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i2));


prodfdiff1_log+=incr1;

prodfdiff1r_log+=incr1r;

prodfdiff1+=(incr1-(1.0/tau)*(-xsta1[i2+1]+NA)*log(1+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2)));

prodfdiff1r+=(incr1r-(1.0/tau)*(-xsta1[i2+1]+NA)*log(1-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i2)));
	}

	//check that the log expression is not negative

	if ((1.0+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2))<0) {

		cout<<" "<<1.0+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2)<<" "<<b<<" log 1\n";
		exit(1);
	}



}
//term which depends on the data xj[]

double at,atr;

for (i3=0;i3<count1+1;++i3){

	at=0;
	atr=0;

	for (i24=xsta1[i3];i24<xsta1[i3+1];++i24) {



		at+=(2*logxj[i24]);
		atr+=(2*logxjr[i24]);

	}



	double incr_at1=at*log(1.0+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i3));

	double incr_at1r=atr*log(1.0-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i3));



	prodfdiff1+=incr_at1;

	prodfdiff1_log+=incr_at1;

	prodfdiff1r+=incr_at1r;

	prodfdiff1r_log+=incr_at1r;

	//check that the log expression is not negative

	if ((1.0+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i3))<=0) {

		cout<<"log 2\n";
		exit(1);
	}

}

//same computation for the right collision point



int xsta2[4992];


int count2=0;

xsta2[0]=minf2;

for (i23=0;i23<L;++i23) {

	if ((xst2[i23]>minf2)&&(xst2[i23]<maxf2)) {

		count2+=1;

		xsta2[count2]=xst2[i23];




	}
}

xsta2[count2+1]=maxf2;



qsort((void*)xsta2,count2+2,sizeof(int),compare_ints);



for (i2=0;i2<count2+1;++i2){

	if (xsta2[i2]<minNA-N1){

		int NA=min(xsta2[i2+1],minNA-N1-1);


double incr2=log(((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))*((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))+(1-b)*(1.0/L)*((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))*sgn2)*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))*(-NA+xsta2[i2]);

double incr2r=log(((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))*((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))-(1-b)*(1.0/L)*((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))*sgn2)*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i2))*(-NA+xsta2[i2]);


prodfdiff2_log+=incr2;

prodfdiff2r_log+=incr2r;

prodfdiff2+=(incr2-(1.0/tau)*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))*(-NA+xsta2[i2]));

prodfdiff2r+=(incr2r-(1.0/tau)*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i2))*(-NA+xsta2[i2]));

	}

	if (xsta2[i2+1]>maxNA-N1){

		int NA=max(xsta2[i2],maxNA-N1+1);


double incr2=log(((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))*((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))+(1-b)*(1.0/L)*((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))*sgn2)*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))*(-xsta2[i2+1]+NA);

double incr2r=log(((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))*((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))-(1-b)*(1.0/L)*((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))*sgn2)*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i2))*(-xsta2[i2+1]+NA);


prodfdiff2_log+=incr2;

prodfdiff2r_log+=incr2r;

prodfdiff2+=(incr2-(1.0/tau)*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))*(-xsta2[i2+1]+NA));

prodfdiff2r+=(incr2r-(1.0/tau)*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i2))*(-xsta2[i2+1]+NA));

	}

	if ((1.0+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))<=0) {

		cout<<b<<" log 3\n";
		exit(1);
	}

}



for (i3=0;i3<count2+1;++i3){

	at=0;

	atr=0;

	for (i24=xsta2[i3]+N1;i24<xsta2[i3+1]+N1;++i24) {

		at+=(2*logxj[i24]);

		atr+=(2*logxjr[i24]);
	}

double incr_at2=at*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i3));

double incr_at2r=atr*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i3));



prodfdiff2+=incr_at2;

prodfdiff2r+=incr_at2r;

prodfdiff2_log+=incr_at2;

prodfdiff2r_log+=incr_at2r;

	if (1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+i3)<=0) {

		cout<<"log 4\n";
		exit(1);
	}


}



//acceptance-rejection step

//likelihood ratio



pow_alpha=(prodfdiff1+prodfdiff2+prodfdiff1r+prodfdiff2r)*(tau/2.0);

//ratio of conditional proposals, priors

pow_alpha+=(pow2-(sig1/6.0)*(2*t21diff[i5]-t23diff[i5]-2*y21+y23)*((1.0/3.0)*(t123sum[i5]-2*y21+y23)+(1.0/3.0)*(t123sum[i5]-2*t21diff[i5]+t23diff[i5])-2*mu1)-(sig2/6.0)*(y23+y21-t21diff[i5]-t23diff[i5])*((1.0/3.0)*(t123sum[i5]+y21+y23)+(1.0/3.0)*(t123sum[i5]+t21diff[i5]+t23diff[i5])-2*mu2)-(sig3/6.0)*(-t21diff[i5]+y21-2*y23+2*t23diff[i5])*((1.0/3.0)*(t123sum[i5]+y21-2*y23)+(1.0/3.0)*(t123sum[i5]+t21diff[i5]-2*t23diff[i5])-2*mu3));




//elimination criteria for flat profiles

if ((((xsty1==0)&&(xsty2==0))||((xsty1==N1)&&(xsty2==N2)))) sk2=1;



if ((pow_alpha>=0)&&(sk2==0)){
	//updating time differences

	t21diff[i5]=y21;
	t23diff[i5]=y23;
	t31diff[i5]=y31;


	temp2[i5]=tempy21;//not used
	temp3[i5]=tempy23;
	temp4[i5]=tempy31;

	if (fi1[i5]+fi2[i5]+fi3[i5]==3) {

		for (int ipr=0;ipr<5;++ipr) obscpr[ipr]=obscpry[ipr];

	}
	else{

		if ((fi1[i5]==0)&&(fi3[i5]==1)){

			obscpr1[0]=obscpr1y[0];

			obscpr1[1]=obscpr1y[1];


		}
		else{

			if ((fi1[i5]==1)&&(fi3[i5]==0)){

				obscpr3[0]=obscpr3y[0];

				obscpr3[1]=obscpr3y[1];

			}
		}

	}

	//acceptance count

	 res2a+=1;

	 //updating profile

	 for (i9=minf1;i9<maxf1;++i9) {

		 F1[i9]=F1[i9]+(1.0/L)*sgn1;


	 }

	 for (i8=minf2+N1;i8<maxf2+N1;++i8) {


		 F1[i8]=F1[i8]+(1.0/L)*sgn2;


	 }

	 //updating square term

	 logvec2-=(prodfdiff1_log+prodfdiff2_log+prodfdiff1r_log+prodfdiff2r_log);

	 //updating collision points

	 xst1[i5]=xsty1;
	 xst2[i5]=xsty2;

	 //updating mean and sd of the proposals

	 diff21m=diff21my;
	 diff21m1=diff21my1;

	 diff21s=diff21sy;

	 diff23m=diff23my;
	 diff23m1=diff23my1;

	 diff23s=diff23sy;

	 diff31m=diff31my;
	 diff31m1=diff31my1;

	 diff31s=diff31sy;



}

else{//same updates but when one has to sample uniform[0,1] variable

	double alpha=exp(pow_alpha);
	 double u=unifrnd(gen);


	 if ((u<=alpha)&&(sk2==0)){


			if (fi1[i5]+fi2[i5]+fi3[i5]==3) {

				for (int ipr=0;ipr<5;++ipr) obscpr[ipr]=obscpry[ipr];

			}
			else{

				if ((fi1[i5]==0)&&(fi3[i5]==1)){

					obscpr1[0]=obscpr1y[0];

					obscpr1[1]=obscpr1y[1];


				}
				else{

					if ((fi1[i5]==1)&&(fi3[i5]==0)){

						obscpr3[0]=obscpr3y[0];

						obscpr3[1]=obscpr3y[1];

					}
				}

			}

		 	t21diff[i5]=y21;
		 	t23diff[i5]=y23;
			t31diff[i5]=y31;

			temp2[i5]=tempy21;
			temp3[i5]=tempy23;
			temp4[i5]=tempy31;



		 	 res2a+=1;

			 for (i9=minf1;i9<maxf1;++i9) {

				 F1[i9]=F1[i9]+(1.0/L)*sgn1;


			 }

			 for (i8=minf2+N1;i8<maxf2+N1;++i8) {


				 F1[i8]=F1[i8]+(1.0/L)*sgn2;


			 }

		 	 logvec2-=(prodfdiff1_log+prodfdiff2_log+prodfdiff1r_log+prodfdiff2r_log);



			 xst1[i5]=xsty1;
			 xst2[i5]=xsty2;


			 diff21m=diff21my;
			 diff21m1=diff21my1;

			 diff21s=diff21sy;

			 diff23m=diff23my;
			 diff23m1=diff23my1;

			 diff23s=diff23sy;

			 diff31m=diff31my;
			 diff31m1=diff31my1;

			 diff31s=diff31sy;




	 }
}





}

//check that the square term is above zero

if (logvec2<0) {

	cout<<" t logvec2\n";

	exit(1);
}

//record every tenth value of the time parameters to a file


if ((i1+1)%10==0){
for (i24=0;i24<L;++i24) {

	f6<<t123sum[i24]<<" "<<t21diff[i24]<<" "<<t23diff[i24]<<"\n";

}
}

//obscuring count for the case when all origins are licensed

cout<<obscpr[0]<<" "<<obscpr[1]<<" "<<obscpr[2]<<" "<<obscpr[3]<<" "<<obscpr[4]<<" obscpr t\n";


//updating licensing indicators

for (i6=0;i6<L;++i6){

int fi1y,fi2y,fi3y;


//normalising factors depending on the identifiability constraint in the current region




double normq1=q1*q2*q3+(1-q1)*q2*q3+q1*q2*(1-q3)+(1-q1)*q2*(1-q3);

double normq2=q1*q2*q3+(1-q1)*q2*q3+q1*q2*(1-q3)+q1*(1-q2)*q3+(1-q1)*q2*(1-q3);

double normq3=q1*q2*q3+(1-q2)*q1*q3+(1-q1)*q2*(1-q3);

double normq4=q1*q2*q3+(1-q2)*q1*q3+(1-q1)*q2*(1-q3)+(1-q1)*q2*q3;

double normq5=q1*q2*q3+(1-q2)*q1*q3+(1-q1)*q2*(1-q3)+q1*q2*(1-q3);

//sampling the proposal

double uf=unifrnd(gen);

if ((fabs(t31diff[i6])>N1+N2)){//|t3-t1|>N1+N2

	if (uf<q1*q2*q3/normq1){

		fi1y=1;
		fi2y=1;
		fi3y=1;
	}
	else{

		if (uf<((1-q1)*q2*q3+q1*q2*q3)/normq1){

			fi1y=0;
			fi2y=1;
			fi3y=1;


		}
		else{

			if (uf<((1-q1)*q2*q3+q1*q2*q3+q1*q2*(1-q3))/normq1){

			fi1y=1;
			fi2y=1;
			fi3y=0;
			}
			else{

				fi1y=0;
				fi2y=1;
				fi3y=0;
			}

		}
	}
}
else{//|t3-t1|<=N1+N2


	if (t21diff[i6]>N1){

		if (t23diff[i6]>N2-2){//t2-t1>N1 and t2-t3>N2-2

			if (uf<q1*q2*q3/normq3){

				fi1y=1;
				fi2y=1;
				fi3y=1;

			}
			else{

				if (uf<((1-q2)*q1*q3+q1*q2*q3)/normq3){

					fi1y=1;
					fi2y=0;
					fi3y=1;
				}
				else{

					fi1y=0;
					fi2y=1;
					fi3y=0;

				}


			}



		}
		else{//t2-t1>N1 but t2-t3<-N2-2

			if (uf<q1*q2*q3/normq4){

				fi1y=1;
				fi2y=1;
				fi3y=1;

			}
			else{

				if (uf<((1-q2)*q1*q3+q1*q2*q3)/normq4){

					fi1y=1;
					fi2y=0;
					fi3y=1;
				}
				else{

					if (uf<((1-q2)*q1*q3+q1*q2*q3+(1-q1)*q2*(1-q3))/normq4){

					fi1y=0;
					fi2y=1;
					fi3y=0;
					}
					else{


						fi1y=0;
						fi2y=1;
						fi3y=1;


					}
				}


			}


		}

	}
	else{



		if (t23diff[i6]>N2-2){//t2-t1<=N1 and t2-t3>N2-2

			if (uf<q1*q2*q3/normq5){

				fi1y=1;
				fi2y=1;
				fi3y=1;

			}
			else{

				if (uf<((1-q2)*q1*q3+q1*q2*q3)/normq5){

					fi1y=1;
					fi2y=0;
					fi3y=1;
				}
				else{

					if (uf<((1-q2)*q1*q3+q1*q2*q3+(1-q1)*q2*(1-q3))/normq5){

					fi1y=0;
					fi2y=1;
					fi3y=0;
					}
					else{


						fi1y=1;
						fi2y=1;
						fi3y=0;


					}
				}


			}



		}

		else{//t2-t1<=N1 and t2-t3<=N2-2


			if (uf<q1*q2*q3/normq2){

				fi1y=1;
				fi2y=1;
				fi3y=1;
			}
			else{

				if (uf<((1-q1)*q2*q3+q1*q2*q3)/normq2){

					fi1y=0;
					fi2y=1;
					fi3y=1;


				}
				else{

					if (uf<((1-q3)*q2*q1+q1*q2*q3+(1-q1)*q2*q3)/normq2){

					fi1y=1;
					fi2y=1;
					fi3y=0;

					}

					else{

						if (uf<((1-q3)*q2*q1+q1*q2*q3+(1-q1)*q2*q3+(1-q2)*q1*q3)/normq2){

						fi1y=1;
						fi2y=0;
						fi3y=1;

						}
						else{

							fi1y=0;
							fi2y=1;
							fi3y=0;
						}
					}
				}
			}


		}
	}




}


//proposal count of occurrence of the regions, all 3 origins are licensed, only right (3)/ left (1) is not licensed


int obscMy,obsc3My,obsc1My;

//proposal count if the region


int obscpry[5],obscprch[5],obscpr3y[2],obscpr1y[2];


if (fi1[i6]*fi2[i6]*fi3[i6]==1){//all 3 origins are currently licensed

	if (fi1y*fi2y*fi3y==0) {//not all 3 origins of the proposal are licensed


		int obscc=obsc_cur(t123sum[i6],t21diff[i6],t23diff[i6],1,1,1,N1,N2);//obscuring status of the current profile

		//computing proposal count of the 'not all 3 origins are licensed' region


		for (int fif=0;fif<5;++fif) obscpry[fif]=obscpr[fif];

		obscpry[obscc]=obscpr[obscc]-1;

		//check it does not go below zero

		if (obscpry[obscc]<0){

			cout<<obscpr[obscc]<<" "<<obscpry[obscc]<<" "<<obscc<<" -1\n";

			exit(1);

		}

		//proposal count of 'all 3 origins are licensed'

		obscMy=obscM-1;


	}

}
else{//currently not all 3 origins are licensed

	if (fi1y*fi2y*fi3y==1) {//all 3 origins of the proposal are licensed


		int obscc=obsc_cur(t123sum[i6],t21diff[i6],t23diff[i6],1,1,1,N1,N2); //obscuring status of the proposal given the proposal

		//computing proposal count of the 'all 3 origins are licensed' region

		for (int fif=0;fif<5;++fif) obscpry[fif]=obscpr[fif];

		obscpry[obscc]=obscpr[obscc]+1;

		//check it does not go below zero


		if (obscpry[obscc]<0){

			cout<<obscpr[obscc]<<" "<<obscpry[obscc]<<" "<<obscc<<" +1\n";

			exit(1);

		}

		//proposal count of 'all 3 origins are licensed'

		obscMy=obscM+1;


	}
}

		if ((fi1[i6]==1)&&(fi3[i6]==0)){//currently only right origin is not licensed

			if ((fi1y!=1)||(fi3y!=0)){//proposed region is different from above



				int obscc=obsc_cur(t123sum[i6],t21diff[i6],t23diff[i6],1,1,0,N1,N2); //obscuring status of the current profile

				obscpr3y[0]=obscpr3[0];

				obscpr3y[1]=obscpr3[1];

				obscpr3y[obscc]=obscpr3[obscc]-1;

				//proposal count of 'only right origin is not licensed'

				obsc3My=obsc3M-1;


			}


		}



		if ((fi1y==1)&&(fi3y==0)){//only right origin of the proposal is not licensed

			if ((fi1[i6]!=1)||(fi3[i6]!=0)){//current region is different from above

				int obscc=obsc_cur(t123sum[i6],t21diff[i6],t23diff[i6],1,1,0,N1,N2); //obscuring status of the current profile given the proposal



				obscpr3y[0]=obscpr3[0];

				obscpr3y[1]=obscpr3[1];

				obscpr3y[obscc]=obscpr3[obscc]+1;

				//proposal count of 'only right origin is not licensed'

				obsc3My=obsc3M+1;



			}


		}


		if ((fi1[i6]==0)&&(fi3[i6]==1)){//currently only left origin is not licensed

			if ((fi1y!=0)||(fi3y!=1)){//proposed region is different from above

				int obscc=obsc_cur(t123sum[i6],t21diff[i6],t23diff[i6],0,1,1,N1,N2);

				obscpr1y[0]=obscpr1[0];

				obscpr1y[1]=obscpr1[1];

				obscpr1y[obscc]=obscpr1[obscc]-1;

				//proposal count of 'only left origin is are licensed'

				obsc1My=obsc1M-1;


			}


		}



		if ((fi1y==0)&&(fi3y==1)){//only left origin of the proposal is not licensed

			if ((fi1[i6]!=0)||(fi3[i6]!=1)){//current region is different from above


				int obscc=obsc_cur(t123sum[i6],t21diff[i6],t23diff[i6],0,1,1,N1,N2);



				obscpr1y[0]=obscpr1[0];

				obscpr1y[1]=obscpr1[1];

				obscpr1y[obscc]=obscpr1[obscc]+1;

				//proposal count of 'only left origin is are licensed'


				obsc1My=obsc1M+1;



			}


		}

		//collision points of the poposal


		int xsty1,xsty2;



xstar(t123sum[i6],t21diff[i6],t23diff[i6],N1,N2,xsty1,xsty2,fi1y,fi2y,fi3y);

//check hat the collision points are feasible

if (((xsty1<0)||(xsty2<0))||((xsty1>N1)||(xsty2>N2))){


	cout<<t123sum[i6]<<" t123 "<<t21diff[i6]<<" t21 "<<t23diff[i6]<<" t23 "<<xsty1<<" xsty1 "<<xsty2<<" xsty2 "<<fi1y<<" f1 "<<fi2y<<" f2 "<<fi3y<<" f3 4\n";

	exit(1);

}

//MH computation step

int minf1=min(xsty1,xst1[i6]);

int maxf1=max(xsty1,xst1[i6]);

int minf2=min(xsty2,xst2[i6]);

int maxf2=max(xsty2,xst2[i6]);





int sgn1=sign(xst1[i6]-xsty1);

int sgn2=sign(xst2[i6]-xsty2);

double Fmin1=0;

double Fmin2=0;





double prodfdiff1=0.0;

double prodfdiff2=0.0;

double prodfdiff1r=0.0;

double prodfdiff2r=0.0;

double prodfdiff1_log=0.0;

double prodfdiff2_log=0.0;

double prodfdiff1r_log=0.0;

double prodfdiff2r_log=0.0;


//identifying collision points affected by the proposal


int xsta1[4992];

xsta1[0]=minf1;

int count1=0;

for (i23=0;i23<L;++i23) {


	if ((xst1[i23]>minf1)&&(xst1[i23]<maxf1)) {

		count1+=1;

		xsta1[count1]=xst1[i23];



	}
}

//identifying parts of the profile affected by these collision points


Fmin1=F1[minf1];

Fmin2=F1[minf2+N1];


xsta1[count1+1]=maxf1;


//sorting the collision points identified above


qsort((void*)xsta1,count1+2,sizeof(int),compare_ints);


//MH step computation

//computing the part depending on the collision points only

//taking into account unsequencable region


for (i2=0;i2<count1+1;++i2){

	if (xsta1[i2]<minNA){

		int NA=min(xsta1[i2+1],minNA-1);

double incr1=log(((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))*((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))+(1-b)*(1.0/L)*((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))*sgn1)*(-NA+xsta1[i2])*log(1+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2));

double incr1r=log(((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))*((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))-(1-b)*(1.0/L)*((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))*sgn1)*(-NA+xsta1[i2])*log(1-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i2));


prodfdiff1_log+=incr1;

prodfdiff1r_log+=incr1r;

prodfdiff1+=(incr1-(1.0/tau)*(-NA+xsta1[i2])*log(1+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2)));

prodfdiff1r+=(incr1r-(1.0/tau)*(-NA+xsta1[i2])*log(1-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i2)));
	}

	if (xsta1[i2+1]>maxNA){

		int NA=max(xsta1[i2],maxNA+1);

double incr1=log(((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))*((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))+(1-b)*(1.0/L)*((1-b)*Fmin1+0.5*b+(1-b)*((double) i2/L))*sgn1)*(-xsta1[i2+1]+NA)*log(1+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2));

double incr1r=log(((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))*((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))-(1-b)*(1.0/L)*((1-b)*(1-Fmin1)+0.5*b-(1-b)*((double) i2/L))*sgn1)*(-xsta1[i2+1]+NA)*log(1-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i2));


prodfdiff1_log+=incr1;

prodfdiff1r_log+=incr1r;

prodfdiff1+=(incr1-(1.0/tau)*(-xsta1[i2+1]+NA)*log(1+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2)));

prodfdiff1r+=(incr1r-(1.0/tau)*(-xsta1[i2+1]+NA)*log(1-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i2)));
	}



	//check that log expression is above zero




	if (1.0+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2)<0) {


		cout<<" "<<1.0+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i2)<<" "<<b<<" log 1\n";

		cout<<prodfdiff1<<" prod\n";

		exit(1);
	}



}

//computing the part depending on the collision points and data xj[]


double at, atr;

for (i3=0;i3<count1+1;++i3){

	at=0;
atr=0;

	for (i24=xsta1[i3];i24<xsta1[i3+1];++i24) {


		at+=(2*logxj[i24]);
		atr+=(2*logxjr[i24]);

	}
	double incr_at1=at*log(1.0+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i3));

	double incr_at1r=atr*log(1.0-(1-b)*((double) sgn1)/(L*((1-b)*(1-Fmin1)+0.5*b)-(1-b)*i3));


	prodfdiff1+=incr_at1;

	prodfdiff1_log+=incr_at1;

	prodfdiff1r+=incr_at1r;

	prodfdiff1r_log+=incr_at1r;


	//check that log expression is above zero


	if ((1.0+(1-b)*((double) sgn1)/(L*((1-b)*Fmin1+0.5*b)+(1-b)*i3))<=0) {

		cout<<"log 2\n";
		exit(1);
	}


}

//the same computations for the right collision point


int xsta2[4992];


int count2=0;

xsta2[0]=minf2;

for (i23=0;i23<L;++i23) {

	if ((xst2[i23]>minf2)&&(xst2[i23]<maxf2)) {

		count2+=1;

		xsta2[count2]=xst2[i23];




	}
}

xsta2[count2+1]=maxf2;


qsort((void*)xsta2,count2+2,sizeof(int),compare_ints);


for (i2=0;i2<count2+1;++i2){



	if (xsta2[i2]<minNA-N1){

		int NA=min(xsta2[i2+1],minNA-N1-1);


double incr2=log(((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))*((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))+(1-b)*(1.0/L)*((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))*sgn2)*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))*(-NA+xsta2[i2]);

double incr2r=log(((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))*((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))-(1-b)*(1.0/L)*((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))*sgn2)*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i2))*(-NA+xsta2[i2]);


prodfdiff2_log+=incr2;

prodfdiff2r_log+=incr2r;

prodfdiff2+=(incr2-(1.0/tau)*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))*(-NA+xsta2[i2]));

prodfdiff2r+=(incr2r-(1.0/tau)*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i2))*(-NA+xsta2[i2]));

	}

	if (xsta2[i2+1]>maxNA-N1){

		int NA=max(xsta2[i2],maxNA-N1+1);


double incr2=log(((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))*((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))+(1-b)*(1.0/L)*((1-b)*Fmin2+0.5*b+(1-b)*(i2/(double) L))*sgn2)*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))*(-xsta2[i2+1]+NA);

double incr2r=log(((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))*((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))-(1-b)*(1.0/L)*((1-b)*(1-Fmin2)+0.5*b-(1-b)*(i2/(double) L))*sgn2)*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i2))*(-xsta2[i2+1]+NA);


prodfdiff2_log+=incr2;

prodfdiff2r_log+=incr2r;

prodfdiff2+=(incr2-(1.0/tau)*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))*(-xsta2[i2+1]+NA));

prodfdiff2r+=(incr2r-(1.0/tau)*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i2))*(-xsta2[i2+1]+NA));

	}


	if ((1.0+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i2))<=0) {

		cout<<b<<" log 3\n";
		exit(1);
	}

}




for (i3=0;i3<count2+1;++i3){

	at=0;

	atr=0;

	for (i24=xsta2[i3]+N1;i24<xsta2[i3+1]+N1;++i24) {

		at+=(2*logxj[i24]);

		atr+=(2*logxjr[i24]);
	}



	double incr_at2=at*log(1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i3));

	double incr_at2r=atr*log(1-(1-b)*((double) sgn2)/(L*((1-b)*(1-Fmin2)+0.5*b)-(1-b)*i3));



	prodfdiff2+=incr_at2;

	prodfdiff2r+=incr_at2r;

	prodfdiff2_log+=incr_at2;

	prodfdiff2r_log+=incr_at2r;



	if (1+(1-b)*((double) sgn2)/(L*((1-b)*Fmin2+0.5*b)+(1-b)*i3)<=0) {

		cout<<"log 4\n";
		exit(1);
	}


}



int sk=0;

//criteria for eliminating flat profiles


if ((((xsty1==0)&&(xsty2==0))||((xsty1==N1)&&(xsty2==N2)))) sk=1;

if (sk==0){

	//acceptance probability

double	pow_alpha=(prodfdiff1+prodfdiff2+prodfdiff1r+prodfdiff2r)*(tau/2.0);


//acceptance-rejection

if (pow_alpha>0){


	//updating the counts of tree regions: all origins are licensed, left/right origin is not licensed


	if ((((fi1[i6]+fi2[i6]+fi3[i6]==3)&&(fi1y+fi2y+fi3y<3))||((fi1[i6]+fi2[i6]+fi3[i6]<3)&&(fi1y+fi2y+fi3y==3)))){

	obscM=obscMy;



	for (int ia=0;ia<5;++ia) obscpr[ia]=obscpry[ia];
	}


	if ((((fi1[i6]==0)&&(fi3[i6]==1))&&((fi1y!=0)||(fi3y!=1)))||(((fi1[i6]!=0)||(fi3[i6]!=1))&&((fi1y==0)&&(fi3y==1)))){

	obsc1M=obsc1My;

	for (int ia1=0;ia1<2;++ia1){

		obscpr1[ia1]=obscpr1y[ia1];

	}
	}


	if ((((fi1[i6]==1)&&(fi3[i6]==0))&&((fi1y!=1)||(fi3y!=0)))||(((fi1[i6]!=1)||(fi3[i6]!=0))&&((fi1y==1)&&(fi3y==0)))){

	obsc3M=obsc3My;

	for (int ia1=0;ia1<2;++ia1){

		obscpr3[ia1]=obscpr3y[ia1];

	}
	}


	//updating the number of times each origin is active


	sumfi3=sumfi3-fi3[i6]+fi3y;

	sumfi1=sumfi1-fi1[i6]+fi1y;

	sumfi2=sumfi2-fi2[i6]+fi2y;

	//updating mean values and sd of the proposals for time differences


	diff21s=diff21s*diff21s*(dp21-1)+t21diff[i6]*t21diff[i6]*(fi2y*fi1y-fi2[i6]*fi1[i6])+dp21*diff21m1*diff21m1;

	diff23s=diff23s*diff23s*(dp23-1)+t23diff[i6]*t23diff[i6]*(fi2y*fi3y-fi3[i6]*fi2[i6])+dp23*diff23m1*diff23m1;

	diff31s=diff31s*diff31s*(dp31-1)+t31diff[i6]*t31diff[i6]*(fi3y*fi1y-fi3[i6]*fi1[i6])+dp31*diff31m1*diff31m1;


	diff21m1=diff21m1*dp21+t21diff[i6]*(fi1y*fi2y-fi1[i6]*fi2[i6]);

	diff23m1=diff23m1*dp23+t23diff[i6]*(fi3y*fi2y-fi3[i6]*fi2[i6]);

	diff31m1=diff31m1*dp31+t31diff[i6]*(fi1y*fi3y-fi1[i6]*fi3[i6]);

	diff21m=diff21m*dp21+temp2[i6]*(fi1y*fi2y-fi1[i6]*fi2[i6]);

	diff23m=diff23m*dp23+temp3[i6]*(fi3y*fi2y-fi3[i6]*fi2[i6]);

	diff31m=diff31m*dp31+temp4[i6]*(fi1y*fi3y-fi1[i6]*fi3[i6]);


	//updating the count of time differences being licensed

	dp21=dp21-fi2[i6]*fi1[i6]+fi2y*fi1y;

	dp23=dp23-fi2[i6]*fi3[i6]+fi2y*fi3y;

	dp31=dp31-fi3[i6]*fi1[i6]+fi3y*fi1y;








	diff21m=diff21m/(double) dp21;

	diff23m=diff23m/(double) dp23;

	diff31m=diff31m/(double) dp31;

	diff21m1=diff21m1/(double) dp21;

	diff23m1=diff23m1/(double) dp23;

	diff31m1=diff31m1/(double) dp31;


	diff21s=sqrt((diff21s-dp21*diff21m1*diff21m1)/(double) (dp21-1));

	diff23s=sqrt((diff23s-dp23*diff23m1*diff23m1)/(double) (dp23-1));

	diff31s=sqrt((diff31s-dp31*diff31m1*diff31m1)/(double) (dp31-1));


	//updating licensing indicators

	fi1[i6]=fi1y;

	fi2[i6]=fi2y;

	fi3[i6]=fi3y;

	//updating profile


	 for (i9=minf1;i9<maxf1;++i9) {


		 F1[i9]=F1[i9]+(1.0/L)*sgn1;


	 }

	 for (i8=minf2+N1;i8<maxf2+N1;++i8) {


		 F1[i8]=F1[i8]+(1.0/L)*sgn2;


	 }

	 //updating the squared term

	 logvec2-=(prodfdiff1_log+prodfdiff2_log+prodfdiff1r_log+prodfdiff2r_log);

	 //updating the collision points


	xst1[i6]=xsty1;

	xst2[i6]=xsty2;



}

else{

	//the same updates in case one has to sampled from uniform[0,1] to accept/reject


	 double u=unifrnd(gen);
	 double alpha=exp(pow_alpha);


	 if (u<=alpha){


			if ((((fi1[i6]+fi2[i6]+fi3[i6]==3)&&(fi1y+fi2y+fi3y<3))||((fi1[i6]+fi2[i6]+fi3[i6]<3)&&(fi1y+fi2y+fi3y==3)))){

			obscM=obscMy;



			for (int ia=0;ia<5;++ia) obscpr[ia]=obscpry[ia];
			}


			if ((((fi1[i6]==0)&&(fi3[i6]==1))&&((fi1y!=0)||(fi3y!=1)))||(((fi1[i6]!=0)||(fi3[i6]!=1))&&((fi1y==0)&&(fi3y==1)))){

			obsc1M=obsc1My;

			for (int ia1=0;ia1<2;++ia1){

				obscpr1[ia1]=obscpr1y[ia1];

			}
			}


			if ((((fi1[i6]==1)&&(fi3[i6]==0))&&((fi1y!=1)||(fi3y!=0)))||(((fi1[i6]!=1)||(fi3[i6]!=0))&&((fi1y==1)&&(fi3y==0)))){

			obsc3M=obsc3My;

			for (int ia1=0;ia1<2;++ia1){

				obscpr3[ia1]=obscpr3y[ia1];

			}
			}

			sumfi3=sumfi3-fi3[i6]+fi3y;

			sumfi1=sumfi1-fi1[i6]+fi1y;

			sumfi2=sumfi2-fi2[i6]+fi2y;



			diff21s=diff21s*diff21s*(dp21-1)+t21diff[i6]*t21diff[i6]*(fi2y*fi1y-fi2[i6]*fi1[i6])+dp21*diff21m1*diff21m1;

			diff23s=diff23s*diff23s*(dp23-1)+t23diff[i6]*t23diff[i6]*(fi2y*fi3y-fi3[i6]*fi2[i6])+dp23*diff23m1*diff23m1;

			diff31s=diff31s*diff31s*(dp31-1)+t31diff[i6]*t31diff[i6]*(fi3y*fi1y-fi3[i6]*fi1[i6])+dp31*diff31m1*diff31m1;


			diff21m1=diff21m1*dp21+t21diff[i6]*(fi1y*fi2y-fi1[i6]*fi2[i6]);

			diff23m1=diff23m1*dp23+t23diff[i6]*(fi3y*fi2y-fi3[i6]*fi2[i6]);

			diff31m1=diff31m1*dp31+t31diff[i6]*(fi1y*fi3y-fi1[i6]*fi3[i6]);

			diff21m=diff21m*dp21+temp2[i6]*(fi1y*fi2y-fi1[i6]*fi2[i6]);

			diff23m=diff23m*dp23+temp3[i6]*(fi3y*fi2y-fi3[i6]*fi2[i6]);

			diff31m=diff31m*dp31+temp4[i6]*(fi1y*fi3y-fi1[i6]*fi3[i6]);

			dp21=dp21-fi2[i6]*fi1[i6]+fi2y*fi1y;

			dp23=dp23-fi2[i6]*fi3[i6]+fi2y*fi3y;

			dp31=dp31-fi3[i6]*fi1[i6]+fi3y*fi1y;





			diff21m=diff21m/(double) dp21;

			diff23m=diff23m/(double) dp23;

			diff31m=diff31m/(double) dp31;

			diff21m1=diff21m1/(double) dp21;

			diff23m1=diff23m1/(double) dp23;

			diff31m1=diff31m1/(double) dp31;


			diff21s=sqrt((diff21s-dp21*diff21m1*diff21m1)/(double) (dp21-1));

			diff23s=sqrt((diff23s-dp23*diff23m1*diff23m1)/(double) (dp23-1));

			diff31s=sqrt((diff31s-dp31*diff31m1*diff31m1)/(double) (dp31-1));



			fi1[i6]=fi1y;

			fi2[i6]=fi2y;

			fi3[i6]=fi3y;

			 for (i9=minf1;i9<maxf1;++i9) {


				 F1[i9]=F1[i9]+(1.0/L)*sgn1;


			 }

			 for (i8=minf2+N1;i8<maxf2+N1;++i8) {


				 F1[i8]=F1[i8]+(1.0/L)*sgn2;


			 }


			 logvec2-=(prodfdiff1_log+prodfdiff2_log+prodfdiff1r_log+prodfdiff2r_log);


			xst1[i6]=xsty1;

			xst2[i6]=xsty2;



	 }

}

//check that the squared term is above zero


if (logvec2<0) {

	cout<<" f logvec2\n";

	exit(1);
}




}



}

//obscuring count for the case 'all three origins are active'


cout<<obscpr[0]<<" "<<obscpr[1]<<" "<<obscpr[2]<<" "<<obscpr[3]<<" "<<obscpr[4]<<" obscpr f\n";


//write licensing indicators into a file every tenth iteration

if ((i1+1)%10==0){

	for (i24=0;i24<L;++i24) f7<<fi1[i24]<<" "<<fi2[i24]<<" "<<fi3[i24]<<"\n";

}



}

f1.close();
f2.close();
f4.close();
f5.close();
f6.close();
f7.close();
f8.close();

//how much time it takes to run it


time=clock()-time;

cout<<"clicks "<<time<<" seconds "<<(float) time/CLOCKS_PER_SEC<<"\n";
}


namespace po = boost::program_options;

using Floats  = std::vector<std::pair<std::string, float>>;
using Ints    = std::vector<std::pair<std::string, int>>;
using Strings = std::vector<std::pair<std::string, std::string>>;

namespace std {
    template <typename V> static inline std::istream& operator>>(std::istream& is, std::pair<std::string, V>& into) {
        char ch;
        while (is >> ch && ch!='=') into.first += ch;
        return is >> into.second;
    }
}



int main(int argc, char** argv)

{
	int i,j;
	double xj[4000], xjr[4000]; //data from the forward xj[] and reverse xjr[] strand of a corresponding size

  po::variables_map vm;
  po::options_description desc("Allowed Options");

  desc.add_options()
  	("help,h", "Help screen")
 //   ("chr,c", po::value<int>()->required(), "Enter the chromosome number")
    ("rat,r", po::value<int>()->default_value(0), "Rat or not")
    ("wt,w", po::value<int>()->default_value(0), "If Rat, WT or not")
    ("left,l", po::value<int>()->required(), "Number of left origin");

  // parse arguments and save them in the variable map (vm)
  po::store(po::parse_command_line(argc, argv, desc), vm);
      if (vm.count("help")){
      std::cout << desc << '\n';
  } else {

  //std::cout << "Hello " << vm["chr"].as<int>() << std::endl;

  int w=vm["wt"].as<int>();
  int rat=vm["rat"].as<int>();
  int l=vm["left"].as<int>();

int N1,N2;
int minNA,maxNA;
	
	//distances between origins

if (rat==0){

		if (l==1014){
//	int 
		N1=510,N2=1957; //ARS1014,1015,1018
		minNA=1098,maxNA=1336; //ARS1014,1015,1018

	} else {
//	int 
		N1=1957,N2=1446; //ARS1015, 1018, 1019
		minNA=1098-510,maxNA=1336-510; //ARS1015,1018,1019
	}

} else {
	if (l==1014){
	//	chromosome 10 data for rat1 example


//	int 
N1=503,N2=1960; //ARS1014,1015,1018

minNA=1108,maxNA=1336; //ARS1014,1015,1018 rat1 data

} else{
//	int 
N1=1960,N2=1446; //ARS1015, 1018, 1019

minNA=1108-503,maxNA=1336-503;
}
}
	int M=4992; //parameter M, number of profile realisations

	//	int it=1000000; //amount of iteration in case of using simulations

//ends of the unsequencable region

//	int minNA=1098,maxNA=1336; //ARS1014,1015,1018
//int minNA=1098-510,maxNA=1336-510; //ARS1015,1018,1019

//	int minNA=1108,maxNA=1336; //ARS1014,1015,1018 rat1 data
//	int minNA=1108-503,maxNA=1336-503; //ARS1015,1018,1019 rat1 data

	ifstream fi;//file for data from the forward strand
	ifstream fir;//file for data from the reverse strand

	ostringstream convert;



	convert<<l;

	string le=convert.str();

	if (rat!=1){
	fi.open("chr10_"+le+"_f.txt");
	fir.open("chr10_"+le+"_r.txt");
} else {
	if (w==0){

		fi.open("chr10_rat1_"+le+"_f.txt");
	fir.open("chr10_rat1_"+le+"_r.txt");
} else {
		fi.open("chr10_wt_"+le+"_f.txt");
	fir.open("chr10_wt_"+le+"_r.txt");

}
}


	for (i=0;i<N1+N2;++i){
		xj[i]=0.0;
		xjr[i]=0.0;

	}

	//	simr(it,N1,N2,xjr,0.7,1.0/sqrt(100)); //simulating reverse strand

	//	sim(-3000.0/3.0,6000.0/3.0,-3000.0/3.0,2000,2000,2000,it,N1,N2,xj,0.05,1.0/sqrt(100),1,1,1); //simulating forward strand

	//writing data into arrays

	for (j=0;j<N1+N2;++j) {


		fi>>xj[j];

		fir>>xjr[j];


	}


	fi.close();
	fir.close();

	MCMC(xj,xjr,M,N1,N2,10000000,minNA,maxNA); //running the algorithm



}


	return 0;
}
