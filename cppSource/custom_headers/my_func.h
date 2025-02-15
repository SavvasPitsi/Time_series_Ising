#define PI 3.1415926535897932384626433832795

#include <math.h>

double random(int max) 																																//Random number generator int threshold (not inclusive)
{
	//double RANDM = ((double) rand() / (RAND_MAX+1))*max;
	//if(floor(RANDM) == max) { RANDM = max - 1; } 
	//RAND_MAX + 1, remove check
	//return RANDM;
	return ((double)rand() / (RAND_MAX + 1.0))*max;
} 
 

double mean(int values[], int size)																													//mean(int)
{
	int sum=0, i;
	double mn=0;
	for(i=0;i<size;i++)
	{
		sum += values[i];
	}
	mn=(1.0*sum)/size;
	return mn;	
}

double mean(double values[], int size)																												//mean(double)
{
	int i;
	double mn=0, sum=0;
	for(i=0;i<size;i++)
	{
		sum += values[i];
	}
	mn=sum/size;
	return mn;
}

double std_dev(int values[], int size)																												//st_dev(int)
{
	double m=mean(values,size),s=0,r=0;
	int i;
	for(i=0;i<size;i++)
	{
		s=s+pow((m-values[i]),2);
	}
	r=sqrt(s/size);
	return r;
}

double std_dev(double values[], int size)																											//st_dev(double)
{
	double m=mean(values,size),s=0,r=0;
	int i;
	for(i=0;i<size;i++)
	{
		s=s+pow((m-values[i]),2);
	}
	r=sqrt(s/size);
	return r;
}
