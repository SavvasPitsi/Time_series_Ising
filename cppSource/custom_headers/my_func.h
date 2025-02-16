#define PI 3.1415926535897932384626433832795

#include <math.h>
#include <random>
#include <time.h>

// double random(int max) 			//simple rand()																													//Random number generator int threshold (not inclusive)
// {
// 	return ((double)rand() / (RAND_MAX + 1.0))*max;
// } 

// double random(int max) { //Mersenne twister algorithm
//     // static std::mt19937 rng(std::random_device{}());
// 	static std::mt19937 rng(time(NULL));
// 	std::uniform_real_distribution<double> dist(0.0, max - 1e-10);  
//     return dist(rng);
// }

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
