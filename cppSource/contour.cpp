#include <iostream> 
#include <stdlib.h>
#include <math.h>
#include <string>
#include <time.h>
#include <thread>
#include <mutex>
#include <matplotlibcpp.h>
#include <my_func.h> //PI, random(max), mean(values,size), std_dev(values,size)

using namespace std;

const int n = 1024, smoothing_times = 5, N = 1000, steps = 25; //resolution: steps x steps
const double temp_step_width = 2 * 0.045 /*0.045*/, magn_step_width = 2.3 * 0.04/*0.04-0.006*/;
double T[steps], B[steps], magnetization[steps][steps], criticalT[steps][steps];
int counter = 0, counter2 = 0, counter3 = 0, total = 0; 

clock_t startTime = clock();

int uu(int v) { return floor(v / steps); }
int aa(int v) { return floor(v % steps); }

bool check(double TT, double BB)																											//check if parameters need smoothing
{
	return (TT < 2.205*exp(-1.5*BB) + 1 / (0.1*BB + 6) - 0.1 - 0.005*exp(BB) || TT < 1.35 - (1.35 / 1.95)*BB) && TT > 2.205*exp(-1.4*BB) + 1 / (BB + 5) - 0.7 + 0.21*BB;
}

void MagCalc(mutex & mtx, mutex & mtx2, mutex & mtx3)
{
	int i = 0, j = 0, s = 0, k = 0, l = 0, b = 0, c = 0, a = 0, u = 0, uptime = 0, flag = 0;
	bool interesting = false;
	double pin[n][n], h[n][n], neigh[4], M[N + 1];
	double S = 0, BB = 0, TT = 0, sum = 0, DE = 0, r = 0, J = 1.0, beta = 1.0, progress = 0;

	srand(time(NULL) + std::hash<std::thread::id>{}(std::this_thread::get_id()));
	
	while (counter < steps*steps && counter2 < steps*steps)
	{
		mtx2.lock();
		a = aa(counter2);
		u = uu(counter2);
		counter2++;
		mtx2.unlock();

		BB = B[a];
		TT = T[u];

		interesting = check(TT, BB);

		mtx.lock();
		uptime = double(clock() - startTime) / (double)CLOCKS_PER_SEC;
		progress = 1.0*counter3 / (1.0*total);
		if (progress > 1) { progress = 1.0; }
		cout << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t \r";
		cout << "Uptime: t = " << floor(uptime) << " sec \t(" << floor(progress*total) << " / " << floor(total) << ", " << floor(100 * progress) << " %) Time remaining: " << floor(1.0*uptime / progress) - floor(uptime) << " sec\t" << "T = " << TT << "\tB = " << BB << "\r";
		mtx.unlock();
		

		for (i = 0; i < n; i++)				                                                       														//Magnetic field
		{

			for (j = 0; j < n; j++)
			{
				h[i][j] = -BB;
				//if (i < n / 2) h[i][j] = -1.0*h[i][j];																										//field topology
			}
		}

		beta = 1.0 / TT;

		if (interesting)
		{
			c = 0;
		}
		else
		{
			c = smoothing_times - 1;
		}

		while (c < smoothing_times)
		{

			for (i = 0; i < n; i++)                                                           																		//Initial lattice
			{

				for (j = 0; j < n; j++)
				{
					//m = floor(random(2));  //T=inf
					//if (m == 2) m = 1;
					//pin[i][j] = 2 * m - 1;

					pin[i][j] = 1;

					//pin[i][j] = 1;
					//if (i < n / 2) pin[i][j] = -1;
				}
			}

			M[0] = 0;
			for (i = 0; i < n; i++)                                                           																//Initial magnetization 
			{
				for (j = 0; j < n; j++)
				{
					M[0] = M[0] + pin[i][j];
				}
			}

			for (s = 0; s < N; s++) 																														//Sweep loop
			{
				M[s + 1] = M[s];																														//Arxikh magnhtish


				for (b = 0; b < n*n; b++)
				{

					k = random(n); l = random(n);


					k = k + 1; if (k == n) { k = 0; } neigh[0] = pin[k][l];
					k = k - 1; if (k == -1) { k = n - 1; } k = k - 1; if (k == -1) { k = n - 1; } neigh[1] = pin[k][l]; k = k + 1; if (k == n) { k = 0; }
					l = l + 1; if (l == n) { l = 0; } neigh[2] = pin[k][l];
					l = l - 1; if (l == -1) { l = n - 1; } l = l - 1; if (l == -1) { l = n - 1; } neigh[3] = pin[k][l]; l = l + 1; if (l == n) { l = 0; }

					for (i = 0; i < 4; i++)
					{
						S = S + neigh[i];
					}

					DE = 2 * J*S*pin[k][l] + 2 * h[k][l] * pin[k][l];

					if (DE > 0)
					{
						r = random(1);

						if (r < exp(-beta * DE))
						{
							pin[k][l] = -1 * pin[k][l];
							M[s + 1] = M[s + 1] + 2 * pin[k][l];
						}
					}

					else
					{
						pin[k][l] = pin[k][l] * (-1);
						M[s + 1] = M[s + 1] + 2 * pin[k][l];
					}
					S = 0;

				}

				if (flag == 0 && M[s + 1] < 0) { criticalT[a][u] = s; flag = 1; break; }																			//measurement of time to flip
			}

			flag = 0;


			for (i = 0; i < N; i++)																															//Sweep counter
			{
				M[i] = M[i + 1] / (n*n);
			}

			sum = sum + mean(M, N);
			c++;
		}
		if (interesting)
		{
			magnetization[a][u] = sum / (1.0*smoothing_times);

			mtx3.lock();
			counter3 = counter3 + smoothing_times;
			mtx3.unlock();
		}
		else
		{
			magnetization[a][u] = sum;

			mtx3.lock();
			counter3++;
			mtx3.unlock();
		}

		//magnetization[a][u] = (N / 2)*(1 + magnetization[a][u]);															//magnetization to flip time conversion

		sum = 0;
		

		mtx2.lock();
		counter++;
		mtx2.unlock();

		mtx.lock();
		uptime = double(clock() - startTime) / (double)CLOCKS_PER_SEC;
		progress = 1.0*counter3 / (1.0*total);
		if (progress > 1) { progress = 1.0; }
		cout << "\t\t\t\t\t\t\t\t\t\t\t\t\t\t \r";
		cout << "Uptime: t = " << floor(uptime) << " sec \t(" << floor(total*progress) << " / " << total << ", " << floor(100. * progress) << " %) Time remaining: " << floor(1.0*uptime / progress) - floor(uptime) << " sec\t" << "T = " << TT << "\tB = " << BB << "\r";
		mtx.unlock();
	}
}




int main()
{
	cout << "This program will calculate the magnetization over the B, T parameter space" << endl;
	cout << "\nInitializing MATLAB engine" << endl;

	// Engine *ep;
	// ep = engOpen("NULL");																															//OPEN ENGINE
	// engEvalString(ep, "clear all;");

	int i = 0, j = 0, m = 0;
	
	mutex mtx, mtx2, mtx3;

	T[0] = 2.205;																																		//initial temperature
	for (i = 0; i<steps; i++)																															//temperatures
	{
		T[i + 1] = T[i] - temp_step_width;
	}
	B[0] = 0;
	for (j = 0; j<steps; j++)
	{
		B[j + 1] = B[j] + magn_step_width;
	}


	for (i = 0; i < steps; i++)
	{
		for (j = 0; j < steps; j++)
		{
			magnetization[i][j] = 10;
			criticalT[i][j] = N;
		}
	}

	for (i = 0; i < steps; i++)
	{
		for (j = 0; j < steps; j++)
		{
			if (check(T[i], B[j]))
			{
				total = total + smoothing_times;
			}
			else
			{
				total++;
			}
		}
	}
	cout << "\n________________________________________________________________" << endl;
	cout << "Lattice size: " << floor(n) << endl;
	cout << "Parameter range: T = " << T[steps - 1] << " to " << T[0] << ", B = " << B[0] << " to " << B[steps - 1] << endl;
	cout << "Total sweeps per parameter set: " << floor(N)  << endl;
	cout << "Smoothing factor: " << floor(smoothing_times) << endl;
	cout << "Parameter space resolution: " << floor(steps) << " x " << floor(steps) << " points" << endl;
	cout << "________________________________________________________________" << endl;
	
	cout << "\nCalculating...\n\n";


	thread t1{ MagCalc, ref(mtx), ref(mtx2), ref(mtx3) };
	thread t2{ MagCalc, ref(mtx), ref(mtx2), ref(mtx3) };
	thread t3{ MagCalc, ref(mtx), ref(mtx2), ref(mtx3) };
	thread t4{ MagCalc, ref(mtx), ref(mtx2), ref(mtx3) };

	t1.join();
	t2.join();
	t3.join();
	t4.join();
	

	cout << "\n\nPlotting... \n\n";

	// mxArray* POINTER6 = mxCreateDoubleMatrix(steps, steps, mxREAL);
	// memcpy((void *)mxGetPr(POINTER6), (void *)criticalT, sizeof(double)*steps*steps);
	// engPutVariable(ep, "FlippingT", POINTER6);

	// mxArray* POINTER7 = mxCreateDoubleMatrix(steps, steps, mxREAL);
	// memcpy((void *)mxGetPr(POINTER7), (void *)magnetization, sizeof(double)*steps*steps);
	// engPutVariable(ep, "Surface", POINTER7);

	// mxArray* POINTER8 = mxCreateDoubleMatrix(1, steps, mxREAL);
	// memcpy((void *)mxGetPr(POINTER8), (void *)T, sizeof(double)*steps);
	// engPutVariable(ep, "T", POINTER8);

	// mxArray* POINTER9 = mxCreateDoubleMatrix(1, steps, mxREAL);
	// memcpy((void *)mxGetPr(POINTER9), (void *)B, sizeof(double)*steps);
	// engPutVariable(ep, "B", POINTER9);



	// engEvalString(ep, "figure(1);contourf(B,T,Surface);ylabel('T');xlabel('B');colorbar;");
	// engEvalString(ep, "title(['Magnetization'],'interpreter','latex');");

	// engEvalString(ep, "figure(2);surf(B,T,Surface);ylabel('T');xlabel('B');colorbar;");
	// engEvalString(ep, "title(['Magnetization'],'interpreter','latex');");

	// engEvalString(ep, "figure(3);contourf(B,T,FlippingT);ylabel('T');xlabel('B');colorbar;");
	// engEvalString(ep, "title(['Flipping time'],'interpreter','latex');");



	cout << "\nDone! \n" << endl;


	cout << "\n--------------------------------\n";
	cout << "Process exited after " << floor(((double)clock() - startTime) / (double)CLOCKS_PER_SEC) << " seconds." << endl;
	system ("pause");
	return 0;
}