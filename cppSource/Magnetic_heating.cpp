#include <iostream> 
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "matplotlibcpp.h"
#include <my_func.h> //PI, random(max), mean(values,size), std_dev(values,size)																		//delete one of the two

using namespace std;
namespace plt = matplotlibcpp;

const int n = 256, N = 1000, steps = 1; //max 15 for figure
int pin[n][n];
double h[n][n], neigh[4], M[N + 1], H[N + 1], M2[N], H2[N], tsweep[N], info[7], B[steps], meanM[steps], heat_cpty[steps], Ni[steps], criticalT[steps];
std::vector<float> vecMat(n * n);
const float* vecPtr = &(vecMat[0]);

void plot_fig(std::vector<double>& v,int figNum = 1, std::string s="Title")
{
    plt::figure(figNum);
    plt::clf();
    plt::plot(v);
    plt::title(s);
    plt::show(false);
    plt::pause(0.1);
}

void plot_fig(const float* p,int figNum = 1, std::string s="Title")
{
    plt::figure(figNum);
    plt::clf();
    plt::imshow(p, n, n, 1);
    plt::title(s);
    plt::show(false);
    plt::pause(0.001);
}


int main()
{
	int i = 0, j = 0, s = 0, k = 0, l = 0, m = 0, u = 0, v = 0, c = 0, count = 0, flag = 0, t_equil = 1, ii = 0, jj = 0;
	int vecI = 0, vecJ = 0;
	double S = 0, DE = 0, r = 0, J = 1.0, step_width = 0.0, T, beta = 0;
	T = 1.90;																																			//initial temperature
	B[0] = 0.1;
	clock_t startTime = clock();


	for (u = 0; u<steps; u++)																															//temperatures
	{
		B[u + 1] = B[u] - step_width;
	}
	//heating 0.4 < B_crit < 0.5
	//cooling 0.009< B_crit <0.010 (128) (256) , 0.01, 0.008 both behaviours

	for (u = 0; u<steps; u++)																															//uncorrelated step calculation
	{
		//Ni[u]=1.0*floor(2*1203*exp(-pow(((T[u]-2.16)/0.1796),2))); //1_Autocorrelation_times
		//if (T[u]<2.35) //2_Autocorrelation_times
		//{
		//	Ni[u] = 1.0*floor(2 * (1600 * exp(-pow(((T[u] - 2.263) / 0.11), 2)) + 125));
		//}
		//else
		//{
		//	Ni[u] = 1.0*floor(2 * (1600 * exp(-pow(((T[u] - 2.263) / 0.11), 2)) + 125)*exp(-(T[u] - 2.35) / 0.05));
		//}
		//if (Ni[u]<125) { Ni[u] = 125.0; }
		//Ni[u] = 125.0;
		Ni[u] = 1;
	}



	// cout << "Initializing MATLAB engine... ";

	srand(time(NULL));
	// Engine *ep;
	// ep = engOpen("NULL");																															//OPEN ENGINE
	// engEvalString(ep, "clear all;");
	// //engEvalString(ep, "close all;");	

	// cout << "\rMATLAB initialized \t \t \n";

    std::cout << "We will use n = " << n << " and go for " << N << " steps\n";

	H[0] = 0;
	for (i = 0; i<n; i++)                                                           																		//Initial hamiltonian
	{

		for (j = 0; j<n; j++)
		{
			k = k + 1; 
            if (k == n) { 
                k = 0; 
            } 
            neigh[0] = pin[k][l];
			k = k - 1; 
            if (k == -1) { 
                k = n - 1; 
            } 
            k = k - 1; 
            if (k == -1) { 
                k = n - 1; 
            } 
            neigh[1] = pin[k][l]; 
            k = k + 1; 
            if (k == n) { 
                k = 0; 
            }
			l = l + 1; 
            if (l == n) { 
                l = 0; 
            } 
            neigh[2] = pin[k][l];
			l = l - 1; 
            if (l == -1) { 
                l = n - 1; 
            } 
            l = l - 1; 
            if (l == -1) { 
                l = n - 1; 
            } 
            neigh[3] = pin[k][l]; 
            l = l + 1; 
            if (l == n) { 
                l = 0; 
            }

			for (s = 0; s<4; s++)
			{
				S = S + neigh[s];
			}

			H[0] = H[0] + pin[i][j] * h[i][j] + J * S*pin[i][j];
			S = 0;
		}
	}



	for (i = 0; i<4; i++)																																//Arxikopoihsh athroismatos geitonwn
	{
		neigh[i] = 0;
	}

	for (i = 0; i<N; i++)																																//Sweep counter
	{
		tsweep[i] = i;
	}

	// mxArray* POINTER1 = mxCreateDoubleMatrix(n, n, mxREAL);

	cout << "Calculating...\n\n";

	//engEvalString(ep, "figure(1);clf;figure(2);clf;figure(3);clf;figure(4);clf;figure(5);clf;figure(6);clf;");
	//engEvalString(ep, "figure(1);clf;figure(2);clf;figure(4);clf;");

	for (u = 0; u<steps; u++)
	{
		beta = 1.0 / T;
		criticalT[u] = N+1;

		for (i = 0; i<n; i++)				                                                       															//Magnetic field
		{
			for (j = 0; j<n; j++)
			{
				// if (i < n / 2) { h[i][j] = -B[u]; }
				// else { h[i][j] = -1 * B[u]; }
                h[i][j] = -B[u];
			}
		}

		for (i = 0; i<n; i++)                                                           																		//Times plegmatos
		{

			for (j = 0; j<n; j++)
			{
				//m = floor(random(2));  //T=inf
				//if (m == 2) m = 1;
				//pin[i][j] = 2 * m - 1;

				pin[i][j] = 1;
				// if (i < 12 && j < 12) { pin[i][j] = -1; }

				//if (i < n / 8 && j < n / 8) { pin[i][j] = -1; }
				//else { pin[i][j] = 1; }
			}
		}

		M[0] = 0;
		for (i = 0; i<n; i++)                                         	                  																	//Ypologismos arxikhs magnhtishs 
		{
			for (j = 0; j<n; j++)
			{
				M[0] = M[0] + pin[i][j];
			}
		}


		cout << "Step " << u + 1 << " of " << steps << " (B = " << B[u] << ", T = " << T << ") \r\n";
		for (s = 0; s < N; s++)																															//Sweep loop
		{
			M[s + 1] = M[s];																															//Arxikh magnhtish
			H[s + 1] = H[s];
			for (v = 0; v<Ni[u]; v++)																													//correlation sweep
			{
				for (ii = 0; ii<n*n; ii++)                                                           															//Times plegmatos
				{
					//c = floor(random((n*n)-2));
					//k = c / n; l = c % n;
					k = random(n); l = random(n);
			

					k = k + 1; if (k == n) { k = 0; } neigh[0] = pin[k][l];
					k = k - 1; if (k == -1) { k = n - 1; } k = k - 1; if (k == -1) { k = n - 1; } neigh[1] = pin[k][l]; k = k + 1; if (k == n) { k = 0; }
					l = l + 1; if (l == n) { l = 0; } neigh[2] = pin[k][l];
					l = l - 1; if (l == -1) { l = n - 1; } l = l - 1; if (l == -1) { l = n - 1; } neigh[3] = pin[k][l]; l = l + 1; if (l == n) { l = 0; }

					for (i = 0; i<4; i++)
					{
						S = S + neigh[i];
					}

					DE = 2 * J*S*pin[k][l] + 2 * h[k][l] * pin[k][l];

					if (DE>0)
					{
						r = random(1);
						if (r == 1) r = 0;

						if (r<exp(-beta * DE))
						{
							pin[k][l] = -1 * pin[k][l];
							M[s + 1] = M[s + 1] + 2 * pin[k][l];
							H[s + 1] = H[s + 1] + DE;
						}
					}

					else
					{
						pin[k][l] = pin[k][l] * (-1);						
						M[s + 1] = M[s + 1] + 2 * pin[k][l];
						H[s + 1] = H[s + 1] + DE;
					}
					S = 0;


				}
				if (count < t_equil)																													//counter for initial equil
				{
					v = v - 1;
					
					if (count == 0) { cout << "Approaching equilibrium: " << t_equil << " sweeps \t Total measurements: " << N << " (Take 1 every " << Ni[u]<< " sweeps)\t\n"; }
					count++;
				}

				for(vecI= 0; vecI<n; vecI++)
				{
					for(vecJ= 0; vecJ<n; vecJ++)
					{
						vecMat.at(vecI * n + vecJ) = 1.0*pin[vecI][vecJ];
					}
				}
				plot_fig(vecPtr, 1, "Lattice N = " + std::to_string(s+1));
				//memcpy((void *)mxGetPr(POINTER1), (void *)pin, sizeof(double)*n*n);
				//engPutVariable(ep, "pin", POINTER1);


				//engEvalString(ep, "figure(10);imagesc(pin);caxis([-1 1]);colormap gray;axis square;");
				//jj++; if (jj > 4) { getchar(); if (jj > 4) { jj = 0; } }
			}
			if (flag == 0 && M[s + 1] < 0) { 
                criticalT[u] = s + 1;
                 flag = 1;
                 std::cout << "Flipped!\n";
                 //break; 
                 }																			//measurement of time to flip
		}
		// if (flag == 0) { u = u - 1; }																													//keep only runs that de-excitate

		count = 0; flag = 0;

		for (i = 0; i<N + 1; i++)																															//Sweep counter
		{
			M[i] = M[i] / (n*n);
		}

		for (i = 0; i<N; i++)																															//removal of the initial measure,ents
		{
			H2[i] = H[i + 1];
			M2[i] = M[i + 1];
		}


		meanM[u] = abs(mean(M2, N));																													//Mean magnetization after equilibrium

		heat_cpty[u] = pow(beta, 2)*pow(std_dev(H2, N), 2) / (n*n);

		info[0] = n * 1.0; info[1] = T; info[2] = 1.0*N; info[3] = 1.0*(u + 1); info[4] = meanM[u]; info[5] = heat_cpty[u]; info[6] = B[u];

		cout << "\nPlotting... \n\n";

		//mxArray* POINTER1 = mxCreateDoubleMatrix(n, n, mxREAL);
		//memcpy((void *)mxGetPr(POINTER1), (void *)pin, sizeof(double)*n*n);
		//engPutVariable(ep, "pin", POINTER1);

		//mxArray* POINTER2 = mxCreateDoubleMatrix(1, N, mxREAL);
		//memcpy((void *)mxGetPr(POINTER2), (void *)M2, sizeof(double)*N);
		//engPutVariable(ep, "Mag", POINTER2);

		//mxArray* POINTER3 = mxCreateDoubleMatrix(1, N, mxREAL);
		//memcpy((void *)mxGetPr(POINTER3), (void *)tsweep, sizeof(double)*N);
		//engPutVariable(ep, "Sweeps", POINTER3);

		//mxArray* POINTER4 = mxCreateDoubleMatrix(1, 7, mxREAL);
		//memcpy((void *)mxGetPr(POINTER4), (void *)info, sizeof(double)*(7));
		//engPutVariable(ep, "Info", POINTER4);

		//mxArray* POINTER5 = mxCreateDoubleMatrix(1, steps, mxREAL);
		//memcpy((void *)mxGetPr(POINTER5), (void *)Ni, sizeof(double)*(steps));
		//engPutVariable(ep, "Ni", POINTER5);

		//mxArray* POINTER6 = mxCreateDoubleMatrix(1, N, mxREAL);
		//memcpy((void *)mxGetPr(POINTER6), (void *)H2, sizeof(double)*N);
		//engPutVariable(ep, "Ham", POINTER6);

		//mxArray* POINTER7 = mxCreateDoubleMatrix(n, n, mxREAL);
		//memcpy((void *)mxGetPr(POINTER7), (void *)h, sizeof(double)*n*n);
		//engPutVariable(ep, "h", POINTER7);

		//
		//engEvalString(ep, "figure(1);subplot(1,2,1);imagesc(pin);caxis([-1 1]);colormap gray;axis square;");
		//engEvalString(ep, "title(['Lattice ',num2str(Info(1)),'x',num2str(Info(1)),', T=',num2str(Info(2)),' B=',num2str(Info(7))],'interpreter','latex');");
		//engEvalString(ep, "xlabel('i');ylabel('j');");
		//engEvalString(ep, "figure(1);subplot(1,2,2);imagesc(h);caxis([-1 1]);axis square;colorbar;");
		//engEvalString(ep, "title(['Applied field',', T=',num2str(Info(2))],'interpreter','latex');");
		//engEvalString(ep, "xlabel('i');ylabel('j');");

		//engEvalString(ep, "figure(2);subplot(3,5,Info(4));plot(Sweeps,Mag);ylim([-1 1]);xlabel('Measurements');ylabel('M(t)');"); //
		//engEvalString(ep, "title(['T=',num2str(Info(2)),' B=',num2str(Info(7))],'interpreter','latex');");

		//engEvalString(ep, "figure(3);subplot(3,5,Info(4));histogram(Mag,20);xlim([-1 1]);xlabel('M');ylabel('Times')");
		//engEvalString(ep, "title(['Sweeps: ',num2str(Info(3)),', T=',num2str(Info(2)),' B=',num2str(Info(7))],'interpreter','latex');");

		//engEvalString(ep, "figure(4);hold on;plot(Info(7),Info(5),'ok');xlabel('B');ylabel('|M|');ylim([0 1])");
		//engEvalString(ep, "title(['Order parameter, Lattice ',num2str(Info(1)),'x',num2str(Info(1))],'interpreter','latex');");

		//engEvalString(ep, "figure(5);subplot(3,5,Info(4));plot(Sweeps,Ham,'-r');xlabel('Measurements');ylabel('H(t)');");
		//engEvalString(ep, "title(['T=',num2str(Info(2)),' B=',num2str(Info(7))],'interpreter','latex');");

		//engEvalString(ep, "figure(6);hold on;plot(Info(7),Info(6),'or');xlabel('B');ylabel('|C|');");
		//engEvalString(ep, "title(['Heat capacity per lattice site'],'interpreter','latex');");

		
	}
	// mxArray* POINTER8 = mxCreateDoubleMatrix(1, steps, mxREAL);
	// memcpy((void *)mxGetPr(POINTER8), (void *)criticalT, sizeof(double)*(steps));
	// engPutVariable(ep, "CriticalT", POINTER8);

	// engEvalString(ep, "figure(7);histogram(CriticalT);");
	// engEvalString(ep, "xlabel('Time to flip (sweeps)');ylabel('Population');title(['Distribution of time to flip for L=',num2str(Info(1)),' '],'interpreter','latex')");


	count = -1;
	cout << "\nDone! \n" << endl;


	cout << "\n--------------------------------\n";
	cout << "Process exited after " << double(clock() - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;


	int n1 = sizeof(M) / sizeof(M[0]);
    std::vector<double> V1(M,M+n1);

	plot_fig(V1, 2, "Magnetization timeseries");




    // Write results in file
    std::cout << "Writing file" << std::endl;
    std::ofstream file("./data/magnetization.csv");

    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Header
    file << "Magnetization\n";

    for (const double& mag : M) {
        file << mag << "\n";
    }

    file.close();
    std::cout << "Data written to magnetization.csv" << std::endl;

	plt::show();

	// system("pause");
	return 0;
}