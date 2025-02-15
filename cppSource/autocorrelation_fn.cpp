#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include <my_func.h> //PI, random(max), mean(values,size), std_dev(values,size), put(DATA,Size1,Size2=1)										//delete one of the two
//using namespace std;

int main()
{
	clock_t startTime = clock();
	int i = 0, j = 0, s = 0, k = 0, l = 0, m = 0, u = 0;
	const int n = 128, N = 30000, steps = 3; //max 15 for figure
	double pin[n][n], h[n][n], neigh[4], pith[4], M[N + 1], H[N + 1], X[N + 1], tsweep[N + 1], info[6], B[steps], meanM[steps], heat_cpty[steps];
	double S = 0, S1 = 0, S2 = 0, DE = 0, r = 0, J = 1.0, step_width = 0.002, T = 1.0, beta = 1.0 / T;
	B[0] = 0.0;																																		//initial temperature
	for (u = 0; u<steps; u++)																															//temperatures
	{
		B[u + 1] = B[u] + step_width;
	}


	// cout << "Initializing MATLAB engine... \n";

	srand(time(NULL));
	// Engine *ep;
	// ep = engOpen("NULL");																															//OPEN ENGINE
	// engEvalString(ep, "clear all;");
	// //engEvalString(ep, "close all;");	


	for (i = 0; i<4; i++)																																//Arxikopoihsh athroismatos geitonwn
	{
		neigh[i] = 0;
	}

	for (i = 0; i<N + 1; i++)																																//Sweep counter
	{
		tsweep[i] = i;
	}

	// engEvalString(ep, "figure(1);clf;figure(2);clf;figure(3);clf;");

	std::cout << "Calculating...\n\n";


	for (u = 0; u<steps; u++)
	{
		for (i = 0; i<n; i++)				                                                       															//Magnetic field
		{
			for (j = 0; j<n; j++)
			{
				h[i][j] = B[u];
				if (i < n / 2) { h[i][j] = -1.0*h[i][j]; }
			}
		}

		for (i = 0; i<n; i++)                                                           																	//Times plegmatos
		{

			for (j = 0; j<n; j++)
			{
				m = floor(random(2));  //T=inf
				if (m == 2) m = 1;
				pin[i][j] = 2 * m - 1;
			}
		}
		M[0] = 0;
		for (i = 0; i<n; i++)                                                           																	//Ypologismos arxikhs magnhtishs 
		{
			for (j = 0; j<n; j++)
			{
				M[0] = M[0] + pin[i][j];
			}
		}
		std::cout << "Step " << u + 1 << " of " << steps << " (B = " << B[u] << ")\t (" << N << " sweeps per step)\r\n";
		for (s = 0; s<(N); s++)																															//Sweep loop
		{
			M[s + 1] = M[s];																															//Arxikh magnhtish

			for (k = 0; k<n; k++)                                                           																//Times plegmatos
			{

				for (l = 0; l<n; l++)
				{
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
							//H[s+1]=H[s+1]+DE;
						}
					}

					else
					{
						pin[k][l] = pin[k][l] * (-1);
						M[s + 1] = M[s + 1] + 2 * pin[k][l];
						//H[s+1]=H[s+1]+DE;
					}
					S = 0;
				}
			}

		}


		for (i = 0; i<N + 1; i++)										//correlations calcolation
		{
			for (j = 0; j<N + 1 - i; j++)
			{
				S = S + M[j] * M[i + j];
				S2 = S2 + M[j];
				S1 = S1 + M[i + j];
			}

			X[i] = (1.0 / (N + 1.0 - i))*(S - S1 * S2*(1.0 / (N + 1.0 - i)));

			S1 = 0; S2 = 0; S = 0;
		}
		for (i = 0; i<N + 1; i++)
		{
			X[i] = X[i] / (X[0]);
		}

		for (i = 0; i<N + 1; i++)																															//magnetization normalization
		{
			M[i] = M[i] / (n*n);
		}

		meanM[u] = mean(M, N + 1);																														//Mean magnetization after equilibrium


		info[0] = n * 1.0; info[1] = B[u]; info[2] = 1.0*N; info[3] = 1.0*(u + 1); info[4] = meanM[u]; info[5] = T;

		std::cout << "\nPlotting... \n\n";

		// mxArray* POINTER1 = mxCreateDoubleMatrix(n, n, mxREAL);
		// memcpy((void *)mxGetPr(POINTER1), (void *)pin, sizeof(double)*n*n);
		// engPutVariable(ep, "pin", POINTER1);

		// mxArray* POINTER2 = mxCreateDoubleMatrix(1, N + 1, mxREAL);
		// memcpy((void *)mxGetPr(POINTER2), (void *)X, sizeof(double)*(N + 1));
		// engPutVariable(ep, "Cor", POINTER2);

		// mxArray* POINTER3 = mxCreateDoubleMatrix(1, N + 1, mxREAL);
		// memcpy((void *)mxGetPr(POINTER3), (void *)tsweep, sizeof(double)*(N + 1));
		// engPutVariable(ep, "Sweeps", POINTER3);

		// mxArray* POINTER4 = mxCreateDoubleMatrix(1, 6, mxREAL);
		// memcpy((void *)mxGetPr(POINTER4), (void *)info, sizeof(double)*(6));
		// engPutVariable(ep, "Info", POINTER4);

		// mxArray* POINTER6 = mxCreateDoubleMatrix(1, N + 1, mxREAL);
		// memcpy((void *)mxGetPr(POINTER6), (void *)M, sizeof(double)*(N + 1));
		// engPutVariable(ep, "Mag", POINTER6);


		// engEvalString(ep, "figure(1);clf;imagesc(pin);caxis([-1 1]);colormap gray;axis square;");
		// engEvalString(ep, "title(['Lattice ',num2str(Info(1)),'x',num2str(Info(1)),', B=',num2str(Info(2)),' T=',num2str(Info(6))],'interpreter','latex');");
		// engEvalString(ep, "xlabel('i');ylabel('j');");

		// engEvalString(ep, "figure(2);subplot(3,5,Info(4));plot(Sweeps,Cor);xlabel('Sweeps');ylabel('X(t)');");
		// engEvalString(ep, "title(['B=',num2str(Info(2)),' T=',num2str(Info(6))],'interpreter','latex');");

		// engEvalString(ep, "figure(3);subplot(3,5,Info(4));plot(Sweeps,Mag);xlabel('Sweeps');ylabel('M(t)');ylim([-1 1])");
		// engEvalString(ep, "title(['B=',num2str(Info(2)),' T=',num2str(Info(6))],'interpreter','latex');");

		//M[0]=M[N]*n*n;																															//begin with the last magnetization at the next step
		//H[0]=H[N];

	}

	std::cout << "\nDone! \n" << std::endl;


	std::cout << "\n--------------------------------\n";
	std::cout << "Process exited after " << double(clock() - startTime) / (double)CLOCKS_PER_SEC << " seconds." << std::endl;


    // Write results in file
    std::cout << "Writing file" << std::endl;
    std::ofstream file("magnetization.csv");

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


	// system("pause");
	return 0;
}