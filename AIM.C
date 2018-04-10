#include "AIM.h"

ParticleScattering::ParticleScattering()
{
	N = 1;
	Np = 1;
	D = 1.0;
	time_dir = 0.0;
	time_dft = 0.0;

	// Allocating complex vectors
	CC = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * 2 );
	XX  = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * 2 );
	YY  = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * 2 );
	Res_DFT = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * 2 );

	plan_XX = fftw_plan_dft_2d(2, 2, CC, XX, FFTW_FORWARD,FFTW_ESTIMATE);
	plan_YY = fftw_plan_dft_2d(2, 2, CC, YY, FFTW_FORWARD,FFTW_ESTIMATE);
	plan_INV = fftw_plan_dft_2d(2, 2, CC, Res_DFT, FFTW_BACKWARD,FFTW_ESTIMATE);
}

ParticleScattering::~ParticleScattering()
{
	return;
}

void ParticleScattering::AllocSize()
{
	T_direct.resize(N);
	T_all_direct.resize(N);
	Ci_FFT.resize(2*N);
	BucketParticle_ndx.resize((N-1)*(N-1));
	Lambda.resize(4*Np);
	Grid_den.resize(N*N);
	Grid_phi.resize(N*N);
	Particle_phi.resize(Np);


	Pattern_Ker.resize(9);


	// Allocating complex vectors
	CC = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * 4*N*N );
	XX  = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * 4*N*N );
	YY  = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * 4*N*N );
	Res_DFT = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * 4*N*N );
}
void ParticleScattering::DeAllocSize()
{
	T_direct.resize(0);
	T_all_direct.resize(0);
	Ci_FFT.resize(0);
	Ci_2_FFT.resize(0);
	BucketParticle_ndx.resize(0);
	Lambda.resize(0);
	Grid_den.resize(0);
	Grid_phi.resize(0);
	Particle_phi.resize(0);
	
	Pattern_Ker.resize(0);

	fftw_destroy_plan(plan_XX);
	fftw_destroy_plan(plan_YY);
	fftw_destroy_plan(plan_INV);

	fftw_free(CC);
	fftw_free(XX);
	fftw_free(YY);
	fftw_free(Res_DFT);

}

void ParticleScattering::ReadParticle()
{

	Particle_X.resize(Np);
	Particle_Y.resize(Np);
	Particle_den.resize(Np);

	int Nb = N-1;
	double dx = D / Nb;

//	Particle_X[0] = 0.0001; Particle_X[1] = 0.9999;
//	Particle_Y[0] = 0.0001; Particle_Y[1] = 0.9999;
//	Particle_X[0] = 0.05*dx; Particle_X[1] = D-0.05*dx; Particle_X[2] = 0.05*dx; Particle_X[3] = 0.95*dx; Particle_X[4] = 0.05*dx; Particle_X[5] = D-0.95*dx;
//	Particle_Y[0] = 0.05*dx; Particle_Y[1] = D-0.05*dx; Particle_Y[2] = D-0.95*dx; Particle_Y[3] = D-0.95*dx; Particle_Y[4] = 0.95*dx; Particle_Y[5] = D-0.95*dx;

//	Particle_X[0] = 0.05; Particle_X[1] = 0.5; Particle_X[2] = 0.95; Particle_X[3] = 0.05; Particle_X[4] = 0.5; Particle_X[5] = 0.95;
//	Particle_X[6] = 0.05; Particle_X[7] = 0.5; Particle_X[8] = 0.95;
//	Particle_Y[0] = 0.05; Particle_Y[1] = 0.05; Particle_Y[2] = 0.05; Particle_Y[3] = 0.5; Particle_Y[4] = 0.5; Particle_Y[5] = 0.5;
//	Particle_Y[6] = 0.95; Particle_Y[7] = 0.95; Particle_Y[8] = 0.95;
//
//	Particle_den[0] = 1.0; Particle_den[1] = 0.4; Particle_den[2] = 0.7; Particle_den[3] = 0.09; Particle_den[4] = 0.26; Particle_den[5] = 0.49;
//	Particle_den[6] = 0.67; Particle_den[7] = 0.87; Particle_den[8] = 0.33;

	char filename[50];
	sprintf(filename, "Particle%d.dat", Np);

	FILE *sample = fopen(filename, "r");

	int ip;
	for (ip = 0; ip < Np; ip++){
		fscanf(sample, "%f ", &Particle_X[ip]);
		fscanf(sample, "%f ", &Particle_Y[ip]);
		fscanf(sample, "%f\n", &Particle_den[ip]);
	}
		fscanf(sample, "%f %f %f\n", &Particle_X[ip], &Particle_Y[ip], &Particle_den[ip]);

	fclose(sample);
}

//Direct Method
void ParticleScattering::ComputeKernel()
{
	int Np = Particle_den.size();
	Kernel_Dir.resize(Np*Np);
	double temp1, temp2;


	for (int i = 0; i < Np; i++) {
		for (int j = 0; j < Np; j++) {
			if (j == i) {
				Kernel_Dir[j+i*Np] = 0.0;
			}
			else {
				temp1 = Particle_X[j] - Particle_X[i];
				temp2 = Particle_Y[j] - Particle_Y[i];

				Kernel_Dir[j+i*Np] = 1.0 / sqrt(temp1*temp1 + temp2*temp2);
			}
		}
	}
}
void ParticleScattering::ComputeDirect()
{
    double alpha = 1.0, beta = 0.0;
    char xx = 'n';
    const int m = N*N, n = N*N, lda = N*N, incx = N*N, incy = N*N;

    double *tmpT= new double[m*n];
    double *tmpX= new double[m];
    double *tmpY= new double[m];

    for (int i = 0; i < m*n; i++)
    	tmpT[i]=T_global_direct[i];

    for (int i = 0; i < m; i++)
    {
    	tmpX[i]=X_global_direct[i];
//    	tmpY[i]=Y_global_direct[i];
    }

    dgemv_(&xx, &m, &n, &alpha, tmpT, &lda, tmpX, &incx, &beta, tmpY, &incy);

    for (int i = 0; i < m; i++)
    	Y_global_direct.push_back(tmpY[i]);

    delete[] tmpT;
    delete[] tmpX;
    delete[] tmpY;
}
void ParticleScattering::ComputeDirect2()
{
	int Np = Particle_den.size();
	Phi_Dir.resize(Np);
	double temp;

	for (int i = 0; i < Np; i++) {
		temp = 0.0;
		for (int j = 0; j < Np; j++)
			temp += Particle_den[j] * Kernel_Dir[j+Np*i];

		Phi_Dir[i] = temp;
	}
}

// AIM Method
void ParticleScattering::FindBucket()
{
	int Nb = N-1;
	double dx = D / Nb;
	int temp;
	double xG, yG;



	for (int ig = 0; ig < Nb*Nb; ig++)
		BucketParticle_ndx[ig].resize(0);

	for (int ip = 0; ip < Np; ip++)
	{
		xG = Particle_X[ip]/dx;
		yG = Particle_Y[ip]/dx;

		if (xG < 0.0) {
			xG = 0.0;
		}
		if (yG < 0.0) {
			yG = 0.0;
		}
		if (xG >= Nb) {
			xG = Nb-1;
		}
		if (yG >= Nb) {
			yG = Nb-1;
		}

		temp = (int) (floor(xG) + Nb * floor(yG));
		BucketParticle_ndx[temp].push_back(ip);

		// Initialization for Phi on particles
		Particle_phi[ip] = 0.0;
	}
}
void ParticleScattering::Compute_W_inv()
{
	W_mu.resize(0);

	double dx = D / (N-1);

	// 1st Row
	W_mu.push_back(1.0);
	W_mu.push_back(-1.0/dx);
	W_mu.push_back(-1.0/dx);
	W_mu.push_back(1.0/(dx*dx));

	// 2nd Row
	W_mu.push_back(0.0);
	W_mu.push_back(1.0/dx);
	W_mu.push_back(0.0);
	W_mu.push_back(-1.0/(dx*dx));

	// 3rd Row
	W_mu.push_back(0.0);
	W_mu.push_back(0.0);
	W_mu.push_back(1.0/dx);
	W_mu.push_back(-1.0/(dx*dx));

	// 4th Row
	W_mu.push_back(0.0);
	W_mu.push_back(0.0);
	W_mu.push_back(0.0);
	W_mu.push_back(1.0/(dx*dx));

}
void ParticleScattering::BucketLambda()
{
	// Compute W inverse
	Compute_W_inv();

	std::vector<double> temp_Den;
	std::vector<double> temp_Q(4, 1.0);
	int iBuck, iNode, temp;
	double dx = D / (N-1), x0, y0, tmp;

	// Initialization for density on grid points
	for (int i = 0; i < N*N; i++)
		Grid_den[i] = 0.0;

	for (int iy = 0; iy < N-1; iy++) {
		for (int ix = 0; ix < N-1; ix++) {

			iBuck = ix + (N-1) * iy;
			iNode = ix + N * iy;
			x0 = ix * dx;
			y0 = iy * dx;

			temp_Den.resize(0);
			for (int itmp = 0; itmp < 4; itmp++)
				temp_Den.push_back(0.0);

			// Loop over particles assigned to a bucket
			for (int ip = 0; ip < BucketParticle_ndx[iBuck].size(); ip++) {

				temp = BucketParticle_ndx[iBuck][ip];

				temp_Q[1] = Particle_X[temp] - x0;
				temp_Q[2] = Particle_Y[temp] - y0;
				temp_Q[3] = temp_Q[1] * temp_Q[2];

				double tmp;

				for (int n = 0; n < 4; n++) {

					// Compute particle lambda
					tmp = 0.0;
					for (int m = 0; m < 4; m++)
						tmp += temp_Q[m] * W_mu[m + 4 * n];

					// Compute density on the grid points locally
					temp_Den[n] += Particle_den[temp] * tmp;

					// Global Lambda
					Lambda[4 * temp + n] = tmp;
				}
			}

			// Density on the grid points globally
			Grid_den[iNode] += temp_Den[0];
			Grid_den[iNode+1] += temp_Den[1];
			Grid_den[iNode+N] += temp_Den[2];
			Grid_den[iNode+N+1] += temp_Den[3];
		}
	}
}

void ParticleScattering::ComputeToeplitzUniques()
{
	double dx = D / (N-1);
	double temp;

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){

			temp = dx * sqrt((double) (1.0*(i*i + j*j)));

			if (temp != 0.0){
				T_direct[i].push_back(1.0/temp);
			}
			else {
				T_direct[i].push_back(0.0);
			}
		}
	}
}

void ParticleScattering::ComputeCirculant4fft()
{
	for (int i = 0; i < N; i++){
		Ci_FFT[i].insert(Ci_FFT[i].end(), T_direct[i].begin(), T_direct[i].end());
		Ci_FFT[i].push_back(0.0);
		for (int j = N-1; j > 0; j--)
			Ci_FFT[i].push_back(T_direct[i][j]);
	}

	for (int i = 0; i < 2*N; i++)
		Ci_FFT[N].push_back(0.0);

	for (int i = N-1; i > 0; i--)
		Ci_FFT[2*N-i].insert(Ci_FFT[2*N-i].end(), Ci_FFT[i].begin(), Ci_FFT[i].end());

	Ci_2_FFT.resize(0);

	for (int i = 0; i < 2*N; i++)
		for (int j = 0; j < 2*N; j++)
			Ci_2_FFT.push_back(Ci_FFT[j][i]);
}
void ParticleScattering::Gen_zeropadded_X4fft()
{
	std::vector<double> temp_M;
	std::vector<double> temp_2M;
	std::vector<double> temp_2M_zp1;
	std::vector<double> temp_2M_zp2;

	for (int i = 0; i < N; i++)
		temp_M.push_back(0.0);

	temp_2M.insert(temp_2M.end(), temp_M.begin(), temp_M.end());
	temp_2M.insert(temp_2M.end(), temp_M.begin(), temp_M.end());

	for (int i = 0; i < N; i++)
	{
		temp_2M_zp1.insert(temp_2M_zp1.end(), Grid_den.begin() + (i*N), Grid_den.begin() + ((i+1)*N));
		temp_2M_zp1.insert(temp_2M_zp1.end(), temp_M.begin(), temp_M.end());

		temp_2M_zp2.insert(temp_2M_zp2.end(), temp_2M.begin(), temp_2M.end());
	}


	temp_2M_zp1.insert(temp_2M_zp1.end(), temp_2M_zp2.begin(), temp_2M_zp2.end());
	X_2_FFT.insert(X_2_FFT.end(), temp_2M_zp1.begin(), temp_2M_zp1.end());
}
void ParticleScattering::Take_DFT_IDFT()
{
	for (unsigned long i = 0; i < 4*N*N; i++)
		CC[i] = Ci_2_FFT[i];

	plan_XX = fftw_plan_dft_2d(2*N, 2*N, CC, XX, FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(plan_XX);

	for (unsigned long i = 0; i < 4*N*N; i++)
		CC[i] = X_2_FFT[i];

	plan_YY = fftw_plan_dft_2d(2*N, 2*N, CC, YY, FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(plan_YY);

	// Multiplying the numbers
	for (unsigned long i = 0; i < 4*N*N; i++)
		CC[i] = XX[i] * YY[i];

	plan_INV = fftw_plan_dft_2d(2*N, 2*N, CC, Res_DFT, FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(plan_INV);
}

void ParticleScattering::MapBack2Particles()
{
	int ns = N*N;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			Grid_phi[j+i*N] = creal(Res_DFT[j+i*2*N]) / (4*ns);

	std::vector<double> temp_Phi;
	int iBuck, iNode, temp;

	// Loop over buckets
	for (int iy = 0; iy < N-1; iy++) {
		for (int ix = 0; ix < N-1; ix++) {

			iBuck = ix + (N-1) * iy;
			iNode = ix + N * iy;

			temp_Phi.resize(0);
			for (int itmp = 0; itmp < 4; itmp++)
				temp_Phi.push_back(0.0);

			// Loop over particles assigned to a bucket
			for (int ip = 0; ip < BucketParticle_ndx[iBuck].size(); ip++) {

				temp = BucketParticle_ndx[iBuck][ip];

				Particle_phi[temp] += Lambda[4*temp] * Grid_phi[iNode];
				Particle_phi[temp] += Lambda[4*temp+1] * Grid_phi[iNode+1];
				Particle_phi[temp] += Lambda[4*temp+2] * Grid_phi[iNode+N];
				Particle_phi[temp] += Lambda[4*temp+3] * Grid_phi[iNode+N+1];
			}
		}
	}
}

void ParticleScattering::CreateTwo_BucketKerPatterns()
{
	double X_pattern[9][4] = {
		{ 0.0, 1.0, 0.0, 1.0 }, { 1.0, 2.0, 1.0, 2.0 }, { 1.0, 2.0, 1.0, 2.0 },
		{ 0.0, 1.0, 0.0, 1.0 }, { -1.0, 0.0, -1.0, 0.0 }, { -1.0, 0.0, -1.0, 0.0 },
		{ -1.0, 0.0, -1.0, 0.0 }, { 0.0, 1.0, 0.0, 1.0 }, { 1.0, 2.0, 1.0, 2.0 } };

	double Y_pattern[9][4] = {
		{ 0.0, 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0, 1.0 }, { 1.0, 1.0, 2.0, 2.0 },
		{ 1.0, 1.0, 2.0, 2.0 }, { 1.0, 1.0, 2.0, 2.0 }, { 0.0, 0.0, 1.0, 1.0 },
		{ -1.0, -1.0, 0.0, 0.0 }, { -1.0, -1.0, 0.0, 0.0 }, { -1.0, -1.0, 0.0, 0.0 } };

	double dx = D / (N-1);
	double temp1, temp2;

	for (int iPat = 0; iPat < 9; iPat++){
		for (int j = 0; j < 4; j++) {
			for (int i = 0; i < 4; i++) {

				temp1 = X_pattern[0][i] - X_pattern[iPat][j];
				temp2 = Y_pattern[0][i] - Y_pattern[iPat][j];

				if (temp1 == 0.0 && temp2 == 0.0) {
					Pattern_Ker[iPat].push_back(0.0);
				}
				else {
					Pattern_Ker[iPat].push_back(1.0 / (dx * sqrt(temp1*temp1 + temp2*temp2)));
				}

			}
		}
	}
}

void ParticleScattering::Assign_NeighKernel(int ID2)
{
	for (int i = 0; i < 16; i++)
		Grid_Ker.push_back(Pattern_Ker[ID2][i]);
}

void ParticleScattering::Two_BucketKernel(int ID1) {

	if (ID1 == 0)
	Assign_NeighKernel(0);
	else if (ID1 == 1)
	Assign_NeighKernel(1);
	else if (ID1 == N+1)
	Assign_NeighKernel(2);
	else if (ID1 == N)
	Assign_NeighKernel(3);
	else if (ID1 == N-1)
	Assign_NeighKernel(4);
	else if (ID1 == -1)
	Assign_NeighKernel(5);
	else if (ID1 == -N-1)
	Assign_NeighKernel(6);
	else if (ID1 == -N)
	Assign_NeighKernel(7);
	else if (ID1 == -N+1)
	Assign_NeighKernel(8);

	// switch (ID1) {
	// case 0:
	// Assign_NeighKernel(0);
	// break;
	// case 1:
	// Assign_NeighKernel(1);
	// break;
	// case N+1:
	// Assign_NeighKernel(2);
	// break;
	// case N:
	// Assign_NeighKernel(3);
	// break;
	// case N-1:
	// Assign_NeighKernel(4);
	// break;
	// case -1:
	// Assign_NeighKernel(5);
	// break;
	// case -N-1:
	// Assign_NeighKernel(6);
	// break;
	// case -N:
	// Assign_NeighKernel(7);
	// break;
	// case -N+1:
	// Assign_NeighKernel(8);
	// break;
	// default:
	// printf("Error");
	// break;
	// }
}

void ParticleScattering::Bucket2Bucket()
{
	int tempC, tempN, nCurr, nNeigh;
	double temp1, temp2;

	// Loop over each Neighboring bucket
	for (int iB = 0; iB < iNeigh.size(); iB++)
	{
		nCurr = BucketParticle_ndx[iNeigh[0]].size();
		nNeigh = BucketParticle_ndx[iNeigh[iB]].size();
		Near_Ker.resize(nCurr * nNeigh);

		// Loop over particles assigned to the current bucket
		for (int ip = 0; ip < nCurr; ip++) {

			tempC = BucketParticle_ndx[iNeigh[0]][ip];

			// Loop over particles assigned to the neighboring bucket
			for (int ipN = 0; ipN < nNeigh; ipN++) {

				tempN = BucketParticle_ndx[iNeigh[iB]][ipN];

				// Computing Direct Kernel for Near Zone
				if (tempC == tempN) {
					Near_Ker[ipN+ip*nNeigh] = 0.0;
				}
				else {
					temp1 = Particle_X[tempC] - Particle_X[tempN];
					temp2 = Particle_Y[tempC] - Particle_Y[tempN];

					Near_Ker[ipN+ip*nNeigh] = 1.0 / sqrt(temp1*temp1 + temp2*temp2);
				}
			}
		}

		Two_BucketKernel(iNeighN[iB]-iNeighN[0]);
		std::vector<double> G_L(4*nNeigh,0.0);

		for (int i = 0; i < 4; i++) {

			// Loop over particles assigned to the neighboring bucket
			for (int ipN = 0; ipN < nNeigh; ipN++) {
				tempN = BucketParticle_ndx[iNeigh[iB]][ipN];

				for (int n = 0; n < 4; n++)
					G_L[nNeigh * i + ipN] += Grid_Ker[4*i + n] * Lambda[n + 4*tempN];
			}
		}

//		std::vector<double> LT_G_L(nCurr * nNeigh,0.0);

		for (int ipC = 0; ipC < nCurr; ipC++) {
			tempC = BucketParticle_ndx[iNeigh[0]][ipC];

			// Loop over particles assigned to the neighboring bucket
			for (int ipN = 0; ipN < nNeigh; ipN++) {

				for (int n = 0; n < 4; n++)
					Near_Ker[nNeigh * ipC + ipN] -= G_L[nCurr*ipC + n] * Lambda[n + 4*tempC];
			}
		}



		for (int ipC = 0; ipC < nCurr; ipC++) {

			tempC = BucketParticle_ndx[iNeigh[0]][ipC];

			for (int ipN = 0; ipN < nNeigh; ipN++) {

				tempN = BucketParticle_ndx[iNeigh[iB]][ipN];

				Particle_phi[tempC] += Near_Ker[nNeigh * ipC + ipN] * Particle_den[tempN];

			}
		}
	}
}

void ParticleScattering::NearZoneCompute()
{

	CreateTwo_BucketKerPatterns();

	int iBuck;
	double dx = D / (N-1);

	int jx, jy;
	// Corner Boundaries
	iNeigh.resize(0); iNeighN.resize(0);
	jx = 0; jy = 0;

	iNeigh.push_back(0); iNeighN.push_back(0);
	iNeigh.push_back(1); iNeighN.push_back(1);
	iNeigh.push_back(N-1); iNeighN.push_back(N);
	iNeigh.push_back(N); iNeighN.push_back(N+1);
	Bucket2Bucket();


	iNeigh.resize(0); iNeighN.resize(0);
	jx = N-2; jy = 0;

	iNeigh.push_back(jx + (N-1) * jy); iNeighN.push_back(iNeigh[0] + jy);
	iNeigh.push_back(jx-1 + (N-1) * jy); iNeighN.push_back(iNeigh[1] + jy);
	iNeigh.push_back(jx + (N-1) * (jy+1)); iNeighN.push_back(iNeigh[2] + jy+1);
	iNeigh.push_back(jx-1 + (N-1) * (jy+1)); iNeighN.push_back(iNeigh[3] + jy+1);
	Bucket2Bucket();


	iNeigh.resize(0); iNeighN.resize(0);
	jx = 0; jy = N-2;

	iNeigh.push_back(jx + (N-1) * jy);  iNeighN.push_back(iNeigh[0] + jy);
	iNeigh.push_back(jx+1 + (N-1) * jy); iNeighN.push_back(iNeigh[1] + jy);
	iNeigh.push_back(jx+1 + (N-1) * (jy-1)); iNeighN.push_back(iNeigh[2] + jy-1);
	iNeigh.push_back(jx + (N-1) * (jy-1)); iNeighN.push_back(iNeigh[3] + jy-1);
	Bucket2Bucket();


	iNeigh.resize(0); iNeighN.resize(0);
	jx = N-2; jy = N-2;

	iNeigh.push_back(jx + (N-1) * jy); iNeighN.push_back(iNeigh[0] + jy);
	iNeigh.push_back(jx-1 + (N-1) * jy); iNeighN.push_back(iNeigh[1] + jy);
	iNeigh.push_back(jx-1 + (N-1) * (jy-1)); iNeighN.push_back(iNeigh[2] + jy-1);
	iNeigh.push_back(jx + (N-1) * (jy-1)); iNeighN.push_back(iNeigh[3] + jy-1);
	Bucket2Bucket();


	int j;
	// Line Boundaries
	for (int i = 1; i < N-2; i++)
	{
		j = 0;

		iNeigh.resize(0); iNeighN.resize(0);

		iNeigh.push_back(i + (N-1) * j); iNeighN.push_back(iNeigh[0] + j);
		iNeigh.push_back(i+1 + (N-1) * j); iNeighN.push_back(iNeigh[1] + j);
		iNeigh.push_back(i-1 + (N-1) * j); iNeighN.push_back(iNeigh[2] + j);
		iNeigh.push_back(i + (N-1) * (j+1)); iNeighN.push_back(iNeigh[3] + j+1);
		iNeigh.push_back(i+1 + (N-1) * (j+1)); iNeighN.push_back(iNeigh[4] + j+1);
		iNeigh.push_back(i-1 + (N-1) * (j+1)); iNeighN.push_back(iNeigh[5] + j+1);
		Bucket2Bucket();


		iNeigh.resize(0); iNeighN.resize(0);

		iNeigh.push_back(j + (N-1) * i); iNeighN.push_back(iNeigh[0] + i);
		iNeigh.push_back(j+1 + (N-1) * i); iNeighN.push_back(iNeigh[1] + i);
		iNeigh.push_back(j + (N-1) * (i+1)); iNeighN.push_back(iNeigh[2] + i+1);
		iNeigh.push_back(j + (N-1) * (i-1)); iNeighN.push_back(iNeigh[3] + i-1);
		iNeigh.push_back(j+1 + (N-1) * (i+1)); iNeighN.push_back(iNeigh[4] + i+1);
		iNeigh.push_back(j+1 + (N-1) * (i-1)); iNeighN.push_back(iNeigh[5] + i-1);
		Bucket2Bucket();


		j = N-2;

		iNeigh.resize(0); iNeighN.resize(0);

		iNeigh.push_back(i + (N-1) * j); iNeighN.push_back(iNeigh[0] + j);
		iNeigh.push_back(i+1 + (N-1) * j); iNeighN.push_back(iNeigh[1] + j);
		iNeigh.push_back(i-1 + (N-1) * j); iNeighN.push_back(iNeigh[2] + j);
		iNeigh.push_back(i + (N-1) * (j-1)); iNeighN.push_back(iNeigh[3] + j-1);
		iNeigh.push_back(i+1 + (N-1) * (j-1)); iNeighN.push_back(iNeigh[4] + j-1);
		iNeigh.push_back(i-1 + (N-1) * (j-1)); iNeighN.push_back(iNeigh[5] + j-1);
		Bucket2Bucket();


		iNeigh.resize(0); iNeighN.resize(0);

		iNeigh.push_back(j + (N-1) * i); iNeighN.push_back(iNeigh[0] + i);
		iNeigh.push_back(j-1 + (N-1) * i);  iNeighN.push_back(iNeigh[1] + i);
		iNeigh.push_back(j + (N-1) * (i+1)); iNeighN.push_back(iNeigh[2] + i+1);
		iNeigh.push_back(j + (N-1) * (i-1)); iNeighN.push_back(iNeigh[3] + i-1);
		iNeigh.push_back(j-1 + (N-1) * (i+1)); iNeighN.push_back(iNeigh[4] + i+1);
		iNeigh.push_back(j-1 + (N-1) * (i-1)); iNeighN.push_back(iNeigh[5] + i-1);
		Bucket2Bucket();

	}


	// Loop over buckets
	for (int iy = 2; iy < N-2; iy++) {
		for (int ix = 2; ix < N-2; ix++) {

			iNeigh.resize(0); iNeighN.resize(0);

			iNeigh.push_back(ix + (N-1) * iy); iNeighN.push_back(iNeigh[0] + iy);
			iNeigh.push_back(ix+1 + (N-1) * iy); iNeighN.push_back(iNeigh[1] + iy);
			iNeigh.push_back(ix-1 + (N-1) * iy); iNeighN.push_back(iNeigh[2] + iy);
			iNeigh.push_back(ix + (N-1) * (iy+1)); iNeighN.push_back(iNeigh[3] + iy+1);
			iNeigh.push_back(ix + (N-1) * (iy-1)); iNeighN.push_back(iNeigh[4] + iy-1);
			iNeigh.push_back(ix+1 + (N-1) * (iy+1)); iNeighN.push_back(iNeigh[5] + iy+1);
			iNeigh.push_back(ix+1 + (N-1) * (iy-1)); iNeighN.push_back(iNeigh[6] + iy-1);
			iNeigh.push_back(ix-1 + (N-1) * (iy+1)); iNeighN.push_back(iNeigh[7] + iy+1);
			iNeigh.push_back(ix-1 + (N-1) * (iy-1)); iNeighN.push_back(iNeigh[8] + iy-1);
			Bucket2Bucket();

		}
	}


//	int i;
//	for (int j = 1; j < N-2; j++)
//	{
//		i = 0;
//		iNeigh.resize(0);
//		iBuck = i + (N-1) * j;
//
//		iNeigh.push_back(i+1 + (N-1) * j);
//		iNeigh.push_back(i + (N-1) * (j+1));
//		iNeigh.push_back(i + (N-1) * (j-1));
//		iNeigh.push_back(i+1 + (N-1) * (j+1));
//		iNeigh.push_back(i+1 + (N-1) * (j-1));
//
//
//		i = N-2;
//		iNeigh.resize(0);
//		iBuck = i + (N-1) * j;
//
//		iNeigh.push_back(i-1 + (N-1) * j);
//		iNeigh.push_back(i + (N-1) * (j+1));
//		iNeigh.push_back(i + (N-1) * (j-1));
//		iNeigh.push_back(i-1 + (N-1) * (j+1));
//		iNeigh.push_back(i-1 + (N-1) * (j-1));
//
//
//	}
}

double ParticleScattering::ErrorEstimate()
{
	double err1 = 0.0, err2 = 0.0, temp1;

	for (int i = 0; i < Np; i++) {
			temp1 = Phi_Dir[i] - Particle_phi[i];
			err1 += temp1 * temp1;
			err2 += Phi_Dir[i] * Phi_Dir[i];
	}

	err1 = sqrt(err1/err2);

	return err1;
}

int main() {


//	for (int iter = 7; iter < 31; iter++) {

		int iter = 33;
		ParticleScattering test;
		test.N = iter;
		test.Np = 10000;
		test.D = 100.0;

		test.AllocSize();
		test.ReadParticle();
		test.ComputeKernel();
		test.ComputeDirect2();

		test.FindBucket();
		test.BucketLambda();

		test.ComputeToeplitzUniques();
		test.ComputeCirculant4fft();
		test.Gen_zeropadded_X4fft();
		test.Take_DFT_IDFT();

		test.MapBack2Particles();

		test.NearZoneCompute();

		double ER = test.ErrorEstimate();

		printf("--------------------------------------\n");
		for (int i = 0; i < test.Np; i++)
			printf("%g  %g\n", test.Phi_Dir[i], test.Particle_phi[i]);

		printf("Size = %d\tError = %g\n", iter, ER);


		test.DeAllocSize();

//	}


	return 0;
}
