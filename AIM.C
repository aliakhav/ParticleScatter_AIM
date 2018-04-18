#include "AIM.h"

ParticleScattering::ParticleScattering()
{
	N = 1;
	Np = 1;
	D = 1.0;

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

				temp1 = X_pattern[0][j] - X_pattern[iPat][i];
				temp2 = Y_pattern[0][j] - Y_pattern[iPat][i];

				if (temp1 == 0.0 && temp2 == 0.0) {
					Pattern_Ker[iPat].push_back(0.0);
				}
				else {
					Pattern_Ker[iPat].push_back(1.0 / (dx * sqrt(temp1*temp1 + temp2*temp2)));
//					Pattern_Ker[iPat].push_back(1.0 / sqrt(temp1*temp1 + temp2*temp2));
				}
			}
		}
	}
}

void ParticleScattering::Assign_NeighKernel(int ID2)
{
	for (int i = 0; i < 16; i++)
		Grid_Ker[i] = Pattern_Ker[ID2][i];
}

void ParticleScattering::Two_BucketKernel(int ID1)
{
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
}

void ParticleScattering::Bucket2Bucket()
{
	int tempC, tempN, nCurr, nNeigh;
	double temp1, temp2;

	nCurr = BucketParticle_ndx[iNeigh[0]].size();

	if (nCurr > 0)
	{
		// Loop over each Neighboring bucket
		for (int iB = 0; iB < iNeigh.size(); iB++)
		{
			nNeigh = BucketParticle_ndx[iNeigh[iB]].size();

			if (nNeigh > 0)
			{
				// Computing Direct Kernel for Near Zone between Current Bucket and Neighboring Bucket
				Near_Ker.resize(nCurr * nNeigh);

				// Loop over particles assigned to the current bucket
				for (int ip = 0; ip < nCurr; ip++)
				{
					tempC = BucketParticle_ndx[iNeigh[0]][ip];

					// Loop over particles assigned to the neighboring bucket
					for (int ipN = 0; ipN < nNeigh; ipN++)
					{
						tempN = BucketParticle_ndx[iNeigh[iB]][ipN];

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

				// Computing Kernel on grid (4*4 matrix) for Near Zone between Current Bucket and Neighboring Bucket
				Two_BucketKernel(iNeighN[iB]-iNeighN[0]);

				// Computing multiplication of Kernel on grid (4*4 matrix) and Lambda matrix for Neighboring
				std::vector<double> G_L(4*nNeigh,0.0);

				// Loop over particles assigned to the neighboring bucket
				for (int ipN = 0; ipN < nNeigh; ipN++)
				{
					tempN = BucketParticle_ndx[iNeigh[iB]][ipN];

					for (int i = 0; i < 4; i++)
						for (int n = 0; n < 4; n++)
							G_L[4 * ipN + i] += Grid_Ker[4*i + n] * Lambda[n + 4*tempN];
				}


				// Computing part of Far zone kernel associated with Current Bucket and the Neighboring Bucket
				// multiplication of Lambda transpose matrix to the pre-computed "Kernel on grid * Neighbor_Lambda"
				// Computing Corrected part of the Global Far Kernel associated with the near zone between current and neighboring buckets
				for (int ipC = 0; ipC < nCurr; ipC++)
				{
					tempC = BucketParticle_ndx[iNeigh[0]][ipC];

					// Loop over particles assigned to the neighboring bucket
					for (int ipN = 0; ipN < nNeigh; ipN++)
						for (int n = 0; n < 4; n++)
							Near_Ker[nNeigh * ipC + ipN] -= Lambda[4*tempC + n] * G_L[n + 4*ipN];
				}


				// Fixing part of the Potential vector affected by incorrect computations due to including near zone interactions inside the global far kernel
				for (int ipC = 0; ipC < nCurr; ipC++)
				{
					tempC = BucketParticle_ndx[iNeigh[0]][ipC];

					for (int ipN = 0; ipN < nNeigh; ipN++)
					{
						tempN = BucketParticle_ndx[iNeigh[iB]][ipN];

						Particle_phi[tempC] += Near_Ker[nNeigh * ipC + ipN] * Particle_den[tempN];
					}
				}
			}
		}
	}
}

void ParticleScattering::NearZoneCompute()
{

	CreateTwo_BucketKerPatterns();
	Grid_Ker.resize(16);

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


	// Loop over inner buckets
	for (int iy = 1; iy < N-2; iy++) {
		for (int ix = 1; ix < N-2; ix++) {

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


	// Error and Accuracy tests

	int SZ1[4] = {500, 2000, 5000, 20000};
	int GPs1[10] = {4, 7, 11, 16, 21, 51, 81, 101, 201, 401};
	int Ns = 4, N_max =10 ;

	char filename[50];

	for (int siz = 0; siz < Ns; Ns++)
	{
		sprintf(filename, "ErrAcc_%d.dat", SZ1[siz]);

		for (int iter = 0; iter < N_max; iter++)
		{
			ParticleScattering ACC;
			ACC.N = GPs1[iter];
			ACC.Np = SZ1[siz];
			ACC.D = 100.0;

			ACC.AllocSize();
			ACC.ReadParticle();
			ACC.ComputeKernel();

			ACC.ComputeDirect();

			ACC.FindBucket();
			ACC.BucketLambda();

			ACC.ComputeToeplitzUniques();
			ACC.ComputeCirculant4fft();
			ACC.Gen_zeropadded_X4fft();

			ACC.Take_DFT_IDFT();
			ACC.MapBack2Particles();
			ACC.NearZoneCompute();

			double ER = ACC.ErrorEstimate();

			FILE* fp = fopen(filename, "a");
			fprintf(fp, "%5d, %3.9g, %3.16g\n", (GPs1[iter]-1)*(GPs1[iter]-1), 1.0/(GPs1[iter]-1), ER);
			fclose(fp);

			ACC.DeAllocSize();

		}
	}

	//////////////////////
	// Complexity tests //
	/////////////////////

	int N_sz =16;

	int SZ2[16] = { 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000,
			200000, 500000, 1000000, 2000000, 5000000, 10000000 };
	int GPs2[16] = { 4, 6, 8, 11, 15, 24, 33, 46, 72, 101, 143, 235, 317, 448,
			708, 1001 };

	timespec before, after;
	timespec time_diff;

	double Time_dir, Time_aim_tot, Time_aim_mv;

	// Timing Direct Method
	sprintf(filename, "Direct.dat");

	for (int iter = 0; iter < 8; iter++)
	{
		ParticleScattering DIR;
		DIR.N = GPs2[iter];
		DIR.Np = SZ2[iter];
		DIR.D = 100.0;

		DIR.AllocSize();
		DIR.ReadParticle();
		DIR.ComputeKernel();

		get_time(&before);
		DIR.ComputeDirect();
		get_time(&after);
		diff(&before,&after,&time_diff);

		// Time in seconds
		Time_dir = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;

		get_time(&before);
		DIR.FindBucket();
		DIR.BucketLambda();

		DIR.ComputeToeplitzUniques();
		DIR.ComputeCirculant4fft();
		DIR.Gen_zeropadded_X4fft();
		get_time(&after);
		diff(&before,&after,&time_diff);

		// Time in seconds
		Time_aim_tot = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;

		get_time(&before);
		DIR.Take_DFT_IDFT();
		DIR.MapBack2Particles();
		DIR.NearZoneCompute();
		get_time(&after);
		diff(&before,&after,&time_diff);

		// Time in seconds
		Time_aim_mv = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;

		Time_aim_tot += Time_aim_mv;


		FILE* fp = fopen(filename, "a");
		fprintf(fp, "%d, %3.16g, %3.9g\n", SZ2[iter], 1.0/(GPs2[iter]-1), Time_dir);
		fclose(fp);


		DIR.DeAllocSize();
	}


	// Timing AIM
	sprintf(filename, "Complexity.dat");

	for (int iter = 0; iter < N_sz; iter++)
	{
		ParticleScattering CPX;
		CPX.N = GPs2[iter];
		CPX.Np = SZ2[iter];
		CPX.D = 100.0;

		CPX.AllocSize();
		CPX.ReadParticle();

		get_time(&before);
		CPX.FindBucket();
		CPX.BucketLambda();

		CPX.ComputeToeplitzUniques();
		CPX.ComputeCirculant4fft();
		CPX.Gen_zeropadded_X4fft();
		get_time(&after);
		diff(&before,&after,&time_diff);

		// Time in seconds
		Time_aim_tot = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;

		get_time(&before);
		CPX.Take_DFT_IDFT();
		CPX.MapBack2Particles();
		CPX.NearZoneCompute();
		get_time(&after);
		diff(&before,&after,&time_diff);

		// Time in seconds
		Time_aim_mv = time_diff.tv_sec + (double)(time_diff.tv_nsec)/1.0e9;

		Time_aim_tot += Time_aim_mv;


		FILE* fp = fopen(filename, "a");
		fprintf(fp, "%d, %3.16g, %3.9g, %3.9g\n", SZ2[iter], 1.0/(GPs2[iter]-1), Time_aim_mv, Time_aim_tot);
		fclose(fp);


		CPX.DeAllocSize();
	}


	return 0;
}
