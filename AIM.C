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
	BucketParticle_ndx.resize(0);
	Lambda.resize(0);
	Grid_den.resize(0);
	Grid_phi.resize(0);
	Particle_phi.resize(0);
	
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

	Particle_X[0] = 0.05; Particle_X[1] = 0.5; Particle_X[2] = 0.95; Particle_X[3] = 0.05; Particle_X[4] = 0.5; Particle_X[5] = 0.95;
	Particle_X[6] = 0.05; Particle_X[7] = 0.5; Particle_X[8] = 0.95;
	Particle_Y[0] = 0.05; Particle_Y[1] = 0.05; Particle_Y[2] = 0.05; Particle_Y[3] = 0.5; Particle_Y[4] = 0.5; Particle_Y[5] = 0.5;
	Particle_Y[6] = 0.95; Particle_Y[7] = 0.95; Particle_Y[8] = 0.95;

	Particle_den[0] = 1.0; Particle_den[1] = 0.4; Particle_den[2] = 0.7; Particle_den[3] = 0.09; Particle_den[4] = 0.26; Particle_den[5] = 0.49;
	Particle_den[6] = 0.67; Particle_den[7] = 0.87; Particle_den[8] = 0.33;
//	char filename[50];
//	sprintf(filename, "Particle%d.dat", Np);
//
//	FILE *sample = fopen(filename, "r");
//
//	int ip;
//	for (ip = 0; ip < Np; ip++){
//		fscanf(sample, "%f ", &Particle_X[ip]);
//		fscanf(sample, "%f ", &Particle_Y[ip]);
//		fscanf(sample, "%f\n", &Particle_den[ip]);
//	}
//		fscanf(sample, "%f %f %f\n", &Particle_X[ip], &Particle_Y[ip], &Particle_den[ip]);

//	fclose(sample);
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
	std::vector<double> temp_Q(4, 1.0);
	int iBuck, iNode, temp;
	double dx = D / (N-1), x0, y0, tmp;

	// Loop over buckets
	for (int iy = 0; iy < N-1; iy++) {
		for (int ix = 0; ix < N-1; ix++) {

			iBuck = ix + (N-1) * iy;
			iNode = ix + N * iy;
			x0 = ix * dx;
			y0 = iy * dx;

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


	for (int iter = 11; iter < 101; iter++) {

		ParticleScattering test;
		test.N = iter;
		test.Np = 9;
		test.D = 1.0;

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
		double ER = test.ErrorEstimate();

		printf("-------------------------\n");
//		for (int i = 0; i < test.Np; i++)
//			printf("%g  %g\n", test.Phi_Dir[i], test.Particle_phi[i]);

		printf("Size = %d\tError = %g\n", iter, ER);


		test.DeAllocSize();

	}


	return 0;
}
