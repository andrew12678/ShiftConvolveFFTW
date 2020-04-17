#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
//#include <time.h>
/*
 * This is the C implemented 'streamlined' approach of convoluting a given probability vector in the Fourier domain
*/
void fftconvPairs(double *p, int *n, double complex *A, double complex *B, double complex *U,  double *result, double *in, double complex *w_vector, double complex *out){
    int H = 4; //Each column belongs to 1 bernoulli prob and starts off with 4 values (p, 1-p, 0, 0)
    int L = *n; //Number of columns (different bernoulli probabilities to start with)
    
    //Setting up fftw so we perform fft on all the columns
    double *hc_to_real; //Double array used for Half Complex -> Real 

    hc_to_real = malloc(sizeof(double) * H);
    fftw_plan my_plan = fftw_plan_dft_r2c_1d(H, in, out, FFTW_ESTIMATE);

    in[2] = 0; //We are always padding the last two entries of the columns with zeroes
    in[3] = 0;
    for(int i = 0; i < L; i++){
        in[0] = p[i];
        in[1] = p[i+L];
        fftw_execute(my_plan);
        A[i*H] = out[0]; //Storing the Fourier domain transformations (complex numbers) in complex array A
        A[i*H+1] = out[1];
        A[i*H+2] = out[2];
        A[i*H+3] = conj(out[1]); //fftw_plan_dft_r2c_1d is not the same as fftw_plan_dft_1d and we need to compute the conjugate for the remaining
        
    }
    
    fftw_destroy_plan(my_plan); //Plans must be destroyed otherwise memory leak, confirmed via Valgrind
    
    
    
    double complex w = 0+1*I; //The imaginary number i
    fftw_plan my_plan2;
    fftw_plan my_plan3;

    while(L!=2){ //Iterate until there are only two columns remaining
        int c; //New number of columns (ceiling)
        int d; //New number of columns (floor)
        if(L % 2 == 0){ //If L even, c and d should be the same
            c = L/2;
            d = L/2;
        }else{ //If L odd, c is the ceiling, d is the floor
            c = L/2 + 1;
            d = L/2;
        }

        for(int i = 0; i < d; i++){
            for(int j = 0; j < H; j++){
                //B is a complex array declared with sufficient memory to handle all potential new dimensions (c||d,H) 
                //and is designated to store the pointwise multiplication of paired columns
                //The memory for B is reused in each iteration.
                B[i*H+j] = A[2*i*H + j] * A[(2*i+1)*H +j]; 
            }
        }
        

        if(d != c){ //If L odd, logic set up by c and d earlier
            for(int i = 0; i < H; i++){
                B[(c-1)*H+i] = A[(L-1)*H+i];
            }
        }
        
        my_plan2 = fftw_plan_r2r_1d(H,hc_to_real,in,FFTW_HC2R, FFTW_ESTIMATE);
        my_plan3 = fftw_plan_dft_1d(H, U, out, FFTW_FORWARD, FFTW_ESTIMATE);

        if(H != 4){ 
          //Dynamically reallocating hc_to_real as it is freed by fftw_destroy_plan
            hc_to_real = realloc(hc_to_real, H * sizeof(double));
        }
        
        //Pre-compute omega
        for(int j = 0; j < H; j++){
          w_vector[j] = cexp(-w*M_PI/H*j);
        }
        
        // Inverse FFT to retreive pmf values into U array
        int im_start;
        for(int i = 0; i < c; i++){
            im_start =  ((H+1)/2) - 1;
            for(int j = 0; j < H; j++){
              //Half complex property: (kth component real part in hc[k] and its imaginary part in hc[n-k] except for k==0 and n/2 for even)
              //Note: H is always even as we are simply pairing
              if(j <= H/2){
                hc_to_real[j] = creal(B[i*H+j]);
              }else{
                hc_to_real[j] = cimag(B[i*H+im_start]);
                im_start -= 1;
              }
              
            }

            fftw_execute_r2r(my_plan2, hc_to_real, in); //Reusing the plan to avoid reallocation
            
            //Executing pointwise multiplication in the same loop for efficiency
            for(int j = 0; j < H; j++){
              U[i*H + j] = in[j] / H * w_vector[j];
            }
            //hc_to_real = malloc(sizeof(double)*H);
        }
        
        for(int i = 0; i < c; i++){
            //Complex --> Complex 1D-FFT of the pointwise multiplication
            fftw_execute_dft(my_plan3, U + i*H, out);
            memcpy(U + i*H, out, sizeof(fftw_complex)*H);
        }
        
        //Advancing odd and even indicies
        for(int i = 0; i < 2*c*H; i++){
            if(i % 2 == 0){
                A[i] = B[i/2];
            }else{
                A[i] = U[i/2];
            }
        }

        //Updating size of columns and number of columns
        H *= 2;
        L = c;
        
        fftw_destroy_plan(my_plan2); //Plans can be destroyed in loop, confirmed via valgrind
        fftw_destroy_plan(my_plan3);
        
    }

    //Final multiplication and inverse FFT to retreive pmf

    for(int i = 0; i < H; i++){
        A[i] = A[i]*A[H+i];
    }

    fftw_plan my_plan4 = fftw_plan_dft_1d(H, A, out, FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_execute(my_plan4);
    for(int i = 0; i < H; i++){
       result[i] = cabs((out[i]/H));
    }
    free(hc_to_real);
    fftw_destroy_plan(my_plan4); //Plans must be destroyed before fftw_cleanup is called
    fftw_cleanup(); //Cleanup otherwise valgrind claims memory is still reachable
    return;
}

