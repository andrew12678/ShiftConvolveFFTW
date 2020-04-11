# ShiftConvolveFFTW

`ShiftConvolve` is a R package which uses exponential shifting and the Fast Fourier Transformations (FFT) to compute the (right) tail of the Poisson Binomial Distribution. 
This package makes use of the `Fastest Fourier Transform in the West` ([FFTW](http://www.fftw.org/)) library to perform the necessary DFT and Inverse DFT computations.

For the `ShiftConvolve` implementation which uses [minFFT](https://github.com/aimukhin/minfft) to perform the Fourier Transformations please go to [https://github.com/andrew12678/ShiftConvolve](https://github.com/andrew12678/ShiftConvolve).

In the interest of speed we all significant computational aspects of our procedure are executed by code written in C and the R code mainly acts as a wrapper around that compiled C code.

## Dependencies

As eluded to above, `FFTW` will be a necessary dependency.

This package was originally compiled in an `Ubuntu 18.04` Linux operating system using `R 3.5.1` and `FFTW 3.3.5` but was later recompiled (without any modification of the configurations) in a `Manjaro 19.02 KDE Plasma` Linux operating system with inbuilt `FFTW 3.3.8` and `R 3.6.3`. 
Since the `FFTW` came along with the `Manjaro` installation, we will only outline installations on the `Ubuntu` system which does not automatically install `FFTW`.

For `Ubuntu` the steps included:

```console
wget http://fftw.org/fftw-3.3.5.tar.gz
tar -xzf fftw-3.3.5.tar.gz
cd fftw-3.3.5
./configure --enable-shared
make
sudo make install
```


## Installation

The most simple installation involves simply cloning this repository and installing `ShiftConvolve` from source. 

```console
git clone https://github.com/andrew12678/ShiftConvolveFFTW.git
cd ShiftConvolve
install.packages('ShiftConvolvePoibin_1.11.0.tar.gz', repos = NULL, type="source")
```

An alternative installation procedure involves cloning the repository, creating a `RStudio` project in the `ShiftConvolvePoibin` folder and then building. 

## Examples

An simple example with the uniform distribution

```R
library(ShiftConvolvePoibin)
set.seed(18)
n = 10000
p = runif(n)
s0 = 5200
shiftpval(p, s0)	# compute the p-value, or right tail at s0
shiftpval(1-p, n-s0)	# compute the p-value, or left tail at s0
```
