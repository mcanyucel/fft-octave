%% Fast Fourier Transform - Basics
%% This example demonstrates the basic FFT use for synthetic waves composed
%% of a combination of sinusoidal functions (each having different frequency, 
%% amplitude, and damping).
%% The effect of different levels of noise is also shown.
%% Mustafa Can Yucel 23.12.2020

clear all;
close all;

# Frequencies
# Sampling frequency in sample per seconds
sps = 100; # Hz or sample/seconds
# Wave frequencies in Hz
f1 = 4;
f2 = 7;
# total signal duration in seconds
td = 30;
# Amplitudes
a1 = 1;
a2 = 2;
# Dampings
c1 = 10;
c2 = 7;

ts = 1/sps; # sample time increment
t = (0:td*sps)*ts; # time axis
n = sps * td; # total number of data points

s = (a1 * sin(2*pi*f1*t) .* exp(-t/c1)) + (a2 * sin(2*pi*f2*t) .* exp(-t/c2));
subplot(2,1,1,"align");
plot(t, s);
title(sprintf("A wave of %d and %d Hz (no noise)", f1, f2));

# FFT
y = fft(s);
freq_spectrum = abs(y); # frequency spectrum (L2 norm / abs )
fft_res = sps / n; # fft resolution
freq_axis = (0:n).*fft_res; # frequency axis
subplot(2,1,2,"align");
plot(freq_axis, freq_spectrum);


# Gaussian noise - Light
wn = randn(size(s));
s_noisy = s + wn;
figure;
subplot(2,1,1,"align");
plot(t, s_noisy);
title(sprintf("A wave of %d and %d Hz (Light Gaussian Noise)", f1, f2));

# FFT
y_noisy = fft(s_noisy);
freq_spectrum_noisy = abs(y_noisy); # frequency spectrum (L2 norm / abs )
subplot(2,1,2,"align");
plot(freq_axis, freq_spectrum_noisy);

# Gaussian noise - Heavy
wn = 5 * randn(size(s));
s_noisy = s + wn;
figure;
subplot(2,1,1,"align");
plot(t, s_noisy);
title(sprintf("A wave of %d and %d Hz (Heavy Gaussian Noise)", f1, f2));

# FFT
y_noisy = fft(s_noisy);
freq_spectrum_noisy = abs(y_noisy); # frequency spectrum (L2 norm / abs )
subplot(2,1,2,"align");
plot(freq_axis, freq_spectrum_noisy);