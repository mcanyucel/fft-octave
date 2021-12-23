%% Fast Fourier Transform - Spectral Leakage
%% As the technical definition, spectral leakage is the creation of additional
%% frequency components in a signal as the result of any operation that is not 
%% linear time-invariant (sampling, for example).
%% It also arises when the recorded signal length is not a multiple of 
%% component periods, in other words, the signal is cut off abruptly.
%% This is almost always the case for real world measurements.
%% This effect can be limited by zero padding and windowing. Windowing allows 
%% the signal to start and end at zero; different window types have different
%% strengths and weaknesses.

%% This example demonstrates spectral leakage caused by the incomplete signal
%% input and the effects of zero padding together with windowing.
%% Mustafa Can Yucel 23.12.2020

clear all;
close all;

# Frequencies
# Sampling frequency in sample per seconds
sps = 100; # Hz or sample/seconds
# Wave frequencies in Hz
f1 = 4;
# Amplitudes
a1 = 1;
ts = 1/sps; # sample time increment
# No Damping


%% EXACT LENGTH
td_exact = 1; # total signal duration in seconds
t_exact = (0:td_exact*sps)*ts; # time axis
n_exact = sps * td_exact; # total number of data points
s_exact = (a1 * sin(2*pi*f1*t_exact));
subplot(4,1,1,"align");
plot(t_exact, s_exact);
title(sprintf("A wave of %d Hz with Exact Duration(no noise)", f1));

# FFT
y_exact = fft(s_exact);
freq_spectrum_exact = abs(y_exact); # frequency spectrum (L2 norm / abs )
fft_res_exact = sps / n_exact; # fft resolution
freq_axis_exact = (0:n_exact).*fft_res_exact; # frequency axis
subplot(4,1,2,"align");
stem(freq_axis_exact, freq_spectrum_exact);

# LEAKY LENGTH
td_leak = 1.35; # total signal duration in seconds
t_leak = (0:td_leak*sps)*ts; # time axis
n_leak = sps * td_leak; # total number of data points
s_leak = (a1 * sin(2*pi*f1*t_leak));
subplot(4,1,3,"align");
plot(t_leak, s_leak);
title(sprintf("A wave of %d Hz with Non-Exact Duration(no noise)", f1));
# FFT
y_leak = fft(s_leak);
freq_spectrum_leak = abs(y_leak); # frequency spectrum (L2 norm / abs )
fft_res_leak = sps / n_leak; # fft resolution
freq_axis_leak = (0:n_leak).*fft_res_leak; # frequency axis
subplot(4,1,4,"align");
stem(freq_axis_leak, freq_spectrum_leak);

%% ZERO PADDING AND WINDOWING THE LEAKY SERIE

## just zero padding
figure;
padding_constant = 3; # we add this times original signal length of zeros
s_long_padded = padarray(s_leak, (size(s_leak) - [0 1]) .* [0 padding_constant], 0, 'post');
td_long_padded = td_leak * (padding_constant+1); # we have original signal plus zeros
t_long_padded = (0:td_long_padded*sps)*ts;
y_long_padded = fft(s_long_padded);
freq_spectrum_long_padded = abs(y_long_padded);
subplot(6,1,1,"align");
plot(t_long_padded, s_long_padded);
n_long_padded = sps * td_long_padded;
title(sprintf("A Leaky wave of %d Hz Zero Padded(%d sps %f seconds)", f1, sps, td_long_padded));
fft_res_long_padded = sps / n_long_padded; # fft resolution
freq_axis_long_padded = (0:n_long_padded).*fft_res_long_padded; # frequency axis
subplot(6,1,2,"align");
stem(freq_axis_long_padded, freq_spectrum_long_padded);

## just windowing (Hann)
window = transpose(hanning(n_leak+1));
s_windowed = window .* s_leak;
subplot(6,1,3,"align");
plot(t_leak, s_windowed);
title(sprintf("A wave of %d Hz with Non-Exact Duration Hann Windowed(no noise)", f1));
# FFT
y_leak = fft(s_windowed);
freq_spectrum_leak = abs(y_leak); # frequency spectrum (L2 norm / abs )
fft_res_leak = sps / n_leak; # fft resolution
freq_axis_leak = (0:n_leak).*fft_res_leak; # frequency axis
subplot(6,1,4,"align");
stem(freq_axis_leak, freq_spectrum_leak);

## First zero padding then windowing
window = transpose(hanning(n_long_padded+1));
s_windowed = window .* s_long_padded;
subplot(6,1,5,"align");
plot(t_long_padded, s_windowed);
title(sprintf("A Leaky wave of %d Hz Zero Padded and Hann Windowed (%d sps %f seconds)", f1, sps, td_long_padded));
# FFT
y_leak = fft(s_windowed);
freq_spectrum_leak = abs(y_leak); # frequency spectrum (L2 norm / abs )
fft_res_leak = sps / n_long_padded; # fft resolution
freq_axis_leak = (0:n_long_padded).*fft_res_leak; # frequency axis
subplot(6,1,6,"align");
stem(freq_axis_leak, freq_spectrum_leak);
