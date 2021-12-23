%% Fast Fourier Transform - Frequency Resolution and Zero Padding
%% This example demonstrates the effect of signal length and zero padding
%% on frequency resolution.
%% Mustafa Can Yucel 23.12.2020

clear all;
close all;
pkg load image;

# Frequencies
# Sampling frequency in sample per seconds
sps = 100; # Hz or sample/seconds
# Wave frequencies in Hz
f1 = 4;
f2 = 7.5;
# Amplitudes
a1 = 1;
a2 = 1;
# Dampings
c1 = 4;
c2 = 4;

ts = 1/sps; # sample time increment

%% CASE 1: SAMPLING DURATION IS NOT ENOUGH
%% When the resolution of original data is not enough to resolve contributing
%% frequencies, zero padding will not achieve pinpoint accuracy results even
%% in case of zero noise, since the information in the original 
%% data is insufficient.

%% In this example, the short signal frequency resolution is 
%% 1/0.5 = 2 Hz, therefore it can distinguish only multiples of 2 Hz in the 
%% signal. Thus, the second peak of 7 Hz (f2) is not visible 
%% in the FFT plot. If f2 was 8 Hz, it would be visible.

# short signal duration in seconds
td_short = 0.5;
t_short = (0:td_short*sps)*ts; # short duration signal time axis
n_short = sps * td_short; # short signal total number of data points
s_short = (a1 * sin(2*pi*f1*t_short) .* exp(-t_short/c1)) + (a2 * sin(2*pi*f2*t_short) .* exp(-t_short/c2));
subplot(6,1,1,"align");
plot(t_short, s_short);
title(sprintf("A short wave of %d and %d Hz (%d sps %f seconds)", f1, f2, sps, td_short));

# FFT
y_short = fft(s_short);
freq_spectrum_short = abs(y_short); # frequency spectrum (L2 norm / abs )
fft_res_short = sps / n_short; # fft resolution
freq_axis_short = (0:n_short).*fft_res_short; # frequency axis
subplot(6,1,2,"align");
plot(freq_axis_short, freq_spectrum_short);

# if f2 was 8 Hz instead of 7
f22 = 8; # Hz
s_short2 = (a1 * sin(2*pi*f1*t_short) .* exp(-t_short/c1)) + (a2 * sin(2*pi*f22*t_short) .* exp(-t_short/c2));
y_short2 = fft(s_short2);
freq_spectrum_short2 = abs(y_short2); # frequency spectrum (L2 norm / abs )
subplot(6,1,3,"align");
plot(t_short, s_short2);
title(sprintf("A short wave of %d and %d Hz (%d sps %f seconds)", f1, f22, sps, td_short));
subplot(6,1,4,"align");
plot(freq_axis_short, freq_spectrum_short2);

% Pad zeroes to the short signal and look again
padding_constant = 3; # we add this times original signal length of zeros
s_short_padded = padarray(s_short, (size(s_short) - [0 1]) .* [0 padding_constant], 0, 'post');
td_short_padded = td_short * (padding_constant+1); # we have original signal plus zeros
t_short_padded = (0:td_short_padded*sps)*ts;
y_short_padded = fft(s_short_padded);
freq_spectrum_short_padded = abs(y_short_padded);
subplot(6,1,5,"align");
plot(t_short_padded, s_short_padded);
n_short_padded = sps * td_short_padded;
title(sprintf("A short wave of %d and %d Hz Zero Padded(%d sps %f seconds)", f1, f2, sps, td_short_padded));
fft_res_short_padded = sps / n_short_padded; # fft resolution
freq_axis_short_padded = (0:n_short_padded).*fft_res_short_padded; # frequency axis
subplot(6,1,6,"align");
plot(freq_axis_short_padded, freq_spectrum_short_padded);

%% CASE 2: SAMPLING DURATION IS ADEQUATE
%% When the resolution of original data is enough to resolve contributing
%% frequencies but it is not a multiple of target frequencies, zero padding
%% will allow to identify them.

%% In this example, the long signal frequency resolution is 
%% 1/3 = 0.091 Hz, therefore it can distinguish f1 (4) and f2(7) frequencies.
%% However, neither 4 nor 7 is a multiple of 0.091, therefore the energy of 
%% these frequencies will leak into neighboring bins, suppressing the peaks.
%% Zero-padding enough values so that padded-frequency resolution will be a 
%% multiple of these frequencies will solve this issue.

# long signal duration in seconds
figure;
td_long = 3;
t_long = (0:td_long*sps)*ts; # long duration signal time axis
n_long = sps * td_long; # long signal total number of data points

s_long = (a1 * sin(2*pi*f1*t_long) .* exp(-t_long/c1)) + (a2 * sin(2*pi*f2*t_long) .* exp(-t_long/c2));
subplot(4,1,1,"align");
plot(t_long, s_long);
title(sprintf("A long wave of %d and %d Hz (%d sps %f seconds)", f1, f2, sps, td_long));

# FFT
y_long = fft(s_long);
freq_spectrum_long = abs(y_long); # frequency spectrum (L2 norm / abs )
fft_res_long = sps / n_long; # fft resolution
freq_axis_long = (0:n_long).*fft_res_long; # frequency axis
subplot(4,1,2,"align");
plot(freq_axis_long, freq_spectrum_long);

% Pad zeroes to long signal and look again
padding_constant = 3; # we add this times original signal length of zeros
s_long_padded = padarray(s_long, (size(s_long) - [0 1]) .* [0 padding_constant], 0, 'post');
td_long_padded = td_long * (padding_constant+1); # we have original signal plus zeros
t_long_padded = (0:td_long_padded*sps)*ts;
y_long_padded = fft(s_long_padded);
freq_spectrum_long_padded = abs(y_long_padded);
subplot(4,1,3,"align");
plot(t_long_padded, s_long_padded);
n_long_padded = sps * td_long_padded;
title(sprintf("A long wave of %d and %d Hz Zero Padded(%d sps %f seconds)", f1, f2, sps, td_long_padded));
fft_res_long_padded = sps / n_long_padded; # fft resolution
freq_axis_long_padded = (0:n_long_padded).*fft_res_long_padded; # frequency axis
subplot(4,1,4,"align");
plot(freq_axis_long_padded, freq_spectrum_long_padded);