%% Fast Fourier Transform - Real World Example

%% This example demonstrates processing of a real world recording, namely
%% a Northridge Earthquake accelerogram.
%% Mustafa Can Yucel 23.12.2021

clear all;
close all;

dataFile = "northridge.txt";
dt = 0.01; # time step in seconds

data = transpose(textread(dataFile));
n_data = size(data)(2); # number of data points

t_original = [0:n_data-1] .* dt;
td_original = n_data * dt;
subplot(2, 1, 1, 'align');
plot(t_original, data);
title(sprintf("Original Data with %d data points (%f seconds)", n_data, td_original));


sps = 1/dt;

# Direct fft with no preprocessing
fs_direct = abs(fft(data));
freq_resolution_direct =  sps / n_data;
freq_axis_direct = (0:n_data-1) .* freq_resolution_direct;
subplot(2, 1, 2, 'align');
stem(freq_axis_direct, fs_direct);
[max_val, idx] = max(fs_direct);
max_freq = idx * freq_resolution_direct;
title(sprintf("FFT with no preprocessing, max frequency is %f Hz", max_freq));

# Zero Padding
figure;
padding_constant = 2; # we add this times original signal length of zeros
data_padded = padarray(data, (size(data)) .* [0 padding_constant], 0, 'post');
td_padded = td_original * (padding_constant+1); # we have original signal plus zeros
n_padded = sps * td_padded;
t_padded = (0:n_padded-1)*dt;
y_padded = fft(data_padded);
freq_spectrum_padded = abs(y_padded);
subplot(6,1,1,"align");
plot(t_padded, data_padded);

title(sprintf("Data Zero Padded(%d sps %f seconds)", sps, td_padded));
fft_res_padded = sps / n_padded; # fft resolution
freq_axis_padded = (0:n_padded-1).*fft_res_padded; # frequency axis
subplot(6,1,2,"align");
plot(freq_axis_padded, freq_spectrum_padded);
[max_val, idx] = max(freq_spectrum_padded);
max_freq = idx * fft_res_padded;
title(sprintf("Max frequency is %f Hz", max_freq));
# Windowing
window = transpose(hanning(n_data));
data_windowed = window .* data;
subplot(6,1,3,"align");
plot(t_original, data_windowed);
title("Data Hann Windowed");
# FFT
y_windowed = fft(data_windowed);
freq_spectrum_windowed = abs(y_windowed); # frequency spectrum (L2 norm / abs )
subplot(6,1,4,"align");
plot(freq_axis_direct, freq_spectrum_windowed);
[max_val, idx] = max(freq_spectrum_windowed);
max_freq = idx * freq_resolution_direct;
title(sprintf("Max frequency is %f Hz", max_freq));

# Both padding and windowing
window = transpose(hanning(n_padded));
data_padded_windowed = window .* data_padded;
subplot(6,1,5,"align");
plot(t_padded, data_padded_windowed);
title("Data Padded and Hann Windowed");
# FFT
y_padded_windowed = fft(data_padded_windowed);
freq_spectrum_padded_windowed = abs(y_padded_windowed); # frequency spectrum (L2 norm / abs )
subplot(6,1,6,"align");
plot(freq_axis_padded, freq_spectrum_padded_windowed);# FFT
[max_val, idx] = max(freq_spectrum_padded_windowed);
max_freq = idx * fft_res_padded;
title(sprintf("Max frequency is %f Hz", max_freq));