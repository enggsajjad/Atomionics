% allanvarianceWilcoxon_comparison_acc.m
% Calculates Allan variance from Wilcoxon accelerometer data and estimates
% interferometer sensitivity.


clear all

close all


%% Importing data

%accel_matrix = importAccelFile('10_5_21_ITE_ ASD.csv');
accel_matrix =  readmatrix('10_5_21_ITE_ ASD.csv');
% The PSD is already calculated in Matt's code and is included in this CSV

freq = accel_matrix(:,1);   % Hz
raw_data_g = accel_matrix(:,2);   % g
psd_g = accel_matrix(:,3);    % g/sqrt(Hz)

raw_data_gal = raw_data_g*980;  % 1 m/s^2 = 100 gal
psd_si = psd_g*9.8; % in m/s^2
psd_gal = psd_g*980;


% figure(925)
% plot(freq,psd_gal)
% xlabel('Frequency (Hz)')
% ylabel('Acceleration spectral density (gal/Hz^{1/2})')

% figure(926)
% loglog(freq,(psd_gal*1000).^2)
% xlabel('Frequency (Hz)')
% ylabel('Acceleration spectral density (mgal^2/Hz)')
% 
% figure(930)
% loglog(freq,psd_si)
% xlabel('Frequency (Hz)')
% ylabel('Acceleration spectral density (m/s^2/Hz^{1/2})')
% title('Measured acceleration spectral density on the minusK')



%% Constants

Omega = 40e3;    % Rabi frequency, Hz        ours is 4e4
tau = 1/(2*Omega);  % pi pulse length, s
TT = 20e-3; % interferometer time T, s      ours is 20e-3 max
wavelength = 780e-9;    % m
keff = 4*pi/wavelength;


%% Transfer function

% freqpoints = 0.1:0.1:1e3;
freqpoints = freq;

% Transfer function
HH = ( 4*2*pi*freqpoints*2*pi*Omega ./ ((2*pi*freqpoints).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*freqpoints*(TT + tau)/2) .* ( cos(2*pi*freqpoints*(TT + tau)/2) + (2*pi*Omega./(2*pi*freqpoints)).*sin(2*pi*freqpoints*TT/2)); % eqn 2.54 Rosi thesis
HH_a_sq = ((keff^2)./((2*pi*freqpoints).^4)) .* abs(HH.^2); 
HH_a_0_sq = ((keff^2)./((2*pi*freqpoints(2)).^4)) .* abs(HH(2).^2);
HH_a_norm = HH_a_sq./HH_a_0_sq;

acc_trans = HH_a_norm.*psd_si;

% Phase noise
S_phi = keff^2 ./ freqpoints.^4 .* psd_si;

% figure(931)
% loglog(freq,acc_trans)
% xlabel('Frequency (Hz)')
% ylabel('Acceleration noise (m/s^2/Hz^{1/2})')
% title('Vibration noise filtering')
% ylim([1e-10 1e-3])



%% Allan variance

Tc = 1; % cycle rep rate, s
fc = 1/Tc;

mm = 10000;   % total number of repetitions
total_time = mm*Tc; % total integration time
tt = 1:1:mm;    % number of measurement cycles    

max_nn = round(max(freq)/fc);
freq_sample_points = min(freq):fc:max(freq);  % Need to sample our measured PSD, so have to pick out the frequencies

freq_index_array = [];

for nn = 1:1:max_nn
    
    freq_index = find(freq > (freq_sample_points(nn) - fc/195) & freq < (freq_sample_points(nn) + fc/195),1);
    freq_index_array = [freq_index_array; freq_index];
    
end

freq_samples = freq(freq_index_array);

S_a = psd_si(freq_index_array);   % acceleration noise density at sampled frequency values

nn = 1:1:max_nn;  % for summing over. fc*nn <= max(freq)


%% Original PSD

TT = 20e-3;

HH_nn = ( 4*2*pi*nn*fc*2*pi*Omega ./ ((2*pi*nn*fc).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*nn*fc*(TT + tau)/2) .* ( cos(2*pi*nn*fc*(TT + tau)/2) + (2*pi*Omega./(2*pi*nn*fc)).*sin(2*pi*nn*fc*TT/2));

HH_nn_sq = abs(HH_nn).^2;
HH_total = HH_nn_sq ./ (2*pi*nn*fc).^4;

av_nn = HH_total.*S_a';
av_sum = sum(av_nn);
av_g = keff^2./(tt*Tc) * av_sum / (keff^2 * TT^4);


%% 10x worse PSD

av_nn_10x = HH_total.*S_a'*10;
av_sum_10x = sum(av_nn_10x);
av_g_10x = keff^2./(tt*Tc) * av_sum_10x / (keff^2 * TT^4);


%% 10x better PSD

av_nn_p1x = HH_total.*S_a'*0.1;
av_sum_p1x = sum(av_nn_p1x);
av_g_p1x = keff^2./(tt*Tc) * av_sum_p1x / (keff^2 * TT^4);


%% Plots

figure(936)
loglog((1:1:mm)*Tc,av_g*100)
hold on
loglog((1:1:mm)*Tc,av_g_10x*100)
loglog((1:1:mm)*Tc,av_g_p1x*100)
hold off
xlabel('Averaging time (s)')
ylabel('Allan variance gal/Hz^{1/2}')
title('Allan variance at different acceleration noise amplitude')
axis tight
legend('Actual acceleration PSD','10x noisier','10x quieter')

% figure(936)
% loglog(freq,psd_si)
% hold on
% loglog(freq_samples_5s,S_a_5s,'x')
% hold off
% xlabel('Frequency (Hz)')
% ylabel('Acceleration spectral density (m/s^2/Hz^{1/2})')
% title('Sampled spectral density points')



