% allanvarianceWilcoxon_comparison.m
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



%% Allan variance - Tc = 5 s

Tc_4s = 4; % cycle rep rate, s
fc_4s = 1/Tc_4s;

mm = 10000;   % total number of repetitions
total_time_4s = mm*Tc_4s; % total integration time
tt_4s = 1:1:mm;    % number of measurement cycles    

max_nn_4s = round(max(freq)/fc_4s);
freq_sample_points_4s = min(freq):fc_4s:max(freq);  % Need to sample our measured PSD, so have to pick out the frequencies

freq_index_array_4s = [];

for nn_4s = 1:1:max_nn_4s
    
    freq_index_4s = find(freq > (freq_sample_points_4s(nn_4s) - fc_4s/35) & freq < (freq_sample_points_4s(nn_4s) + fc_4s/35),1);
    freq_index_array_4s = [freq_index_array_4s; freq_index_4s];
    
end

freq_samples_4s = freq(freq_index_array_4s);

nn_4s = 1:1:max_nn_4s;  % for summing over. fc*nn <= max(freq)

HH_nn_4s = ( 4*2*pi*nn_4s*fc_4s*2*pi*Omega ./ ((2*pi*nn_4s*fc_4s).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*nn_4s*fc_4s*(TT + tau)/2) .* ( cos(2*pi*nn_4s*fc_4s*(TT + tau)/2) + (2*pi*Omega./(2*pi*nn_4s*fc_4s)).*sin(2*pi*nn_4s*fc_4s*TT/2));

HH_nn_sq_4s = abs(HH_nn_4s).^2;
HH_total_4s = HH_nn_sq_4s ./ (2*pi*nn_4s*fc_4s).^4;

S_a_4s = psd_si(freq_index_array_4s);   % acceleration noise density at sampled frequency values

av_nn_4s = HH_total_4s.*S_a_4s';
av_sum_4s = sum(av_nn_4s);
av_g_4s = keff^2./(tt_4s*Tc_4s) * av_sum_4s / (keff^2 * TT^4);


%% Allan variance - Tc = 2

Tc_2s = 2; % cycle rep rate, s
fc_2s = 1/Tc_2s;

mm = 10000;   % total number of repetitions
total_time_2s = mm*Tc_2s; % total integration time
tt_2s = 1:1:mm;    % number of measurement cycles    

max_nn_2s = round(max(freq)/fc_2s);
freq_sample_points_2s = min(freq):fc_2s:max(freq);  % Need to sample our measured PSD, so have to pick out the frequencies

freq_index_array_2s = [];

for nn_2s = 1:1:max_nn_2s
    
    freq_index_2s = find(freq > (freq_sample_points_2s(nn_2s) - fc_2s/95) & freq < (freq_sample_points_2s(nn_2s) + fc_2s/95),1);
    freq_index_array_2s = [freq_index_array_2s; freq_index_2s];
    
end

freq_samples_2s = freq(freq_index_array_2s);

nn_2s = 1:1:max_nn_2s;  % for summing over. fc*nn <= max(freq)

HH_nn_2s = ( 4*2*pi*nn_2s*fc_2s*2*pi*Omega ./ ((2*pi*nn_2s*fc_2s).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*nn_2s*fc_2s*(TT + tau)/2) .* ( cos(2*pi*nn_2s*fc_2s*(TT + tau)/2) + (2*pi*Omega./(2*pi*nn_2s*fc_2s)).*sin(2*pi*nn_2s*fc_2s*TT/2));

HH_nn_sq_2s = abs(HH_nn_2s).^2;
HH_total_2s = HH_nn_sq_2s ./ (2*pi*nn_2s*fc_2s).^4;

S_a_2s = psd_si(freq_index_array_2s);   % acceleration noise density at sampled frequency values

av_nn_2s = HH_total_2s.*S_a_2s';
av_sum_2s = sum(av_nn_2s);
av_g_2s = keff^2./(tt_2s*Tc_2s) * av_sum_2s / (keff^2 * TT^4);


%% Allan variance - Tc = 1

Tc_1s = 1; % cycle rep rate, s
fc_1s = 1/Tc_1s;

mm = 10000;   % total number of repetitions
total_time_1s = mm*Tc_1s; % total integration time
tt_1s = 1:1:mm;    % number of measurement cycles    

max_nn_1s = round(max(freq)/fc_1s);
freq_sample_points_1s = min(freq):fc_1s:max(freq);  % Need to sample our measured PSD, so have to pick out the frequencies

freq_index_array_1s = [];

for nn_1s = 1:1:max_nn_1s
    
    freq_index_1s = find(freq > (freq_sample_points_1s(nn_1s) - fc_1s/195) & freq < (freq_sample_points_1s(nn_1s) + fc_1s/195),1);
    freq_index_array_1s = [freq_index_array_1s; freq_index_1s];
    
end

freq_samples_1s = freq(freq_index_array_1s);

nn_1s = 1:1:max_nn_1s;  % for summing over. fc*nn <= max(freq)

HH_nn_1s = ( 4*2*pi*nn_1s*fc_1s*2*pi*Omega ./ ((2*pi*nn_1s*fc_1s).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*nn_1s*fc_1s*(TT + tau)/2) .* ( cos(2*pi*nn_1s*fc_1s*(TT + tau)/2) + (2*pi*Omega./(2*pi*nn_1s*fc_1s)).*sin(2*pi*nn_1s*fc_1s*TT/2));

HH_nn_sq_1s = abs(HH_nn_1s).^2;
HH_total_1s = HH_nn_sq_1s ./ (2*pi*nn_1s*fc_1s).^4;

S_a_1s = psd_si(freq_index_array_1s);   % acceleration noise density at sampled frequency values

av_nn_1s = HH_total_1s.*S_a_1s';
av_sum_1s = sum(av_nn_1s);
av_g_1s = keff^2./(tt_1s*Tc_1s) * av_sum_1s / (keff^2 * TT^4);


%% Allan variance - Tc = 0.5 s

Tc_500ms = 0.5; % cycle rep rate, s
fc_500ms = 1/Tc_500ms;

mm = 10000;   % total number of repetitions
total_time_500ms = mm*Tc_500ms; % total integration time
tt_500ms = 1:1:mm;    % number of measurement cycles    

max_nn_500ms = round(max(freq)/fc_500ms);
freq_sample_points_500ms = min(freq):fc_500ms:max(freq);  % Need to sample our measured PSD, so have to pick out the frequencies

freq_index_array_500ms = [];

for nn_500ms = 1:1:max_nn_500ms
    
    freq_index_500ms = find(freq > (freq_sample_points_500ms(nn_500ms) - fc_500ms/215) & freq < (freq_sample_points_500ms(nn_500ms) + fc_500ms/215),1);
    freq_index_array_500ms = [freq_index_array_500ms; freq_index_500ms];
    
end

freq_samples_500ms = freq(freq_index_array_500ms);

nn_500ms = 1:1:max_nn_500ms;  % for summing over. fc*nn <= max(freq)

HH_nn_500ms = ( 4*2*pi*nn_500ms*fc_500ms*2*pi*Omega ./ ((2*pi*nn_500ms*fc_500ms).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*nn_500ms*fc_500ms*(TT + tau)/2) .* ( cos(2*pi*nn_500ms*fc_500ms*(TT + tau)/2) + (2*pi*Omega./(2*pi*nn_500ms*fc_500ms)).*sin(2*pi*nn_500ms*fc_500ms*TT/2));

HH_nn_sq_500ms = abs(HH_nn_500ms).^2;
HH_total_500ms = HH_nn_sq_500ms ./ (2*pi*nn_500ms*fc_500ms).^4;

S_a_500ms = psd_si(freq_index_array_500ms);   % acceleration noise density at sampled frequency values

av_nn_500ms = HH_total_500ms.*S_a_500ms';
av_sum_500ms = sum(av_nn_500ms);
av_g_500ms = keff^2./(tt_500ms*Tc_500ms) * av_sum_500ms / (keff^2 * TT^4);


%% Allan variance - Tc = 0.25 s

Tc_250ms = 0.25; % cycle rep rate, s
fc_250ms = 1/Tc_250ms;

mm = 10000;   % total number of repetitions
total_time_250ms = mm*Tc_250ms; % total integration time
tt_250ms = 1:1:mm;    % number of measurement cycles    

max_nn_250ms = round(max(freq)/fc_250ms);
freq_sample_points_250ms = min(freq):fc_250ms:max(freq);  % Need to sample our measured PSD, so have to pick out the frequencies

freq_index_array_250ms = [];

for nn_250ms = 1:1:max_nn_250ms
    
    freq_index_250ms = find(freq > (freq_sample_points_250ms(nn_250ms) - fc_250ms/215) & freq < (freq_sample_points_250ms(nn_250ms) + fc_250ms/215),1);
    freq_index_array_250ms = [freq_index_array_250ms; freq_index_250ms];
    
end

freq_samples_250ms = freq(freq_index_array_250ms);

nn_250ms = 1:1:max_nn_250ms;  % for summing over. fc*nn <= max(freq)

HH_nn_250ms = ( 4*2*pi*nn_250ms*fc_250ms*2*pi*Omega ./ ((2*pi*nn_250ms*fc_250ms).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*nn_250ms*fc_250ms*(TT + tau)/2) .* ( cos(2*pi*nn_250ms*fc_250ms*(TT + tau)/2) + (2*pi*Omega./(2*pi*nn_250ms*fc_250ms)).*sin(2*pi*nn_250ms*fc_250ms*TT/2));

HH_nn_sq_250ms = abs(HH_nn_250ms).^2;
HH_total_250ms = HH_nn_sq_250ms ./ (2*pi*nn_250ms*fc_250ms).^4;

S_a_250ms = psd_si(freq_index_array_250ms);   % acceleration noise density at sampled frequency values

av_nn_250ms = HH_total_250ms.*S_a_250ms';
av_sum_250ms = sum(av_nn_250ms);
av_g_250ms = keff^2./(tt_250ms*Tc_250ms) * av_sum_250ms / (keff^2 * TT^4);


%% Plots

figure(935)
loglog((1:1:mm)*Tc_4s,av_g_4s*100)
hold on
loglog((1:1:mm)*Tc_2s,av_g_2s*100)
loglog((1:1:mm)*Tc_1s,av_g_1s*100)
loglog((1:1:mm)*Tc_500ms,av_g_500ms*100)
loglog((1:1:mm)*Tc_250ms,av_g_250ms*100)
hold off
xlabel('Averaging time (s)')
ylabel('Allan variance gal/Hz^{1/2}')
title('Allan variance for difference cycle times')
axis tight
legend('T_c = 4 s','T_c = 2 s','T_c = 1 s','T_c = 0.5 s','T_c = 0.25 s')

% figure(936)
% loglog(freq,psd_si)
% hold on
% loglog(freq_samples_5s,S_a_5s,'x')
% hold off
% xlabel('Frequency (Hz)')
% ylabel('Acceleration spectral density (m/s^2/Hz^{1/2})')
% title('Sampled spectral density points')



