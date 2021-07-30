% allanvarianceWilcoxon.m
% Calculates Allan variance from Wilcoxon accelerometer data and estimates
% interferometer sensitivity.


clear all


%% Importing data

%accel_matrix = importAccelFile('10_5_21_ITE_ ASD.csv');
accel_matrix =  readmatrix('10_5_21_ITE_ ASD.csv');
% The PSD is already calculated in Matt's code and is included in this CSV

freq = accel_matrix(:,1);   % Hz
raw_data_g = accel_matrix(:,2);   % g
psd_g = accel_matrix(:,3);    % g/sqrt(Hz)

% Converting to different units for convenience
raw_data_gal = raw_data_g*980;  % 1 m/s^2 = 100 gal
psd_si = psd_g*9.8; % in m/s^2
psd_gal = psd_g*980;    % in gal


% figure(925)
% plot(freq,psd_gal)
% xlabel('Frequency (Hz)')
% ylabel('Acceleration spectral density (gal/Hz^{1/2})')

figure(926)
loglog(freq,psd_gal)
xlabel('Frequency (Hz)')
ylabel('Acceleration spectral density (gal/Hz^{1/2})')

figure(930)
loglog(freq,psd_si)
xlabel('Frequency (Hz)')
ylabel('Acceleration spectral density (m/s^2/Hz^{1/2})')
title('Measured acceleration spectral density on the minusK')



%% Constants

Omega = 4e4;    % Rabi frequency, Hz        ours is 4e4
tau = 1/(2*Omega);  % pi pulse length, s
TT = 20e-3; % interferometer time T, s      ours is 20e-3 max
wavelength = 780e-9;    % m
keff = 4*pi/wavelength;


%% Transfer function

% freqpoints = 0.1:0.1:1e3; % use this to generate a purely theoretical
% curve (does not work with the Allan variance calculation, bu you can view
% the transfer function at arbitrary frequencies

freqpoints = freq;  % takes the frequency points from the accelerometer data for easy matching with the transfer function
Tc = 6; % cycle rep rate, s
fc = 1/Tc;

% Transfer function
% Direct from equation 2.54 from Rosi's thesis
HH = ( 4*2*pi*freqpoints*2*pi*Omega ./ ((2*pi*freqpoints).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*freqpoints*(TT + tau)/2) .* ( cos(2*pi*freqpoints*(TT + tau)/2) + (2*pi*Omega./(2*pi*freqpoints)).*sin(2*pi*freqpoints*TT/2)); % eqn 2.54 Rosi thesis
HH_a_sq = ((keff^2)./((2*pi*freqpoints).^4)) .* abs(HH.^2); % transfer function squared
HH_a_0_sq = ((keff^2)./((2*pi*freqpoints(2)).^4)) .* abs(HH(2).^2); % H(0) is undefined, so taking the smallest non-zero frequency and normalising to that
HH_a_norm = HH_a_sq./HH_a_0_sq; % Normalised transfer function

acc_trans = HH_a_norm.*psd_si;  % Applying the transfer function to the acceleration profile
% This should show how the interferometer itself filters acceleration noise

% Phase noise
S_phi = keff^2 ./ freqpoints.^4 .* psd_si;  % I don't use this for anything, can ignore

figure(931)
loglog(freq,acc_trans)
xlabel('Frequency (Hz)')
ylabel('Acceleration noise (m/s^2/Hz^{1/2})')
title('Vibration noise filtering')
ylim([1e-10 1e-3])



%% Allan variance

mm = 10000;   % total number of repetitions of the experiment
total_time = mm*Tc; % total integration time
tt = 1:1:mm;    % number of measurement cycles    

% Sampling the transfer function and acceleration PSD.
% This follows equation 2.68 from Rosi.
% For the acceleration PSD data, we have to sample the spectrum in steps of
% fc. In the code below, I set up the frequency sample points that I want,
% then use the "find" function to grab the closest frequency to the one I
% want, since the accelerometer data won't have the exact frequency.

max_nn = round(max(freq)/fc);   % Total number of frequency points
freq_sample_points = min(freq):fc:max(freq);  % Need to sample our measured PSD, so have to pick out the frequencies

freq_index_array = [];

% Loop through each multple of fc and grab the index of the nearest
% frequency in the data.
% We are looking for frequencies that fall in a window of +/- fc/30. In
% case of multiple results, we just grab the first result that's found.
% If Tc is reduced, then we must decrease the window size. Eg. Tc = 1 s
% needs +/-fc/195 in order to retrieve sensible results.

for nn = 1:1:max_nn
    
    freq_index = find(freq > (freq_sample_points(nn) - fc/30) & freq < (freq_sample_points(nn) + fc/30),1);
    freq_index_array = [freq_index_array; freq_index];
    
end

freq_samples = freq(freq_index_array);  % for plotting later

nn = 1:1:max_nn;  % for summing over. fc*nn <= max(freq)

% Now we use Rosi equation 2.68:
HH_nn = ( 4*2*pi*nn*fc*2*pi*Omega ./ ((2*pi*nn*fc).^2 - (2*pi*Omega)^2) ) .* sin(2*pi*nn*fc*(TT + tau)/2) .* ( cos(2*pi*nn*fc*(TT + tau)/2) + (2*pi*Omega./(2*pi*nn*fc)).*sin(2*pi*nn*fc*TT/2));

HH_nn_sq = abs(HH_nn).^2;
HH_total = HH_nn_sq ./ (2*pi*nn*fc).^4;

S_a = psd_si(freq_index_array);   % acceleration noise density at sampled frequency values

av_nn = HH_total.*S_a';
av_sum = sum(av_nn);

% Allan variance for acceleration, as per Freier thesis equation 2.16:
av_g = keff^2./(tt*Tc) * av_sum / (keff^2 * TT^4);


figure(927)
loglog(freqpoints,HH.^2)
xlabel('Frequency (Hz)')
ylabel('H(2{\pi}f)^2')
ylim([0.001 20])

figure(928)
loglog(freqpoints,HH_a_norm)
xlabel('Frequency (Hz)')
ylabel('H_a(2{\pi}f)^2 (normalised)')
title('Normalised transfer function')
ylim([1e-7 1])

figure(929)
loglog((1:1:mm)*Tc,av_g)
xlabel('Averaging time (s)')
ylabel('Allan variance m/s^2/Hz^{1/2}')

% The important plot!
figure(934)
loglog((1:1:mm)*Tc,av_g*100)
xlabel('Averaging time (s)')
ylabel('Allan variance (gal/Hz^{1/2})')
title('T = 20 ms, T_c = 6 s')
axis tight

% Useful plot to check if the acceleration PSD sampling is working as expected
figure(933)
loglog(freq,psd_si)
hold on
loglog(freq_samples,S_a,'x')
hold off
xlabel('Frequency (Hz)')
ylabel('Acceleration spectral density (m/s^2/Hz^{1/2})')
title('Sampled spectral density points')



