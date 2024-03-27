% Read file
cd 'C:\Users\chph1\OneDrive\Documents\GitHub\EEGals_Firmware\matlab_test_environment'
results_file = "results.csv";
data_files = dir("data");

data = readtable("data/" + data_files(3).name);
data2=data;
Ntotal = size(data2, 1);

%% Convert voltage values to actual voltage values in uV
gain=24;

for col_index=2:9
    for row_index=1:Ntotal
        data.(col_index)(row_index) = data2.(col_index)(row_index)*(4500000/gain/(2^23-1));
    end
end

writetable(data, 'd2.csv');
edata = readtable('d2.csv');

% Sampling frequency and number of samples per channel 

fs = 250;

% Period of the signal 
t = 1/fs;

% Calculate the time duration of the signal 
T = Ntotal/fs; % approximately 60 seconds 

% Desired number of samples in an epoch 
No_Desired_Pts_Per_Epoch = 2048; % 2 to the power of 11

% Calculate the number of time points per epoch where each epoch is 8.192s
% seconds
Epoch_Length = No_Desired_Pts_Per_Epoch/fs; % 8.192 s
No_Time_Pts_Per_Epoch = fs * Epoch_Length; % 2048 pts

% Calculate the number of epochs within the time frame 
Total_Epochs = fix(T/Epoch_Length); % approximately (60 s/ 8.192s)

N = 0;
i = 0;

%Define frequency bands in Hz (may/may not need all of them)
alpha_band = [8 12];
beta_band = [13 30];
delta_band = [1 3];
gamma_band = [30 100]; 
theta_band = [4 7];

% Initialize sum of relative powers of each EEG band to zero (used for
% computing average value later)
sum_alpha_rel = 0;
sum_beta_rel = 0;
sum_delta_rel = 0;
sum_gamma_rel = 0;
sum_theta_rel = 0;


% Calculates the frequency bin width (frequency resolution)
% which is sampling rate(250 Hz) divided by samples or fft length (2048)per
% epoch
frequency_bin_width = fs / 2048;

% Calculates the indices of the fft that correspond to a certain frequency
% bin which contains the frequncy ranges of interest

% We need to search for the nearest frequency bin 

ind_alpha1 = alpha_band(1) / frequency_bin_width;
ind_alpha2 = alpha_band(2) / frequency_bin_width;

ind_alpha_low = floor(ind_alpha1);
ind_alpha_high = ceil(ind_alpha2);

ind_beta1 = beta_band(1) / frequency_bin_width;
ind_beta2 = beta_band(2) / frequency_bin_width;

ind_beta_low = floor(ind_beta1);
ind_beta_high = ceil(ind_beta1);

ind_delta1 = delta_band(1) / frequency_bin_width;
ind_delta2 = delta_band(2) / frequency_bin_width;

ind_delta_low = floor(ind_delta1);
ind_delta_high = ceil(ind_delta2);

ind_gamma1 = gamma_band(1) / frequency_bin_width;
ind_gamma2 = gamma_band(2) / frequency_bin_width;

ind_gamma_low = floor(ind_gamma1);
ind_gamma_high = ceil(ind_gamma2);

ind_theta1 = theta_band(1) / frequency_bin_width;
ind_theta2 = theta_band(2) / frequency_bin_width;

ind_theta_low = floor(ind_theta1);
ind_theta_high = ceil(ind_theta2);


% Loop through each channel 
for col_index=2:9
   
    writematrix(["Channel",col_index-1], results_file, "WriteMode", "append");

    end_point = No_Time_Pts_Per_Epoch;
    start_point = 1;
    
    % sum of periodograms for all epochs for one channel. So col are the
    % channel and row are the number of values in an epoch
    sum_periodograms = zeros(No_Time_Pts_Per_Epoch/2 , 1);
    
    i = 1;
    i2 = 2;
    % For each channel, loop through the number of epochs 
    for epoch_index=1:Total_Epochs
      
      
      voltageData = edata{start_point:end_point,col_index};
      
      % To reduce the spike at frequency 0
      voltageData = voltageData - mean(voltageData);
      
      N = length(voltageData);
      
      % Note the fft produces complex numbers because it is the best way to
      % contain information about both the phase and amplitude of a signal 
      Y = fft(voltageData);
      writematrix(abs(Y), results_file, "WriteMode", "append");

      % We are looking at one-sided spectrum so we want half of the frequencies. 
      % We also don't require the signal at frequency at 0 since it is just
      % the DC component
      Y_one_side = Y(1:N/2);
         
      Y_one_side(2:end) = 2*Y_one_side(2:end);
      
      % This line produces the frequency bin starting with 0 Hz all the way
      % to the frequency bin N/2-1 Hz. We will not get the N/2 Hz frequency
      % bin
      f_one_side = (0:N/2-1) * fs / N; 
%{      
      figure(col_index*100 + i);
      i = i+2;
      plot(f_one_side, abs(Y_one_side));
      title('FFT of Voltage Data');
      xlabel('Frequency (Hz)');
      ylabel('Magnitude');
      grid on;
%}      
    
      % Taking the absolute value squared is the same as multiplying a complex
      % vector by its conjugate which is used as a step in a series
      % of steps to get the power spectrum 
      psdx = (1/N*fs) * abs(Y_one_side).^2; 
  
      % We need to multiply by 2 because we are only looking at one side of
      % the spectrum so we need to double the power. We don't need the
      % first part of psdx because the DC component does not need doubling
      psdx(2:end-1) = 2*psdx(2:end-1);
      
%      figure(col_index*100 + i2);
      i2 = i2+2;
      
      % Plot the PSD for the frequency 
%      plot(f_one_side, pow2db(psdx));
      
      % Sum the periodograms/psd for each epoch and divide by the number of epochs.
      sum_periodograms = sum_periodograms + psdx;
%{      
      title('Periodogram of Voltage Data');
      xlabel('Frequency (Hz)');
      ylabel('Power (dB)');
      grid on;
%}      
      start_point = end_point;
      end_point = end_point + No_Time_Pts_Per_Epoch;
      
    end 
    
    % Average all the periodograms across the epochs (Welch's Method)
    pwelch = sum_periodograms/Total_Epochs;
%{
    % Plot the pwelch for each channel 
    figure(col_index*100);
    f_one_side = (0:N/2-1) * fs / N;
    plot(f_one_side, pow2db(pwelch));
    
    title('PWelch of Voltage Data');
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    grid on;
%}
    % Band powers of the different frequency bands
    alpha_bp = sum(pwelch(ind_alpha1:ind_alpha2));
    beta_bp = sum(pwelch(ind_beta1:ind_beta2));
    delta_bp = sum(pwelch(ind_delta1:ind_delta2));
    gamma_bp = sum(pwelch(ind_gamma1:ind_gamma2));
    theta_bp = sum(pwelch(ind_theta1:ind_theta2));
    
    % Total Power 
    total_power = sum(pwelch);
    
    % Relative Power (Band Power/Total Power)
    alpha_rel = alpha_bp / total_power;
    beta_rel = beta_bp / total_power;
    delta_rel = delta_bp / total_power;
    gamma_rel = gamma_bp / total_power;
    theta_rel = theta_bp / total_power;
   
    writematrix([alpha_rel, beta_rel, delta_rel, gamma_rel, theta_rel], results_file, "WriteMode", "append");
end












