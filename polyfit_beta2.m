lambda_range = 6e-7:4e-8:1.6e-6;
beta2_range = [5.51e-26 4.99e-26 4.6e-26 4.35e-26 3.9e-26 3.52e-26 3.32e-26 2.95e-26 ...
    2.71e-26 2.33e-26 2.02e-26 1.7e-26 1.5e-26 1.15e-26 0.95e-26 0.5e-26 0.2e-26 -0.1e-26 ...
    -0.5e-26 -1.1e-26 -1.5e-26 -2.0e-26 -2.35e-26 -2.7e-26 -3.0e-26 -3.5e-26];

lambda_range_extended = 6e-7:4e-8:1.8e-6;

p = polyfit(lambda_range, beta2_range, 2);
beta2_range_extended = polyval(p,lambda_range_extended);

%Constants and variables
h = 6.6261*10^-34;
c = 3*10^8;
nsp = 1.58;
G = 16;
ns = 5;
gamma = 1.2*10^-3;
L = 80;
alpha = 0.0461;
Ptx_dbm = -3;
Nch = 100;
Rs = 5*10^10;
Brx = 5*10^10;

%Linearizing and calculating Leff
Leff = (1-exp(-alpha*L))/alpha;
linear_G = 10^(G/10);
Ptx = 10^-3*10^(Ptx_dbm/10);

SER_simulated = zeros(size(lambda_range_extended));
SER_theoretical = zeros(size(lambda_range_extended));
BER_theoretical = zeros(size(lambda_range_extended));

for lambda_idx = 1:length(lambda_range_extended)
    current_lambda = lambda_range_extended(lambda_idx);
    current_beta2 = beta2_range_extended(lambda_idx);
    current_f = c/current_lambda;
    %Calculating PSD of noise
    Sase = h*current_f*nsp*(linear_G-1);
    Snli = ((2/3)^3*gamma^2*(Leff*10^3)*Ptx^3)*(log(pi^2*abs(current_beta2)*(Leff*10^3)*(Nch*Rs)^2)/(pi*abs(current_beta2)*Rs^3));
    %Calculating SNR
    SNR = Ptx/((Sase+Snli)*ns*Brx);
    
    % QPSK BER Simulation and Theoretical Calculation in AWGN
    
    N = 1e7;                  % Number of bits
    M = 64;
    SNR_dB = 10*log10(SNR);              % Signal-to-Noise Ratio in dB
    
    % Generate random bits
    bits = randi([0 M-1], N, 1);
    
    % Modulate bits using QPSK
    modulated_signal = qammod(bits, M);
    
    % Pass the modulated signal through an AWGN channel
    received_signal = awgn(modulated_signal, SNR_dB, 'measured');
    
    % Demodulate the received signal
    received_bits = qamdemod(received_signal, M);
    
    % Calculate the simulated BER
    SER_simulated(lambda_idx) = sum(bits ~= received_bits)/N;
    
    % Calculate the theoretical BER for QPSK in AWGN
    BER_theoretical(lambda_idx) = (2/log2(M))*(1-1/sqrt(M))*erfc(sqrt(1.5*SNR/(M-1)));
    SER_theoretical(lambda_idx) = 1-(1-BER_theoretical(lambda_idx))^log2(M);

    fprintf('SNR = %d dB, Simulated BER = %e, Theoretical BER = %e\n', ...
            SNR_dB, SER_simulated(lambda_idx), SER_theoretical(lambda_idx));
end
plot(lambda_range, beta2_range, 'r', lambda_range_extended, beta2_range_extended, 'b')
fprintf('\nSimulation complete. Plotting BER curve.\n');
figure;
semilogy(lambda_range_extended, SER_simulated, 'bo-', 'DisplayName', 'Simulated 64-QAM');
hold on;
semilogy(lambda_range_extended, SER_theoretical, 'rx-', 'DisplayName', 'Theoretical 64-QAM');
hold on;
semilogy(lambda_range_extended, BER_theoretical, 'm*-', 'DisplayName', 'Theoretical 64-QAM BER');
xlabel('wavelength');
ylabel('Sample Error Rate (SER)');
title('64-QAM Modulation over AWGN Channel');
grid on;
ylim([1e-8 1]);
legend('show');
% Display results