%Constants and variables
h = 6.6261*10^-34;
c = 3*10^8;
lambda = 1.55e-6;
f = c/lambda;
nsp = 1.58;
G = 16;
ns_range = 0:2:20;
gamma = 1.2*10^-3;
L = 80;
alpha = 0.0461;
Ptx_dbm = -3;
beta2 = -2.17*10^-26;
Nch = 100;
Rs = 5*10^10;
Brx = 5*10^10;

%Linearizing and calculating Leff
Leff = (1-exp(-alpha*L))/alpha;
linear_G = 10^(G/10);
Ptx = 10^-3*10^(Ptx_dbm/10);

SER_simulated = zeros(size(ns_range));
SER_theoretical = zeros(size(ns_range));
BER_theoretical = zeros(size(ns_range));

for ns_idx = 1:length(ns_range)
    current_ns = ns_range(ns_idx);
    %Calculating PSD of noise
    Sase = h*f*nsp*(linear_G-1);
    Snli = ((2/3)^3*gamma^2*(Leff*10^3)*Ptx^3)*(log(pi^2*abs(beta2)*(Leff*10^3)*(Nch*Rs)^2)/(pi*abs(beta2)*Rs^3));
    %Calculating SNR
    SNR = Ptx/((Sase+Snli)*current_ns*Brx);
    
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
    SER_simulated(ns_idx) = sum(bits ~= received_bits)/N;
    
    % Calculate the theoretical BER for QPSK in AWGN
    BER_theoretical(ns_idx) = (2/log2(M))*(1-1/sqrt(M))*erfc(sqrt(1.5*SNR/(M-1)));
    SER_theoretical(ns_idx) = 1-(1-BER_theoretical(ns_idx))^log2(M);

    fprintf('SNR = %d dB, Simulated BER = %e, Theoretical BER = %e\n', ...
            SNR_dB, SER_simulated(ns_idx), SER_theoretical(ns_idx));
end
fprintf('\nSimulation complete. Plotting BER curve.\n');
figure;
semilogy(ns_range, SER_simulated, 'bo-', 'DisplayName', 'Simulated 64-QAM');
hold on;
semilogy(ns_range, SER_theoretical, 'rx-', 'DisplayName', 'Theoretical 64-QAM');
xlabel('number of spans');
ylabel('Symbol Error Rate (SER)');
title('64-QAM Modulation over AWGN Channel');
grid on;
ylim([1e-8 1]);
legend('show');
% Display results