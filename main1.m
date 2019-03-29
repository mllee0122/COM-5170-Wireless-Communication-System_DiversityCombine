% Wireless Communication Network
% Homework4_20181224
% 107064522

clc;
clear;

%%%%%%%%%%%%%%%%%%%%
%%%%% Rayleigh %%%%%
%%%%%%%%%%%%%%%%%%%%

%% Transmission Signal & Channel Condition
no_samples = 3e5;   % Number of bits to be transmitted
SNR = [1 3 5 7 9];  % Es/N0 in dB scale
L = [1 2 3 4];   % Diversity branches
I_data = 2*(rand(5, no_samples) > 0.5)-1;   % Inphase bipolar seq.
Q_data = 2*(rand(5, no_samples) > 0.5)-1;   % Quadrature bipolar seq.
Tx = I_data+j*Q_data;   % QPSK signal
P_tx = sum(abs(Tx).^2)/length(Tx);  % QPSK signal power

% Noise
for i = 1:length(SNR)
    SNR_linear  = 10^(-SNR(i)/10);  % SNR in linear scale
    AWGN(i, :) = (1/sqrt(2))*(sqrt(SNR_linear)* normrnd(0,1,1,no_samples)+j*sqrt(SNR_linear)* normrnd(0,1,1,no_samples));   % Var = 0.5
end

% Channel fading gain
for i = 1:4
    G(:,:,i) = (1/sqrt(2))*(normrnd(0, 1, 5, no_samples) + j*normrnd(0, 1, 5, no_samples));   
    max_G(:,:,i) = max(G, [], 3);   % SC-> Choose the maximum gain
end

%% SC
% Receiver
for i = 1:length(L)
    Rx_SC = Tx.*abs(max_G(:,:,i))+AWGN; % Y = g*X+N
    err = sum(max(-(0.5*(I_data.*sign(real(Rx_SC))+1)-1),...
        -(0.5*(Q_data.*sign(imag(Rx_SC))+1)-1)), 2);
    Pe_SC(:, i) = err/no_samples;
end

% plot
figure;
semilogy(SNR, Pe_SC);
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');	
ylabel('P_e');
title ('Selective Combining with Rayleigh Signals');

%% MRC
% Receiver
Rx_MRC = 0;
for i = 1:length(L)
    Rx_MRC = (Tx.*G(:,:,i)+AWGN).*conj(G(:,:,i)) + Rx_MRC; % Y = g*X+N
    err_MRC = sum(max(-(0.5*(I_data.*sign(real(Rx_MRC))+1)-1),...
        -(0.5*(Q_data.*sign(imag(Rx_MRC))+1)-1)), 2);
    Pe_MRC(:, i) = err_MRC/no_samples;
end

% plot
figure;
semilogy(SNR, Pe_MRC);
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');	
ylabel('P_e');
title ('Maximal Ratio Combining with Rayleigh Signals');

%% EGC
% Receiver
Rx_EGC = 0;
for i = 1:length(L)
    Rx_EGC = (Tx.*G(:,:,i)+AWGN).*exp(-j*angle(G(:,:,i))) + Rx_EGC; % Y = g*X+N
    err_EGC = sum(max(-(0.5*(I_data.*sign(real(Rx_EGC))+1)-1),...
        -(0.5*(Q_data.*sign(imag(Rx_EGC))+1)-1)), 2);
    Pe_EGC(:, i) = err_EGC/no_samples;
end

% plot
figure;
semilogy(SNR, Pe_EGC);
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');	
ylabel('P_e');
title ('Equal Gain Combining with Rayleigh Signals');

%% DC
Tx_dc = sign(normrnd(0,1,5,no_samples));
Rx_dc = 0;

sum_G = zeros(5, no_samples, length(L));
sum_G(:,:,1) = G(:,:,1);
for i = 1:(length(L)-1)
    sum_G(:,:,i+1) = sum_G(:,:,i) + G(:,:,i+1);
end

for i = 1:length(L)
    for k = 1:length(SNR)
        SNR_linear  = 10^(-SNR(k)/10);  % SNR in linear scale
        AWGN(k, :) = (1/sqrt(2))*(sqrt(SNR_linear)* normrnd(0,1,1,no_samples)+j*sqrt(SNR_linear)* normrnd(0,1,1,no_samples));   % Var = 0.5
    end
    Rx_dc = (Tx_dc.*G(:,:,i)+AWGN) + Rx_dc; % Y = g*X+N
    Rx_DC = Rx_dc.*exp(-j*angle(sum_G(:,:,i)));
    err_DC = sum(-(((Tx_dc.*sign(real(Rx_DC)) + 1)/2)-1), 2);
    Pe_DC(:,i) = err_DC/no_samples;  % Error Probability
end

% plot
figure;
semilogy(SNR, Pe_DC);
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');	
ylabel('P_e');
title ('Direct Combining with Rayleigh Signals');