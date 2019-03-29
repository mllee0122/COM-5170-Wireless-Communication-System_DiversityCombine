% Wireless Communication Network
% Homework4_20181224
% 107064522

clc;
clear;
%close all;

%%%%%%%%%%%%%%%%%%
%%%%% Ricean %%%%%
%%%%%%%%%%%%%%%%%%

%% Transmission Signal & Channel Condition
no_samples = 3e5;   % Number of bits to be transmitted
SNR = [1 3 5 7 9];  % Es/N0 in dB scale
L = [1 2 3 4];   % Diversity branches
I_data = 2*(rand(5, no_samples) > 0.5)-1;   % Inphase bipolar seq.
Q_data = 2*(rand(5, no_samples) > 0.5)-1;   % Quadrature bipolar seq.
Tx = I_data+j*Q_data;   % QPSK signal
P_tx = sum(abs(Tx).^2)/length(Tx);  % QPSK signal power

% Ricean fading
K = 1;
var = (2*(K+1))^-1;
mu_I=sqrt(K/(K+1))*sin(pi/4);
mu_Q=sqrt(K/(K+1))*cos(pi/4);


% Noise
for i = 1:length(SNR)
    SNR_linear  = 10^(-SNR(i)/10);  % SNR in linear scale
    AWGN(i, :) = (1/sqrt(2))*(sqrt(SNR_linear)* normrnd(0,1,1,no_samples)+j*sqrt(SNR_linear)* normrnd(0,1,1,no_samples));   % Var = 0.5
end

% Channel fading gain
for i = 1:length(L)
    G(:,:,i) = normrnd(mu_Q, sqrt(var), 5, no_samples) + j*normrnd(mu_I, sqrt(var), 5, no_samples);
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
ylim([10^-4 1]);
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');	
ylabel('P_e');
title ('Selective Combining with Ricean Signals');

%% MRC
% Receiver
%%%%Rx_MRC = zeros(5, 3e5, length(L));
Rx_MRC(:,:,1) = (Tx.*G(:,:,1)+AWGN).*conj(G(:,:,1)) ; % Y = g*X+N
for i = 1:length(L)-1
    for k = 1:length(SNR)
        SNR_linear  = 10^(-SNR(k)/10);  % SNR in linear scale
        AWGN(k, :) = (1/sqrt(2))*(sqrt(SNR_linear)* normrnd(0,1,1,no_samples)+j*sqrt(SNR_linear)* normrnd(0,1,1,no_samples));   % Var = 0.5
    end
    Rx_MRC(:,:,i+1) = (Tx.*G(:,:,i+1)+AWGN).*conj(G(:,:,i+1)) + Rx_MRC(:,:,i); % Y = g*X+N
end

for i = 1:length(L)
    for k = 1:length(SNR)
        Rx_re = sign(real(Rx_MRC(:,:,i)));
        Rx_im = sign(imag(Rx_MRC(:,:,i)));
        err_I(k, i) = (no_samples-sum(I_data(k,:) == Rx_re(k,:)))/no_samples;
        err_Q(k, i) = (no_samples-sum(Q_data(k,:) == Rx_im(k,:)))/no_samples;
        Pe_MRC(k, i) = mean([err_I(k,i) err_Q(k,i)]);
    end
end

% plot
figure;
semilogy(SNR, Pe_MRC);
ylim([10^-4 1]);
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');	
ylabel('P_e');
title ('Maximal Ratio Combining with Ricean Signals');

%% EGC
% Receiver
Rx_EGC(:,:,1) = (Tx.*G(:,:,1)+AWGN).*exp(-j*angle(G(:,:,1))) ; % Y = g*X+N
for i = 1:length(L)-1
    for k = 1:length(SNR)
        SNR_linear  = 10^(-SNR(k)/10);  % SNR in linear scale
        AWGN(k, :) = (1/sqrt(2))*(sqrt(SNR_linear)* normrnd(0,1,1,no_samples)+j*sqrt(SNR_linear)* normrnd(0,1,1,no_samples));   % Var = 0.5
    end
    Rx_EGC(:,:,i+1) = (Tx.*G(:,:,i+1)+AWGN).*exp(-j*angle(G(:,:,i+1))) + Rx_EGC(:,:,i); % Y = g*X+N
end

for i = 1:length(L)
    for k = 1:length(SNR)
        Rx_re = sign(real(Rx_EGC(:,:,i)));
        Rx_im = sign(imag(Rx_EGC(:,:,i)));
        err_I(k, i) = (no_samples-sum(I_data(k,:) == Rx_re(k,:)))/no_samples;
        err_Q(k, i) = (no_samples-sum(Q_data(k,:) == Rx_im(k,:)))/no_samples;
        Pe_EGC(k, i) = mean([err_I(k,i) err_Q(k,i)]);
    end
end



% plot
figure;
semilogy(SNR, Pe_EGC);
ylim([10^-4 1]);
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');	
ylabel('P_e');
title ('Equal Gain Combining with Ricean Signals');

%% DC
% Noise
for i = 1:length(SNR)
	SNR_linear  = 10^(-SNR(i)/10);  % SNR in linear scale
    AWGN(i,:) = normrnd(0,sqrt(SNR_linear/2),1,no_samples) +...
        j*normrnd(0,sqrt(SNR_linear/2),1,no_samples);
end

% Channel fading gain
for i = 1:length(L)
    G(:,:,i) = normrnd(mu_Q, sqrt(var), 5, no_samples) + j*normrnd(mu_I, sqrt(var), 5, no_samples);
end

Tx = I_data + j*Q_data;
Rx_DC(:,:,1) = (AWGN + Tx.*G(:,:,1));
G1(:,:,1) = (G(:,:,1));
for i = 1:length(L)-1
    for k = 1:length(SNR)
        SNR_linear  = 10^(-SNR(k)/10);  % SNR in linear scale
        AWGN(k,:) = normrnd(0,sqrt(SNR_linear/2),1,no_samples) +...
            j*normrnd(0,sqrt(SNR_linear/2),1,no_samples);
    end
    Rx_DC(:,:,i+1) = Rx_DC(:,:,i)+(AWGN + Tx.*G(:,:,i+1));
    G1(:,:,i+1)=G1(:,:,i)+(G(:,:,i+1));
end
for i=1:length(L)
    Rx_DC(:,:,i) = Rx_DC(:,:,i).*exp(-j*angle(G1(:,:,i)));
end

% Demodulation
for i = 1:length(L)
	for k = 1:length(SNR)
		Rx_I=sign(real(Rx_DC(:,:,i)));
		Rx_Q=sign(imag(Rx_DC(:,:,i)));
		err_I(k,i)=(no_samples-sum(I_data(k,:)==Rx_I(k,:)))/no_samples;
		err_Q(k,i)=(no_samples-sum(Q_data(k,:)==Rx_Q(k,:)))/no_samples;
        Pe_DC(k,i)=mean([err_I(k,i) err_Q(k,i)]);
	end
end

% plot
figure;
semilogy(SNR, Pe_DC);
ylim([10^-4 1]);
legend('L = 1','L = 2','L = 3','L = 4');
xlabel('SNR (dB)');	
ylabel('P_e');
title ('Direct Combining with Ricean Signals');