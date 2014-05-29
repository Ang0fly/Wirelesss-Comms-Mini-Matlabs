%% MIMO through Rayleigh Channel
%  Mark Bryk and the Traubs
%  ECE 408 - Wireless Comms
%  5/1/14
%  MIMO
% 

%%
tic
clc, clear, close all

%% Define various simulation parameters

M = 2;                                % alphabet size
k = log2(M);                          % bits per symbol

EbNo = 0:2.5:25;                 % bit energy/noise pwr spect density (dB)
EsNo = EbNo + 10*log10(k);       % symb energy/noise pwr spect density (dB)
TsymTsmp = 1;                    % Ratio of symbol period to sample period
snr = EsNo - 10*log10(TsymTsmp); % SNR (dB)

%% Part 1
% BPSK through flat-fading rayleigh channel w/ ZF Equalizer

Nt = 3;                     % Transmitters
Nr = 3;                     % Receivers
N = 1e5;                    % messsage length per stream in symbols

x = randi([0 M-1],Nt*N,1);

% Modulate message
modulated = reshape(qammod(x,M,0,'gray'),Nt,1,N);  % Nt x 1 x N

% Transmit through flat-fading Rayleigh Chan. & AWGN chan.
H = sqrt(1/2)*(randn(Nr,Nt,N)+1j*randn(Nr,Nt,N));
Hx=sum(H.*repmat(permute(modulated,[2 1 3]),[Nr,1,1]),2);
n = sqrt(1/2)*(randn(Nr,1,N)+1j*randn(Nr,1,N));
transmitted = zeros(Nr,length(snr),N);

for i=1:length(snr)
    transmitted(:,i,:) = Hx+10^(-snr(i)/20)*n;
end

% Form ZF Equalization Matrix W
W = zeros(Nt,Nr,N);
for i=1:N
    W(:,:,i) = (H(:,:,i)'*H(:,:,i))^-1*H(:,:,i)';
end

% Equalize using W
y = zeros(Nt,length(snr),N);
for i=1:length(snr)
    y(:,i,:) = sum(W.*repmat(permute(transmitted(:,i,:),[2 1 3]),[Nt,1,1]),2);
end

% Demodulate
demodulated = qamdemod(reshape(permute(y,[1 3 2]),Nt*N,[]),M,0,'gray'); 

% Compute BER
[~,ber1] = biterr(demodulated,x); 

%% PART 2 
% BPSK through flat-fading rayleigh channel w/ MMSE Equalizer

Nt = 3;
Nr = 3;
N = 1e5;

x = randi([0 M-1],Nt*N,1);

% Modulate message
modulated = reshape(qammod(x,M,0,'gray'),Nt,1,N);  % Nt x 1 x N

% Transmit through flat-fading Rayleigh Chan. & AWGN chan.
H = sqrt(1/2)*(randn(Nr,Nt,N)+1j*randn(Nr,Nt,N));
Hx=sum(H.*repmat(permute(modulated,[2 1 3]),[Nr,1,1]),2);
n = sqrt(1/2)*(randn(Nr,1,N)+1j*randn(Nr,1,N));
transmitted = zeros(Nr,length(snr),N);

for i=1:length(snr)
    transmitted(:,i,:) = Hx+10^(-snr(i)/20)*n;
end

% Form MMSE Equalization Matrix W and Equalize
W = zeros(Nt,Nr,N);
y = zeros(Nt,length(snr),N);
for i=1:length(snr)
    N0 = 10^(-snr(i)/10);
    for j=1:N
        W(:,:,j) = (H(:,:,j)'*H(:,:,j)+N0*eye(Nt))^-1*H(:,:,j)';
    end
    y(:,i,:) = sum(W.*repmat(permute(transmitted(:,i,:),[2 1 3]),[Nt,1,1]),2);
end

% Demodulate
demodulated = qamdemod(reshape(permute(y,[1 3 2]),Nt*N,[]),M,0,'gray'); 

% Compute BER
[~,ber2] = biterr(demodulated,x); 

%% Plot results

figure('Name','BER for MIMO BPSK through a Rayleigh Channel');
semilogy(EbNo,ber1,'-o',EbNo,ber2,'-x');
title('BER for MIMO BPSK through  a Rayleigh Channel');
xlabel('Eb/No (dB)'); ylabel('Bit Error Rate');
legend('ZF Equalizer','MMSE Equalizer');

toc