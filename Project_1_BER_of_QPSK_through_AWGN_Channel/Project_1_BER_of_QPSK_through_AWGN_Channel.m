% Elie Weintraub 
% ECE408 - Wireless Comms.
% Project 1: BER of QPSK (4-QAM) through an AWGN Channel
% 1/25/14

clc,clear all,close all

%% PART 1 - UNENCODED QPSK
%% Define various parameters
M = 4;                           % alphabet size
k = log2(M);                     % bits per symbol

EbNo = 0:10;                     % bit energy/noise pwr spect density (dB)
EsNo = EbNo + 10*log10(k);       % symb energy/noise pwr spect density (dB)

TsymTsmp = 1;                    % Ratio of symbol period to sample period
snr = EsNo - 10*log10(TsymTsmp); % SNR (dB)

%% Simulate Unencoded QPSK
%Define simulation parameters
msg_length = 1e6;               % messsage length in symbols
%Generate random message
x = randi([0 M-1],msg_length,1);
%Modulate message
modulated = qammod(x,M,0,'gray');
%Transmit through AWGN channel
transmitted = zeros(length(modulated),length(snr));
for i=1:length(snr)
   transmitted(:,i) = awgn(modulated,snr(i),'measured');
end
%Demodulate message
demodulated = qamdemod(transmitted,M,0,'gray');
%Compute BER
[~,ber] = biterr(demodulated,x); 
%% Compute theoretical BER
theoretical_ber = berawgn(EbNo,'qam',M);
%% Plot results
figure('Name','BER for Unencoded QPSK through an AWGN Channel');
semilogy(EbNo,theoretical_ber,EbNo,ber,'x');
title('BER for Unencoded QPSK through an AWGN Channel');
xlabel('Eb/No (dB)'); ylabel('Bit Error Rate');
legend('theoretical BER','empirical BER');

%% PART 2 - ENCODED QPSK
%% Define various parameters
M = 4;                           % alphabet size
clen = 3;                        % trellis constraint length   
trel = poly2trellis(clen,[6,7]); % trellis struct for convolutional encoder
code_rate = 1/2;                 % information bits per coded bits
bps = log2(M);                   % bits per symbol
k = bps*code_rate;               % informaion bits per symbol

EbNo = 0:10;                     % bit energy/noise pwr spect density (dB)
EsNo = EbNo + 10*log10(k);       % symb energy/noise pwr spect density (dB)

TsymTsmp = 1;                    % Ratio of symbol period to sample period
snr = EsNo - 10*log10(TsymTsmp); % SNR (dB)

%% Simulate Encoded QPSK
%Define simulation parameters
msg_length = bps*1e6;            % messsage length in bits
%Generate random message
x = randi([0 1],msg_length,1);
%Encode message
encoded = convenc(x,trel);
%Convert bits to symbols
enc_symbols = bi2de(reshape(encoded,bps,[]).','left-msb');
%Modulate message
modulated = qammod(enc_symbols,M,0,'gray');
%Transmit through AWGN channel
transmitted = zeros(length(modulated),length(snr));
for i=1:length(snr)
   transmitted(:,i) = awgn(modulated,snr(i),'measured');
end
%Demodulate message
demodulated = qamdemod(transmitted,M,0,'gray');
%Convert symbols to bits
demod_bits = reshape(de2bi(demodulated,bps,'left-msb').',[],length(snr));
%Decode message
tblen = 5*clen;        % traceback length
decoded = zeros(length(x),length(snr));
for i=1:length(snr)
    decoded(:,i) = vitdec(demod_bits(:,i),trel,tblen,'trunc','hard');
end
%Compute BER
[~,ber] = biterr(decoded,x); 
%% Compute theoretical BER
dspec = distspec(trel,4); % distance spectrum of the convolutional code
theoretical_ber = bercoding(EbNo,'conv','hard',code_rate,dspec);
%% Plot results
figure('Name','BER for Encoded QPSK through an AWGN Channel');
semilogy(EbNo,theoretical_ber,EbNo,ber,'x');
title('BER for Encoded QPSK through an AWGN Channel');
xlabel('Eb/No (dB)'); ylabel('Bit Error Rate');
legend('theoretical BER','empirical BER');