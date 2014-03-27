% Elie Weintraub 
% ECE408 - Wireless Comms.
% Project 2: BER of 4-QAM through a Rayleigh Channel
% 2/23/14

clc,clear all,close all

%% PART 1 - QAM Through a Flat-Fading Rayleigh Channel
%% Define various parameters
M = 4;                           % alphabet size
k = log2(M);                     % bits per symbol

EbNo = 0:10;                     % bit energy/noise pwr spect density (dB)
EsNo = EbNo + 10*log10(k);       % symb energy/noise pwr spect density (dB)

TsymTsmp = 1;                    % Ratio of symbol period to sample period
snr = EsNo - 10*log10(TsymTsmp); % SNR (dB)

Ts = 1e-4;                       % Sample time of Rayleigh chan. input sig.
fd = 100;                        % Maximum Doppler shift
r_chan = rayleighchan(Ts,fd);    % Rayleigh channel object
r_chan.StorePathGains = 1;
%% Simulate QAM through flat-fading rayleigh channel
%Define message length
msg_length = 1e6;               % messsage length in symbols
%Generate random message
x = randi([0 M-1],msg_length,1);
%Modulate message
modulated = qammod(x,M,0,'gray');
%Transmit through flat-fading Rayleigh Channel and AWGN channel
filtered=filter(r_chan,modulated);
transmitted = zeros(length(modulated),length(snr));
for i=1:length(snr)
   transmitted(:,i) = awgn(filtered,snr(i),'measured')./r_chan.PathGains;
end
%Demodulate message
demodulated = qamdemod(transmitted,M,0,'gray');
%Compute BER
[~,ber] = biterr(demodulated,x); 
%% Compute theoretical BER
theoretical_ber = berfading(EbNo,'qam',M,1);
%% Plot results
figure('Name','BER for QAM through a Flat-Fading Rayleigh Channel');
semilogy(EbNo,theoretical_ber,EbNo,ber,'x');
title('BER for QAM through a Flat-Fading Rayleigh Channel');
xlabel('Eb/No (dB)'); ylabel('Bit Error Rate');
legend('theoretical BER','empirical BER');

%% PART 2 - QAM Through a Frequency Selective Fading Rayleigh Channel
%% Define various parameters
M = 4;                           % alphabet size
k = log2(M);                     % bits per symbol

EbNo = 0:10;                     % bit energy/noise pwr spect density (dB)
EsNo = EbNo + 10*log10(k);       % symb energy/noise pwr spect density (dB)

TsymTsmp = 1;                    % Ratio of symbol period to sample period
snr = EsNo - 10*log10(TsymTsmp); % SNR (dB)

Ts = 1e-4;                            % Sample time 
fd = 0;                               % Max. doppler shift(0=>quasi-static)
tau=[0 2*Ts  3*Ts];                   % path delays(sec)
pdb = [0 -6 -9];                      % average path dains (dB)
r_chan = rayleighchan(Ts,fd,tau,pdb); % Rayleigh channel object
eq = lineareq(7,lms(0.01));           % linear equalizer object
trainLen = 500;                       % num symbols used to train equalizer
eq.SigConst = qammod(0:M-1,M,0,'gray');
%% Plot frequency response of rayleigh channel
figure('Name','Frequency Response of Frequency Selective Rayleigh Channel')
chan_resp = r_chan.ChannelFilter.TapGains.Values;
freqz(chan_resp);
%% Simulate QAM through frequency selective fading rayleigh channel
%Define message length
msg_length = 1e4;               % messsage length in symbols
%Generate random message
x = randi([0 M-1],msg_length,1);
%Modulate message
modulated = qammod(x,M,0,'gray');
%Transmit through freq. selective fading Rayleigh Channel and AWGN channel
filtered=filter(r_chan,modulated);
transmitted = zeros(length(modulated),length(snr));
equalized = zeros(length(modulated),length(snr));
for i=1:length(snr)
   transmitted(:,i) = awgn(filtered,snr(i),'measured');
   [sym_est,equalized(:,i)] = equalize(eq,transmitted(:,i),...
                                       modulated(1:trainLen)); % equalize
end
%Demodulate message
demodulated_uneq = qamdemod(transmitted,M,0,'gray');
demodulated_eq = qamdemod(equalized,M,0,'gray');
%Compute BER
[~,ber_uneq] = biterr(demodulated_uneq(trainLen+1:end,:),x(trainLen+1:end)); 
[~,ber_eq] = biterr(demodulated_eq(trainLen+1:end,:),x(trainLen+1:end)); 
%% Plot results
figure('Name',...
      'BER for QAM through a Frequency Selective Fading Rayleigh Channel');
semilogy(EbNo,ber_uneq,'-o',EbNo,ber_eq,'-x');
title('BER for QAM through a Frequency Selective Fading Rayleigh Channel');
xlabel('Eb/No (dB)'); ylabel('Bit Error Rate');
legend('unequalized','equalized');
