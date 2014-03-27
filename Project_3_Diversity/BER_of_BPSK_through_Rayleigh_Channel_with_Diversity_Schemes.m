% Elie Weintraub 
% ECE408 - Wireless Comms.
% Project 3: BER of BPSK through Rayleigh using various Diversity Schemes 
% 3/25/14

clc,clear all,close all
%% Define various simulation parameters

%Modulation Scheme
M = 4;                                % alphabet size
k = log2(M);                          % bits per symbol
PSK_Mod = comm.PSKModulator(M,0);    
PSK_Demod = comm.PSKDemodulator(M,0);

%EbNo,EsNo,SNR
EbNo = 0:2:20;                   % bit energy/noise pwr spect density (dB)
EsNo = EbNo + 10*log10(k);       % symb energy/noise pwr spect density (dB)
TsymTsmp = 1;                    % Ratio of symbol period to sample period
snr = EsNo - 10*log10(TsymTsmp); % SNR (dB)

%Rayleigh Channels
Ts = 1e-4;                       % Sample time of Rayleigh chan. input sig.
fd = 100;                        % Maximum Doppler shift
r_chan1 = rayleighchan(Ts,fd);    % Rayleigh channel object
r_chan1.StorePathGains = 1;
r_chan2 = rayleighchan(Ts,fd);    % Rayleigh channel object
r_chan2.StorePathGains = 1;
r_chan3 = rayleighchan(Ts,fd);    % Rayleigh channel object
r_chan3.StorePathGains = 1;
r_chan4 = rayleighchan(Ts,fd);    % Rayleigh channel object
r_chan4.StorePathGains = 1;

%% PART 1 - BPSK through flat-fading rayleigh channel w/o diversity
%% Simulate BPSK through flat-fading rayleigh channel w/o diversity

%Define message length
N = 1e6;               % messsage length in symbols
%Generate random message
x = randi([0 M-1],N,1);
%Modulate message
modulated = step(PSK_Mod,x);
%Transmit through flat-fading Rayleigh Chan. & AWGN chan. Then demodulate.
filtered=filter(r_chan1,modulated);
transmitted = zeros(length(modulated),length(snr));
demodulated = zeros(length(modulated),length(snr));
for i=1:length(snr)
   transmitted(:,i) = awgn(filtered,snr(i),'measured')./r_chan1.PathGains;
   demodulated(:,i) = step(PSK_Demod,transmitted(:,i)); 
end
%Compute BER
[~,ber1] = biterr(demodulated,x); 

%% PART 2 - BPSK through flat-fading rayleigh channel w/ MRRC  (2 Rx)
%% Simulate BPSK through flat-fading rayleigh channel w/ MRRC (2 Rx)
%Define message length
N = 1e6;               % messsage length in symbols
%Generate random message
x = randi([0 M-1],N,1);
%Modulate message
modulated = step(PSK_Mod,x);
%Transmit through flat-fading Rayleigh channels 
filtered1 = filter(r_chan1,modulated);
filtered2 = filter(r_chan2,modulated);
h = [r_chan1.PathGains,r_chan2.PathGains]; %path gains for each branch
filtered = [filtered1,filtered2]; %rayleigh filtered signal for each branch
% Transmit through AWGN chan. Then MRRC combine and demodulate.
combined = zeros(length(modulated),length(snr));
demodulated = zeros(length(modulated),length(snr));
for i=1:length(snr)
   n_filtered = awgn(filtered,snr(i),'measured'); 
   combined(:,i) = sum(conj(h).*n_filtered,2)./sum(h.*conj(h),2);
   demodulated(:,i) = step(PSK_Demod,combined(:,i)); 
end
%Compute BER
[~,ber2] = biterr(demodulated,x);

%% PART 3 - BPSK through flat-fading rayleigh channel w/ MRRC  (4 Rx)
%% Simulate BPSK through flat-fading rayleigh channel w/ MRRC (4 Rx)
%Define message length
N = 1e6;               % messsage length in symbols
%Generate random message
x = randi([0 M-1],N,1);
%Modulate message
modulated = step(PSK_Mod,x);
%Transmit through flat-fading Rayleigh channels 
filtered1 = filter(r_chan1,modulated);
filtered2 = filter(r_chan2,modulated);
filtered3 = filter(r_chan3,modulated);
filtered4 = filter(r_chan4,modulated);
h = [r_chan1.PathGains,r_chan2.PathGains,...
     r_chan3.PathGains,r_chan4.PathGains]; %path gains for each branch
filtered = [filtered1,filtered2,...
            filtered3,filtered4]; %rayleigh filtered signal for each branch
% Transmit through AWGN chan. Then MRRC combine and demodulate.
combined = zeros(length(modulated),length(snr));
demodulated = zeros(length(modulated),length(snr));
for i=1:length(snr)
   n_filtered = awgn(filtered,snr(i),'measured'); 
   combined(:,i) = sum(conj(h).*n_filtered,2)./sum(h.*conj(h),2);
   demodulated(:,i) = step(PSK_Demod,combined(:,i)); 
end
%Compute BER
[~,ber3] = biterr(demodulated,x);

%% PART 4 - BPSK through rayleigh channel w/ Alamouti (2 Tx, 1 Rx)
%% Simulate BPSK through rayleigh channel w/ Alamouti (2 Tx, 1 Rx)
%Define message length
N = 1e6;               % messsage length in symbols
%Generate random message
x = randi([0 M-1],N,1);
%Modulate message
modulated = step(PSK_Mod,x);
%Encode using Alamouti space-time block coding
encoded=zeros(N,2);
s1=modulated(1:2:end); s2=modulated(2:2:end);
encoded(1:2:end,:)=sqrt(0.5)*[s1,s2];
encoded(2:2:end,:)=sqrt(0.5)*[-conj(s2),conj(s1)];
%Transmit through Rayleigh channels with constant responses over 2 symbols 
h = sqrt(0.5)*kron((randn(N/2,2) + 1j*randn(N/2,2)),[1;1]); 
filtered = h.*encoded;
% Transmit through AWGN chan. and demodulate.
combined = zeros(length(modulated),length(snr));
decoded = zeros(length(modulated),length(snr));
demodulated = zeros(length(modulated),length(snr));
for i=1:length(snr)
   %Add AWGN  
   n_filtered = awgn(filtered,snr(i),'measured'); 
   %Combine the branches "in the air"
   combined(:,i) = sum(n_filtered,2);
   %Alamouti decode
   y1 = combined(1:2:end,i); y2 = combined(2:2:end,i);
   y = [kron(y1,[1;1]), kron(conj(y2),[1;1])]; 
   h_trans = zeros(N,2); 
   h1=h(1:2:end,1); h2=h(1:2:end,2);
   h_trans(1:2:end,:) = [conj(h1), h2];
   h_trans(2:2:end,:) = [conj(h2), -h1];
   decoded(:,i) = sum(h_trans.*y,2)./sum(h_trans.*conj(h_trans),2);
   %Demodulate
   demodulated(:,i) = step(PSK_Demod,decoded(:,i)); 
end
%Compute BER
[~,ber4] = biterr(demodulated,x);

%% PART 5 - BPSK through rayleigh channel w/ Alamouti (2 Tx, 2 Rx)
%% Simulate BPSK through rayleigh channel w/ Alamouti (2 Tx, 2 Rx)
%Define message length
N = 1e6;               % messsage length in symbols
%Generate random message
x = randi([0 M-1],N,1);
%Modulate message
modulated = step(PSK_Mod,x);
%Encode using Alamouti space-time block coding
encoded=zeros(N,4);
s1=modulated(1:2:end); s2=modulated(2:2:end);
encoded(1:2:end,:)=sqrt(0.5)*[s1,s2, s1,s2];
encoded(2:2:end,:)=sqrt(0.5)*[-conj(s2),conj(s1), -conj(s2),conj(s1)];
%Transmit through Rayleigh channels with constant responses over 2 symbols 
h = sqrt(0.5)*kron((randn(N/2,4) + 1j*randn(N/2,4)),[1;1]); 
filtered = h.*encoded;
% Transmit through AWGN chan. and demodulate.
combined = zeros(length(modulated),2,length(snr));
decoded = zeros(length(modulated),length(snr));
demodulated = zeros(length(modulated),length(snr));
for i=1:length(snr)
   %Add AWGN  
   n_filtered = awgn(filtered,snr(i),'measured'); 
   %Combine the branches "in the air"
   combined(:,:,i) = [sum(n_filtered(:,1:2),2), sum(n_filtered(:,3:4),2)];
   %Alamouti decode
   y1_1 = combined(1:2:end,1,i); y2_1 = combined(2:2:end,1,i);
   y1_2 = combined(1:2:end,2,i); y2_2 = combined(2:2:end,2,i);
   y = [kron(y1_1,[1;1]),kron(conj(y2_1),[1;1]),...
        kron(y1_2,[1;1]),kron(conj(y2_2),[1;1])]; 
   h_trans = zeros(N,4); 
   h1=h(1:2:end,1); h2=h(1:2:end,2); h3=h(1:2:end,3); h4=h(1:2:end,4);
   h_trans(1:2:end,:) = [conj(h1),h2,  conj(h3),h4];
   h_trans(2:2:end,:) = [conj(h2),-h1, conj(h4),-h3];
   decoded(:,i) = sum(h_trans.*y,2)./sum(h_trans.*conj(h_trans),2);
   %Demodulate
   demodulated(:,i) = step(PSK_Demod,decoded(:,i)); 
end
%Compute BER
[~,ber5] = biterr(demodulated,x);

%% PART 6 - Plot simulation results
%% Plot results
figure('Name','BER for BPSK through  a Rayleigh Channel');
semilogy(EbNo,ber1,'-o',EbNo,ber2,'-v',EbNo,ber3,'-s',...
         EbNo,ber4,'-d',EbNo,ber5,'-^');
title('BER for BPSK through  a Rayleigh Channel'); grid on;
xlabel('Eb/No (dB)'); ylabel('Bit Error Rate');
legend('no diversity (1 Tx, 1 Rx)','MRRC (1 Tx, 2 Rx)',...
       'MRRC (1 Tx, 4 Rx)', 'Alamouti (2 Tx, 1 Rx)',...
       'Alamouti (2 Tx, 2 Rx)');
   