%  Mark Bryk, Elie Weintraub, Hillel Weintraub
%  ECE 408 - Wireless Comms
%  5/1/14
% OFDM
clear

%% Constants
M = 2; N = 100;
nFFT = 64;
nTap = 10;
% BPS = 52; %bits/sym
% SNR = 0:15;
% k = log2(M);
% EsNo = SNR; EbNo = EsNo - 10*log10(k);
EbN0dB = [0:2:40]; % bit to noise ratio
EsN0dB = EbN0dB + 10*log10(64/80);

%% Data

% Message:
x = randi([0 M-1], 1, nFFT*N);
x_modulated = qammod(x, M, [], 'gray');
ipMod = reshape(x_modulated, nFFT, N).'; % grouping into multiple symbolsa

% Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
%xF = [zeros(N,6) ipMod(:,[1:BPS/2]) zeros(N,1) ipMod(:,[BPS/2+1:BPS]) zeros(N,5)] ;

% Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
xt = ifft(fftshift(ipMod.')).';

% Appending cylic prefix
xt = [xt(:,[49:64]) xt];

ht = 1/sqrt(2)*1/sqrt(nTap)*(randn(N,nTap) + j*randn(N,nTap));
hF = fftshift(fft(ht,64,2));

ht=rayleighchan(1/10000,0, [0 2e-4 3e-4 4e-4], [0 -6 -9 -12]);
ht.StorePathGains = 1;

% convolution of each symbol with the random channel
for jj = 1:N
  %xht(jj,:) = conv(ht(jj,:),xt(jj,:));
  xht(jj,:) = filter(ht,xt(jj,:));
  hFreq(jj,:) = ht.PathGains;
end
xt = xht;

% Concatenating multiple symbols to form a long vector
xt = reshape(xt.',1,N*(80+nTap-1));

%% Calaculations
for ii = 1:length(EbN0dB)

   % Gaussian noise of unit variance, 0 mean
   nt = 1/sqrt(2)*[randn(1,N*(80+nTap-1)) + j*randn(1,N*(80+nTap-1))];

   % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
   yt = sqrt(80/64)*xt + 10^(-EsN0dB(ii)/20)*nt;

   % Receiver
   yt = reshape(yt.',80+nTap-1,N).'; % formatting the received vector into symbols
   yt = yt(:,[17:80]); % removing cyclic prefix

   N0 = 10^(-EbN0dB(ii)/10); % noise power
   % converting to frequency domain
   yF = fftshift(fft(yt.')).';
   yF_MMSE = yF./(hF + N0);
   yF_ZF = yF./(hF);
   
   
   %yMod_MMSE = yF_MMSE(:,[6+[1:BPS/2] 7+[BPS/2+1:BPS] ]); 
   %yMod_ZF = yF_ZF(:,[6+[1:BPS/2] 7+[BPS/2+1:BPS] ]); 
    yMod_MMSE = yF_MMSE;
    yMod_ZF = yF_ZF;
	yPredicted_ZF = qamdemod(yMod_ZF, M, [], 'gray');
    yPredicted_ZF = reshape(yPredicted_ZF.',nFFT*N,1).';
    
    yPredicted_MMSE = qamdemod(yMod_MMSE, M, [], 'gray');
    yPredicted_MMSE = reshape(yPredicted_MMSE.',nFFT*N,1).';
      
    [n ber_sim_ZF(ii)] = biterr(x, yPredicted_ZF);
    [n ber_sim_MMSE(ii)] = biterr(x, yPredicted_MMSE);

end


%% Figure
figure
semilogy(EbN0dB, ber_sim_ZF, EbN0dB, ber_sim_MMSE);
title('OFDM: Compare ZF and MMSE Equalizers');
xlabel('Eb/No, dB'); ylabel('BER');
legend('ZF','MMSE');
