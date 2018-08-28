clc;
clear all;
close all;

NFFT = 2^12;
Fs = 200e6;
F = linspace(-1,1,NFFT)*Fs/2;
rate = Fs/2;
fig_iter = 1;

%-------------------- PN sequence generation --------------------
num_samps = 1023;
init = randi(2,1,10) - 1;
while init == 0
    init = randi(2,1,10) - 1;
end
poly = [10 7 0];
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = 2*pnSequence() - 1;

figure(fig_iter);
fig_iter = fig_iter+1;
plot(xcorr(PN_Seq)/num_samps);
title('Auto-correlation of PN sequence');
tx_buf = 0;


% while 1
% --------------- QPSK Data generation -------------------
spread_factor = num_samps;
qpsk_nsamps = ceil(num_samps/spread_factor);
for i = 1:10
X_Input = randi([0 3],qpsk_nsamps,1);
QPSKModulatorObject = comm.QPSKModulator('BitInput',false);
QPSKOutput(i) = step(QPSKModulatorObject,X_Input);
QPSKOutput_spread = reshape(transpose(repmat(QPSKOutput(i),1,spread_factor)),num_samps,1);


% Multiplication of PN and QPSK
Y((i-1)*num_samps+1:i*num_samps) = PN_Seq .* QPSKOutput_spread;
end

t_data = (0:length(Y)-1)/rate;

% Raised Cosine Transmit Filter
Nsym = 6;           % Filter span in symbol durations
beta = 0.1;         % Roll-off factor
sampsPerSym = 8;    % Upsampling factor

rctFilt = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Normal', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'OutputSamplesPerSymbol', sampsPerSym);

y_rct = rctFilt([transpose(Y)]);
% y_rct = reshape(y_rct,length(y_rct)*sampsPerSym,1);
t_rct = (0:length(y_rct)*sampsPerSym-1)/Fs;
% figure(fig_iter);
% fig_iter = fig_iter+1;
% stem(t_data(1:10),real(Y(1:10)), 'k'); hold on;
% plot(t_rct(1:10*sampsPerSym),real(y_rct(1+(Nsym/(2*rate))*Fs:10*sampsPerSym+(Nsym/(2*rate))*Fs)), '-');

% Passing signal through channel
chan = ricianchan(1/Fs,0,[rand(1,1),rand(1,1),rand(1,1),rand(1,1),rand(1,1)],[0,10e-9*randi([1 10]),20e-9*randi([1 10]),30e-9*randi([1 10]),40e-9*randi([1 10])],[1,1,1,1,1]);
y_filt = filter(chan,y_rct);

% Y = cat(2,zeros(1,10),Y);
figure(fig_iter);
fig_iter = fig_iter+1;
subplot(211);
plot(F,10*log10(abs(fftshift(fft(QPSKOutput_spread,NFFT)/NFFT)).^2));
xlabel('Frequency characteristics of QPSK input');


subplot(212);
plot(F,10*log10(abs(fftshift(fft(Y,NFFT)/NFFT)).^2));
xlabel('Frequency characteristics of signal to be transmitted');

% Raised Cosine Receive Filter

Nsym = 6;           % Filter span in symbol durations
beta = 0.1;         % Roll-off factor
sampsPerSym = 8;    % Downsampling factor

rcrFilt = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Normal', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'InputSamplesPerSymbol', sampsPerSym, ...
  'DecimationFactor', 8 );

y_rcr = rcrFilt([y_filt]);
y_rcr = awgn(y_rcr,-10);

%------------ Correlating with PN sequence at receiver ------------
iter = 0;
while iter <= length(Y)-length(PN_Seq)
    temp = transpose(y_rcr(1+iter:iter+length(PN_Seq)))*PN_Seq;
    out(iter+1) = temp/length(PN_Seq);
    iter = iter+1;
% pause(0.1);
end
figure(fig_iter);
fig_iter = fig_iter+1;
plot(abs(out));
title('Complex correlation at receiver');
% 
% pause(1);
% end