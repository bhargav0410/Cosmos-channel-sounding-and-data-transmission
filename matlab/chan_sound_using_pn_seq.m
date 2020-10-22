clc;
clear all;
close all;

NFFT = 2^10;% Used for displying frequency domain stats
Fs = 200e6;% Used during simulation of up and down conversion
F = linspace(-1,1,NFFT)*Fs/2;% Used for displying frequency domain stats
rate = Fs/2; % Sampling rate used for simulation

% -------------------- PN sequence generation ---------------------------
% The code below generates PN Sequence using Linear Feedback Shift Register
% method.

num_samps = 1023; % PN Sequence length
init = randi(2,1,10) - 1; % Initial values for the register 
while sum(init) == 0
    init = randi(2,1,10) - 1;
end
poly = [10 7 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = 2*pnSequence() - 1;


Y = PN_Seq;
% -----------------------------------------------------------------------
% Following code simulates the following processes
% Up-conversion - Transmission through channel - Reception - Down-converison
% The up-conversion and down-conversion processes are simulated using the
% transmit and received raised cosine filters respectively.

% Raised Cosine Transmit Filter

Nsym = 6;           % Filter span in symbol durations
beta = 0.1;         % Roll-off factor
sampsPerSym = 8;    % Upsampling factor

rctFilt = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Normal', ...
  'RolloffFactor',          beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'OutputSamplesPerSymbol', sampsPerSym);

y_rct = rctFilt([Y]);
% y_rct = reshape(y_rct,length(y_rct)*sampsPerSym,1);
t_rct = (0:length(y_rct)*sampsPerSym-1)/Fs;
% figure(fig_iter);
% fig_iter = fig_iter+1;
% stem(t_data(1:10),real(Y(1:10)), 'k'); hold on;
% plot(t_rct(1:10*sampsPerSym),real(y_rct(1+(Nsym/(2*rate))*Fs:10*sampsPerSym+(Nsym/(2*rate))*Fs)), '-');

% Passing signal through channel
chan = ricianchan(1/Fs,0,[rand(1,1),rand(1,1),rand(1,1),rand(1,1),rand(1,1)],[0,10e-9*randi([1 10]),20e-9*randi([1 10]),30e-9*randi([1 10]),40e-9*randi([1 10])],[1,1,1,1,1]);
y_filt = filter(chan,y_rct);

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

% Adding noise
y_rcr = awgn(y_rcr,-10);

%------------ Correlating with PN sequence at receiver ------------
iter = 0;
thres_iter = 0; thres_flag = 1;
while iter <= length(Y)-length(PN_Seq)
    temp = transpose(y_rcr(1+iter:iter+length(PN_Seq)))*PN_Seq;
    out(iter+1) = temp/(PN_Seq'*PN_Seq);
    if abs(out(iter+1)) > 0.3 && thres_flag == 1
        thres_iter = iter+1;
        thres_flag = 0;
    end
    iter = iter+1;
% pause(0.1);
end
figure;
plot(abs(out(Nsym:end)));
title('Complex correlation at receiver');

% -----------------------------------------------------------------------

PN_fft = fft(PN_Seq,NFFT); Y_fft = fft(y_rcr,NFFT);
H_t = pinv(PN_Seq)*y_rcr;

figure; plot(abs(H_t));
Y = y_rcr*H_t;