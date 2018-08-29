clc;
clear all;
close all;

% Parameters for simulation

NFFT = 2^10;% Used for displying frequency domain stats
Fs = 200e6;% Used during simulation of up and down conversion
F = linspace(-1,1,NFFT)*Fs/2;% Used for displying frequency domain stats
rate = Fs/2; % Sampling rate used for simulation

% -------------------- PN sequence generation ---------------------------
% The code below generates PN Sequence using Linear Feedback Shift Register
% method.

num_samps = 1023; % PN Sequence length
init = randi(2,1,10) - 1; % Initial values for the register 
while init == 0
    init = randi(2,1,10) - 1;
end
poly = [10 7 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = 2*pnSequence() - 1;

% PLotting auto-correlation of PN Sequence
figure;
plot(xcorr(PN_Seq)/num_samps);
title('Auto-correlation of PN sequence');
tx_buf = 0;
% -----------------------------------------------------------------------

% while 1
% --------------- QPSK Data generation ----------------------------------
% The code below generates QPSK sequence and then spreads the sequence
% according to the required spread factor.
% Here, spread_factor = chip rate/sampling rate, where, sampling rate is
% the rate of actual transmission and chip rate is the rate of PN sequence.
% Spread factor is the number of PN sequence bits to be XORed with the input
% QAM sample.

spread_factor = num_samps; % Decides spread factor
qpsk_nsamps = ceil(num_samps/spread_factor);
% Generating 10 symbols and spreading them according to spread factor
for i = 1:100
    if i == 1
        QPSKOutput(i) = 1;
        QPSKOutput_spread = reshape(transpose(repmat(QPSKOutput(i),1,spread_factor)),num_samps,1);
    else
        X_Input = randi([0 3],qpsk_nsamps,1);
        QPSKModulatorObject = comm.QPSKModulator('BitInput',false);
        QPSKOutput(i) = step(QPSKModulatorObject,X_Input);
        QPSKOutput_spread = reshape(transpose(repmat(QPSKOutput(i),1,spread_factor)),num_samps,1);
    end
    % Multiplication of PN and QPSK
    Y((i-1)*num_samps+1:i*num_samps) = PN_Seq .* QPSKOutput_spread;
end

t_data = (0:length(Y)-1)/rate;

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

y_rct = rctFilt([transpose(Y)]);
% y_rct = reshape(y_rct,length(y_rct)*sampsPerSym,1);
t_rct = (0:length(y_rct)*sampsPerSym-1)/Fs;
% figure(fig_iter);
% fig_iter = fig_iter+1;
% stem(t_data(1:10),real(Y(1:10)), 'k'); hold on;
% plot(t_rct(1:10*sampsPerSym),real(y_rct(1+(Nsym/(2*rate))*Fs:10*sampsPerSym+(Nsym/(2*rate))*Fs)), '-');

% Passing signal through channel
chan = ricianchan(1/Fs,0,[rand(1,1),rand(1,1),rand(1,1),rand(1,1),rand(1,1)],[0,10e-9*randi([1 10]),20e-9*randi([1 10]),30e-9*randi([1 10]),40e-9*randi([1 10])],[1,1,1,1,1]);

%%%% If MIMO channel is to be used %%%%
h = comm.MIMOChannel;
h.SampleRate = rate;
h.SpatialCorrelation = false; % Independent channels
h.NumTransmitAntennas = 1;
h.NumReceiveAntennas = 1;
h.FadingDistribution = 'Rayleigh';
% h.DirectPathDopplerShift = Doppshift;
% h.MaximumDopplerShift = Doppshift;
% h.KFactor = 5;
h.PathDelays = [0,1,2,3]*10e-9; 
h.NormalizePathGains = false;
h.AveragePathGains = [0,-0.9,-4.9,-8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y_filt = filter(chan,y_rct);
y_filt = step(h,y_rct);

% Plotting frequency characteristics of the up-converted signal passed
% through the channel.
figure;
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

% Adding noise
y_rcr = awgn(y_rcr,6);

% -----------------------------------------------------------------------

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
% 
% pause(1);
% end

% -----------------------------------------------------------------------



% Finding received signal

chan_est = transpose(out(1:1023));

thres = 0.2;

for i = 1:length(chan_est)
    if abs(chan_est(i)) < thres
        chan_est(i) = 0;
    end
end

chan_est_thres = chan_est(find(abs(chan_est)>=thres));
chan_est_thres_delay = find(abs(chan_est)>=thres);

figure; plot(abs(chan_est));

y_reshaped = reshape(y_rcr(length(PN_Seq)+1:end),spread_factor,99);
to_be_decoded = y_rcr(spread_factor+1:end);
for i = 1:size(y_reshaped,2)
    Y_f = fft(y_reshaped(:,i))./fft(chan_est,size(y_reshaped,1));
    y_copy = zeros(spread_factor,length(chan_est_thres_delay));
    for j = 1:length(chan_est_thres_delay)
        y_copy(1:end-chan_est_thres_delay(j)+1,j) = conj(chan_est_thres(j))*y_reshaped(chan_est_thres_delay(j):end,i);
    end
    Y_t = sum(y_copy,2)./sum(abs(chan_est_thres).^2);
    Y_final(i) = (PN_Seq'*Y_t)/length(PN_Seq);
end
scatterplot(QPSKOutput);
scatterplot(Y_final);
