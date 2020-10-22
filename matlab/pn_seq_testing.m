clc;
clear all;
close all;

NFFT = 1024;% Used for displying frequency domain stats
Fs = 10e6;% Used during simulation of up and down conversion
F = linspace(-1,1,NFFT-1)*Fs/2;% Used for displying frequency domain stats
rate = Fs/2; % Sampling rate used for simulation

% -------------------- PN sequence generation ---------------------------
% The code below generates PN Sequence using Linear Feedback Shift Register
% method.


num_samps = NFFT-1; % PN Sequence length
init = randi(2,1,10) - 1; % Initial values for the register 
while sum(init) == 0
    init = randi(2,1,10) - 1;
end

PN_Seq = [];
poly = [10 7 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [10 3 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [10 4 3 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [10 5 2 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

h = comm.MIMOChannel;
h.SampleRate = Fs;
h.SpatialCorrelation = false; % Independent channels
h.NumTransmitAntennas = 1;
h.NumReceiveAntennas = 1;
h.FadingDistribution = 'Rician';
h.PathDelays = [0,1,2,3]*10e-8;
h.NormalizePathGains = true;
h.AveragePathGains = [0,-0.9,-4.9,-8];

ofdm_mod = comm.OFDMModulator;
ofdm_mod.FFTLength = num_samps + 1;
ofdm_mod.NumGuardBandCarriers = [0;0];
ofdm_mod.InsertDCNull = 1;
ofdm_mod.CyclicPrefixLength = 72;
ofdm_mod.NumSymbols = 500;

ofdm_mod_p = comm.OFDMModulator;
ofdm_mod_p.FFTLength = num_samps + 1;
ofdm_mod_p.NumGuardBandCarriers = [0;0];
ofdm_mod_p.InsertDCNull = 1;
ofdm_mod_p.CyclicPrefixLength = 72;
ofdm_mod_p.NumSymbols = 1;

ofdm_demod = comm.OFDMDemodulator;
ofdm_demod.FFTLength = num_samps + 1;
ofdm_demod.NumGuardBandCarriers = [0;0];
ofdm_demod.RemoveDCCarrier = 1;
ofdm_demod.CyclicPrefixLength = 72;
ofdm_demod.NumSymbols = 500;

for type = [1]
for mod_size = [2 4 16 64 256]
in_len = randi([0 mod_size - 1],(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
input = reshape(qammod(in_len,mod_size, 'UnitAveragePower', true),num_samps,ofdm_mod.NumSymbols);

ofdm_out = ofdm_mod(input);
iter = 1;
for snr = 0:50
    rng(snr);
if type == 1 | type == 3
        in_seq(:,1) = [PN_Seq(num_samps - ofdm_demod.CyclicPrefixLength + 1:end,1); ...
            PN_Seq(:,1);ofdm_out];
    in_seq = in_seq/max(abs(in_seq));
%     y = awgn(h(in_seq),snr,'measured');
end
if type == 2
    pil = ofdm_mod_p(ones(ofdm_mod_p.FFTLength - 1,1));
    in_seq = [pil;ofdm_out];
    in_seq = in_seq/max(abs(in_seq));
    %y = awgn(h(in_seq),snr,'measured');
end

if type == 3
    y = [h(in_seq)];
    y = [y(1:length([PN_Seq(num_samps - ofdm_demod.CyclicPrefixLength + 1:end);PN_Seq])); ...
        awgn(y(length([PN_Seq(num_samps - ofdm_demod.CyclicPrefixLength + 1:end);PN_Seq])+1:end),snr,'measured')];
else
    y = awgn(h(in_seq),snr,'measured');
end
% y_test = h([1; zeros(num_samps-1,1)]);

% plot(abs(fft(xcorr(PN_Seq,PN_Seq) + xcorr(PN_Seq,PN_Seq_2) + xcorr(PN_Seq,PN_Seq_3),1024))/(num_samps+1));
if type == 1 | type == 3
    for i = 1:ofdm_mod.CyclicPrefixLength
        C_mat(i,:) = circshift(PN_Seq(:,1),i-1);
    end

    H_est = C_mat*y(ofdm_mod.CyclicPrefixLength + 1:ofdm_mod.CyclicPrefixLength + num_samps)/num_samps;
    thres = 0.1;
    for i = 1:length(H_est)
        if abs(H_est(i)) < thres
            H_est(i) = 0;
        end
    end
    H_fft = fftshift(fft(([H_est(:);zeros(num_samps - ofdm_mod.CyclicPrefixLength+1,1)]),ofdm_demod.FFTLength));
    H_fft(ofdm_mod.FFTLength/2+1) = [];
end
if type == 2
    H_fft = fftshift(fft(y(ofdm_mod.CyclicPrefixLength+1:ofdm_mod.CyclicPrefixLength + ofdm_mod.FFTLength),ofdm_demod.FFTLength));
    H_fft(ofdm_mod.FFTLength/2+1) = [];
end
% figure; plot(abs(H_est));
% figure; plot(F + 5.4e9,abs(H_fft));
% figure; plot(F + 5.4e9,angle(H_fft));

if type == 1 | type == 3
    ofdm_demod_out = ofdm_demod(y(ofdm_mod.FFTLength + ofdm_mod.CyclicPrefixLength:end));
end
if type == 2
    ofdm_demod_out = ofdm_demod(y(ofdm_mod.FFTLength + ofdm_mod.CyclicPrefixLength+1:end));
end
for i = 1:ofdm_mod.NumSymbols
    rec_sym(:,i) = ofdm_demod_out(:,i).*conj(H_fft)./abs(H_fft).^2;
end
rec_for_qam_demod = reshape(rec_sym,(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
% figure; plot(real(rec_for_qam_demod),imag(rec_for_qam_demod), '.');% axis([-1.5 1.5 -1.5 1.5]);

rec_final = qamdemod(rec_for_qam_demod,mod_size, 'UnitAveragePower', true);
ser(iter) = sum(abs(rec_final - in_len) > 0)/length(in_len);
per(iter) = 0;
for num = 1:ofdm_mod.NumSymbols
    if sum(abs(rec_final((num - 1)*num_samps + 1:num*num_samps) - in_len((num - 1)*num_samps + 1:num*num_samps))) > 0
        per(iter) = per(iter) + 1;
    end
end
per(iter) = per(iter)/ofdm_mod.NumSymbols;
iter = iter+1;

end
if type == 1
    figure(10); semilogy(ser, '-'); hold on;
    figure(11); semilogy(per, '-'); hold on;
end
if type == 2
    figure(10); semilogy(ser, '-.'); hold on;
    figure(11); semilogy(per, '-.'); hold on;
end
if type == 3
    figure(10); semilogy(ser, '--'); hold on;
    figure(11); semilogy(per, '--'); hold on;
end
end
end

