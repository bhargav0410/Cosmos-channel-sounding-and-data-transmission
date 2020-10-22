clc;
clear all;
close all;

NFFT = 2^11;% Used for displying frequency domain stats
cp_len = 72;
Fs = 200e6;% Used during simulation of up and down conversion
F = linspace(-1,1,NFFT-1)*Fs/2;% Used for displying frequency domain stats
rate = Fs/2; % Sampling rate used for simulation

num_samps = 2047;

init = randi(2,1,11) - 1; % Initial values for the register 
while sum(init) == 0
    init = randi(2,1,11) - 1;
end

PN_Seq = [];
poly = [11 2 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 4 2 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 5 3 2 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 6 2 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 6 5 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 6 5 2 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 6 5 4 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 6 5 4 3 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 7 3 2 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 7 6 4 2 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 7 4 2 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 7 6 5 4 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 7 4 3 2 1 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 7 5 3 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 7 6 5 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

poly = [11 7 6 5 4 2 0]; % Seed polynoimal
pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',num_samps,'InitialConditions',init);
PN_Seq = [PN_Seq,2*pnSequence() - 1];

% PN_Seq = repmat(transpose(exp(-1i*((pi*1*(0:num_samps-1).*((0:num_samps-1) + rem(num_samps,2)))./(num_samps)))),1,16);

tx_chan = 16;rx_chan = 1;
h = comm.MIMOChannel;
h.SampleRate = Fs;
h.SpatialCorrelation = false; % Independent channels
h.NumTransmitAntennas = tx_chan;
h.NumReceiveAntennas = rx_chan;
h.FadingDistribution = 'Rician';
%     h.PathDelays = [0,1,2,3]*10e-8;
%     h.NormalizePathGains = true;
%     h.AveragePathGains = [0,-0.9,-4.9,-8];
% h.FadingDistribution = 'Rayleigh';
h.PathDelays = [0,1,2,3,4,5,6,7,8]*2*1e-7;
h.NormalizePathGains = true;
 h.AveragePathGains = transpose(-randi([0 29],9,1)-rand(9,1));%[0,-0.9,-4.9,-8,-10,-13,-9,-3,-25];
%     h.Visualization = 'Impulse response';

test_imp = [1;zeros(num_samps-1,1)];
test = zeros(num_samps*tx_chan,tx_chan);
for i = 1:tx_chan
    test(num_samps*(i-1) + (1:num_samps),i) = test_imp;
end
imp_resp = h(test(1:end,:));
figure; plot(abs(imp_resp(:,1)));

sim_tx = 16;

wh_mat = hadamard(sim_tx);
wh_mat_rep = repelem(wh_mat,(num_samps+1)/sim_tx,1);

in_seq = zeros(num_samps*(tx_chan/sim_tx),tx_chan);
for i = 0:sim_tx:tx_chan-1
    for tx = 1:sim_tx
        in_seq(num_samps*(i/sim_tx) + (1:num_samps),i + tx) = PN_Seq(1:num_samps,tx);%.*wh_mat_rep(1:num_samps,tx);
    end
end
in_seq_rep = repmat(in_seq,sim_tx,1);

for i = 0:sim_tx:tx_chan-1
    for tx = 1:sim_tx
        for j = 1:sim_tx
            in_seq_rep(num_samps*sim_tx*(i/sim_tx) + (j-1)*num_samps + (1:num_samps),i + tx) ...
                = in_seq_rep(num_samps*sim_tx*(i/sim_tx) + (j-1)*num_samps + (1:num_samps),i + tx)*wh_mat(j,tx);
        end
    end
end

y = awgn(h(in_seq_rep),10);

for j = 1:sim_tx
    y(num_samps*(j-1) + (1:num_samps)) = y(num_samps*(j-1) + (1:num_samps))*wh_mat(j,2);
end
y_avg = mean(reshape(y,num_samps,length(y)/num_samps),2);
for i = 0:num_samps-1
    C_mat(:,i+1) = circshift(PN_Seq(:,1),i);
end
ccorr_out = ifft(fft(y_avg,2048).*conj(fft(PN_Seq(:,2),2048)),2048);

figure; plot(abs(xcorr(y_avg(:,1),PN_Seq(1:num_samps,2)))/num_samps);
figure; plot(abs(ccorr_out)/num_samps);



