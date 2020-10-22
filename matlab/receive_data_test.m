clc;
clear all;
close all;

FC = 5400e6;
FS = 10e6;
NFFT = 1024;
cp_size = 64;
range = FC + (FS/2 * linspace(-1,1,NFFT - 1));
fid = fopen('D:\Cosmos-channel-sounding-and-data-transmission\ch_dumped_0_binary','r');
rec = fread(fid, 'float');
rec_comp = rec(1:2:end) + 1i*rec(2:2:end);
plot(real(rec_comp));

fid = fopen('D:\Cosmos-channel-sounding-and-data-transmission\global_out.dat','r');
chan = fread(fid, 'float');
chan_comp = chan(1:2:end) + 1i*chan(2:2:end);
figure; plot(10*log10(abs(chan_comp)));

% fid2 = fopen('D:\Cosmos-channel-sounding-and-data-transmission\PNSeq_255_MaxLenSeq.dat','r');
% pn_seq = fread(fid2, 'float');
% pn_comp = pn_seq(1:2:end) + 1i*pn_seq(2:2:end);

% PN_Seq = [];
% init = [0 0 0 0 0 0 0 0 0 1];
% poly = [10 7 0]; % Seed polynomial
% pnSequence = comm.PNSequence('Polynomial',poly,'SamplesPerFrame',1023,'InitialConditions',init);
% PN_Seq = [PN_Seq,2*pnSequence() - 1];
% zc_seq(:,1) = exp(-1i*((pi*1*(0:255-1).*((0:255-1) + 1))./(255)));


% for i = 1:length(rec_comp)-length(zc_seq)+1
%     
%     corr_out(i) = sum(rec_comp(i:(i-1)+length(zc_seq)).*conj(circshift(zc_seq,45)))/length(zc_seq);
% end

% for i = 1:length(PN_Seq)
%     corr_out(i) = max(abs(xcorr(rec_comp, circshift(PN_Seq,i-1)))/1023);
% end
% [~,maxarg] = max(corr_out);

% figure; plot(abs(xcorr(rec_comp, circshift(PN_Seq,maxarg-1)))/1023);

fid = fopen('D:\Cosmos-channel-sounding-and-data-transmission\PN_Seq_used.dat','r');
PN_Seq = fread(fid, 'int');

for i = 1:length(rec_comp)-length(PN_Seq)+1
    peak_out(i) = sum(rec_comp(i:(i-1)+length(PN_Seq)).*PN_Seq)/255;
end
start_peak = find(abs(peak_out) >= 0.01);
signal = rec_comp(start_peak(1) + 255 :start_peak(1) + ((NFFT + cp_size)*2*1 + 254));
figure; plot(real(signal));

signal_fft_1_temp = [];
signal_fft_1 = [];
for i = 1:1
    signal_fft_1_temp(:,i) = fftshift(fft(signal((NFFT + cp_size)*(i-1)*2 + cp_size + 1:(NFFT + cp_size)*(i-1)*2 + cp_size+ NFFT)));
end
signal_fft_1 = mean(signal_fft_1_temp,2);
signal_fft_1(NFFT/2 + 1) = [];
figure; plot(range,10*log10(abs(signal_fft_1)));

signal_fft_2_temp = [];
signal_fft_2 = [];
for i = 1
    chan_fft_2_temp(:,i) = fftshift(fft(signal((NFFT + cp_size)*((i-1)*2 + 1) + cp_size + 1:(NFFT + cp_size)*((i-1)*2 + 1) + cp_size+ NFFT)));
end
for i = 1:1
    signal_fft_2_temp(:,i) = fftshift(fft(signal((NFFT + cp_size)*((i-1)*2 + 1) + cp_size + 1:(NFFT + cp_size)*((i-1)*2 + 1) + cp_size+ NFFT)));
end
% signal_fft_2 = mean(signal_fft_2_temp,2);
signal_fft_2_temp(NFFT/2 + 1,:) = [];
signal_fft_2_temp = signal_fft_2_temp./repmat(chan_comp(1*(NFFT-1) + (1:NFFT-1),1),1,1);
figure; scatter(real(signal_fft_2_temp(:,1)), imag(signal_fft_2_temp(:,1)));
% figure; plot(range,10*log10(abs(signal_fft_2)));

fclose all;