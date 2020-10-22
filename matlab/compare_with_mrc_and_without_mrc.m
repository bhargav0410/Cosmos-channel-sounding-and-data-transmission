clc;
clear all;
close all;

for mod_size = [2 4 16 64 256]
Fs = 20e6;
NFFT = 2048;
tx_chan = 1;
fading_dist = 'Rayleigh';
DIR = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
        num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(tx_chan),'ants_',fading_dist,'_ls_with_interp_ber.txt'];
DIR2 = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
        num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(tx_chan),'ants_',fading_dist,'_single_pn_seq_ber.txt'];
% DIR3 = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
%         num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_1ants_',fading_dist,'_dist_mrc.txt'];
    
in = dlmread(DIR);
in_mrc = dlmread(DIR2);
% in_1ant = dlmread(DIR3);

figure;
for i = 1:length(in)/51
    semilogy(in(51*(i-1)+1:51*i,1), in(51*(i-1)+1:51*i,2)); hold on;
%     semilogy(in_mrc(51*(i-1)+2:51*i,1), in_mrc(51*(i-1)+2:51*i,2), '-.'); hold on;
end
for i = 1:length(in_mrc)/51
%     semilogy(in(51*(i-1)+2:51*i,1), in(51*(i-1)+2:51*i,2)); hold on;
    semilogy(in_mrc(51*(i-1)+1:51*i,1), in_mrc(51*(i-1)+1:51*i,2), '-.'); hold on;
end
% for i = 1:length(in_1ant)/51
% %     semilogy(in(51*(i-1)+2:51*i,1), in(51*(i-1)+2:51*i,2)); hold on;
%     semilogy(in_1ant(51*(i-1)+2:61*i,1), in_1ant(51*(i-1)+2:51*i,2), '--'); hold on;
% end
grid on;
end