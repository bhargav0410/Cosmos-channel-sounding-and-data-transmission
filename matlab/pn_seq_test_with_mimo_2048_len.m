clc;
clear all;
close all;

NFFT = 2^11;% Used for displying frequency domain stats
cp_len = 72;
Fs = 20e6;% Used during simulation of up and down conversion
F = linspace(-1,1,NFFT-1)*Fs/2;% Used for displying frequency domain stats
rate = Fs/2; % Sampling rate used for simulation
for mod_size = [2]
    num_ofdm_syms = ceil(1e3/((log2(mod_size))*(NFFT-1)));
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

for rx_chan = [1]
for tx_chan = [1]

    h = comm.MIMOChannel;
    h.SampleRate = Fs;
    h.SpatialCorrelation = false; % Independent channels
    h.NumTransmitAntennas = tx_chan;
    h.NumReceiveAntennas = rx_chan;
%     h.FadingDistribution = 'Rician';
%     h.PathDelays = [0,1,2,3]*10e-8;
%     h.NormalizePathGains = true;
%     h.AveragePathGains = [0,-0.9,-4.9,-8];
    h.FadingDistribution = 'Rayleigh';
    h.PathDelays = [0,1,2,3,4,5,6,7,8]*2*1e-7;
    h.NormalizePathGains = true;
     h.AveragePathGains = transpose(-randi([0 29],9,1)-rand(9,1));%[0,-0.9,-4.9,-8,-10,-13,-9,-3,-25];
%     h.Visualization = 'Frequency response';

    ofdm_mod = comm.OFDMModulator;
    ofdm_mod.FFTLength = NFFT;
    ofdm_mod.NumGuardBandCarriers = [0;0];
    ofdm_mod.InsertDCNull = 1;
    ofdm_mod.CyclicPrefixLength = cp_len;
    ofdm_mod.NumSymbols = num_ofdm_syms;

    ofdm_mod_p = comm.OFDMModulator;
    ofdm_mod_p.FFTLength = NFFT;
    ofdm_mod_p.NumGuardBandCarriers = [0;0];
    ofdm_mod_p.InsertDCNull = 1;
    ofdm_mod_p.CyclicPrefixLength = cp_len;
    ofdm_mod_p.NumSymbols = 1;

    ofdm_demod = comm.OFDMDemodulator;
    ofdm_demod.FFTLength = NFFT;
    ofdm_demod.NumGuardBandCarriers = [0;0];
    ofdm_demod.RemoveDCCarrier = 1;
    ofdm_demod.CyclicPrefixLength = cp_len;
    ofdm_demod.NumSymbols = num_ofdm_syms;

%     DIR = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
%         num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(rx_chan),'ants_',h.FadingDistribution,'_single_pn_seq_ber.txt'];
%     file1 = fopen(DIR,'w');
%     DIR = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
%         num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(rx_chan),'ants_',h.FadingDistribution,'_single_pn_seq_mse.txt'];
%     file2 = fopen(DIR,'w');
%     DIR = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
%         num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(rx_chan),'ants_',h.FadingDistribution,'_single_pn_seq_mse_td.txt'];
%     file3 = fopen(DIR,'w');
    
    mod_train = qammod([0:mod_size-1], mod_size, 'UnitAveragePower', true);
    
for sim_tx = [1 2 4 8 16]
%     corr_val = [];
%     corr_thres = [];
%     corr_val = abs(corr(PN_Seq(:,1:sim_tx))) - eye(sim_tx);
%     corr_thres = sum(corr_val);
    
    for sim_rx = [1]
    
in_len = randi([0 mod_size - 1],(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
input = reshape(qammod(in_len,mod_size, 'UnitAveragePower', true),NFFT - 1,ofdm_mod.NumSymbols);


        % Testing sequence for pure channel coefficients
        for i = 1:tx_chan
            pn_test = zeros(num_samps,tx_chan);
            for j = 1:tx_chan
                if i == j
                    pn_test(:,i) = [1;zeros(num_samps-1,1)];%PN_Seq(:,1);
                end
            end
            mse_test(:,i) = h(pn_test);
        end
        % for i = 1:ofdm_mod.CyclicPrefixLength
        %     C_mat(i,:) = circshift(PN_Seq(:,1),i-1);
        % end
        % H_test = C_mat*mse_test./num_samps;
        for t = 1:tx_chan
        %     for i = 1:length(H_test(:,t))
        %         if abs(H_test(i,t)) <= 0.01;%sqrt((num_samps - 1)/(num_samps^3))
        %             H_test(i) = 0;
        %         end
        %     end
            H_fft_t = fftshift(fft(([mse_test(:,t);zeros(num_samps - ofdm_mod.CyclicPrefixLength+1,1)]),ofdm_demod.FFTLength));
            H_fft_t(ofdm_mod.FFTLength/2) = [];
            H_compare(:,t) = H_fft_t;
        end

iter = 1;
for snr = 0:50
    rng(snr);
    curr_pts = 1:sim_tx:NFFT-1;
    in_seq = zeros((ofdm_mod.FFTLength + ofdm_mod.CyclicPrefixLength)*((h.NumTransmitAntennas/sim_tx)), ...
        h.NumTransmitAntennas);
    
    for total_tx = 0:sim_tx:h.NumTransmitAntennas-1
        pil = [];
        for tx = 1:sim_tx
%             pil(:,tx) = [PN_Seq(num_samps - cp_len-1:num_samps-1,tx);PN_Seq(:,tx)];
%             pil(:,tx) = [circshift(PN_Seq(num_samps - cp_len-1:num_samps-1,1),(tx-1)*floor(num_samps/sim_tx)); ...
%                 circshift(PN_Seq(:,1),(tx-1)*floor(num_samps/sim_tx))];
            pil(:,tx) = [zeros(length(num_samps - cp_len-1:num_samps-1),1); ...
                circshift(PN_Seq(:,1),(tx-1)*floor(num_samps/sim_tx))];
            in_seq((num_samps + cp_len+1)*(total_tx)/sim_tx + 1:(num_samps + cp_len + 1)*((total_tx)/sim_tx + 1),tx + total_tx) = [pil(:,tx)];
%             in_seq(:,tx + total_tx - 1) = in_seq(:,tx + total_tx - 1)/max(abs(in_seq(:,tx + total_tx - 1)));
        end
    end
    %in_seq = [in_seq;repmat(ofdm_out,1,h.NumTransmitAntennas)];
    
    
    
            
                
    y = awgn(h(in_seq),snr,'measured');
% y_test = h([1; zeros(NFFT - 1-1,1)]);

% plot(abs(fft(xcorr(PN_Seq,PN_Seq) + xcorr(PN_Seq,PN_Seq_2) + xcorr(PN_Seq,PN_Seq_3),1024))/(NFFT - 1+1));
    H_final = zeros(NFFT-1,tx_chan,rx_chan);
    H_interp = [];
    for rx = 1:rx_chan
        for total_tx = 0:sim_tx:tx_chan-1
            
            if sim_tx == 1
                for i = 1:ofdm_mod.CyclicPrefixLength
                    C_mat(i,:) = circshift(PN_Seq(:,1),i-1);
                end
%                 for i = 1:num_samps
%                     C_mat_total(i,:) = circshift(PN_Seq(:,1),i-1);
%                 end
%                 del_sub = C_mat_total'*C_mat_total./num_samps;
                
%                 del_sub =  ((C_mat*C_mat')./(num_samps)) - eye(cp_len);
%                 del_sub = (1/num_samps)*[eye(cp_len),zeros(cp_len,num_samps-cp_len)]*(C_mat_total*C_mat_total')*[eye(cp_len),zeros(cp_len,num_samps-cp_len)]';

                H_est(:,total_tx+1) = C_mat*y((num_samps + cp_len+1)*(total_tx)/sim_tx + cp_len+2: ...
                    (num_samps + cp_len + 1)*((total_tx)/sim_tx + 1),rx)/num_samps;
%                 thres = 0.1;
%                 H_est(:,total_tx+1) = inv(del_sub(1:cp_len,1:cp_len))*H_est(:,total_tx+1);
                for i = 1:length(H_est(:,1+total_tx))
                    if abs(H_est(i,1+total_tx)) <= 0%sqrt((10^(-snr/10))/(num_samps) + (num_samps - 1)/(num_samps^3))
                        H_est(i,1+total_tx) = 0;
                    end
                end
                H_fft = (fft(([H_est(:,1+total_tx);zeros(num_samps - ofdm_mod.CyclicPrefixLength+1,1)]),ofdm_demod.FFTLength));
                H_fft(ofdm_mod.FFTLength/2+1) = [];
                H_interp(:,1+total_tx,rx) = H_fft;
            else
                for tx = 1:sim_tx
                    for i = 1:ofdm_mod.CyclicPrefixLength
                        C_mat(i,:) = circshift(PN_Seq(:,1),(tx-1)*floor(num_samps/sim_tx)+(i-1));
                    end
%                     for i = 1:num_samps
%                         C_mat_total(i,:) = circshift(PN_Seq(:,1),i-1);
%                     end
%                     del_sub =  ((C_mat*C_mat')./(num_samps));
%                     del_sub = (1/num_samps)*[eye(cp_len),zeros(cp_len,num_samps-cp_len)]*(C_mat_total*C_mat_total')*[eye(cp_len),zeros(cp_len,num_samps-cp_len)]';


                    H_est(:,tx+total_tx) = C_mat*y((num_samps + cp_len+1)*(total_tx)/sim_tx + cp_len+2: ...
                        (num_samps + cp_len+1)*((total_tx)/sim_tx + 1),rx)/num_samps;
%                     thres = 0.1;
%                     H_est = inv(del_sub)*H_est;
                    for i = 1:length(H_est(:,tx+total_tx))
                        if abs(H_est(i,tx+total_tx)) < 0%sqrt((10^(-snr/10))/(num_samps) + (num_samps - 1)/(num_samps^3))
                            H_est(i,tx+total_tx) = 0;
                        end
                    end
                    H_fft = (fft(([H_est(:,tx+total_tx);zeros(num_samps - ofdm_mod.CyclicPrefixLength+1,1)]),ofdm_demod.FFTLength));
                    H_fft(ofdm_mod.FFTLength/2+1) = [];
                    H_interp(:,tx+total_tx,rx) = H_fft;

                end
%                 H_interp(ofdm_mod.FFTLength/2 + 1,:) = [];
            end
%             H_final(:,total_tx+(1:sim_tx),rx) = H_interp;
        end
    end
    H_final = H_interp;
    
    % Channel estimation part done
    % Performing MRT
    pilot = ones(NFFT-1,1);
%     for j = 1:NFFT-1
%         H_temp(:,:) = H_final(j,:,:);
%         H_inv(j,:,:) = transpose(H_temp(:,:))'*inv(transpose(H_temp(:,:))*transpose(H_temp(:,:))');
%     end
    input_ = [];
    ofdm_out = [];
    pilot_out = [];
    for tx = 1:tx_chan
        for syms = 1:ofdm_mod.NumSymbols
            input_(:,syms) = conj(H_final(:,tx)).*input(:,syms)./sum(abs(H_final(:,:)).^2,2);
%             input_(:,syms) = H_inv(:,1,tx).*input(:,syms);
        end
%         pilot_out(:,tx) = [PN_Seq(NFFT - cp_len-1:NFFT-1,tx);PN_Seq(:,tx)];
        pilot_out(:,tx) = ofdm_mod_p(H_final(:,tx).*pilot./sum(abs(H_final(:,:)).^2,2));
        ofdm_out(:,tx) = [pilot_out(:,tx);ofdm_mod(input_)];
%         ofdm_out(:,tx) = ofdm_out(:,tx)/max(abs(ofdm_out(:,tx)));
    end
    y = awgn(h([ofdm_out]),snr,'measured');
    
%     y = y./max(abs(y));
        
    err = ones(NFFT-1,1);
    H_fft = [];
    ofdm_demod_out = [];
    rec_sym = [];
    for rx = 1:rx_chan
%         H_fft(:,rx) = fftshift(fft(y(cp_len+1:NFFT+cp_len,rx),NFFT));
        
        ofdm_demod_out(:,rx) = reshape(ofdm_demod(y(NFFT+cp_len+1:end,rx)),(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
    end
%     H_fft(ofdm_mod.FFTLength/2+1,:) = [];
    for i = 1:ofdm_mod.NumSymbols
%         rec_sym(:,i) = sum(ofdm_demod_out((i-1)*(NFFT-1)+1:(NFFT-1)*i,:).*conj(H_fft),2)./sum(abs(H_fft).^2,2);
        rec_sym(:,i) = ofdm_demod_out((i-1)*(NFFT-1)+1:(NFFT-1)*i,1);%./err;
%         for j = 1:NFFT-1
%             [~,min_arg] = min(rec_sym(j,i) - mod_train);
%             err(j) = rec_sym(j,i)/mod_train(min_arg);
%         end
    end
    
    rec_for_qam_demod = [];
    rec_final = [];
    out_bin_val = [];
    in_bin_val = [];
    in_gray_val = [];
    out_gray_val = [];
    rec_for_qam_demod = reshape(rec_sym,(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
    % figure; plot(real(rec_for_qam_demod),imag(rec_for_qam_demod), '.');% axis([-1.5 1.5 -1.5 1.5]);

    rec_final = qamdemod(rec_for_qam_demod,mod_size, 'UnitAveragePower', true);
    ser(iter) = sum(abs(rec_final - in_len) > 0)/length(in_len);
    out_bin_val = fliplr(de2bi(rec_final));
    in_bin_val = fliplr(de2bi(in_len));
    for i = 1:length(rec_final)
        in_gray_val(i,1) = in_bin_val(i,1);
        out_gray_val(i,1) = out_bin_val(i,1);
       for j = 2:size(out_bin_val,2)
           in_gray_val(i,j) = mod((in_bin_val(i,j-1)+in_bin_val(i,j)),2);
           out_gray_val(i,j) = mod((out_bin_val(i,j-1)+out_bin_val(i,j)),2);
       end
    end
    ber(iter) = sum(sum(out_gray_val ~= in_gray_val))/(length(rec_final)*log2(mod_size));
    % Mean Square Error of channel estimation
    mse(iter) = mean(mean(abs(H_final - H_compare).^2));
    mse_td(iter) = mean(mean((abs(ifft(H_final,num_samps) - mse_test).^2)));
    
%     fprintf(file1,'%f %f\n',[snr ser(iter)]);
%     fprintf(file2,'%f %f\n',[snr mse(iter)]);
%     fprintf(file3,'%f %f\n',[snr mse_td(iter)]);
    
    iter = iter+1;
end


    figure(mod_size); semilogy(0:50,ber, '-'); hold on;
    title([num2str(mod_size),'-QAM, ',num2str(Fs/1e6),' MHz BW, ',num2str(NFFT),' sub-carriers, and ',num2str(tx_chan),' ants']);
    grid on;
    xlabel('SNR (dB)');
    ylabel('BER');
    legend('1 ant','2 ants','4 ants','8 ants','16 ants','32 ants','64 ants');
    
    figure(mod_size+1); semilogy(0:50,mse_td, '-'); hold on;
    title([num2str(mod_size),'-QAM, ',num2str(Fs/1e6),' MHz BW, ',num2str(NFFT),' sub-carriers, and ',num2str(tx_chan),' ants']);
    grid on;
    xlabel('SNR (dB)');
    ylabel('MSE');
    legend('1 ant','2 ants','4 ants','8 ants','16 ants','32 ants','64 ants');
    end
    
end
   fclose('all'); 
end
end
end
% legend('1','2','4','8','16','32');