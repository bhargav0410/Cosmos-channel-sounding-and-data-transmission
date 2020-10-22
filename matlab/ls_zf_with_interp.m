clc;
clear all;
close all;

for rx_chan = [64]
for tx_chan = [4]
NFFT = 2048;% Used for displying frequency domain stats
cp_len = 72;
Fs = 200e6;% Used during simulation of up and down conversion
F = linspace(-1,1,NFFT-1)*Fs/2;% Used for displying frequency domain stats
rate = Fs/2; % Sampling rate used for simulation
num_ofdm_syms = ceil(1e3/(NFFT-1));


    h = comm.MIMOChannel;
    h.SampleRate = Fs;
    h.SpatialCorrelation = false; % Independent channels
    h.NumTransmitAntennas = tx_chan;
    h.NumReceiveAntennas = rx_chan;
    h.FadingDistribution = 'Rayleigh';
    h.PathDelays = [0,1,2,3]*10e-8;
    h.NormalizePathGains = true;
    h.AveragePathGains = [0,-0.9,-4.9,-8];

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

for mod_size = [2 4 16 64 256 1024]
%     DIR = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
%         num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(rx_chan),'ants_',h.FadingDistribution,'_dist_mrc.txt'];
%     file1 = fopen(DIR,'w');
    
    mod_train = qammod([0:mod_size-1], mod_size, 'UnitAveragePower', true);
    
for sim_tx = [1 2 4]
    for sim_rx = [1]
    
in_len = randi([0 mod_size - 1],(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
input = reshape(qammod(in_len,mod_size, 'UnitAveragePower', true),NFFT - 1,ofdm_mod.NumSymbols);

iter = 1;
for snr = 0:50
    rng(snr);
    curr_pts = 1:sim_tx:NFFT-1;
    in_seq = zeros((ofdm_mod.FFTLength + ofdm_mod.CyclicPrefixLength)*((h.NumTransmitAntennas/sim_tx)), ...
        h.NumTransmitAntennas);
    
    for total_tx = 1:sim_tx:h.NumTransmitAntennas
        pil = [];
        for tx = 1:sim_tx
            pil_idx = zeros(ofdm_mod.FFTLength - 1,1);
            for idx = curr_pts
                if idx + tx - 1 <= NFFT - 1
                    pil_idx(idx + tx - 1) = 1;
                end
            end
            pil(:,tx) = [ofdm_mod_p(pil_idx)];
            in_seq((NFFT + cp_len)*(total_tx - 1)/sim_tx + 1:(NFFT + cp_len)*((total_tx-1)/sim_tx + 1),tx + total_tx - 1) = [pil(:,tx)];
%             in_seq(:,tx + total_tx - 1) = in_seq(:,tx + total_tx - 1)/max(abs(in_seq(:,tx + total_tx - 1)));
        end
    end
    %in_seq = [in_seq;repmat(ofdm_out,1,h.NumTransmitAntennas)];

    %y = awgn(h(in_seq),snr,'measured');

    y = awgn(h(in_seq),snr,'measured');
% y_test = h([1; zeros(NFFT - 1-1,1)]);

% plot(abs(fft(xcorr(PN_Seq,PN_Seq) + xcorr(PN_Seq,PN_Seq_2) + xcorr(PN_Seq,PN_Seq_3),1024))/(NFFT - 1+1));
    H_final = zeros(NFFT-1,tx_chan,rx_chan);
    for rx = 1:rx_chan
        for total_tx = 1:sim_tx:tx_chan
            H_fft = fftshift(fft(y(((total_tx - 1)/sim_tx)*(cp_len + NFFT) + (cp_len+1:cp_len + NFFT),rx),NFFT));
            H_fft(ofdm_mod.FFTLength/2+1) = [];
            H_interp = [];
            if sim_tx == 1

                H_interp = H_fft;
            else
                for tx = 1:sim_tx
                    curr_pts = tx:sim_tx:NFFT-1;
                    H_new_fft = H_fft(tx:sim_tx:end);

                    new_pts = 1:NFFT;
                    H_interp(:,tx) = transpose(interp1(curr_pts,H_new_fft,new_pts));
                    tail_iter = 1;
                    head_iter = tx;
                    for tail = NFFT-sim_tx+1:NFFT
                        H_interp(tail,tx) = (H_interp(NFFT-sim_tx,tx) - H_interp(NFFT-sim_tx-1,tx))*(tail_iter) + H_interp(NFFT-sim_tx - 1,tx);
                        tail_iter = tail_iter + 1;
                    end
                    for head = 1:tx-1
                        H_interp(head,tx) = (H_interp(2*tx,tx) - H_interp(tx,tx))*(head_iter) + H_interp(tx,tx);
                        head_iter = head_iter - 1;
                    end
                end
                H_interp(ofdm_mod.FFTLength/2+1,:) = [];
            end
            H_final(:,total_tx+(0:sim_tx-1),rx) = H_interp;
        end
    end
    
    % Channel estimation part done
    % Performing MRT
    pilot = ones(NFFT-1,1);
    H_fft = [];
    for i = 1:ofdm_mod.NumSymbols
        for j = 1:NFFT-1
            H_fft(:,:) = H_final(j,:,:);
            H_inv = pinv(transpose(H_fft));
            tx_sym((i-1)*(NFFT-1)+j,i) = transpose(repmat(input((i-1)*(NFFT-1)+j,i),1,tx_chan));
        end
%         for syms = 1:ofdm_mod.NumSymbols
%             input_(:,syms) = sum(conj(H_final(:,tx,1)) .* input(:,syms),2)./sum(abs(H_final(:,:,1)).^2,2);
%         end
%         pilot_out(:,tx) = ofdm_mod_p(sum(conj(H_final(:,tx,1)) .* pilot,2)./sum(abs(H_final(:,:,1)).^2,2));
        pilot_out(:,tx) = ofdm_mod_p(pilot);
        ofdm_out(:,tx) = [pilot_out(:,tx);ofdm_mod(input)];
%         ofdm_out(:,tx) = ofdm_out(:,tx)/max(abs(ofdm_out(:,tx)));
    end
    y = awgn(h([ofdm_out]),snr,'measured');
        
    err = zeros(NFFT-1,1);
    H_fft = [];
    for rx = 1:rx_chan
%         H_fft(:,rx) = fftshift(fft(y(cp_len+1:NFFT+cp_len,rx),NFFT));
        
        ofdm_demod_out(:,rx) = reshape(ofdm_demod(y(NFFT+cp_len+1:end,rx)),(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
    end
%     H_fft(ofdm_mod.FFTLength/2+1,:) = [];
    for i = 1:ofdm_mod.NumSymbols
        for j = 1:NFFT-1
            H_fft(:,:) = H_final(j,:,:);
            H_inv = pinv(transpose(H_fft));
            rec_sym((i-1)*(NFFT-1)+j,i) = mean(H_inv * transpose(ofdm_demod_out((i-1)*(NFFT-1)+j,:)));
        end
%         rec_sym(:,i) = sum(ofdm_demod_out((i-1)*(NFFT-1)+1:(NFFT-1)*i,:).*conj(H_fft),2)./sum(abs(H_fft).^2,2);
%         rec_sym(:,i) = ofdm_demod_out((i-1)*(NFFT-1)+1:(NFFT-1)*i,1) - err;
%         for j = 1:NFFT-1
%             err(j) = min(rec_sym(j,i) - mod_train);
%         end
    end

    rec_for_qam_demod = reshape(rec_sym,(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
    % figure; plot(real(rec_for_qam_demod),imag(rec_for_qam_demod), '.');% axis([-1.5 1.5 -1.5 1.5]);

    rec_final = qamdemod(rec_for_qam_demod,mod_size, 'UnitAveragePower', true);
    ser(iter) = sum(abs(rec_final - in_len) > 0)/length(in_len);
%     per(iter) = 0;
%     for num = 1:ofdm_mod.NumSymbols
%         if sum(abs(rec_final((num - 1)*(NFFT - 1) + 1:num*(NFFT - 1)) - in_len((num - 1)*(NFFT - 1) + 1:num*(NFFT - 1)))) > 0
%             per(iter) = per(iter) + 1;
%         end
%     end
%     per(iter) = per(iter)/ofdm_mod.NumSymbols;
    
%     fprintf(file1,'%f %f\n',[snr ser(iter)]);
    iter = iter+1;
end


    figure(mod_size); semilogy(0:50,ser, '-'); hold on;
    title([num2str(mod_size),'-QAM, ',num2str(Fs/1e6),' MHz BW, ',num2str(NFFT),' sub-carriers, and ',num2str(tx_chan),' ants']);
    grid on;
    xlabel('SNR (dB)');
    ylabel('SER');
    legend('1 ant','2 ants','4 ants','8 ants','16 ants','32 ants','64 ants');
%     figure(11); semilogy(per, '-.'); hold on;
    end
end
   fclose('all'); 
end
end
end
% legend('1','2','4','8','16','32');