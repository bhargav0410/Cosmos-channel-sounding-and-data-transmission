clc;
clear all;
close all;

for rx_chan = [1]
for tx_chan = [64]
NFFT = 2^11;% Used for displying frequency domain stats
cp_len = 72;
Fs = 20e6;% Used during simulation of up and down conversion
F = linspace(-1,1,NFFT-1)*Fs/2;% Used for displying frequency domain stats
rate = Fs/2; % Sampling rate used for simulation
for mod_size = [2]
    num_ofdm_syms = ceil(1e3/((log2(mod_size))*(NFFT-1)));


    h = comm.MIMOChannel;
    h.SampleRate = Fs;
    h.SpatialCorrelation = false; % Independent channels
    h.NumTransmitAntennas = tx_chan;
    h.NumReceiveAntennas = rx_chan;
    h.FadingDistribution = 'Rayleigh';
    h.PathDelays = [0,1,2,3,4,5,6,7,8]*4*10e-8;
    h.NormalizePathGains = true;
    h.AveragePathGains = transpose(-randi([0 29],9,1)-rand(9,1));%[0,-0.9,-4.9,-8,-10,-13,-9,-3,-25];

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
%         num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(rx_chan),'ants_',h.FadingDistribution,'_ls_with_interp_ber.txt'];
%     file1 = fopen(DIR,'w');
%     DIR = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
%         num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(rx_chan),'ants_',h.FadingDistribution,'_ls_with_interp_mse.txt'];
%     file2 = fopen(DIR,'w');
%     DIR = ['D:\Cosmos-channel-sounding-and-data-transmission\',num2str(mod_size),'qam_', ...
%         num2str(Fs/1e6),'mhz_',num2str(NFFT),'subcar_',num2str(rx_chan),'ants_',h.FadingDistribution,'_ls_with_interp_mse_td.txt'];
%     file3 = fopen(DIR,'w');
    
    mod_train = qammod([0:mod_size-1], mod_size, 'UnitAveragePower', true);
    
for sim_tx = [1 2 4 8 16]
    for sim_rx = [1]
    
         % Testing sequence for pure channel coefficients
        for i = 1:tx_chan
            pn_test = zeros(NFFT,tx_chan);
            for j = 1:tx_chan
                if i == j
                    pn_test(:,i) = [1;zeros(NFFT-1,1)];%PN_Seq(:,1);
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
            H_fft = fftshift(fft(([mse_test(:,t);zeros(NFFT - cp_len + 1,1)]),ofdm_demod.FFTLength));
            H_fft(ofdm_mod.FFTLength/2) = [];
            H_compare(:,t) = H_fft;
        end
        
        
        
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
%     y = y./max(abs(y));
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
%                 H_interp(NFFT/2 + 1) = [];
            else
                for tx = 1:sim_tx
                    curr_pts = tx:sim_tx:NFFT-1;
                    H_new_fft = H_fft(tx:sim_tx:end);
%                     H_interp(:,tx) = repelem(H_new_fft,sim_tx);
                    new_pts = 1:NFFT;
                    H_interp(:,tx) = transpose(interp1(curr_pts,H_new_fft,new_pts,'spline'));
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
%     for j = 1:NFFT-1
%         H_inv(j,:) = transpose(H_final(j,:))'*inv(transpose(H_final(j,:))*transpose(H_final(j,:))');
%     end
    input_ = [];
    ofdm_out = [];
    pilot_out = [];
    for tx = 1:tx_chan
        for syms = 1:ofdm_mod.NumSymbols
            input_(:,syms) = sum(conj(H_final(:,tx,1)) .* input(:,syms),2)./sum(abs(H_final(:,:,1)).^2,2);
%             input_(:,syms) = input(:,syms).*H_inv(:,tx);
        end
        pilot_out(:,tx) = ofdm_mod_p(sum(conj(H_final(:,tx,1)) .* pilot,2)./sum(abs(H_final(:,:,1)).^2,2));
%         pilot_out(:,tx) = ofdm_mod_p(pilot);
        ofdm_out(:,tx) = [pilot_out(:,tx);ofdm_mod(input_)];
%         ofdm_out(:,tx) = ofdm_out(:,tx)/max(abs(ofdm_out(:,tx)));
    end
    y = awgn(h([ofdm_out]),snr,'measured');
%     y = y./max(abs(y));
        
    err = zeros(NFFT-1,1);
%     H_fft = [];
ofdm_demod_out = [];
    for rx = 1:rx_chan
%         H_fft(:,rx) = fftshift(fft(y(cp_len+1:NFFT+cp_len,rx),NFFT));
        
        ofdm_demod_out(:,rx) = reshape(ofdm_demod(y(NFFT+cp_len+1:end,rx)),(ofdm_mod.FFTLength - 1)*ofdm_mod.NumSymbols,1);
    end
%     H_fft(ofdm_mod.FFTLength/2+1,:) = [];
rec_sym = [];
    for i = 1:ofdm_mod.NumSymbols
%         rec_sym(:,i) = sum(ofdm_demod_out((i-1)*(NFFT-1)+1:(NFFT-1)*i,:).*conj(H_fft),2)./sum(abs(H_fft).^2,2);
        rec_sym(:,i) = ofdm_demod_out((i-1)*(NFFT-1)+1:(NFFT-1)*i,1) - err;
%         for j = 1:NFFT-1
%             err(j) = min(rec_sym(j,i) - mod_train);
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
    mse(iter) = mean(mean(abs(H_final - H_compare).^2));
    mse_td(iter) = mean(mean(abs(ifft(H_final,NFFT) - mse_test).^2));
    
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