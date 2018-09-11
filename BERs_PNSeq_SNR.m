clc;
clear;
close all;

% Parameters for simulation

NFFT = 2^10; % Used for displying frequency domain stats
FS = 200e6; % Used during simulation of up and down conversion
F = linspace(-1,1,NFFT) * FS/2; % Used for displying frequency domain stats
RATE = FS/2; % Sampling rate used for simulation 

UPSAMPLING_FACTOR = FS / RATE; 
DOWNSAMPLING_FACTOR = FS / RATE; 
DECIMATION_FACTOR = FS / RATE;

NUMBER_OF_SYMBOLS = 1000;
NUMBER_OF_REF_SYMBOLS = 1;


% Raised Cosine Transmit Filter

FILTER_SPAN = 6;           % Filter span in symbol durations
BETA = 0.1;         % Roll-off factor

rctFilt = comm.RaisedCosineTransmitFilter(...
    'Shape',                  'Normal', ...
    'RolloffFactor',          BETA, ...
    'FilterSpanInSymbols',    FILTER_SPAN, ...
    'OutputSamplesPerSymbol', UPSAMPLING_FACTOR);

% Raised Cosine Receive Filter

FILTER_SPAN = 6;           % Filter span in symbol durations
BETA = 0.1;         % Roll-off factor

rcrFilt = comm.RaisedCosineReceiveFilter(...
    'Shape',                  'Normal', ...
    'RolloffFactor',          BETA, ...
    'FilterSpanInSymbols',    FILTER_SPAN, ...
    'InputSamplesPerSymbol', DOWNSAMPLING_FACTOR, ...
    'DecimationFactor', DECIMATION_FACTOR );


% -------------------- PN sequence generation ---------------------------
% The code below generates PN Sequence using Linear Feedback Shift Register
% method. 

PN_SEQ_LENGTHS = [15, 31, 63, 127, 255];

SNR_VALUES = -20:1:20; 

ber_values = zeros(length(PN_SEQ_LENGTHS), length(SNR_VALUES));

h = comm.MIMOChannel;
h.SampleRate = FS;
h.SpatialCorrelation = false; % Independent channels
h.NumTransmitAntennas = 1;
h.NumReceiveAntennas = 1;
h.FadingDistribution = 'Rayleigh';
h.PathDelays = [0,1,2,3]*10e-9;
h.NormalizePathGains = false;
h.AveragePathGains = [0,-0.9,-4.9,-8];


% --------------- QPSK Data generation ----------------------------------
X_Input = randi([0 3],NUMBER_OF_SYMBOLS,1);
QPSKOutput = qammod(X_Input, 4);

for z = 1:length(PN_SEQ_LENGTHS)
    
    PN_SEQ_LEN = PN_SEQ_LENGTHS(z); % PN Sequence length 
    disp(['PN_SEQ_LEN = ',int2str(PN_SEQ_LEN)])
    INIT = randi(2,1,nextpow2(PN_SEQ_LEN)) - 1; % Initial values for the register
    while INIT == 0
        INIT = randi(2,1,nextpow2(PN_SEQ_LEN)) - 1;
    end 
    % seed polynomial is stored in POLY
    if PN_SEQ_LEN == 15
        POLY = [4 3 0]; 
    elseif PN_SEQ_LEN == 31
        POLY = [5 3 0];
    elseif PN_SEQ_LEN == 63
        POLY = [6 5 0];
    elseif PN_SEQ_LEN == 127
        POLY = [7 6 0];
    elseif PN_SEQ_LEN == 255
        POLY = [8 6 5 4 0];
    end
    pnSequence = comm.PNSequence('Polynomial',POLY,'SamplesPerFrame',PN_SEQ_LEN,'InitialConditions',INIT);
    PN_Seq = 2*pnSequence() - 1;
    
    % -----------------------------------------------------------------------
    
    
    SPREAD_FACTOR = ceil(PN_SEQ_LEN/2); % Decides spread factor
    % Decides the number of zeros to be added at the end of packet since
    % the PN sequence langth and the spreaded symbol length does not always
    % match
    if mod(NUMBER_OF_SYMBOLS*SPREAD_FACTOR,PN_SEQ_LEN) == 0
        TAIL_ZEROS = 0;
    else
        TAIL_ZEROS = PN_SEQ_LEN - mod(NUMBER_OF_SYMBOLS*SPREAD_FACTOR,PN_SEQ_LEN);
    end
    % Spreading QPSK symbols according to spread factor
        % Multiplication of PN and QPSK
        
        ChanEstOutput(1:NUMBER_OF_REF_SYMBOLS) = 1;
        ChanEst_spread = repmat(ChanEstOutput,PN_SEQ_LEN,1);
        QPSKOutput_spread = repmat(QPSKOutput,1,SPREAD_FACTOR);
        QPSKOutput_spread = transpose(QPSKOutput_spread);
        QPSKOutput_spread = [QPSKOutput_spread(:);zeros(TAIL_ZEROS,1)];

        Y = repmat(PN_Seq,length(QPSKOutput_spread)/PN_SEQ_LEN + NUMBER_OF_REF_SYMBOLS,1) .* [ChanEst_spread;QPSKOutput_spread];
%     end
    
    t_data = (0:length(Y)-1)/RATE;
 
    y_rct = rctFilt(Y);
    t_rct = (0:length(y_rct)*UPSAMPLING_FACTOR-1)/FS;
    
    ber_values_for_pn = [];
    
    y_filt = step(h,y_rct);
    
    y_rcr = rcrFilt([y_filt]);
    
    % Get BERs for different SNR values 
    for index = 1:length(SNR_VALUES)
        
        snr = SNR_VALUES(index); 
        disp(['snr = ',num2str(snr)])
        
        % Add noise
        y_rcr_noise = awgn(y_rcr,snr);
        
        % -----------------------------------------------------------------------
        
        %------------ Correlating with PN sequence at receiver ------------
        iter = 0;
        while iter <= length(Y)-length(PN_Seq)
            temp = transpose(y_rcr_noise(1+iter:iter+length(PN_Seq)))*PN_Seq;
            out(iter+1) = temp/(PN_Seq'*PN_Seq);
            iter = iter+1;
        end
        
        % Finding received signal
        
        chan_est = transpose(out(1:PN_SEQ_LEN));
        
        THRESHOLD_FOR_CORR_PEAKS = 0.5;
        
        for i = 1:length(chan_est)
            if abs(chan_est(i)) < THRESHOLD_FOR_CORR_PEAKS
                chan_est(i) = 0;
            end
        end
        chan_est_thres_delay = find(abs(chan_est)>=THRESHOLD_FOR_CORR_PEAKS);
        chan_est_thres = chan_est(chan_est_thres_delay);
        
%         y_reshaped = reshape(y_rcr_noise(length(PN_Seq)+1:end),SPREAD_FACTOR,NUMBER_OF_SYMBOLS-1);
        y_mod_data = y_rcr_noise(PN_SEQ_LEN*NUMBER_OF_REF_SYMBOLS+1:end);
%         y_mod_data_reshape = reshape()
%         for i = 1:size(y_reshaped,2)
            y_copy = zeros(length(y_mod_data),length(chan_est_thres_delay));
            for j = 1:length(chan_est_thres_delay)
                y_copy(1:end+(-chan_est_thres_delay(j)+1),j) = conj(chan_est_thres(j))*y_mod_data(chan_est_thres_delay(j):end);
            end
            Y_t = sum(y_copy,2)./sum(abs(chan_est_thres).^2);
            y_pn_mult = (repmat(PN_Seq,length(QPSKOutput_spread)/PN_SEQ_LEN,1).*Y_t);
            y_reshaped = transpose(reshape(y_pn_mult(1:end-TAIL_ZEROS),SPREAD_FACTOR,NUMBER_OF_SYMBOLS));
            Y_final = sum(y_reshaped,2)/SPREAD_FACTOR;
%         end
%         scatterplot(QPSKOutput);
%         scatterplot(Y_final);
        
        output_symbols = qamdemod(Y_final, 4);
        ber = sum(output_symbols ~= X_Input) / length(X_Input);
        
        ber_values_for_pn = [ber_values_for_pn ber];
    end  
    ber_values(z,:) = ber_values_for_pn;
end

% dlmwrite('BER_DIFF_PN_LENGTHS.txt',ber_values);
[X,Y] = meshgrid(1:1:length(PN_SEQ_LENGTHS), -20:1:20);
surf(Y, X, transpose(ber_values));
yticklabels(PN_SEQ_LENGTHS);
% savefig('BER_DIFF_PN_LENGTHS.fig');


