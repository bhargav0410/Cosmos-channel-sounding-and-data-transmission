clc;
clear all;
close all;

Fs = 100e6;
pn_len = 255;
thres = 1;
delay_prev = 0;
delay_now = 0;
multi = 0;
num_chan = 8;

while 1
for i = 1:num_chan
    str = ['corr_ch_',char(num2str(i-1)),'_binary'];
    fid = fopen(str);
    rec(:,i) = fread(fid,pn_len*2,'int32');
    fclose('all');
end

for i = 1:num_chan
%    if i == 1
        [maxval, maxarg] = max(abs(double(rec(:,i))/double(intmax)));
        for j = 1:length(rec(:,i))
        if abs(double(rec(j,i))/double(intmax))/maxval > 0.05
            thres = max(j,1);
            break;
        else
            continue;
        end
        end
%    end
    rec(:,i) = circshift(rec(:,i),length(rec(:,i)) - thres + 3);
    figure(1);
    subplot(round(num_chan/4),round(num_chan/2),i);
    %plot((0:49)*1/Fs, 10*log10(abs(double(rec(1:50,i))/double(intmax))));
    plot((0:49)*1/Fs, abs(double(rec(1:50,i))/double(intmax)));
    axis([0 49*1/Fs 0 0.1]);
    title_str = ['Channel ',char(num2str(i-1))];
    title(title_str);
    xlabel('Time Delay (in \mus)');
    ylabel('Power Delay Profile');
%    plot_mat(:,i*5) = (abs(double(rec(1:100,i))/double(intmax)));
end

%fin_plot = mesh(((0:length(plot_mat)-1)*1/Fs)*1e6,(1:5*i)/5,plot_mat'); hold on;
%set(gca,'xdir','reverse');
% rotate(fin_plot,[0 0 1],90);
%xlabel('Time Delay (in \mus)');
%ylabel('Channel');
%zlabel('Power Delay Profile');
%axis([0 1 1 16 0 1]);
pause(0.5);
%hold off;
end
    
