clc;
clear all;
close all;

one_server_10_syms = [16.7028	32.7935	67.3108	132.853	274.376
10.7245	16.4022	29.2937	54.7454	104.548
5.95846	11.0207	19.2771	31.7841	57.7622
3.84777	7.62413	14.2242	24.2066	42.9735
2.78708	5.41192	11.5815	20.4644	34.5055
2.48489	4.38483	8.67923	17.8141	29.4529
];

two_server_10_syms = [9.09215	18.3436	36.0178	71.07	142.487
6.61584	8.58864	15.7513	30.9576	58.7259
3.52348	6.27721	11.3373	20.14	34.4448
2.41877	4.1522	7.65029	15.309	26.9795
1.87191	3.45541	6.43579	12.7673	22.7664
3.38779	3.79387	7.8917	11.6674	22.7941
];

four_server_10_syms = [4.64579	9.79735	18.2539	37.6213	73.5145
4.00555	4.93082	7.20293	12.9342	31.2191
10.4074	10.9151	12.7247	16.8517	27.5552
1.57344	2.76717	4.80775	8.70521	18.146
2.23725	2.78196	4.63908	8.1414	15.0749
3.99815	4.46723	7.79328	10.4164	15.2795
];

one_server_100_syms = [175.412	348.212	678.269	1334.42	2770.99
69.6183	131.681	219.937	399.614	806.595
40.0653	71.4343	123.63	231.285	466.651
29.8383	52.2445	94.2441	170.395	337.849
24.544	42.0383	74.3639	137.394	271.003
22.2956	35.8421	63.0749	117.798	225.302
];

two_server_100_syms = [88.5869	173.869	348.507	671.573	1348.39
40.2584	71.7422	118.092	201.45	412.7
24.2403	40.5504	72.7513	145.134	265.243
18.9405	30.0897	56.9702	100.155	194
16.4877	25.3197	43.0654	76.8043	163.35
16.7646	21.9239	41.6253	68.3439	138.646
];

four_server_100_syms = [44.659	86.3253	170.474	344.118	697.381
21.5075	38.9625	65.9285	124.159	214.091
21.0413	29.6031	48.2161	84.3668	151.134
12.0711	19.5648	31.3025	55.6089	114.497
9.21265	17.4188	28.7632	45.2597	86.4483
10.4009	16.5716	26.513	44.1969	87.6958
];


figure;
NFFT = 2.^(6:10);
semilogy(NFFT,one_server_10_syms*1e-3,'*-'); hold on;
semilogy(NFFT, ((NFFT + 64).*16.*10 + 255)./10e6, 'o--');
semilogy(NFFT, ((NFFT + 64).*16.*10 + 255)./20e6, 'x-.');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./5e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./2e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./1e6, '--');
legend('Processing time (1 cores/server)','Processing time (4 cores/server)','Processing time (8 cores/server)','Processing time (12 cores/server)', ...
    'Processing time (16 cores/server)','Processing time (20 cores/server)','Propagation time (10 MHz)','Propagation time (20 MHz)');
grid on;
ylabel('Time');
xlabel('FFT size');

figure;
NFFT = 2.^(6:10);
semilogy(NFFT,two_server_10_syms*1e-3, '*-'); hold on;
semilogy(NFFT, ((NFFT + 64).*16.*10 + 255)./10e6, 'o--');
semilogy(NFFT, ((NFFT + 64).*16.*10 + 255)./20e6, 'x-.');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./5e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./2e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./1e6, '--');
legend('Processing time (1 cores/server)','Processing time (4 cores/server)','Processing time (8 cores/server)','Processing time (12 cores/server)', ...
    'Processing time (16 cores/server)','Processing time (20 cores/server)','Propagation time (10 MHz)','Propagation time (20 MHz)');
grid on;
ylabel('Time');
xlabel('FFT size');

figure;
NFFT = 2.^(6:10);
semilogy(NFFT,four_server_10_syms*1e-3, '*-'); hold on;
semilogy(NFFT, ((NFFT + 64).*16.*10 + 255)./10e6, 'o--');
semilogy(NFFT, ((NFFT + 64).*16.*10 + 255)./20e6, 'x-.');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./5e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./2e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./1e6, '--');
legend('Processing time (1 cores/server)','Processing time (4 cores/server)','Processing time (8 cores/server)','Processing time (12 cores/server)', ...
    'Processing time (16 cores/server)','Processing time (20 cores/server)','Propagation time (10 MHz)','Propagation time (20 MHz)');
grid on;
ylabel('Time');
xlabel('FFT size');

figure;
NFFT = 2.^(6:10);
semilogy(NFFT,one_server_100_syms*1e-3,'*-'); hold on;
semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./10e6, 'o--');
semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./20e6, 'x-.');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./5e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./2e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./1e6, '--');
legend('Processing time (1 cores/server)','Processing time (4 cores/server)','Processing time (8 cores/server)','Processing time (12 cores/server)', ...
    'Processing time (16 cores/server)','Processing time (20 cores/server)','Propagation time (10 MHz)','Propagation time (20 MHz)');
grid on;
ylabel('Time');
xlabel('FFT size');

figure;
NFFT = 2.^(6:10);
semilogy(NFFT,two_server_100_syms*1e-3, '*-'); hold on;
semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./10e6, 'o--');
semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./20e6, 'x-.');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./5e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./2e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./1e6, '--');
legend('Processing time (1 cores/server)','Processing time (4 cores/server)','Processing time (8 cores/server)','Processing time (12 cores/server)', ...
    'Processing time (16 cores/server)','Processing time (20 cores/server)','Propagation time (10 MHz)','Propagation time (20 MHz)');
grid on;
ylabel('Time');
xlabel('FFT size');

figure;
NFFT = 2.^(6:10);
semilogy(NFFT,four_server_100_syms*1e-3, '*-'); hold on;
semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./10e6, 'o--');
semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./20e6, 'x-.');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./5e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./2e6, '--');
% semilogy(NFFT, ((NFFT + 64).*16.*100 + 255)./1e6, '--');
legend('Processing time (1 cores/server)','Processing time (4 cores/server)','Processing time (8 cores/server)','Processing time (12 cores/server)', ...
    'Processing time (16 cores/server)','Processing time (20 cores/server)','Propagation time (10 MHz)','Propagation time (20 MHz)');
grid on;
ylabel('Time');
xlabel('FFT size');