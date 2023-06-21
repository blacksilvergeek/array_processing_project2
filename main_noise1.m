% function STOI = evaculate_stoi(test_ratio, noise_type, mode, bf)
noise_type = '1';
mode = 'Rn unknown';
bf = 'MVDR';

test_ratios = 0:4:20;
STOIs = zeros(1, length(test_ratios));
for i = 1:length(test_ratios)
    test_ratio = test_ratios(i);
    [STOIs(i)] = evaculate_stoi(test_ratio, noise_type, mode, bf);
end
STOI_1 = STOIs;
% plot
figure;
plot(test_ratios, STOIs, 'LineWidth', 2);
xlabel('Noise ratio (dB)');
ylabel('STOI');
title('STOI vs. Test ratio');
grid on;


%% Data generation
load('impulse_responses')
[s_target, fs] = audioread('clean_speech.wav');
s_target=speechnormalize(s_target);

% Noise generation
[bbn, fs] = audioread('babble_noise.wav');
[s2, fs] = audioread('clean_speech_2.wav');
[afn, fs] = audioread('aritificial_nonstat_noise.wav');
[ssn, fs] = audioread('Speech_shaped_noise.wav');

% Trim the audio data
lenght_min = min([size(s_target, 1), size(s2, 1), size(bbn, 1), size(afn, 1), size(ssn, 1)]);
s_target = s_target(1:lenght_min);
bbn = bbn(1:lenght_min);
s2 = s2(1:lenght_min);
afn = afn(1:lenght_min);
ssn = ssn(1:lenght_min);


% Convolution
s_1 = conv(s_target, h_target(1,:));
s_2 = conv(s_target, h_target(2,:));
s_3 = conv(s_target, h_target(3,:));
s_4 = conv(s_target, h_target(4,:));
s = [s_1, s_2, s_3, s_4];
s_max = max(abs(s_target));

bbn = speechnormalize(bbn)*s_max;
bbn_1 = conv(bbn, h_inter2(1,:));
bbn_2 = conv(bbn, h_inter2(2,:));
bbn_3 = conv(bbn, h_inter2(3,:));
bbn_4 = conv(bbn, h_inter2(4,:));
bbn_conv = [bbn_1, bbn_2, bbn_3, bbn_4];


s2 = speechnormalize(s2)*s_max;
cs_1 = conv(s2, h_inter1(1,:));
cs_2 = conv(s2, h_inter1(2,:));
cs_3 = conv(s2, h_inter1(3,:));
cs_4 = conv(s2, h_inter1(4,:));
cs_conv = [cs_1, cs_2, cs_3, cs_4];

afn = speechnormalize(afn)*s_max;
afn_1 = conv(afn, h_inter3(1,:));
afn_2 = conv(afn, h_inter3(2,:));
afn_3 = conv(afn, h_inter3(3,:));
afn_4 = conv(afn, h_inter3(4,:));
afn_conv = [afn_1, afn_2, afn_3, afn_4];

ssn = speechnormalize(ssn)*s_max;
ssn_1 = conv(ssn, h_inter4(1,:));
ssn_2 = conv(ssn, h_inter4(2,:));
ssn_3 = conv(ssn, h_inter4(3,:));
ssn_4 = conv(ssn, h_inter4(4,:));
ssn_conv = [ssn_1, ssn_2, ssn_3, ssn_4];



ratios = [ratio1; ratio2; ratio3; ratio4];
% convert SNR to linear scale
ratio_lin = 10.^(ratios/10);
n = cs_conv / ratio_lin(1) + bbn_conv / ratio_lin(2) + afn_conv / ratio_lin(3) + ssn_conv / ratio_lin(4);
x = s + n;



% n=bbn_conv + cs_conv + afn_conv + ssn_conv;
% n = n*1; %???
% x = s+n;
x = speechnormalize(x);

% sound(x(:,1),fs)
% sound(s_est,fs)
% clear sound
%% Shortern 
Len = 160*500;
x = x(20000:20000+Len-1, :);
s = s(20000:20000+Len-1, :);
s_target=s_target(20000:20000+Len-1,:);
n = n(20000:20000+Len-1,:);

%plot
figure
plot(s_target)
title('s_{target}')
xlabel('Time')
ylabel('Amplitude')
figure
plot(s(:,1))
title('Convolutional s_{target}')
xlabel('Time')
ylabel('Amplitude')
figure
plot(n(:,1))
title('noise')
xlabel('Time')
ylabel('Amplitude')
figure
plot(x(:,1))
title('noisy convolutional signal')
xlabel('time')
ylabel('Amplitude')


%% STFT
LenS = length(s);
M = 4;

fs = 16000; % Sampling frequency
WinLen_sec=0.02; % (Second)
Winlen_sample = WinLen_sec*fs; % Window length
overlap = 0.5;
overlap_sample = overlap*Winlen_sample;
NumFFT = Winlen_sample;
NumFrame = (LenS - overlap_sample)/(Winlen_sample-overlap_sample);
win = hamming(Winlen_sample, 'periodic');
xx = stft(x, fs, 'Window', win, 'OverlapLength', overlap_sample, 'FFTLength',NumFFT);
ss = stft(s_target, fs, 'Window', win, 'OverlapLength', overlap_sample, 'FFTLength', NumFFT);
nn = stft(n, fs, 'Window', win, 'OverlapLength', overlap_sample, 'FFTLength', NumFFT);

%% Beamformers
a = fft(transpose(h_target), NumFFT);
L = size(xx, 2); % Number of frames
ss_est = zeros(NumFFT, L);

SS = zeros(NumFFT, L);
NN = zeros(NumFFT, L);


% estimate Rn
for k = 1:NumFFT
    Rn_all{k} = zeros(M, M);
    for l = 1:L
        nnn=squeeze(nn(k,l,:));
        Rn_all{k} = Rn_all{k}+nnn*nnn'/L;
    end
end

% s_1 = conv(s_target, h_target(1,:));

ss_1 = stft(s_1, fs, 'Window', win, 'OverlapLength', overlap_sample, 'FFTLength', NumFFT);

% estimate sigma_s
for k = 1:NumFFT
    sigma_s_all{k} = 0;
    for l = 1:L
        sss=squeeze(ss_1(k,l,:));
        sigma_s_all{k} = sigma_s_all{k}+sss*sss'/L;
    end
end


for k = 1:NumFFT
    Rx_ini = 0;
    for l = 1:L
        aa = transpose(a(k,:));
        xxx = squeeze(xx(k,l,:));
        alpha = 1-1e-6;
        %alpha=1;
%         Rx = (l-1)/l*Rx_ini+1/l*(xxx*xxx');
        Rx = alpha*Rx_ini+(1-alpha)*(xxx*xxx');
        Rx = Rx+1e-10*eye(4);
        Rx_ini = Rx; % Update Rx
        %% a is known
        w_delay = pinv(aa);
        %% a is known, Rn is unknown
        % if mode == "a known"
        %     Rn = Rn_all{k};
        %     w_mvdr = inv(Rn)*aa/(aa'*inv(Rn)*aa+eps);
        % end
        %% a is known, Rn is known
        % if mode == "all known"
        %     nnn=squeeze(nn(k,l,:));
        %     Rn = nnn*nnn';
        %     Rn = Rn+1e-5*eye(4);
        %     w_mvdr = inv(Rn)*aa/(aa'*inv(Rn)*aa+eps);
        % end
        if mode == "Rn known"
            nnn=squeeze(nn(k,l,:));
            Rn = nnn*nnn';
            Rn = Rn+1e-5*eye(4);
            [Rs_est, a_est]=GEVD(Rx, Rn);
            w_mvdr=inv(Rn)*a_est/(a_est'*inv(Rn)*a_est+eps);
            w_single_wiener = sigma_s_all{k}/(sigma_s_all{k}+inv(a_est'*inv(Rn)*a_est));
            w_wiener = w_single_wiener*w_mvdr;
            w_delay = (a_est)/(a_est'*a_est+eps);
        end

        %% a is unknown, Rn is known
        %GEVD
        if mode == "Rn unknown"
            nnn=squeeze(nn(k,l,:));
            Rn = Rn_all{k};
            Rn = Rn+1e-5*eye(4);
            [Rs_est, a_est]=GEVD(Rx, Rn);
            w_mvdr=inv(Rn)*a_est/(a_est'*inv(Rn)*a_est+eps);
            w_single_wiener = sigma_s_all{k}/(sigma_s_all{k}+inv(a_est'*inv(Rn)*a_est)+eps);
            w_wiener = w_single_wiener*w_mvdr;
            w_delay = (a_est)/(a_est'*a_est+eps);
        end

        
        %% est_s
        if bf=="wiener"
            ss_est(k,l)=w_wiener'*xxx;
        elseif bf=="mvdr"
            ss_est(k,l)=w_mvdr'*xxx;
        elseif bf=="delay and sum"
            ss_est(k,l)=w_delay'*xxx;
        end
        % ss_est(k,l)=w_mvdr'*xxx;

        %% performance
        SS(k, l)=(abs(ss_est(k,l))).^2;
        nnn=squeeze(nn(k,l,:));
        NN(k,l)=(abs(w_mvdr'*nnn)).^2;
    end
end
OutputSNR_mean=mean(mean(ss))/mean((mean(NN)));

%% ISTFT
s_est = real(istft(ss_est, fs, 'Window', win, 'OverlapLength', overlap_sample, 'FFTLength', NumFFT));

%% Plot
s_est = speechnormalize(s_est);
figure
plot(s_est)
title('MVDR_estimated_signal')
xlabel('Time')
ylabel('Amplitude')

figure
plot(s_target)
hold on 
plot(s_est)
title('s_{target} and MVDR estimated signal')
legend('s_{target}', 's_{est}')
mse(s_target, s_est)
%% STOI
STOI = stoi(s_target, s_est, fs)







