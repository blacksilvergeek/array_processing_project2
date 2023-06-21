% function STOI = evaculate_stoi(test_ratio, noise_type, mode, bf)
noise_type = '4';
mode = 'Rn unknown';
bf = 'MVDR';

test_ratios = 0:4:20;
STOIs = zeros(1, length(test_ratios));
for i = 1:length(test_ratios)
    test_ratio = test_ratios(i);
    [STOIs(i)] = evaculate_stoi(test_ratio, noise_type, mode, bf);
end
STOI_4 = STOIs;
% plot
figure;
plot(test_ratios, STOIs, 'LineWidth', 2);
xlabel('Noise ratio (dB)');
ylabel('STOI');
title('STOI vs. Test ratio');
grid on;


