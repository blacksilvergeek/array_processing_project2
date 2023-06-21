%% figure;
test_ratios = 0:4:20;
load('STOI_1.mat')
load('STOI_2.mat')
load('STOI_3.mat')
load('STOI_4.mat')
load('STOI_all.mat')
plot(test_ratios, STOI_1, 'LineWidth', 2);
hold on 
plot(test_ratios, STOI_2, 'LineWidth', 2);
hold on 
plot(test_ratios, STOI_3, 'LineWidth', 2);
hold on
plot(test_ratios, STOI_4, 'LineWidth', 2);
hold on 
plot(test_ratios, STOI_all, 'LineWidth', 2);
legend('Clean speech 2', 'Babble noise', 'Artificial nonstationary noise','Stationary speech shaped noise','all noise')

xlabel('Noise range ratio (dB)');
ylabel('STOI');
title('STOI vs. Test ratio');
grid on;