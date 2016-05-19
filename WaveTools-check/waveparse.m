% Reads in free surface elevations and true spectra created by WaveTools
[f, S] = read_data('waves3D.csv',5);
true_spectrum = csvread('true_spectrum.csv');
size(f)
figure()
hold on
plot(f, S)
plot(true_spectrum(:,1), true_spectrum(:,2))
xlabel('f(Hz)');
ylabel('S(f)');
legend('Calculated Spectrum', 'True Spectrum');
title('Control')
print('spectra_t', '-dpng')







