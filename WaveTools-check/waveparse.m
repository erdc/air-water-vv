% Reads in free surface elevations and true spectra created by WaveTools
[f, S] = read_data('waves3D.csv',2);
true_spectrum = csvread('true_spectrum.csv');
size(f)
size(S)
figure()
hold on
smooth_fac = floor(0.02*length(f));
S = smoothdata(S,smooth_fac);
if length(f)  ~= length(S)
    f = [f;0];
end
plot(f, S)
plot(true_spectrum(:,1), true_spectrum(:,2))
xlabel('f(Hz)');
ylabel('S(f)');
legend('Calculated Spectrum', 'True Spectrum');
title('Control')
print('spectra_t', '-dpng')







