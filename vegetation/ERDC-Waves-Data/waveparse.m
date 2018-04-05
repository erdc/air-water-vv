function y = waveparse
%close all;

h = 0.533;
T = 1.5;
H0_ref = 0.111;

dt = 0.01;

a_cds = importdata('anderson_calcs.mat');

[f, S, H0, Hrms] = read_data('column_gauge3D_control_53.3_1.1.csv', dt, T);
idx = find_trial(a_cds, h, T, H0_ref);

figure()
hold on
plot(f, S(:,1))
plot(a_cds{idx}.norm.freq, a_cds{idx}.norm.spec{2}(:,1) / 10^4)
legend('WT', 'Anderson')
title('Control')


% 
% H0
%  Hrms
% Hrms2
% 
% % figure()
% % hold on
% % plot(fn, Sn)
% % title('N = 400')
% 
% figure()
% hold on
% plot(f2, S2)
% title('2D')

% figure()
% title('Side by side, gauge 5')
% hold on 
% plot(f, S(:,1))
% plot(fn, Sn(:,1))

end

%------------------------------------------------------
% FUNCTION: read_data
%   Takes file of gauge data and returns spectral data and height.
%   Takes peaks and troughs of each wave to determine root mean square
%   height.
%
% Inputs: 
%   - filename: [string], name of file.
%   - dt: [double], Desired timestep interval (will interpolate data to this time
%   interval)
%   - T: [double], dominant period of waves.
%
% Outputs:
%   - f: [double array], array of spectral frequencies.
%   - S: [double array, size(f) x num_gauges], spectral data.
%   - H0: [double array, 1 x num_gauges], calculated zero moment wave height.
%   - Hrms: [double array, 1 x num_gauges], root mean square wave height.
%------------------------------------------------------

function [f, S, H0, Hrms] = read_data(filename, dt, T)
% Variables, change if necessary.
A = csvread(filename,1,0);
time = A(:,1);
data = 1.05 - A(:,2:end);

peak_dist = 0.3;
time_shift = 20;


[n, bin] = histc(time, unique(time));
multiple = find(n > 1);

for i=1:length(multiple)
    index = find(ismember(bin, multiple()));
    first = index(1);
    last = index(end);
    data(first,:) = mean( data(first:end, : ) );
    data(first+1:end, :) = [];
    time(first+1:end) = [];
end

Hrms = zeros(length(data(1,:)), 1);
H0 = zeros(length(data(1,:)), 1);

for i = 1:length(data(1,:))

    if i == 1
        
        [f, df, S] = fft_gauge(time, data(:,1), dt);
        tmp = S;
    else
        [~, ~, tmp] = fft_gauge(time, data(:,i), dt);
        S = [S, tmp];
    end
    H0(i) = z_mom_wh(tmp, df);
    demeaned = data(:,i) - mean(data(:,i));
    [pks_high, idx_high] = findpeaks(demeaned(time>time_shift), time(time>time_shift), 'MinPeakDistance', T-peak_dist);
%    

    pks_low = zeros(length(pks_high), 1);
    ind_low = zeros(length(pks_high), 1);
    
    for j = 1:(length(idx_high)-1)
        idx_tmp_low = find( time == idx_high(j) );
        idx_tmp_hgh = find( time == idx_high(j+1) );
        pks_low(j) = min(demeaned(idx_tmp_low:idx_tmp_hgh));
        ind_low(j) = time( demeaned == pks_low(j) );
    end
    
    pks_low(end) = min(demeaned(find( time == idx_high(end)):end));
    if pks_low(end) > 0
       pks_low = pks_low(1:end-1);
       ind_low = ind_low(1:end-1);
       pks_high = pks_high(1:end-1);
    else
       ind_low(end) = time( demeaned == pks_low(end) );
    end
%
%    hold on
%    findpeaks(demeaned, time, 'MinPeakDistance', T-peak_dist);
%    plot(ind_low, pks_low, 'o')
    
    H_temp = pks_high - pks_low;

    Hrms(i) = sqrt( 1/length(pks_high) * sum(H_temp.^2));
end

% [pks_low, idx_low] = findpeaks(-1*demeaned, time, 'MinPeakDistance', T-0.25);
% hold on
% findpeaks(demeaned, time, 'MinPeakDistance', T-0.25)
% plot(idx_low, -1*pks_low, 'o')
end

%------------------------------------------------------
% FUNCTION fft_gauge
% Takes data about wave heights over time, returns FFT decomposition.
% INPUT:
%   - time: array, time each data point was collected.
%   - data_col: height data for one individual gauge over all time.
%   - dt: the desired time interval, used for interpolation.
% OUTPUT:
%   - f: spectral frequency.
%   - df: increment between two calculated frequencies.
%   - P1: contribution of each frequency, as calculated through FFT.
% ------------------------------------------------------

function [f, df, P1] = fft_gauge(time, data_col, dt)
    %given that times are uneven, we need to smooth them.
    tt = time(1):dt:time(end);
    
    if mod(length(tt),2) ~= 0
        tt = tt(1:end-1);
    end
    
    yy = spline(time, data_col, tt);
    L = length(tt);
%    plot(tt,yy,'.')

    yy = yy - mean(yy);
    
    [P1, f] = powerspec(yy',length(yy), dt, 1);
    P1 = P1*2*pi;
    f = f/(2*pi);
    P1( f > 2 ) = [];
    f ( f > 2 ) = [];
    P1 = smoothdata(P1, floor(0.008*length(P1)));
     
%     Y = fft(yy);
%     P2 = abs(Y/L);
%     P1 = P2(1:L/2+1);
%     P1(2:end-1) = 2*P1(2:end-1);
%     f = (1/dt)*(0:L/2)/L;
%     
%     P1( f > 2 ) = [];
%     f ( f > 2 ) = [];
%     
%     P1 = P1 + [0, 0, P1(1:end-2)] + [0, P1(1:end-1)] + [P1(2:end), 0] + [P1(3:end), 0, 0];
%     P1 = P1'/5;
%     
    df = f(2) - f(1);
end

function wave_animate(t,d,gl)
    
    padding = 1;
    max_range=max(d(:))+padding;
    min_range=min(d(:))-padding;
    h = plot(NaN, NaN,'o');
    axis([0 40 min_range max_range])
    
    prev = t(1);
    
    for ii = 2:length(t)
    pause((t(ii) - prev)/10)
    prev = t(ii);
    set(h, 'XData', gl, 'YData', d(ii,:));
%    title(strcat('Heights at ', num2str(t(ii)), ' seconds'))
    drawnow 
    end
    
end

function wave_gif(t,d,gl,filename)
 
figure(1)

    padding = 1;
    max_range=max(d(:))+padding;
    min_range=min(d(:))-padding;
%    h = plot(NaN, NaN,'o');
%    axis([0 40 min_range max_range])

% for n = 1:length(t)
%       plot(gl,d(n,:),'o')
%       axis([0 40 min_range max_range])
%       xlabel('Gauge Location (m)')
%       ylabel('Height of water (m)')
%       title(strcat('Waves at Time :', num2str(t(n),'%1.2f'), ' Seconds'))
%       drawnow
%       frame = getframe(1);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if n == 1;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.01);
%       end
% end
end

function h = z_mom_wh(S, df)
    h = 4*sqrt(sum(S*df));
end

function idx = find_specs(calc, lsh, h0h, hlp)
    idx = zeros(length(calc), 1);
    ls = 0.415;
    tol = 0.02;
    for i = 1:length(calc)
        lshi = abs(ls / calc{i}.h - lsh);
        h0hi = abs(calc{i}.H0 / (100*(calc{i}.h)) - h0h);
        hlpi = abs(calc{i}.h / calc{i}.Lp - hlp);
        
        [lshi, h0hi, hlpi];
        
        if lshi < tol && h0hi < tol && hlpi < tol
            idx(i) = i;
        elseif lshi < tol && h0h == -1 && hlp == -1
            idx(i) = i;
        elseif h0hi < tol && lsh == -1 && hlp == -1
            idx(i) = i;
        elseif hlpi < tol && hlp == -1 && lsh == -1
            idx(i) = i;
        end 
    end
    
    idx = idx(idx > 0);
end

function idx = find_trial(calc, h, T, H0)
    idx = zeros(length(calc), 1);
    tol = 0.01;
    
    
    for i = 1:length(calc)
        h_i = calc{i}.h;
        h_T = calc{i}.T;
        h_H0 = abs(calc{i}.H0 - H0);
        if h_H0 < tol && h_i == h && h_T == T
            idx(i) = i;
        end 
    end
    
    idx = idx(idx > 0);
end



