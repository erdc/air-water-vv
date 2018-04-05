function y = check_anderson
clc
close all
gauges = [1; 5; 7; 10; 13];
ls = 0.415;
smooth_factor = 2;

data_path = strcat(pwd, '/data/');
filebase = textread(strcat(data_path, 'roots.txt'), '%s');

calc = cell( length(filebase), 1);
%fid=fopen('cd.txt','w');
%fprintf(fid, '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s\n', ...
%    ['H', 'T', 'Lp', 'H0', 'Cd_high', 'Cdr2_high', 'Cd_med', 'Cdr2_med']);

for i = 1:length(filebase)
    calc{i} = get_calcs(data_path, filebase{i}, smooth_factor); 
    Cd_high = calc{i}.high.cd;
    Cd_med = calc{i}.med.cd;
    Cdr2_high = calc{i}.high.cdr2;
    Cdr2_med = calc{i}.med.cdr2;
    h = calc{i}.h;
    T = calc{i}.T;
    Lp = calc{i}.Lp;
    H0 = calc{i}.H0;
    
    H0rms_med = calc{i}.med.H0rms;
    H0rms_high = calc{i}.high.H0rms;
    
%    fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',...
%        [h, T, Lp, H0, 400, H0rms_high, Cd_high, Cdr2_high]);
%    fprintf(fid, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',...
%        [h, T, Lp, H0, 200, H0rms_med, Cd_med, Cdr2_med]);
end
%fclose(fid);

filename='anderson_calcs.mat';
save(filename, 'calc')

idx = find_specs(calc, 0.78, 0.21, 0.18);
spec = calc{idx}.norm.spec{2} / 10^4;
freq = calc{idx}.norm.freq;

hold on
plot(freq, spec(:, 5:end))
title('Anderson data, control')

figure()
hold on
spec = calc{idx}.norm.spec{2} / 10^4;
freq = calc{idx}.norm.freq;
plot(freq, spec(:,1))
title('Anderson data, N = 400')

%---------------
%    PLOT 4
%---------------
% fig4a = figure();
% hold on
% fig4b = figure();
% hold on
% idx = find_specs(calc, 1.36, 0.37, 0.11);
% 
% %SET MIN/MAX SOME OTHER TIME.
% for i = 1:length(gauges)
%     spec_high = calc{idx(1)}.high.spec;
%     spec_low = calc{idx(1)}.norm.spec;
%     
%     figure(fig4a)
%     plot(spec_low.freq, spec_low.spec{2}(:,gauges(i)))
%     
%     figure(fig4b)
%     plot(spec_high.freq, spec_high.spec{2}(:,gauges(i)))
% end

%---------------
%    PLOT 5
%---------------
% 
% idx = find_specs(calc, 0.91, -1, -1);
% fig5b = figure();
% hold on
% for i = 1:length(idx)
%     x = ( calc{idx(i)}.H0 / 100 ) / calc{idx(i)}.h;
%     y_norm = calc{idx(i)}.norm.spec.k{2};
%     y_med = calc{idx(i)}.med.spec.k{2};
%     y_high = calc{idx(i)}.high.spec.k{2};
%     plot(x, y_norm, 'o')
%     plot(x, y_med, 'o')
%     plot(x, y_high, 'o')
% end

% ---------------
%   PLOT 9
% ---------------
fig9a = figure();
fig9b = figure();
fig10a = figure();
fig10b = figure();
for i = 1:length(calc)
    lsh = ls/calc{i}.h;
    ccode = '';
    if abs( lsh - 0.78 ) < 0.1
        ccode = 'b*';
    elseif abs( lsh - 0.91 ) < 0.1
        ccode = 'ks';
    elseif abs( lsh - 1.36 ) < 0.1
        ccode = 'ro';
    end

    figure(fig9a)
    hold on
    xlabel('Re')
    ylabel('C_D')
    title('Figure 9a')
    plot(calc{i}.med.Re, calc{i}.med.cd,ccode)
    plot(calc{i}.high.Re, calc{i}.high.cd,ccode)

    figure(fig9b)
    hold on
    xlabel('KC')
    ylabel('C_D')
    title('Figure 9b')
    plot(calc{i}.med.KC, calc{i}.med.cd,ccode)
    plot(calc{i}.high.KC, calc{i}.high.cd,ccode)

    figure(fig10a)
    hold on
    xlabel('Q_{Re}')
    ylabel('C_D')
    title('Figure 10a')
    plot(calc{i}.med.Re/(lsh^1.5), calc{i}.med.cd,ccode)
    plot(calc{i}.high.Re/(lsh^1.5), calc{i}.high.cd,ccode)

    figure(fig10b)
    hold on
    xlabel('Q_{KC}')
    ylabel('C_D')
    title('Figure 10b')
    plot(calc{i}.med.KC/(lsh^1.5), calc{i}.med.cd,ccode)
    plot(calc{i}.high.KC/(lsh^1.5), calc{i}.high.cd,ccode)
end


% spec_plot(s, freq, sdata, gauges, i, spec_files{i});
% spec_plot(us, freq, smoothed, gauges, i, spec_files{i});
% h_plot(h, H_unsmoothed, xlocs, hdata, i, h_files{i});
% h_plot(hs, H_smoothed, xlocs, hdata, i, h_files{i});
% d = linspace(xlocs(veg_start_gauge) - xlocs(5), xlocs(end) - xlocs(5));
% hold on
% figure(kfig)
% plot(xlocs(veg_start_gauge:end) - xlocs(5), H_unsmoothed(veg_start_gauge:end)/H_unsmoothed(5), 'o')
% plot(d, exp(-k_unsmoothed*d))
% hold on
% plot(freq, dE_s(:,3))
% plot(freq, log10(E_u(:,5)))
% plot(freq, log10(E_u(:,13)))



end

function y = get_calcs(data_path, root, smooth_factor)

hstr = root(1:2);

if strcmp(hstr, '12')
    h = 0.305;
elseif strcmp(hstr, '18')
    h = 0.457;
elseif strcmp(hstr, '21')
    h = 0.533;
else
    fprintf(strcat('No such file \t', root, '\n') )
    return
end

tstr = root(6:7);
if strcmp(tstr,'12')
    T = 1.25;
elseif strcmp(tstr,'15')
    T = 1.5;
elseif strcmp(tstr,'17')
    T = 1.75;
elseif strcmp(tstr,'20')
    T = 2.0;
elseif strcmp(tstr,'22')
    T = 2.25;
else
    fprintf(strcat('No such file \t', root, '\n') )
    return
end

%-------------------------------------------------------
s_opts= {'_spec','_med_spec','_high_spec'};
h_opts = {'_H','_med_H','_high_H'};
s_files = strcat(data_path, root, s_opts, '.txt');
h_files = strcat(data_path, root, h_opts, '.txt');
gauge_lon=[6.1,6.4,7.0,26.0,26.9,27.4,27.9,28.5,29.5,31,32.7,34.4,36.2];
veg_start_gauge=5;
ls=0.415;
b=0.0064;
inc=5;

%table_data = read_table('./data/table.csv')

%s = figure();
%us = figure();
%h = figure();
%hs = figure();
% kfig = figure();
%-------------------------------------------------------

y = struct('name',root,'h',h, 'T',T,...
    'Lp',0,'H0',[0 0 0],'norm',{{}},'med',{{}},'high',{{}});
y.Lp = wavelen_calc(y.h,y.T,10,9.81);

for i = 1:length(s_files)

    xlocs = gauge_lon';
    
    h_skip = 0;
    s_skip = 0;
    if ~exist(s_files{i}, 'file')
        warningMessage = sprintf('spec file does not exist:\n%s\n%s', ...
            s_files{i}, 'Cannot calculate k, epsilon, leaving N/A');
        warning(warningMessage);
        s_skip = 1;
    else
        [freq, sdata, smoothed] = read_data_spec(s_files{i}, smooth_factor);
    end
    
    if ~exist(h_files{i}, 'file')
        warningMessage = sprintf('spec file does not exist:\n%s\n%s', ...
            h_files{i}, 'Cannot calculate Cd');
        warning(warningMessage);
        h_skip = 1;
    else
        [~, hdata] = read_data_height(h_files{i});
    end
    
    if ~s_skip
        df = freq(2) - freq(1);
        H_u = h_calc(sdata, df);
        H_s = h_calc(smoothed, df);
        [k_u, ~, kur2] = k_calc(H_u/H_u(inc), xlocs, veg_start_gauge);
        [k_s, ~, ksr2] = k_calc(H_s/H_s(inc), xlocs, veg_start_gauge);
        [~,E_u, dE_u] = a_calc(sdata,df);
        [~,E_s, dE_s] = a_calc(smoothed,df);
        y.H0(i) = H_u(veg_start_gauge);
    else
        freq = [];
        df = NaN;
        H_u = [];
        H_s = [];
        k_u = NaN;
        kur2 = NaN;
        k_s = NaN;
        ksr2 = NaN;
        E_u = [];
        dE_u = [];
        E_s = [];
        dE_s = [];
        y.H0(i) = NaN;
    end
    
    if ~h_skip
        if i == 1
            data_struct = struct('freq',freq,'H',{{H_u; H_s}},...
              'k',{{k_u; k_s}}, 'kr2', {{kur2; ksr2}}, ... 
              'E', {{[E_u, dE_u];[E_s,dE_s]}}, ...
              'spec', {{sdata; smoothed}}, ...
              'cd',NaN,'cdr2',NaN,'u_c',NaN,...
              'Re',NaN,'KC',NaN,'H0rms',hdata(veg_start_gauge,3));            
            y.norm = data_struct;
            continue
        elseif i == 2
            N = 200;
        elseif i == 3
            N= 400;
        end
        [cd, cdr2, u_c, Re, KC] = cd_calc(hdata(:,3), N, xlocs, h, T, ls, b, veg_start_gauge);
    end

    data_struct = struct('freq',freq,'H',{{H_u; H_s}},...
        'k',{{k_u; k_s}}, 'kr2', {{kur2; ksr2}}, ... 
        'E', {{[E_u, dE_u];[E_s,dE_s]}}, ...
        'spec', {{sdata; smoothed}}, ...
        'cd',cd,'cdr2',cdr2,'u_c',u_c,...
        'Re',Re,'KC',KC,'H0rms',hdata(veg_start_gauge,3));
   
    if i == 2
       y.med = data_struct; 
    elseif i == 3
       y.high = data_struct;
    end
end

y.H0 = nanmean(y.H0) / 100;

end

function [freq, sdata, smoothed] = read_data_spec(filename, sf)

%PLOT CALCULATED SPECTRA (frequencies under 2 Hz)
M = dlmread(filename);
M( M(:,1) > 2, : ) = [];
freq = M(:,1);
sdata = M(:,2:end);

%smooth this data
smoothed = zeros(size(sdata));
for j = (1+sf):length(sdata)-sf
    smoothed(j,:) = mean(sdata(j-sf:j+sf,:));
end

end

function [xlocs, hdata] = read_data_height(filename)

N = dlmread(filename);
xlocs = N(:,1);
hdata = N(:,2:end);

end

function spec_plot(fig, freq, data, gauges, index, tit)
figure(fig)
subplot(2,2,index)
hold on

for i = 1:length(gauges)
    plot(freq, data(:,gauges(i)))
    axis([0 2 0 50])
end

legend('WG1', 'WG5', 'WG7', 'WG10', 'WG13')
t=title(tit);
set(t,'Interpreter','none')
xlabel('Frequency (Hz)')
ylabel('S(f) (m^2*s)')
end

function H = h_calc(data, df)
H = 4*sqrt(sum(data*df, 1));
H = H';
end

function h_plot(fig, H, xlocs, hdata, index, tit)
figure(fig)
subplot(2,2,index)
hold on
plot(xlocs, H,'x')
plot(xlocs,hdata(:,2),'o')
legend('Calculated', 'Text File','Location','SouthWest')
t=title(strcat('Zero Moment Wave Height Comparisons',tit));
set(t,'Interpreter','none')
xlabel('Gauge Location (m)')
ylabel('Wave height (cm)')
axis([0 50 0 14])
end

function [k, res, rsq] = k_calc(h_ratio, xlocs, shift)
%We have to shift xlocs, because this is taken with x = 0 at gauge 5.

xlocs = xlocs(shift:end) - xlocs(5);
fun = @(k)exp(-k*xlocs) - h_ratio(shift:end);
x0 = 0.1;
opts1 = optimset('display','off');
[k, res] = lsqnonlin(fun,x0,[],[],opts1);

yresid = fun(k);
SSresid = sum(yresid.^2);
SStotal = (length(h_ratio(shift:end)) - 1) * var( h_ratio(shift:end) );
rsq = 1 - SSresid / SStotal;
end

function [a, E, dE] = a_calc(data, df)
a = sqrt(2*data*df);
E = 0.5*1*9.81*a.^2;
E5 = E(:,5);
dE = [(E(:,7)-E5)./E5, (E(:,10) - E5)./E5, (E(:,13)-E5)./E5];
end

function [L] = wavelen_calc(d, T, n, g)
%
% function wavelength(d,T,n,g);
%
% n = # of iterations
%
% d = depth, this can be an array
% T = period
% n = number of iterations
% g = 32.2 or 9.81

Leck = (g*(T^2)*0.5/pi)*sqrt(tanh(4*pi*pi*d/(T*T*g)));  % 1984 SPM, p.2-7
L1 = Leck;

for k=1:n
    L2 = (g*(T^2)*0.5/pi)*tanh(2*pi*d./L1);   % check if it's right
    L1 = L2;                                 % redo
end

L = L2;
k=2*pi./L;

ko = find(d<=0);
L(ko) = 0;
end

function [cd, cdr2, u_c, Re, KC] = cd_calc(hrms, N, xlocs, h, T, ls, b, shift)
xlocs = xlocs(shift:end)-xlocs(5);
k = 2*pi/wavelen_calc(h,T,10,9.81);
H0=hrms(5);
ratio = hrms(shift:end)/H0;

H0 = H0/100;

if ls/h > 1
    s = 1;
else
    s=ls/h;
end

term1=(1/(3*sqrt(pi)))*b*N*H0*k;
term2=(sinh(k*s*h)^3+3*sinh(k*s*h))/((sinh(2*k*h)+2*k*h)*sinh(k*h));

beta=term1*term2;

fun = @(cd) 1./(1 + beta*xlocs*cd) - ratio;
opts1 = optimset('display','off');
[cd, res] = lsqnonlin(fun, 2,[],[],opts1);

cdres = fun(k);
SSresid = sum(cdres.^2);
SStotal = (length(ratio) - 1) * var( ratio );
cdr2 = 1 - SSresid / SStotal;

u_c = H0/2 * (2*pi/T) * (cosh(k*s*h)/sinh(k*h));

d=linspace(xlocs(1), xlocs(end));

nu = 10^-6;
Re = u_c * b / nu;
KC = u_c * T / b;

%hold on
%plot(xlocs, ratio, 'o')
%plot(d, 1./(1+beta*cd*d))

end

function [B, K, c] = cdfit_calc(calc)

for i=1:length(calc)
    
end

end

function idx = find_specs(calc, lsh, h0h, hlp)
    idx = zeros(length(calc), 1);
    ls = 0.415;
    tol = 0.02;
    for i = 1:length(calc)
        lshi = abs(ls / calc{i}.h - lsh);
        h0hi = abs(calc{i}.H0 / ((calc{i}.h)) - h0h);
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