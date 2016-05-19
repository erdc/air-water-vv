
%------------------------------------------------------
% FUNCTION: read_data
%   Takes file of gauge data and returns spectral data.
%
% Inputs: 
%   - filename: [string], name of file.
%
% Outputs:
%   - f: [double array], array of spectral frequencies.
%   - S: [double array, size(f) x num_gauges], spectral data.
%------------------------------------------------------

function [f, S] = read_data(filename)
A = csvread(filename,1,0);
time = A(:,1);
data = A(:,2);
[f, df, S] = fft_gauge(time,data);
end
