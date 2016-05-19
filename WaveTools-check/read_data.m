
%------------------------------------------------------
% FUNCTION: read_data
%   Takes file of gauge data and returns spectral data.
%
% Inputs: 
%   - filename: [string], name of file.
%   - avg_over: number of values on each side over with to averag
%
% Outputs:
%   - f: [double array], array of spectral frequencies.
%   - S: [double array, size(f) x num_gauges], spectral data.
%------------------------------------------------------

function [f, S] = read_data(filename, avg_over)
A = csvread(filename,1,0);
time = A(:,1);
data = A(:,2);
[f, df, S] = fft_gauge(time,data, avg_over);
end
