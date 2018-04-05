function [ data2 ] = smoothdata(data, N )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
P1=data;
P1new = data;
    for I=1:N
        P1new = P1new + [zeros(I,1); P1(1:end-I)] + [P1(I+1:end) ;zeros(I,1)];
    end
    data2 = P1new/(N*2+1);
end

