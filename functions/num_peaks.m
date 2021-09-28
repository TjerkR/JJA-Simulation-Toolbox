function [num] = num_peaks(Ic_data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

M = max(Ic_data);
peaks = findpeaks(Ic_data);
peaks = peaks(peaks > M/4);

num = numel(peaks);

end

