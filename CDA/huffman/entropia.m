function [h] = entropia(cadena)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p = hist(cadena, 0:255);
p = p / length(cadena);
p = p(p > 0);
h = -sum(p.*log2(p));
end

