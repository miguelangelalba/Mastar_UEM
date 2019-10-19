function [h] = entropia(vector)
cadena = vector(:);
p = hist(cadena,[0:255]);
%figure, subplot(p);
p = p /length(cadena);
p = p(p>0);
h = -sum(p.*log2(p));
end