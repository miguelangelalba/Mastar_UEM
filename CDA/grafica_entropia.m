clear all, close all, clc
[lena, MAP_lena] = imread ('lena512.bmp');
lena = imresize(lena,[64 64]);
cadena = lena(:);

D =2
eje_x = 0:0.1:1
%eje_x = 0:0.1:255
for i = 1:length(eje_x)
    p = eje_x(i);
    q =1-p;
    if (p > 0)||(q > 0)
        x =[p q];
        h(i)= -sum(x.*log2(x));
    end
end
plot(eje_x,h)

for i = 1:length(h)
	% weight of this source symbol
	% W(i) = i;
	% label for this source symbol
	L{i} = sprintf('I am node #%d',h(i));
end
C = huffman(D, h, L, 'screen', 'dot')
