clear all, close all, clc
[lena, MAP_lena] = imread ('lena512.bmp');
lena = imresize(lena,[64 64]);
cadena = lena(:);
p = hist(cadena,[0:255]);
p = p /length(cadena);
p = p(p>0);


D =2
%eje_x = 0:0.1:1
%eje_x = 0:1:255
for i = 1:length(cadena)
    h(i)= -sum(p.*log2(p));
end
plot(cadena,h)

%h = -sum(p.*log2(p));

for i = 1:length(h)
	% weight of this source symbol
	% W(i) = i;
	% label for this source symbol
	L{i} = sprintf('I am node #%d',h(i));
end
C = huffman(D, h, L, 'screen', 'dot')
