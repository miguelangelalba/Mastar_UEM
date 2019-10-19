% Demonstrates how to use huffman.m

clear all;
home;

% Number of source symbols
N = 4;

W = zeros(N, 1);
W = [0.5; 0.25; 0.125; 0.125];
for i = 1:length(W)
	% weight of this source symbol
	% W(i) = i;
	% label for this source symbol
	L{i} = sprintf('I am node #%d',W(i));
end

% size of output code alphabet (e.g. D = 2 for binary)
D = 2;

% generate Huffman code and display code tree on the screen
% C = huffman(D, W, L, 'screen')

% generate Huffman code and display code tree on the screen and using
% Graphviz DOT
C = huffman(D, W, L, 'screen', 'dot')
