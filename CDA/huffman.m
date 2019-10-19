function C = huffman(D,W,varargin)
% HUFFMAN   Huffman encoder.
%    C = HUFFMAN(D,W) generates a static minimum-variance Huffman tree and
%    corresponding codebook C for the source symbols with nonnegative
%    weights given by vector W, using a D-ary output code alphabet,
%    e.g. D = 2 for a binary output code alphabet.
%    W is an N-vector of nonnegative source symbol weights, where N is the
%    number of source symbols.
%    The generated Huffman code C minimizes the weighted codeword length,
%    e.g. if the weight is the probability of the source symbol, then C
%    minimizes the expected codeword length.
%    C is an N-vector of strings (i.e. cell array) giving the
%    codewords for each of the N source symbols in W.
%
%    C = HUFFMAN(D,W,L,'screen') additionally displays a simple text
%    rendition of the resulting Huffman code tree, with the source symbol
%    labels given by L. L is an N-vector of strings (i.e. cell array)
%    containing the labels for each of the N source symbols.
%
%    C = HUFFMAN(D,W,L,'dot') additionally creates a DOT file describing the
%    resulting Huffman code tree, and runs Graphviz DOT to create the
%    corresponding PNG image. The output files are automatically
%    time-stamped. The source symbol labels are given by L, which is an
%    N-vector of strings (i.e. cell array) containing the labels for each
%    of the N source symbols. Please ensure that Graphviz DOT is accessible
%    from the current directory (e.g. by adding it to the path).
%
%    Multiple output can be requested using a single command, e.g.
%    C = HUFFMAN(D,W,L,'dot','screen')


% VALIDATE INPUT ARGUMENTS

% check vector W describing the nonnegative source symbol weights
if ~isvector(W)
	error('Input argument W must be a vector.');
end

if min(W) < 0
	error('Input argument W must contain nonnegative entries.');
end

% number of source symbols
N = length(W);

% check code alphabet size D
if ~isscalar(D)
	error('Input argument D must be a scalar.');
end

codeAlphabet = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ';
if D > length(codeAlphabet)
	error(['Not enough literal characters to accommodate requested code alphabet; '...
		'please modify the source code manually to add more characters.']);
end

% source symbol labels
if ~isempty(varargin)
	L = varargin{1};

	if ~isvector(L)
		error('Input argument L must be a vector.');
	end

	if  length(L) ~= N
		error('Input arguments W and L must be vectors of the same length.');
	end
end

% output types
outputScreen = false;
outputDot = false;

for i = 2:length(varargin)
	output = varargin{i};

	if strcmp(output,'screen')
		outputScreen = true;
	elseif strcmp(output,'dot')
		outputDot = true;
	else
		error(['Input argument ''' output ''' is not a recognized output type.']);
	end
end


% BUILD HUFFMAN CODE TREE

% create a tree node for each of the N source symbols
for i = 1:N
	H(i).weight = W(i);  % weight
	H(i).children = [];  % this node has no children
	H(i).height = 0;     % this subtree has height 0
	H(i).exposed = true; % this node is exposed
	H(i).dummy = false;  % this node is not a dummy
end

% check if dummy nodes are needed
k = floor((N-1) / (D-1));
nDummyNodes = 0;
if (N-1) - k*(D-1) ~= 0
	nDummyNodes = (k+1) * (D-1) - (N-1);

	for i = 1:nDummyNodes
		H(N+i).weight = 0;     % weight
		H(N+i).children = [];  % this node has no children
		H(N+i).height = 0;     % this subtree has height 0
		H(N+i).exposed = true; % this node is exposed
		H(N+i).dummy = true;   % this node is a dummy
	end
end

% iteratively merge D subtrees at a time, until only 1 subtree is left
while sum([H.exposed] == true) > 1

	% sort exposed nodes in ascending order of weight
	exposedIndices = find([H.exposed] == true);
	allWeights = [H.weight];
	allHeights = [H.height];
	[Y,I] = sortrows([allWeights(exposedIndices)' allHeights(exposedIndices)'], [1 2]);
	sortedIndices = exposedIndices(I);

	% combine first D nodes into a new node
	combinedIndices = sortedIndices(1:D);
	for i = 1:D
		H(combinedIndices(i)).exposed = false;
	end

	ii = length(H) + 1;
	H(ii).weight = sum(allWeights(combinedIndices));  % sum of weights
	H(ii).children = sortedIndices(1:D);  % this node has D children
	H(ii).height = max(allHeights(combinedIndices)) + 1; % height of new subtree
	H(ii).exposed = true; % this node is exposed
	H(ii).dummy = false;  % this node is not a dummy
end

% PARSE CODE TREE TO GENERATE CODEWORDS

rootIndex = find([H.exposed] == true);
if isscalar(rootIndex) == false
	error('Internal error building Huffman tree.');
end

% initialize return value, i.e. codewords
for i = 1:N
	C{i} = '';
end

if outputDot
	% open DOT output file for writing
	filename = ['dot.' datestr(now,'yyyymmdd.HHMMSS')];
	fid = fopen([filename '.txt'],'w');

	fprintf(fid,'digraph Huffman {\nrankdir=TB');
end

% parse the code tree, beginning from the root
parseCodeTree(rootIndex,'');

if outputScreen
	% display statistics
	fprintf('\n\nSTATISTICS:');
	fprintf('\n   Number of source symbols N = %d', N);

	if nDummyNodes > 0
		fprintf('\n   Number of dummy nodes added = %d', nDummyNodes);
	end

	tw = 0;
	for i = 1:N
		tw = tw + W(i);
	end
	fprintf('\n   Total weight = sum{i=1:N} W_i = %g', tw);

	wcl = 0;
	for i = 1:N
		wcl = wcl + length(C{i}) * W(i);
	end
	fprintf('\n   (A) Weighted codeword length = sum{i=1:N} length(C_i) * W_i = %g', wcl);

	wll = 0;
	for i = 1:N
		wll = wll + length(L{i}) * W(i);
	end
	fprintf('\n   (B) Weighted label length    = sum{i=1:N} length(L_i) * W_i = %g', wll);

	fprintf('\n   Ratio (A) / (B) = %g', wcl / wll);

	fprintf('\n\n');
end

% call Graphviz dot if requested
if outputDot
	fprintf(fid,'\n}\n');
	fclose(fid);
	eval(['!dot -Tpng ' filename '.txt -o ' filename '.png']);

	% display image using default external program
	eval(['!' filename '.png']);

	% display image in MATLAB
	%[X, map] = imread([filename '.png']);
	%figure(1);
	%image(X);
	%colormap(map);
	%axis off;
	%axis image;
end


% NESTED FUNCTIONS

	function parseCodeTree(nodeIndex,nodeCodeword)
		% Parse nodes in the code tree recursively.

		node = H(nodeIndex);

		if outputScreen
			% create padding string
			padding = [''];
			for l = 1:length(nodeCodeword);
				padding = [padding '   '];
			end
		end

		if outputScreen
			if nodeIndex == rootIndex
				rootIndicator = ' (ROOT)';
			else
				rootIndicator = '';
			end
		end

		if isempty(node.children)
			if node.dummy == true
				% this is a dummy leaf node
				if outputScreen
					fprintf('\n%s[''%s'',%g]%s (DUMMY)',...
						padding, nodeCodeword, node.weight,rootIndicator);
				end

				if outputDot
					fprintf(fid,'\nNODEID%s [shape=record, label="''%s''\\n[%g]\\n(DUMMY)%s", style=filled, fillcolor=lightgrey]',...
						nodeCodeword, nodeCodeword, node.weight, rootIndicator);
				end
			else
				% this is a source symbol (non-dummy leaf node)
				C{nodeIndex} = nodeCodeword;

				if outputScreen
					fprintf('\n%s[''%s'',%g] ''%s''%s (SOURCE)',...
						padding, nodeCodeword, node.weight, L{nodeIndex}, rootIndicator);
				end

				if outputDot
					fprintf(fid,'\nNODEID%s [shape=record, label="''%s''\\n[%g]\\n''%s''\\n(SOURCE)%s", style=filled, fillcolor=salmon]',...
						nodeCodeword, nodeCodeword, node.weight, L{nodeIndex}, rootIndicator );
				end
			end
		else
			% this node is not a leaf
			if outputScreen
				fprintf('\n%s[''%s'',%g]%s',...
					padding, nodeCodeword, node.weight, rootIndicator);
			end

			if outputDot
				fprintf(fid,'\nNODEID%s [shape=record, label="''%s''\\n[%g]\\n%s", style=filled, fillcolor=white]',...
					nodeCodeword, nodeCodeword, node.weight, rootIndicator);
			end

			% recurse into child nodes
			for l = 1:length(node.children)
				childCodeword = [nodeCodeword, codeAlphabet(l)];

				if outputDot
					if H(node.children(l)).dummy == true
						% child is a dummy node
						fprintf(fid,'\nNODEID%s -> NODEID%s [sametail, label="%s", style=dotted]',...
							nodeCodeword, childCodeword, codeAlphabet(l));
					else
						% child is not a dummy node
						fprintf(fid,'\nNODEID%s -> NODEID%s [sametail, label="%s"]',...
							nodeCodeword, childCodeword, codeAlphabet(l));
					end
				end

				parseCodeTree(node.children(l), childCodeword);
			end
		end
	end

end
