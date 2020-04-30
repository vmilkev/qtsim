function [G] = GetSubGN(...
                        genes,...
                        beta,...
                        type,...
                        rngState...
                        )

% Generate the number 'n' of simplex networks of different types determined by 'type'.

% G is a cell vector, elements of which are the particular graph's adjacency matrices
% 'genes' is the number of genes in subnetwork;
% 'beta' parameter of subnetwork;
% 'type' = 'd1', 'd2', 'd3'.
% G - adjacency matrix of generated network.

szG = size(genes);
szB = size(beta);

if ( szG(1,1) ~= szB(1,1) )
    error('the number of elements in vectors of genes in network and beta parameters should be equal \n');
end

n = szG(1,1);

G = cell(n,1);

for i = 1:n
    if type == 'd1'
        G{i,1} = GNd1( genes(i,1), 1, beta(i,1), rngState );
    elseif type == 'd2'
        G{i,1} = GNd2( genes(i,1), 1, beta(i,1), rngState );
    elseif type == 'd3'
        G{i,1} = GNd3( genes(i,1), 1, beta(i,1), rngState );
    end
end

end