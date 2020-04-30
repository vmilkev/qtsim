
function [GNA] = GetGNA(...
                        archType,...
                        genesNum_ntw1,...
                        beta_ntw1,...
                        genesNum_ntw2,...
                        beta_ntw2,...
                        genesNum_ntw3,...
                        beta_ntw3,...
                        mergingStrength,...
                        rngState...
                        )

% archType - parameter of string type; determines the type of GN Architecture; should be either 'downstream' or 'upstream';
% genesNum_ntw1 - column vector of number of genes (integer type) for each of basal networks of type 1;
% genesNum_ntw2 - column vector of number of genes (integer type) for each of basal networks of type 2;
% genesNum_ntw3 - column vector of number of genes (integer type) for each of basal networks of type 3;
% beta_ntw1 - column vector of network configuration parameter (float type) for each of basal networks of type 1,
%             beta = [0 10];
% beta_ntw2 - column vector of network configuration parameter (float type) for each of basal networks of type 2,
%             beta = [0 10];
% beta_ntw3 - column vector of network configuration parameter (float type) for each of basal networks of type 3,
%             beta = [0 10];
% mergingStrength - parameter which determines the number of contacts betweeng merging basal networks, mergingStrength = [0 inf], if
%                   mergeStrength = 0 - merging networks will be completely disconnected;
% createReport - logical parameter, if createReport = TRUE - a network statistics will be generated;
% viewGraphics - parameter of string type, if viewGraphics = 'on'('off') - a network statistics will be (not) printed on screen;
%                requires createReport = TRUE;
% printResToFile - logical parameter, if printResToFile = TRUE - a network statistics will be additionally saved to *.dat files;
%                  requires createReport = TRUE;
% GNA - adjacency matrix of GN Architecture.


%---------------------------------------------------------------------
% START of the main calculation stage
%---------------------------------------------------------------------

% check the minimum number of genes in networks
if ( min(genesNum_ntw1) < 2 )
    error('the minimum number of genes in the subnetwork of type 1 should be >= 2 \n');
elseif ( min(genesNum_ntw2) < 3 )
    error('the minimum number of genes in the subnetwork of type 1 should be >= 3 \n');
elseif ( min(genesNum_ntw3) < 4 )
    error('the minimum number of genes in the subnetwork of type 1 should be >= 4 \n');
end

% get expected total number of genes in resulting network
totalGenes = sum(genesNum_ntw1) + sum(genesNum_ntw2) + sum(genesNum_ntw3);

% get the number of basal networls of each types
sz_type1 = size(genesNum_ntw1);
sz_type2 = size(genesNum_ntw2);
sz_type3 = size(genesNum_ntw3);

% Generate basal networks of different types
% here we return ADJACENCY MATRICES (not matlab graph objects)
Gsub_type1 = GetSubGN( genesNum_ntw1, beta_ntw1, 'd1', rngState );
Gsub_type2 = GetSubGN( genesNum_ntw2, beta_ntw2, 'd2', rngState );
Gsub_type3 = GetSubGN( genesNum_ntw3, beta_ntw3, 'd3', rngState );

% define empty cell container for the generated networks
Gsub = {};

% fill-in the cell container
Gsub( 1:sz_type1(1,1), 1 ) = Gsub_type1;
Gsub( sz_type1(1,1)+1:sz_type1(1,1)+sz_type2(1,1), 1 ) = Gsub_type2;
Gsub( sz_type1(1,1)+sz_type2(1,1)+1:sz_type1(1,1)+sz_type2(1,1)+sz_type3(1,1), 1 ) = Gsub_type3;

% get total number of generated basan networks
nwSz = size(Gsub);

% we don't need these data any more
clearvars Gsub_type1 Gsub_type2 Gsub_type3;
clearvars genesNum_ntw1 genesNum_ntw2 genesNum_ntw3;
clearvars beta_ntw1 beta_ntw2 beta_ntw3;
clearvars sz_type1 sz_type2 sz_type3;

% merging the generated basal networks
% working only with adjacency matrices

if ( nwSz(1,1) > 0 ) % check whether the network container is empty
    
    % proceed with different cases (in order to effectively use of MergeGN() function
    % if the number of generated networks are small)
    
    if (nwSz(1,1) == 1)
        G = Gsub{1,1};
    elseif (nwSz(1,1) == 2)
        sz1 = size( Gsub{1,1} );
        sz2 = size( Gsub{2,1} );
        n = min( sz1(1,1), sz2(1,1) );
        G = MergeGN(Gsub{1,1},Gsub{2,1}, [], mergingStrength*n, rngState);
    elseif (nwSz(1,1) == 3)
        sz1 = size( Gsub{1,1} );
        sz2 = size( Gsub{2,1} );
        sz3 = size( Gsub{3,1} );
        n = min( sz1(1,1), sz2(1,1) );
        n = min( n, sz3(1,1) );
        G = MergeGN(Gsub{1,1},Gsub{2,1}, Gsub{3,1}, mergingStrength*n, rngState);
    else
        G = Gsub{1,1};
        for i = 2:nwSz(1,1)
            G2 = Gsub{i,1};
            sz1 = size( G );
            sz2 = size( G2 );
            n = min( sz1(1,1), sz2(1,1) );
            G = MergeGN(G,G2, [], mergingStrength*n, rngState);
        end
    end
    
else
    %fprintf('the number of subnetworks should be > 2 \n');
    error('the number of subnetworks should be >= 1 \n');
end

% transform G to a directed graph (make directed gene architecture)

% 1) get genes connectivity in the network (graph degree)
sz = size( G );
dG = zeros( sz(1,1), 1);
for i = 1:sz(1,1)
    dG(i,1) = numel( find (G(i,:)>0) );
end

% 2) sort the graph nodes according to their degree (lowest -first)
[D,I] = sort(dG);

% 3) renumber the adjacency matrix rows and columns according to the vector I
G = G(:,I);
G = G(I,:);

% generate a directed architecture of two different architectures using the
% renumbered adjacency matrix G
if ( strcmp(archType, 'downstream') )
    % downstream architecture: node hubs have higher outgoing edge degrees
    GNA = tril(G);
    
    for l = 1:sz(1,1)
        for p = 1:l
            if ( GNA(l,p) > 0 )
                GNA(l,p) = 1;
            end
        end
    end
    
else
    % upstream architecture: node hubs have higher incoming edge degrees
    GNA = triu(G);
    
    for l = 1:sz(1,1)
        for p = l:sz(1,1)
            if ( GNA(l,p) > 0 )
                GNA(l,p) = 1;
            end
        end
    end
    
end

% we don't need these vectors any more
clearvars D I sz sz1 sz2 sz3 G2 n;

% clear the rest of the data
clearvars G Gsub nwSz dG;

end

