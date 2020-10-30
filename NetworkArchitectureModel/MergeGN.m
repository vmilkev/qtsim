function [G] = MergeGN(...
                       G1,...
                       G2,...
                       G3,...
                       strength,...
                       rngState...
                       )

% Function to merge at least two (max three) networks.
% G(i=1,2,3) - adjacency matrix of undirected network (graph) to merge;
% strength - the parameter which determines the merge strength of combined
%            graphs; its range is [0, +inf) and if strength = 0 merged
%            graphs are completely disconnected, the number of newly created edges
%            between the merging graphs are proportional to strength
%            parameter.
% G - adjacency matrix of merged network (graph).

rng(rngState);

% calculate sizes of adjacency matrices
sz = size(G1);
n1 = sz(1,1);
sz = size(G2);
n2 = sz(1,1);
sz = size(G3);
n3 = sz(1,1);

n = n1 + n2 + n3;

% allocate memory for adjacency matrix of resulting merged graph G
% zeros(n,n) consumes memory but works faster; alternative is to use sparse(n,n)
A = zeros( n,n );

% combining all the adjacency matrices into one

% first graph
for i = 1:n1
    for j = 1:i
        A(i,j) = G1(i,j);
    end
end

% second graph
for i = n1+1:n1+n2
    for j = n1+1:i
        A(i,j) = G2(i-n1,j-n1);
    end
end

% third graph (if not empty)
for i = n1+n2+1:n
    for j = n1+n2+1:i
        A(i,j) = G3(i-n1-n2,j-n1-n2);
    end
end

% merging loop, createing new edges (only) between the merging networks;
% here, a random merge is implemented (nodes degrees are not taken into account)
for i = 1:floor(strength)
    % merge Graph_1 & Graph_2
    if ( (n1 && n2) ~= 0 )
        num1 = randi([1 n1]);
        num2 = randi([n1+1 n1+n2]);
        ri = randi([num1 n1]);    
        rj = randi([num2 n1+n2]);
        A(ri,rj) = 1;
        A(rj,ri) = 1;
    end
    % merge Graph_1 & Graph_3
    if ( (n1 && n3) ~= 0 )
        num1 = randi([1 n1]);
        num2 = randi([n1+n2+1 n]);
        ri = randi([num1 n1]);    
        rj = randi([num2 n]);
        A(ri,rj) = 1;
        A(rj,ri) = 1;
    end
    % merge Graph_2 & Graph_3
    if ( (n2 && n3) ~= 0 )
        num1 = randi([n1+1 n1+n2]);
        num2 = randi([n1+n2+1 n]);
        ri = randi([num1 n1+n2]);    
        rj = randi([num2 n]);
        A(ri,rj) = 1;
        A(rj,ri) = 1;
    end
end

% making a final merged network in terms of adjacency matrix
G = tril(A) + tril(A)' - diag( diag( tril(A) ) );

end

