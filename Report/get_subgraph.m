function [ labels ] = get_subgraph( G, node )

% allocate memory for sub-graph adjacency matrix
N = zeros( size(G,1) );

% first list of col/row is for target (input) node
nodes = node;

N = subgraph( G, N, nodes );

nodes_list = 0;
nodes_index = 0;
nodes_list2 = 0;
nodes_index2 = 0;
for i = 1:size(N,1)
    k = find(N(:,i));
    p = find(N(i,:));
    if ( ~isempty(k) )
        nodes_index = nodes_index + 1;
        nodes_list(nodes_index,1) = i;
    end
    if ( ~isempty(p) )
        nodes_index2 = nodes_index2 + 1;
        nodes_list2(nodes_index2,1) = i;
    end
end

red_nodes = unique( cat(1,nodes_list,nodes_list2) );

%N = N(red_nodes,red_nodes);
labels = red_nodes;

end

function N = subgraph( G, N, node )

    if ( isempty(node) )
        return;
    end
    
    nodes = node';
    
    for j = 1:size(nodes,1)
        
        col = nodes(j,1);
        row = G(:,col);
        N( :,col ) = row;
        new_nodes = find( row );
        
        N = subgraph( G, N, new_nodes );
        
    end

end