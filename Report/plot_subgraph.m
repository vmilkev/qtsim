function plot_subgraph( G, node, orient, which_plot, view_fig, print_fig, name_fig )

if ( isempty(name_fig) )
    name_fig.full.title = ' (selected) Genes Sub-Network within ';
    name_fig.full.file = 'Gene_Sub-Network_Within_Whole_Architecture.jpg';
    name_fig.red.title = ' (selected) Genes Sub-Network Architecture';
    name_fig.red.file = 'Gene_Sub-Network_Architecture.jpg';
end

if ( strcmp( orient, 'out' ) )
    [ lb ] = get_subgraph( G', node );
elseif ( strcmp( orient, 'in' ) )
    [ lb ] = get_subgraph( G, node );
elseif ( strcmp( orient, 'both' ) )
    [ lb1 ] = get_subgraph( G', node );
    [ lb2 ] = get_subgraph( G, node );
    lb = [ lb1; lb2 ];
    lb = unique( lb( lb ~= 0 ) );
else
    runtime_log( 0, 'ERROR: not correct option in plot_subgraph() function!' );
    return;
end

totalGenes = size(G,1);

fig = figure('Name', 'Gene Sub-Network Architecture', 'visible', view_fig, 'PaperPositionMode', 'auto',...
        'units','normalized','position',[0.2,0.2,0.6,0.7]);

set(gca, 'XColor', 'none', 'YColor', 'none');

if ( lb ~= 0 )
    
    gr_full = digraph(G);

    gr_red = digraph( G(lb,lb) );
    nLabels = cellstr( string( lb ) )';

    eWeight = ones( size(gr_full.Edges.EndNodes,1),1 )*0.5;
    weight_coef = 2.0;

    eColor = zeros( size(gr_full.Edges.EndNodes,1),3 );
    eColor(:,3) = 1.0; % make default color BLUE

    % node size and color
    nWeight = ones( size(G,1),1 )*3.0;
    nWeight(lb,1) = 5;
    nColor = zeros( size(G,1),3 );
    nColor(:,3) = 1.0; % make default color BLUE
    nColor(lb,1) = 1.0; % highlight subgraph by red nodes
    nColor(lb,3) = 0.0;


    for i = 1:size(lb,1)
       node = lb(i,1);
       k = find( gr_full.Edges.EndNodes(:,1) == node );
       if ( ~isempty(k) )
          for j = 1:size(k,1)
              if ( any( lb == gr_full.Edges.EndNodes( k(j,1),2 ) ) )
                eWeight( k(j,1) ) = weight_coef;
                eColor( k(j,1),3 ) = 0.0;
                eColor( k(j,1),1 ) = 1.0; % red color
              end
          end
       end
    end

    gr_full.Edges.Weight = eWeight;
    gr_full.Edges.LWidths = gr_full.Edges.Weight;
    gr_full.Edges.EdgeColors = eColor;
    gr_full.Nodes.NodeData = nWeight;
    gr_full.Nodes.NodeColor = nColor;

    if ( strcmp( which_plot,'full' ) )
        p_full = plot( gr_full,'MarkerSize',gr_full.Nodes.NodeData, 'LineWidth',gr_full.Edges.LWidths );
        p_full.EdgeColor = gr_full.Edges.EdgeColors;
        p_full.NodeColor = gr_full.Nodes.NodeColor;
        p_full.ArrowSize = 12;
        title([num2str( size(G(lb,lb),1) ), name_fig.full.title, num2str(totalGenes), ' Genes Architecture'],'Interpreter','latex');
        
    if ( strcmp( print_fig, 'on' ) )
        print(gcf, '-djpeg', name_fig.full.file);
    end

    end

    if ( strcmp( which_plot,'reduced' ) )
        
        di = indegree(gr_red);
        do = outdegree(gr_red);
        dg = di+do;
        nSizes = 2*sqrt( dg - min( dg ) + 0.2 );
        nColors = dg;

        p_red = plot( gr_red, 'NodeLabel',nLabels, 'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.9 );
        p_red.ArrowSize = 10;
%         p_red.EdgeColor = 'r';
%         p_red.NodeColor = 'black';
        title([num2str( size(G(lb,lb),1) ), name_fig.red.title],'Interpreter','latex');
        cb = colorbar;
        set( get(cb,'label'), 'string', 'Nodes Degree' );
       
    if ( strcmp( print_fig, 'on' ) )
        name = strcat(name_fig.red.file, '.jpeg');
        print(gcf, '-djpeg', name);
    end

    end

    
    
end

    if ( strcmp( print_fig, 'on' ) )
        close(fig);
    end

end