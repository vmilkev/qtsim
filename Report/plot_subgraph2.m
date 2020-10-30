function plot_subgraph2( G, node, orient, which_plot, view_fig, print_fig, name_fig, numSubplots )

if ( isempty(name_fig) )
    name_fig.full.title = ' (selected) Genes Sub-Network within ';
    name_fig.full.file = 'Gene_Sub-Network_Within_Whole_Architecture.jpg';
    name_fig.red.title = ' (incoming)Sub-Network';
    name_fig.red.file = 'Gene_(incoming)Sub-Network.jpg';
end

gNumAll = size(node,2);

% we allow only 'numSubplots' sub-plots on one fig.
figsNum = ceil( gNumAll/numSubplots );

iPlot = 1;

for iFig = 1:figsNum
    
    fName = strcat('part',num2str(iFig));
    fName = strcat('_',fName);
    if ( strcmp( which_plot,'full' ) )
        name_fig2 = strcat(name_fig.full.file,fName);
    elseif ( strcmp( which_plot,'reduced' ) )
        name_fig2 = strcat(name_fig.red.file,fName);
    end   
    name = strcat(name_fig2, '.jpeg');

    
    gNum = numSubplots;
    
    halfG_1 = ceil( sqrt(gNum) );
    halfG_2 = ceil( gNum/halfG_1 );
    
    fig = figure('Name', 'Gene Sub-Network Architecture', 'visible', view_fig, 'PaperPositionMode', 'auto',...
        'units','normalized','position',[0.2,0.2,0.6,0.7]);
    
    
    iPlot2 = 1;
    for i = 1:halfG_1
        for j = 1:halfG_2
            if ( iPlot <= gNumAll )
                subplot( halfG_1,halfG_2, iPlot2 );
                
                %-----------------------------------------------------------
                if ( strcmp( orient, 'out' ) )
                    [ lb ] = get_subgraph( G', node(1,iPlot) );
                elseif ( strcmp( orient, 'in' ) )
                    [ lb ] = get_subgraph( G, node(1,iPlot) );
                elseif ( strcmp( orient, 'both' ) )
                    [ lb1 ] = get_subgraph( G', node(1,iPlot) );
                    [ lb2 ] = get_subgraph( G, node(1,iPlot) );
                    lb = [ lb1; lb2 ];
                    lb = unique( lb( lb ~= 0 ) );
                else
                    runtime_log( 0, 'ERROR: not correct option in plot_subgraph() function!' );
                    return;
                end
                
                totalGenes = size(G,1);
                
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
                        node1 = lb(i,1);
                        k = find( gr_full.Edges.EndNodes(:,1) == node1 );
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
                        title(['Gene no. ', num2str( node(1,iPlot) ), name_fig.full.title],'Interpreter','latex');
                    end
                    
                    if ( strcmp( which_plot,'reduced' ) )
                        
                        di = indegree(gr_red);
                        do = outdegree(gr_red);
                        dg = di+do;
                        %nSizes = 2*sqrt( dg - min( dg ) + 0.2 );
                        nSizes = (dg./dg) * mean( dg ) * 2;
                        nColors = dg;
                        %p_red = plot( gr_red, 'NodeLabel',nLabels, 'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.9 );
                        p_red = plot( gr_red, 'NodeLabel',nLabels, 'MarkerSize',nSizes,'EdgeAlpha',0.9 );
                                                
                        %p_red = plot(gr_red,'MarkerSize',nSizes,'EdgeAlpha',0.9);
                        layout(p_red,'layered','AssignLayers','alap');
                        highlight(p_red,find(lb == node(1,iPlot)),'NodeColor','r','MarkerSize',mean( dg ) * 3);

                        p_red.ArrowSize = 11;
                        title(['Gene no. ', num2str( node(1,iPlot) ), name_fig.red.title],'Interpreter','latex');
                    end
                    
                    set(gca, 'XColor', 'none', 'YColor', 'none');
                    
                end
                iPlot2 = iPlot2 + 1;
                iPlot = iPlot + 1;
            end
        end
    end
    
    if ( strcmp( print_fig, 'on' ) )
        print(gcf, '-djpeg', name);
    end
    
end

if ( strcmp( print_fig, 'on' ) )
    close(fig);
end

end