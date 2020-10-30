function plot_fullgraph( G, cnodes, view_fig, print_fig, file_name )

if nargin<=3
    file_name = 'Network_architecture.jpeg';
end

totalGenes = size(G,1);

multipl = 2.0;

if (totalGenes > 100)
    multipl = 1.0;
elseif (totalGenes > 500)
    multipl = 0.1;
elseif (totalGenes > 1000)
    multipl = 0.0001;
end

fig = figure('Name', 'Gene Network Architecture', 'visible', view_fig, 'PaperPositionMode', 'auto',...
    'units','normalized','position',[0.2,0.2,0.6,0.7]);
fig.Color = 'w';

dirG = digraph( G );

di = indegree(dirG);
do = outdegree(dirG);
dg = di+do;

nSizes = multipl*sqrt( dg - min( dg ) + 0.2 ) * 2;
nColors = dg;

if (totalGenes < 1000)
    %p = plot(dirG,  'Layout', 'Layered','MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.9);
    p = plot(dirG,'MarkerSize',nSizes,'EdgeAlpha',0.9);
    layout(p,'layered','AssignLayers','alap','Sinks',cnodes);
    highlight(p,cnodes,'NodeColor','r','MarkerSize',mean( dg ) * 4);
else
    p = plot(dirG,'MarkerSize',nSizes,'NodeCData',nColors,'EdgeAlpha',0.9);
end
p.ArrowSize = 12;
title([num2str(totalGenes), ' Gene Network Architecture'],'Interpreter','latex');
%cb = colorbar;
%set( get(cb,'label'), 'string', 'Nodes Degree' );
set(gca, 'XColor', 'none', 'YColor', 'none');

if ( strcmp( print_fig, 'on' ) )
    %print(gcf, '-dpdf', 'Gene_Network_Architecture.pdf');
    name = strcat(file_name, '.jpeg');
    print(gcf, '-djpeg', name);
    close(fig);
end

end