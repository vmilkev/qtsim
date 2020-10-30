function [ core ] = plot_products( product_files, name_fig, config, numSubplots )

%     config_view.view_fig = 'on';
%     config_view.print_fig = 'on';
%     name_fig = 'Core_products_distribution.jpeg';
%     numSubplots = 6;

% get products from file
P = read_products( product_files );

binNum = ceil(size(P,1)/5);

gNumAll = size(P,2);

core = zeros(gNumAll,1);
for cg = 1:gNumAll
    core(cg,1) = P(1,cg);
end

% we allow only 'numSubplots' sub-plots on one fig.
figsNum = ceil( gNumAll/numSubplots );

iPlot = 1;

for iFig = 1:figsNum
    
    fName = strcat('part',num2str(iFig));
    fName = strcat('_',fName);
    name_fig2 = strcat(name_fig,fName);
    name = strcat(name_fig2, '.jpeg');
    
    gNum = numSubplots;
    
    fig = figure('Name', 'Gene Sub-Network Architecture', 'visible', config.view_fig, 'PaperPositionMode', 'auto',...
        'units','normalized','position',[0.2,0.2,0.6,0.7]);
    
    halfG_1 = ceil( sqrt(gNum) );
    halfG_2 = ceil( gNum/halfG_1 );
    
    iPlot2 = 1;
    for i = 1:halfG_1
        for j = 1:halfG_2
            if ( iPlot <= gNumAll )
                subplot( halfG_1,halfG_2, iPlot2 );
                histogram( P(3:end,iPlot), binNum ,'Normalization','probability' );
                ylabel('probability', 'FontSize',14);
                xlabel('gene product, [molecules]', 'FontSize',14);
                t = title(['gene no. ', num2str( P(1,iPlot) )],'Interpreter','latex');
                t.FontSize = 14;
                set(gca,'FontSize',14)
                iPlot2 = iPlot2 + 1;
                iPlot = iPlot + 1;
            end
        end
    end
    
    if ( strcmp( config.print_fig, 'on' ) )
        print(gcf, '-djpeg', name);
        close(fig);
    end
    
end

end