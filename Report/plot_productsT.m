function plot_productsT( res_folder, id, genes, product, norm, name_fig, config, numSubplots )

gNumAll = size(id,1);


% we allow only 'numSubplots' sub-plots on one fig.
figsNum = ceil( gNumAll/numSubplots );

iPlot = 1;

for iFig = 1:figsNum
    
    fName = strcat('part',num2str(iFig));
    fName = strcat('_',fName);
    name_fig2 = strcat(name_fig,fName);
    name = strcat(name_fig2, '.jpeg');

    
    gNum = numSubplots;
    
    fig = figure('Name', 'Genes Products Dynamics', 'visible', config.view_fig, 'PaperPositionMode', 'auto',...
        'units','normalized','position',[0.2,0.2,0.6,0.7]);
    
    halfG_1 = ceil( sqrt(gNum) );
    halfG_2 = ceil( gNum/halfG_1 );
    
    iPlot2 = 1;
    for i = 1:halfG_1
        for j = 1:halfG_2
            if ( iPlot <= gNumAll )
                                
                sp = subplot( halfG_1,halfG_2, iPlot2 );
                spW = sp.Position(3);
                spX = sp.Position(1);
                %sp.Position(3) = width * 0.9;
                hold on;
                               
                %----------------------------------------------------------------------
                % get products files
                folder = strcat('/', num2str( id(iPlot,1) ));
                folder = strcat( res_folder, folder );
                t_name = strcat(folder, '/time.txt');
                r_name = strcat(folder, '/rna.txt');
                p_name = strcat(folder, '/p.txt');
                
                t = read_products( t_name );
                
                if ( strcmp( product, 'rna' ) )
                    p = read_products( r_name );
                    pLabel = 'mRNA';
                elseif ( strcmp( product, 'prot' ) )
                    p = read_products( p_name );
                    pLabel = 'Protein';
                else
                    runtime_log( 0, 'ERROR: not correct option in plot_productsT() function!' );
                    return;
                end
                %----------------------------------------------------------------------
                
                if ( strcmp(norm,'on') )
                max_p = max(p(3:end,genes(:)));

                h = plot( t./max(t), p(3:end,genes)./max(max_p), '.-' );
                box on;
                %lgd = legend({num2str(genes)}, 'Location', 'Best');
                lgd = legend({num2str(genes)});
                lX = lgd.Position(1);
                lW = lgd.Position(3);
                sp.Position(3) = spW - lW * 0.6;
                lgd.Position(1) = spX + spW - lW * 0.4;                
                %lgd.Position(1) = xPos * 1.1;
                legend('boxoff');
                title(lgd,'Genes:');
                ylabel(pLabel);
                xlabel('Time');
                elseif (strcmp(norm,'off'))
                h = plot( t, p(3:end,genes), '.-' );
                box on;
                %lgd = legend({num2str(genes)}, 'Location', 'Best');
                lgd = legend({num2str(genes)});
                lX = lgd.Position(1);
                lW = lgd.Position(3);
                sp.Position(3) = spW - lW * 0.6;
                lgd.Position(1) = spX + spW - lW * 0.4;                
                legend('boxoff');
                title(lgd,'Genes:');
                pLabel = 'mRNA, [1/min]';
                ylabel(pLabel);
                xlabel('Time, [min]');
                else
                    runtime_log( 0, 'ERROR: not correct option in plot_productsT() function => normalization option!' );
                    return;
                end                    
                title(['Id no. ', num2str( id(iPlot,1) )],'Interpreter','latex');
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