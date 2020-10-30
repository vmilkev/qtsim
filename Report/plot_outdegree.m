function plot_outdegree( G, view_fig, print_fig, name_fig  )

    dirG = digraph( G );
    
    fig = figure('Name', 'out-Connectivity Distribution', 'visible', view_fig, 'PaperPositionMode', 'auto',...
        'units','normalized','position',[0.2,0.2,0.6,0.7]);
    
    do = outdegree(dirG);
    [xo, fo] = degree_distr(do, 'out_degree_distr.dat', false);
    
    [ao,bo,rsquareo,xPredo,yPredo] = fit_outdegree(xo, fo);
        
    mrkSz = 9;
    line = 2;
    col = 'k';
    
    loglog(xo,fo, 'o','MarkerSize',mrkSz, 'MarkerEdgeColor','y','MarkerFaceColor',col);
    hold on;
    loglog(xPredo,yPredo, col, 'LineWidth', line);
    lgd=legend('data','fit','Location', 'Best' );
    legend('boxoff');
    title(lgd,{ ['power-law fit: ${y} = {a} \cdot x^{b}$'] ['a = ', num2str(ao), '; b = ', num2str(bo)] ['${R}^2$ = ', num2str(rsquareo)] },'Interpreter','latex');
    ylabel('TF','Interpreter','latex');
    xlabel('regulated genes per TF','Interpreter','latex');
    title('out-Connectivity distribution (log/log scale)','Interpreter','latex');
    ylim([min(fo) inf]);
    
    if ( strcmp( print_fig, 'on' ) )
        %print(gcf, '-dpdf', 'Gene_Network_Architecture.pdf');
        name = strcat(name_fig, '.jpeg');
        print(gcf, '-djpeg', name);
        close(fig);
    end
    
end