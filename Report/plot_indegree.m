function plot_indegree( G, view_fig, print_fig, name_fig  )

    dirG = digraph( G );
    
    fig = figure('Name', 'in-Connectivity Distribution', 'visible', view_fig, 'PaperPositionMode', 'auto',...
        'units','normalized','position',[0.2,0.2,0.6,0.7]);
    
    di = indegree(dirG);
    
    [xi, fi] = degree_distr(di, 'in_degree_distr.dat', false);
    
    [ai,bi,rsquarei,xPredi,yPredi] = fit_indegree(xi, fi);
        
    mrkSz = 9;
    line = 2;
    col = 'k';
    
    loglog(xi,fi, 'o','MarkerSize',mrkSz, 'MarkerEdgeColor','y','MarkerFaceColor',col);
    hold on;
    loglog(xPredi,yPredi, col, 'LineWidth', line);
    lgd=legend('data','fit','Location', 'Best' );
    legend('boxoff');
    title(lgd,{ ['exponential fit: ${y} = {a} \cdot {exp}({b} \cdot x)$'] ['a = ', num2str(ai), '; b = ', num2str(bi)] [ '${R}^2$ = ', num2str(rsquarei)] },'Interpreter','latex');
    ylabel('genes','Interpreter','latex');
    xlabel('TF per gene','Interpreter','latex');
    title('in-Connectivity distribution (log/log scale)','Interpreter','latex');
    ylim([min(fi) inf]);
        
    if ( strcmp( print_fig, 'on' ) )
        %print(gcf, '-dpdf', 'Gene_Network_Architecture.pdf');
        name = strcat(name_fig, '.jpeg');
        print(gcf, '-djpeg', name);
        close(fig);
    end
    %----------------------------------------------------------------------------------------------------------------------
   
end