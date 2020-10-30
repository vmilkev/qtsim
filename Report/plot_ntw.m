function plot_ntw( ntwName, config, name_fig, subplots  )

% Input parameters.
%
% ntwName - network file name;
%
% config.view = 'all' | 'sub' | 'degree';
%               'all' - draw full network architecture;
%               'sub' - draw a sub-network for particular gene(s);
%               'sub2' - draw a sub-network for particular gene(s) as sub-plots;
%               'degree' - draw in|out degree distributions;
%
% the rest of 'config' structure is elevant only if config.view == 'sub':
%
%   config.node = [ genes_ID(s) ];
%
%   config.orient = 'in' | 'out' | 'both';
%               'in' - sub-network with in-coming connections;
%               'out' - sub-network with out-coming connections;
%               'both' - full sub-network (with in- & out connections);
%
%   config.which_plot = 'full' | 'reduced';
%               'full' - drawing a highlited sub-network within a full genomic network;
%               'reduced' - drawing a sub-network as separate figure;
%
%   config.view_fig = 'on' | 'off';
%
%   config.print_fig = 'on' | 'off';

if nargin<=2
    name_fig = [];
    subplots = [];
end
    G = read_ntw( ntwName );
    
    if ( strcmp( config.view, 'all' ) )
        
        plot_fullgraph( G, config.node, config.view_fig, config.print_fig, name_fig );
        
    elseif ( strcmp( config.view, 'sub' ) )
        
        plot_subgraph( G, config.node, config.orient, config.which_plot, config.view_fig, config.print_fig, name_fig);
        
    elseif ( strcmp( config.view, 'sub2' ) )
        
        plot_subgraph2( G, config.node, config.orient, config.which_plot, config.view_fig, config.print_fig, name_fig, subplots );
        
    elseif ( strcmp( config.view, 'indegree' ) )
        
        plot_indegree( G, config.view_fig, config.print_fig, name_fig  );
        
    elseif ( strcmp( config.view, 'outdegree' ) )
        
        plot_outdegree( G, config.view_fig, config.print_fig, name_fig  );
        
    else
        runtime_log( 0, 'ERROR: not correct option in plot_ntw() function!' );
        return;
        
    end
        
end