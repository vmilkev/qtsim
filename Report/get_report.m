function get_report(datapath, inparam )

resname1 = strcat( datapath, '/');
resname = strcat( resname1, 'Report' );

status = mkdir( resname );
if ( ~status )
    runtime_log( 0, 'ERROR: cannot create Report directory!' );
    return;
end

config_view.view_fig = 'off';
config_view.print_fig = 'on';

% DEFAULT option
if ( inparam.report_def )
    
    %--------------------------------------------------------------------------
    
    network = strcat(datapath,'/N_snpID.gntw');
    networkA = strcat(datapath,'/A_snpID.gntw');
    networkR = strcat(datapath,'/R_snpID.gntw');
    %--------------------------------------------------------------------------
    % display full network
    config_view.view = 'all';
    plot_ntw( network, config_view, strcat(resname,'/Complete_network')  );
    plot_ntw( networkA, config_view, strcat(resname,'/Activators_sub-network')  );
    plot_ntw( networkR, config_view, strcat(resname,'/Repressors_sub-network')  );
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % display network degree distribution
    config_view.view = 'indegree';
    plot_ntw( network, config_view, strcat(resname,'/in-Connectivity_Distribution') );
    
    config_view.view = 'outdegree';
    plot_ntw( network, config_view, strcat(resname,'/out-Connectivity_Distribution') );
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % display products distribution
    pFiles = FindParFile2( datapath, 'products' );
    for iP = 1:size(pFiles,1)
        product{1,iP} =  strcat(resname1,pFiles{iP,1}) ;
        %product{1,iP} = pFiles{iP,1};
    end
    core = plot_products( product, strcat(resname,'/Core_products_distribution'), config_view, 20 );
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % display sub-network for core genes
    nameFig.red.title = ' Sub-Network';
    nameFig.red.file = strcat(resname,'/Core_Genes_Sub-Network');
    
    config_view.view = 'sub2';
    config_view.which_plot = 'reduced'; % network representation mode
    config_view.node = core'; % selected genes
    config_view.orient = 'both'; % type of connections
    
    plot_ntw( network, config_view, nameFig, min( 20,size(core,1) )  );
    
    %--------------------------------------------------------------------------
    
    %--------------------------------------------------------------------------
    % display sub-network for core genes
    nameFig.full.title = ' Sub-Network';
    nameFig.full.file = strcat(resname,'/Core_Genes_Sub-Network_in_Whole_Architecture');
    
    config_view.view = 'sub2';
    config_view.which_plot = 'full'; % network representation mode
    config_view.node = core'; % selected genes
    config_view.orient = 'both'; % type of connections
    
    plot_ntw( network, config_view, nameFig, min( 4,size(core,1) )  );
    %--------------------------------------------------------------------------
end

% GENE option
if ( inparam.report_genS )
    
    % display sub-network for particular genes
    nameFig.red.title = ' Sub-Network';
    nameFig.red.file = strcat(resname,'/Selected_Genes_Sub-Network');
    nameFig.full.title = ' Gene Sub-Network within ';
    nameFig.full.file = strcat(resname,'/Selected_Genes_Sub-Network_Within_Whole_Architecture');
    
    config_view.view = 'sub2';
    config_view.which_plot = 'reduced'; % network representation mode
    config_view.node = inparam.report_genV; % selected genes
    config_view.orient = 'both'; % type of connections
    plot_ntw( network, config_view, nameFig, min( 16,size(config_view.node,2) )  );
    
    config_view.which_plot = 'full'; % network representation mode
    plot_ntw( network, config_view, nameFig, min( 4,size(config_view.node,2) )  );
    
end

% ID option
if ( inparam.report_idS )
    
    id = inparam.report_idVi';
    gen = inparam.report_idVg';
    res_folder = {datapath};
    prod_type = 'rna';
    
    normalize = 'on';
    
    plot_productsT( res_folder, id, gen, prod_type, normalize, strcat(resname,'/selected_products_dynamics'), config_view, min( 16,size(id,1) ) );
    
end

end
