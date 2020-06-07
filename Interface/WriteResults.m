function WriteResults( rslt_path, rngState, trait, products, coreMap, ncoreMap, genot_id, dynSol, genes, runAll, grange )

if ~exist('grange','var')
    grange.start = 0;
end

% rows: diff. genotypes; col.: (additive) trait for corresponding genotype
res_trait = zeros( size(genot_id,2), 1+1);

% rows: gene's products; cols: (1) snp ID, (2) corresponding node ID,
% (3)...(end) diff. genotypes (in consecutive order)
res_prod = zeros( size(products,1), size(genot_id,2)+2 );

res_trait(:,1) = genot_id(1,:)';
res_trait(:,2) = trait(:,1);

node = cell2mat( keys( coreMap ) )';
snp = cell2mat( values( coreMap ) )';
res_prod(:,1) = snp(:,1);
res_prod(:,2) = node(:,1);
res_prod(:,3:size(products,2)+2) = products(:,1:end);

if ( ~rslt_path )
    res_folder = strcat( 'traitsim_results-', date );
else
    res_folder = rslt_path;
end

if ( ~MakeFolder(res_folder) )
    return;
end

% if (grange.start)
%     iter_name1 = num2str( grange.start );
%     iter_name2 = num2str( grange.end );
%     name0 = strcat( iter_name1,'-' );
%     name0 = strcat( name0,iter_name2 );
%     name0 = strcat( name0,'.txt' );
%     name0 = strcat( '/traits_',name0 );
%     fName = strcat( res_folder,name0);
% else
    fName = strcat( res_folder, '/traits.txt' );
% end

dlmwrite( fName,res_trait,'-append','Delimiter',' ' );

if (grange.start)
    iter_name1 = num2str( grange.start );
    iter_name2 = num2str( grange.end );
    name0 = strcat( iter_name1,'-' );
    name0 = strcat( name0,iter_name2 );
    name0 = strcat( name0,'.txt' );
    name0 = strcat( '/products_',name0 );
    fName = strcat( res_folder,name0);
else
    fName = strcat( res_folder, '/products.txt' );
end
dlmwrite( fName,res_prod,'Delimiter',' ' ); 

fName = strcat( res_folder, '/randomstate.txt' );
fileID = fopen(fName,'w');
fprintf(fileID,'%d\n',rngState.Seed);
fclose(fileID);

if ( runAll )
    loop_sz = size(genot_id,2);
else
    loop_sz = 1;
end

for igt = 1:loop_sz
    solX = dynSol(igt).x;
    Y_rna = zeros( (size(dynSol(igt).y,1)-1)/2, size(dynSol(igt).y,2)+2 );
    Y_prt = zeros( (size(dynSol(igt).y,1)-1)/2, size(dynSol(igt).y,2)+2 );
    
    % core genes
    k = keys( coreMap );
    for ik = 1:size(k,2)
        ir = k{ik}; % node mRNA
        ip = k{ik} + genes; % node Protein
        Y_rna(ik,1) = coreMap(ir);
        Y_prt(ik,1) = coreMap(ir);
        Y_rna(ik,2) = ir;
        Y_prt(ik,2) = ir;
        Y_rna(ik,3:end) = dynSol(igt).y(ir,:);
        Y_prt(ik,3:end) = dynSol(igt).y(ip,:);
    end
    
    % non-core genes
    k2 = keys( ncoreMap );
    for ik2 = 1:size(k2,2)
        ir = k2{ik2}; % node mRNA
        ip = k2{ik2} + genes; % node Protein
        Y_rna(ik+ik2,1) = ncoreMap(ir);
        Y_prt(ik+ik2,1) = ncoreMap(ir);
        Y_rna(ik+ik2,2) = ir;
        Y_prt(ik+ik2,2) = ir;
        Y_rna(ik+ik2,3:end) = dynSol(igt).y(ir,:);
        Y_prt(ik+ik2,3:end) = dynSol(igt).y(ip,:);
    end
    
    name_i = num2str( genot_id(1,igt) );
    fname = strcat( res_folder,'/');
    if ( runAll )
        fname = strcat( fname,name_i);
    else
        fname = strcat( fname, 'ref_genotype' );
    end

    if ( ~MakeFolder(fname) )
        return;
    end
    
    name1 = '/time.txt';
    fileName = strcat( fname, name1 );
    dlmwrite( fileName,solX,'Delimiter',' ' );

    name1 = '/rna.txt';
    fileName = strcat( fname, name1 );
    dlmwrite( fileName,Y_rna,'Delimiter',' ' );

    name1 = '/p.txt';
    fileName = strcat( fname, name1 );
    dlmwrite( fileName,Y_prt,'Delimiter',' ' );

end

end