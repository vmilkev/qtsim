function [ cM, ncM, cW ] = GetCoreNonCoreMap( coreFileName, ntwGenes, ntwArch )

% read core IDs from file
[ core_snp, core_snp_num, core_w ] = getCoreID(coreFileName);

if ( core_snp_num > ntwGenes )
    error('the number of genes in the used network isdoc dbstack lower than the number of core genes in the "core.dat" file \n');
end

[ n_core_snp, ncore_snp_num ] = getNonCoreID( core_snp, ntwGenes );

% map core IDs to the Graph (Network) nodes

% (1) create a vector of Graph nodes correspondent to core genes (or SNPs)

if ( strcmp(ntwArch, 'downstream') )
    % The Graph nodes starts from core genes (for "Downstream architecture" !?)
    % The last nodes (with highest IDs) are always TF.
    for ci = 1:core_snp_num
        coreNodes(ci) = ci;
    end
    % (2) the rest of nodes are TFs
    cIndex = 0;
    for ci = core_snp_num+1:ntwGenes
        cIndex = cIndex + 1;
        n_coreNodes(cIndex) = ci;
    end
else
    range_val = 1;
    for ci = ntwGenes:-1:ntwGenes-core_snp_num+1
        coreNodes(range_val) = ci;
        range_val = range_val + 1;
    end
    % (2) the rest of nodes are TFs
    cIndex = 0;
    for ci = ntwGenes-core_snp_num:-1:1
        cIndex = cIndex + 1;
        n_coreNodes(cIndex) = ci;
    end
end
cM = containers.Map(coreNodes,core_snp);
cW = containers.Map(core_snp,core_w);
ncM = containers.Map(n_coreNodes,n_core_snp);

end
%---------------------------------------------------------------


%---------------------------------------------------------------
% Local Functions
%---------------------------------------------------------------
function [ noncID, noncoreNum ] = getNonCoreID( core_id, ntwGenes )

index = 0;
for i6 = 1:ntwGenes
    if ( isempty( find( core_id == i6 ) ) )
        index = index + 1;
        noncID(index) = i6;
    end
end

sz = size(noncID);
noncoreNum = sz(1,2);

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function [ coreID, coreNum, coreW ] = getCoreID(dataFileName)

delimiter = ' ';

formatSpec = '%s%s%[^\n\r]';

fileID = fopen(dataFileName,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);

fclose(fileID);

param = table(dataArray{1:end-1});
paramSz = size(param);

for i5 = 1:paramSz(1,1)
    coreID(i5) = str2double( table2array( param(i5,1) ) );
    coreW(i5) = str2double( table2array( param(i5,2) ) );
end

coreNum = paramSz(1,1);

end
%---------------------------------------------------------------