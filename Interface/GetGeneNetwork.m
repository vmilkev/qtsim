function [ N, A, R, P, s ] = GetGeneNetwork(...
                                            coreFileName,...
                                            rslt_path,...
                                            dataStructName,...
                                            snpNum,...
                                            reproduce,...
                                            rngState,...
                                            savedNetwork,...
                                            savedActivators,...
                                            savedRepressors...
                                           )

% set the state of random number generator
if (reproduce)
    rng(rngState);
    s = rng;
else
    rng('shuffle');
    s = rng;
end

if ~exist('savedNetwork','var')
    savedNetwork = 0;
end

% % read parameters data from *.grm file
% delimiter = {''};
% 
% formatSpec = '%s%[^\n\r]';
% 
% fileID = fopen(dataFileName,'r');
% 
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
% 
% fclose(fileID);
% 
% param = table(dataArray{1:end-1});
% 
% paramSz = size(param);

gnum1 = 0.0;
gnum2 = 0.0;
gnum3 = 0.0;

P = setNtwParamStruct( dataStructName, snpNum );

% for i4 = 2:2:paramSz(1,1)
%     switch i4
%         case 2
%             P.type = table2array( param(i4,1) );
%             if ( ~strcmp(P.type,'downstream') )
%                 if ( ~strcmp(P.type,'upstream') )
%                     error('not correct NAME for the network ARCHITECTURE \n');
%                 end
%             end
%         case 4
%             gnum1 = floor( snpNum * str2double( table2array( param(i4,1) ) ) );
%             if (gnum1 == 0.0)
%                 P.gnum_ntw1 = [];
%             else
%                 P.gnum_ntw1 = gnum1;
%             end
%         case 6
%             gnum2 = floor( snpNum * str2double( table2array( param(i4,1) ) ) );
%             if (gnum2 == 0.0)
%                 P.gnum_ntw2 = [];
%             else
%                 P.gnum_ntw2 = gnum2;
%             end
%         case 8
%             gnum3 = floor( snpNum * str2double( table2array( param(i4,1) ) ) );
%             if (gnum3 == 0.0)
%                 P.gnum_ntw3 = [];
%             else
%                 P.gnum_ntw3 = gnum3;
%             end
%         case 10
%             if (gnum1 == 0.0)
%                 P.beta_ntw1 = [];
%             else
%                 P.beta_ntw1 = str2double( table2array( param(i4,1) ) );
%             end
%         case 12
%             if (gnum2 == 0.0)
%                 P.beta_ntw2 = [];
%             else
%                 P.beta_ntw2 = str2double( table2array( param(i4,1) ) );
%             end
%         case 14
%             if (gnum3 == 0.0)
%                 P.beta_ntw3 = [];
%             else
%                 P.beta_ntw3 = str2double( table2array( param(i4,1) ) );
%             end
%         case 16
%             P.strength = str2double( table2array( param(i4,1) ) );
%         case 18
%             mode = table2array( param(i4,1) );
%             if ( strcmp(mode,'strict') )
%                 P.ar_mode = 'strict';
%             elseif ( strcmp(mode,'mix') )
%                 P.ar_mode = 'mix';
%             elseif ( strcmp(mode,'tf') )
%                 P.ar_mode = 'tf';
%             else
%                 error('not correct NAME for ACTIVATOR|REPRESSOR SETTING of the network \n');
%             end
%         case 20
%             P.ar_val = str2double( table2array( param(i4,1) ) );
%             if (P.ar_val < 0.0 || P.ar_val > 1.0)
%                 error('not correct VALUE for ACTIVATOR|REPRESSOR SETTING of the network \n');
%             end
%     end
% end

% Get the number of core genes 'core_num'.
% If core_num < num. of nodes in N with no out-connections,
% such connections will be created in order to avoid a blind
% non-core genes.

core_num = getCoreNum( coreFileName );

N = [];

if (savedNetwork)
    N = ReadNetworkFromFile( savedNetwork );
    oldNtwSize = size(N);
    if (oldNtwSize(1,1) < snpNum)
        error('the number of genes in saved and reused network does not correspond to current number of SNPs >= 1 \n');
    end
    runtime_log( 0, 'Custom NETWORK file will be used in simulation.' );
else
    N = GetGNA( ...
        P.type,...
        P.gnum_ntw1,...
        P.beta_ntw1,...
        P.gnum_ntw2,...
        P.beta_ntw2,...
        P.gnum_ntw3,...
        P.beta_ntw3,...
        P.strength,...
        s... % state of RNG
        );
    
    if ( size(N,1) > core_num )
        N = rmvBlindGenes( core_num, N );
    end
    
end

if ( isempty(N) )
    runtime_log( 0, 'ERROR: gene network is empty!' );
    return;
end

A = [];
R = [];

if (savedActivators & savedRepressors)
    A = ReadNetworkFromFile( savedActivators );
    oldNtwSize = size(A);
    if (oldNtwSize(1,1) < snpNum)
        error('the number of genes in saved and reused activators file does not correspond to the current number of SNPs >= 1 \n');
    end
    R = ReadNetworkFromFile( savedRepressors );
    oldNtwSize = size(R);
    if (oldNtwSize(1,1) < snpNum)
        error('the number of genes in saved and reused activators file does not correspond to the current number of SNPs >= 1 \n');
    end
    runtime_log( 0, 'Custom ACTIVATORS and REPRESSORS network files will be used in simulation.' );
    if (~savedNetwork)
        runtime_log( 0, 'WARNING! The generated NETWORK might not be consistent with custom ACTIVATORS and REPRESSORS networks.' );
        runtime_log( 0, '........ Check if A & R are the sub-sets of N!' );
    end
elseif (savedActivators & ~savedRepressors)
    A = ReadNetworkFromFile( savedActivators );
    oldNtwSize = size(A);
    if (oldNtwSize(1,1) < snpNum)
        error('the number of genes in saved and reused activators file does not correspond to the current number of SNPs >= 1 \n');
    end
    [ ~, R ] = SetRegulators( N, P.ar_mode, P.ar_val, s );
    runtime_log( 0, 'Custom ACTIVATORS network file will be used in simulation, REPRESSORS network files will be generated.' );
    if (~savedNetwork)
        runtime_log( 0, 'WARNING! The generated NETWORK might not be consistent with custom ACTIVATORS network.' );
        runtime_log( 0, '........ Check if A is the sub-set of N!' );
    end
elseif (~savedActivators & savedRepressors)
    R = ReadNetworkFromFile( savedRepressors );
    oldNtwSize = size(R);
    if (oldNtwSize(1,1) < snpNum)
        error('the number of genes in saved and reused activators file does not correspond to the current number of SNPs >= 1 \n');
    end
    [ A, ~ ] = SetRegulators( N, P.ar_mode, P.ar_val, s );
    runtime_log( 0, 'Custom REPRESSORS network file will be used in simulation, ACTIVATORS network files will be generated.' );
    if (~savedNetwork)
        runtime_log( 0, 'WARNING! The generated NETWORK might not be consistent with custom REPRESSORS network.' );
        runtime_log( 0, '........ Check if R is the sub-set of N!' );
    end
else
    [ A, R ] = SetRegulators( N, P.ar_mode, P.ar_val, s );
    runtime_log( 0, 'Both ACTIVATORS and REPRESSORS network files will be generated.' );
end

if ( isempty(A) )
    runtime_log( 0, 'ERROR: activators network is empty!' );
    return;
end

if ( isempty(R) )
    runtime_log( 0, 'ERROR: repressors network is empty!' );
    return;
end

ntwSize = size(N);

P.genes = ntwSize(1,1);

% write N, A, R to file

if ( ~rslt_path )
    res_folder = strcat( 'mescot_results-', date );
else
    res_folder = rslt_path;
end

if ( ~MakeFolder(res_folder) )
    return;
end

fName = strcat( res_folder, '/N.gntw' );
[row, col, v] = find(N);
row( size(row,1)+1, 1 ) = ntwSize(1,1);
col( size(col,1)+1, 1 ) = ntwSize(1,1);
v( size(v,1)+1, 1 ) = 0;
dlmwrite(fName,[row col v], 'delimiter', '\t');

fName = strcat( res_folder, '/A.gntw' );
[row, col, v] = find(A);
row( size(row,1)+1, 1 ) = ntwSize(1,1);
col( size(col,1)+1, 1 ) = ntwSize(1,1);
v( size(v,1)+1, 1 ) = 0;
dlmwrite(fName,[row col v], 'delimiter', '\t');

fName = strcat( res_folder, '/R.gntw' );
[row, col, v] = find(R);
row( size(row,1)+1, 1 ) = ntwSize(1,1);
col( size(col,1)+1, 1 ) = ntwSize(1,1);
v( size(v,1)+1, 1 ) = 0;
dlmwrite(fName,[row col v], 'delimiter', '\t');

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function coreNum = getCoreNum( coreFileName )

delimiter = ' ';

formatSpec = '%s%s%[^\n\r]';

fileID = fopen(coreFileName,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);

fclose(fileID);

param = table(dataArray{1:end-1});
paramSz = size(param);

coreNum = paramSz(1,1);

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function net = rmvBlindGenes( coreGenes, net )
    
    genes = size(net,1);
    
    if (genes <= coreGenes)
        return;
    end
    
    for i = coreGenes+1:genes
        blind = ~any( net(i,:) );
        if ( blind )
            j = 1 + (i-2)*rand(1,1);
            net( i,floor(j) ) = 1;
%         else
%             return;
        end
    end
    
end
%---------------------------------------------------------------

%---------------------------------------------------------------
function [ P ] = setNtwParamStruct( ntwpar,snpNum )

P.type = ntwpar.arch;
if ( ~strcmp(P.type,'downstream') )
    if ( ~strcmp(P.type,'upstream') )
        error('not correct NAME for the network ARCHITECTURE \n');
    end
end

gnum1 = floor( snpNum * ntwpar.n1 );
if (gnum1 == 0.0)
    P.gnum_ntw1 = [];
else
    P.gnum_ntw1 = gnum1;
end

gnum2 = floor( snpNum * ntwpar.n2 );
if (gnum2 == 0.0)
    P.gnum_ntw2 = [];
else
    P.gnum_ntw2 = gnum2;
end

gnum3 = floor( snpNum * ntwpar.n3 );
if (gnum3 == 0.0)
    P.gnum_ntw3 = [];
else
    P.gnum_ntw3 = gnum3;
end

if (gnum1 == 0.0)
    P.beta_ntw1 = [];
else
    P.beta_ntw1 = ntwpar.b1;
end

if (gnum2 == 0.0)
    P.beta_ntw2 = [];
else
    P.beta_ntw2 = ntwpar.b2;
end

if (gnum3 == 0.0)
    P.beta_ntw3 = [];
else
    P.beta_ntw3 = ntwpar.b3;
end

P.strength = ntwpar.mrg;

mode = ntwpar.mode;
if ( strcmp(mode,'strict') )
    P.ar_mode = 'strict';
elseif ( strcmp(mode,'mix') )
    P.ar_mode = 'mix';
elseif ( strcmp(mode,'tf') )
    P.ar_mode = 'tf';
else
    error('not correct NAME for ACTIVATOR|REPRESSOR SETTING of the network \n');
end

P.ar_val = ntwpar.ratio;
if (P.ar_val < 0.0 || P.ar_val > 1.0)
    error('not correct VALUE for ACTIVATOR|REPRESSOR SETTING of the network \n');
end

end
