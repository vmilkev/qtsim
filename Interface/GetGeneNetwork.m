function [ N, A, R, P, s ] = GetGeneNetwork( dataFileName, snpNum, reproduce, rngState, savedNetwork )

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

% read parameters data from *.grm file
delimiter = {''};

formatSpec = '%s%[^\n\r]';

fileID = fopen(dataFileName,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);

fclose(fileID);

param = table(dataArray{1:end-1});

paramSz = size(param);

gnum1 = 0.0;
gnum2 = 0.0;
gnum3 = 0.0;
for i4 = 2:2:paramSz(1,1)
    switch i4
        case 2
            P.type = table2array( param(i4,1) );
            if ( ~strcmp(P.type,'downstream') )
                if ( ~strcmp(P.type,'upstream') )
                    error('not correct NAME for the network ARCHITECTURE \n');
                end
            end
        case 4
            gnum1 = floor( snpNum * str2double( table2array( param(i4,1) ) ) );
            if (gnum1 == 0.0)
                P.gnum_ntw1 = [];
            else
                P.gnum_ntw1 = gnum1;
            end
        case 6
            gnum2 = floor( snpNum * str2double( table2array( param(i4,1) ) ) );
            if (gnum2 == 0.0)
                P.gnum_ntw2 = [];
            else
                P.gnum_ntw2 = gnum2;
            end
        case 8
            gnum3 = floor( snpNum * str2double( table2array( param(i4,1) ) ) );
            if (gnum3 == 0.0)
                P.gnum_ntw3 = [];
            else
                P.gnum_ntw3 = gnum3;
            end
        case 10
            if (gnum1 == 0.0)
                P.beta_ntw1 = [];
            else
                P.beta_ntw1 = str2double( table2array( param(i4,1) ) );
            end
        case 12
            if (gnum2 == 0.0)
                P.beta_ntw2 = [];
            else
                P.beta_ntw2 = str2double( table2array( param(i4,1) ) );
            end
        case 14
            if (gnum3 == 0.0)
                P.beta_ntw3 = [];
            else
                P.beta_ntw3 = str2double( table2array( param(i4,1) ) );
            end
        case 16
            P.strength = str2double( table2array( param(i4,1) ) );
        case 18
            mode = table2array( param(i4,1) );
            if ( strcmp(mode,'strict') )
                P.ar_mode = 'strict';
            elseif ( strcmp(mode,'mix') )
                P.ar_mode = 'mix';
            elseif ( strcmp(mode,'tf') )
                P.ar_mode = 'tf';
            else
                error('not correct NAME for ACTIVATOR|REPRESSOR SETTING of the network \n');
            end
        case 20
            P.ar_val = str2double( table2array( param(i4,1) ) );
            if (P.ar_val < 0.0 || P.ar_val > 1.0)
                error('not correct VALUE for ACTIVATOR|REPRESSOR SETTING of the network \n');
            end
    end
end

N = [];

if (savedNetwork)
    N = ReadNetworkFromFile( savedNetwork );
    oldNtwSize = size(N);
    if (oldNtwSize(1,1) < snpNum)
        error('the number of genes in saved and reused network does not correspond to current number of SNPs >= 1 \n');
    end
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
end

if ( isempty(N) )
    runtime_log( 0, 'ERROR: gene network is empty!' );
    return;
end

A = [];
R = [];

[ A, R ] = SetRegulators( N, P.ar_mode, P.ar_val, s );

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

res_folder = strcat( 'traitsim_results-', date );

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
