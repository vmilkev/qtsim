function RunTraitSimulator( simSetUpFile )

runtime_log( 1, 'STARTING TraitSim' );

if ( isempty(simSetUpFile) )
    fParam = FindParFile( pwd, '.tsim' );
else
    fParam = simSetUpFile;
end

if ( isempty(fParam) )
    runtime_log( 0, 'ERROR: program cannot find "*.tsim" file!' );
    return;
else
    runtime_log( 0, 'TraitSim parameters file name..........:', fParam );
end

%par = ReadTraitSimParam( fParam );
par = fReadSimPar( fParam );

if ( isempty(par.allele) )
    runtime_log( 0, 'ERROR: there is no SNPs file name in "*.tsim" file, otherwise there are problems while opening file!' );
    return;
end
if ( isempty(par.core) )
    runtime_log( 0, 'ERROR: there is no Core SNPs file name in "*.tsim" file, otherwise there are problems while opening file!' );
    return;
end
if ( isempty(par.param) )
    runtime_log( 0, 'ERROR: there is no model parameters file name in "*.tsim" file, otherwise there are problems while opening file!' );
    return;
end
% if ( isempty(par.network) )
%     runtime_log( 0, 'ERROR: there is no network parameters file name in "*.tsim" file, otherwise there are problems while opening file!' );
%     return;
% end
if ( isempty(par.scip) )
    par.scip = 0.2;
end
if ( isempty(par.all) )
    par.all = 1;
end
if ( isempty(par.reproduce) )
    par.reproduce = 0;
end
if ( isempty(par.rng) )
    par.rng = 0;
end

file_allele = par.allele;
file_core = par.core;
% file_model_param = par.model;
% file_ntw_param = par.network;
file_param = par.param;
skip_sol = par.scip;
run_all = par.all;
reproduce = par.reproduce;
rngState = par.rng;

if (reproduce)
    rng(rngState);
    rng_state = rng;
else
    rng('shuffle');
    rng_state = rng;
end

runtime_log( 0, 'Alleles file name......................:', file_allele );
runtime_log( 0, 'Core alleles file name.................:', file_core );
% runtime_log( 0, 'Model parameters file name.............:', file_model_param );
% runtime_log( 0, 'Network parameters file name...........:', file_ntw_param );
runtime_log( 0, 'Parameters file name...................:', file_param );
runtime_log( 0, 'Amount of data to remove from solution.:', num2str(skip_sol) );
runtime_log( 0, 'Simulate ALL genotypes.................:', num2str(run_all) );
runtime_log( 0, 'Reproduce simulation...................:', num2str(reproduce) );
runtime_log( 0, 'State of random number generator.......:', num2str(rngState) );

if ( isempty(par.savedfileN) )
    fSavedNetwork = 0;
else
    fSavedNetwork = par.savedfileN;
    runtime_log( 0, 'Reused Network file name...............:', fSavedNetwork );
end

if ( isempty(par.savedfileA) )
    fSavedActivators = 0;
else
    fSavedActivators = par.savedfileA;
    runtime_log( 0, 'Reused Activators file name............:', fSavedActivators );
end

if ( isempty(par.savedfileR) )
    fSavedRepressors = 0;
else
    fSavedRepressors = par.savedfileR;
    runtime_log( 0, 'Reused Repressors file name............:', fSavedRepressors );
end

if ( isempty(par.savedpath) )
    fSavedPath = 0;
else
    fSavedPath = par.savedpath;
    runtime_log( 0, 'Output results path....................:', fSavedPath );
end

% if ( strcmp( fSavedNetwork, 'foo' ) )
%     fSavedNetwork = 0;
% else
%     runtime_log( 0, 'Reused Network file name...............:', fSavedNetwork );
% end
% 
% if ( strcmp( fSavedActivators, 'foo' ) )
%     fSavedActivators = 0;
% else
%     runtime_log( 0, 'Reused Activators file name............:', fSavedActivators );
% end
% 
% if ( strcmp( fSavedRepressors, 'foo' ) )
%     fSavedRepressors = 0;
% else
%     runtime_log( 0, 'Reused Repressors file name............:', fSavedRepressors );
% end
% 
% if ( strcmp( fSavedPath, 'foo' ) )
%     fSavedPath = 0;
% else
%     runtime_log( 0, 'Output results path....................:', fSavedPath );
% end

[ parN, parM ] = fReadParam( file_param );

if ( isempty(parN.arch) )
    parN.arch = 'downstream';
end
if ( isempty(parN.n1) )
    parN.n1 = 0.0;
end
if ( isempty(parN.n2) )
    parN.n2 = 0.0;
end
if ( isempty(parN.n3) )
    parN.n3 = 0.0;
end
if ( isempty(parN.b1) )
    parN.b1 = 0.1;
end
if ( isempty(parN.b2) )
    parN.b2 = 0.1;
end
if ( isempty(parN.b3) )
    parN.b3 = 0.1;
end
if ( isempty(parN.mrg) )
    parN.mrg = 1;
end
if ( isempty(parN.mode) )
    parN.mode = 'strict';
end
if ( isempty(parN.ratio) )
    parN.ratio = 0.5;
end

if ( isempty(parM.raterna) )
    parM.raterna = 4.0;
end
if ( isempty(parM.ratep) )
    parM.ratep = 8.0;
end
if ( isempty(parM.degrna) )
    parM.degrna = 0.2;
end
if ( isempty(parM.degp) )
    parM.degp = 0.6;
end
if ( isempty(parM.ebind) )
    parM.ebind = 9.0;
end
if ( isempty(parM.kbind) )
    parM.kbind = 8100;
end
if ( isempty(parM.kact) )
    parM.kact = 500000;
end
if ( isempty(parM.eact) )
    parM.eact = 6.0;
end
if ( isempty(parM.krep) )
    parM.krep = 500;
end
if ( isempty(parM.erep) )
    parM.erep = 7.0;
end
if ( isempty(parM.gendiff) )
    parM.gendiff = 0.0;
end
if ( isempty(parM.snpdiff) )
    parM.snpdiff = 0.25;
end
if ( isempty(parM.tmax) )
    parM.tmax = 500;
end
if ( isempty(parM.tdil) )
    parM.tdil = 0;
end
if ( isempty(parM.stoch) )
    parM.stoch = 0;
end

runtime_log( 0, 'Type of network architecture...........:', parN.arch );
runtime_log( 0, 'Proportion of genes in network type 1..:', num2str(parN.n1) );
runtime_log( 0, 'Proportion of genes in network type 2..:', num2str(parN.n2) );
runtime_log( 0, 'Proportion of genes in network type 3..:', num2str(parN.n3) );
runtime_log( 0, 'Parameter Beta for network type 1......:', num2str(parN.b1) );
runtime_log( 0, 'Parameter Beta for network type 2......:', num2str(parN.b2) );
runtime_log( 0, 'Parameter Beta for network type 3......:', num2str(parN.b3) );
runtime_log( 0, 'Networks merging parameter.............:', num2str(parN.mrg) );
runtime_log( 0, 'Mode of Activator|Repressor setting....:', parN.mode );
runtime_log( 0, 'Proportion of Activator|Repressor......:', num2str(parN.ratio) );
runtime_log( 0, 'mRNA transcription rate................:', num2str(parM.raterna) );
runtime_log( 0, 'Protein translation rate...............:', num2str(parM.ratep) );
runtime_log( 0, 'mRNA degradation rate..................:', num2str(parM.degrna) );
runtime_log( 0, 'Protein degradation rate...............:', num2str(parM.degp) );
runtime_log( 0, 'Polymerase_II complex binding energy...:', num2str(parM.ebind) );
runtime_log( 0, 'Polymerase_II complex binding constant.:', num2str(parM.kbind) );
runtime_log( 0, 'Activator complex binding constant.....:', num2str(parM.kact) );
runtime_log( 0, 'Activator complex binding energy.......:', num2str(parM.eact) );
runtime_log( 0, 'Repressor complex binding constant.....:', num2str(parM.krep) );
runtime_log( 0, 'Repressor complex binding energy.......:', num2str(parM.erep) );
runtime_log( 0, 'Difference among genes in network......:', num2str(parM.gendiff) );
runtime_log( 0, 'Strength of SNP polymorphism...........:', num2str(parM.snpdiff) );
runtime_log( 0, 'Max time of model dynamix solution.....:', num2str(parM.tmax) );
runtime_log( 0, 'Time delay due to molecules diffusion..:', num2str(parM.tdil) );
runtime_log( 0, 'Stochastic perturbation................:', num2str(parM.stoch) );

[ snpMatr, indID ] = ImportSNPtoMatrix(file_allele);

if ( isempty(snpMatr) )
    runtime_log( 0, 'ERROR: program cannot read snp data from allele file or its content is empty, otherwise there is some access problem!' );
    return;
end

runtime_log( 0, 'Number of alleles......................:', num2str( size(snpMatr,2) ) );
runtime_log( 0, 'Number of genotypes....................:', num2str( size(snpMatr,1) ) );

for icontent = 1:min(10,size(snpMatr,1))
    runtime_log( 0, 'Alleles content (first 20 SNPs)........:', num2str( snpMatr(icontent, 1:min(20,size(snpMatr,2))) ) );
end

if ( isempty(indID) )
    runtime_log( 0, 'ERROR: snp IDs matrix returned empty!' );
    return;
end

for icontent = 1:min(10,size(indID,1))
    runtime_log( 0, 'SNP IDs content (first 20 IDs).........:', num2str( indID(1, 1:min(20,size(indID,2))) ) );
end

snpSz = size(snpMatr);
snp_num = snpSz(1,2);
ind_num = snpSz(1,1);

runtime_log( 0, 'MAKING REFERENCE GENOTYPE' );

%[ refGenotype, refParam, rngSeed1 ] = GetRefGenotype( snpMatr, file_model_param, 0.8, true, rng_state, indID );
[ refGenotype, refParam, rngSeed1 ] = GetRefGenotype( snpMatr, parM, 0.8, true, rng_state, indID );

runtime_log( 0, 'MAKING GENOMIC REGULATORY NETWORK' );

%[ N, A, R, N_par, rngSeed2 ] = GetGeneNetwork( file_core, fSavedPath, file_ntw_param, snp_num, true, rng_state, fSavedNetwork, fSavedActivators, fSavedRepressors );
[ N, A, R, N_par, rngSeed2 ] = GetGeneNetwork( file_core, fSavedPath, parN, snp_num, true, rng_state, fSavedNetwork, fSavedActivators, fSavedRepressors );

if ( isempty(N) || isempty(A) || isempty(R) )
    runtime_log( 0, 'STOP: due to network creation problem!' );
    return;
end

ntw_genes = N_par.genes;
ntw_arch = N_par.type;

if ( ntw_genes == 0 )
    runtime_log( 0, 'ERROR: number of genes in the network is 0!' );
    return;
end

runtime_log( 0, 'Number of genes in the network.........:', num2str( ntw_genes ) );
runtime_log( 0, 'Type of network architecture...........:', ntw_arch );

runtime_log( 0, 'MAKING ASSOCIATION BETWEEN SNP ( == values ) AND NETWORK NODES ( == keys )' );

[ core_map, ncore_map, core_weit, rec_map ] = GetCoreNonCoreMap( file_core, ntw_genes, ntw_arch );

WriteRecodedNetworks( fSavedPath, rec_map, N, A, R );

runtime_log( 0, 'SOLVING FOR REFERENCE GENOTYPE' );

% define structure which contains dynamic solutions
for isol = 1:size(refGenotype.genotypes_id,2)
    
    dynsol(isol).x = [];
    dynsol(isol).y = [];
    
end

Xr = [];
Yr = [];

[ Xr, Yr, rngSeed3, r_weit ] = SolveGRM( parM.snpdiff, refGenotype, refParam, refGenotype.alleles, A, R, core_map, ncore_map, core_weit, true, rng_state );

if ( isempty(Xr) || isempty(Yr) )
    runtime_log( 0, 'STOP: due to GRM solution problem -> solution is empty!' );
    return;
end

[ r_trait, r_vtrait ] = CalcTrait( Xr, Yr, core_map, r_weit, ntw_genes, skip_sol, []);

if ( ~run_all )
    dynsol(1).x = Xr;
    dynsol(1).y = Yr;
    runtime_log( 0, 'WRITING CORE TRAITS AND CORE PROTUCTS TO A FILE' );
    WriteResults( fSavedPath, rng_state, r_trait, r_vtrait, core_map, ncore_map, refGenotype.genotypes_id, dynsol, ntw_genes, run_all );
    runtime_log( 0, 'SIMULATION COMPLETED SUCCESSFULLY' );
    return;
end

runtime_log( 0, 'SOLVING FOR EACH GENOTYPE' );

jobsize = 1000;
blocks = ceil( ind_num/jobsize );

for par_i = 1:blocks
    
    % Explanation for outer/inner 'FOR/PARFOR' organisation (but not for inner/outer 'PARFOR/FOR').
    % The role of outer 'FOR' loop is to unload the memory during computation. 
    % In order to inner 'PARFOR' loop be efficient, the 'blocks' variable should be
    % small, just enough to keep memory consumption constant during computations.
    % If the jobsize = ind_num, the case become trivial where the outer 'FOR' loop
    % losts its role and inner 'PARFOR' loop takes full load with maximum
    % memory consumption.
    
    ini_i = (par_i-1) * jobsize + 1;
    end_i = par_i * jobsize;
    
    if (end_i > ind_num)
        end_i = ind_num;
    end
    
    %a_trait = zeros(ind_num,1);
    
    runtime_log( 0, 'Solving for the next block of genotypes:', num2str( end_i-ini_i+1 ) );
    
    %parfor ig = 1:ind_num
    istart = 1;
    iend = end_i-ini_i+1;
    %ig = ini_i:end_i
    parfor ig = istart:iend
        
        [ Xi, Yi, rngSeed4, cweit ] = SolveGRM( parM.snpdiff, refGenotype, refParam, snpMatr(ig+ini_i-1,:), A, R, core_map, ncore_map, core_weit, true, rng_state );
        dynsol(ig).x = Xi;
        dynsol(ig).y = Yi;
        [ a_trait(ig,1), vtrait(:,ig) ] = CalcTrait( Xi, Yi, core_map, cweit, ntw_genes, skip_sol, r_vtrait);

    end

    runtime_log( 0, 'WRITING CORE TRAITS AND CORE PROTUCTS TO A FILE' );
    
    gRange.start = ini_i;
    gRange.end = end_i;
    
    WriteResults( fSavedPath, rng_state, a_trait, vtrait, core_map, ncore_map, refGenotype.genotypes_id(1,ini_i:end_i), dynsol, ntw_genes, run_all, gRange );
    
    %WriteResults( rng_state, a_trait, vtrait, core_map, ncore_map, refGenotype.genotypes_id, dynsol, ntw_genes );

end

runtime_log( 0, 'SIMULATION COMPLETED SUCCESSFULLY' );

end


%---------------------------------------------------------------
% Local Functions
%---------------------------------------------------------------
function [ p ] = ReadTraitSimParam( filename )

delimiter = ',';

formatSpec = '%s%*s%*s%[^\n\r]';

fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);

fclose(fileID);

parameters = [dataArray{1:end-1}];

clearvars filename delimiter formatSpec fileID dataArray ans;

for i = 1:size(parameters,1)
    switch i
        case 2
            p.allele = char( parameters(i,1) );
        case 4
            p.core = char( parameters(i,1) );
        case 6
            p.model = char( parameters(i,1) );
        case 8
            p.network = char( parameters(i,1) );
        case 10
            p.scip = str2double( char( parameters(i,1) ) );
        case 12
            p.all = str2double( char( parameters(i,1) ) );
        case 14
            p.reproduce = str2double( char( parameters(i,1) ) );
        case 16
            p.rng = str2double( char( parameters(i,1) ) );
        case 18
            p.savedfileN = char( parameters(i,1) );
        case 20
            p.savedfileA = char( parameters(i,1) );
        case 22
            p.savedfileR = char( parameters(i,1) );
        case 24
            p.savedpath = char( parameters(i,1) );
    end
end

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function [ p ] = fReadSimPar( filename )

fid = fopen( filename );

tline = '0';
msg = '0';
i = 1;
info = {};
pos = 0;

while ( tline ~= -1 )
    tline = fgetl(fid);
    
    while ( isempty(tline) )
        tline = fgetl(fid);
    end
    
    if ( ~isempty(tline) && tline(1,1) == '$' )
        
        pos = ftell(fid);
        msg = fgetl(fid);
        
        while ( isempty(msg) || msg(1,1) == '#' )
            pos = ftell(fid);
            msg = fgetl(fid);
        end
        if ( msg(1,1) == '$' )
            fseek(fid,pos,'bof');
        else
            if (msg(1,1) ~= -1)
                info{i,2} = msg;
            end
        end
        info{i,1} = tline(2:end);
    end
    i = i + 1;
end

p.allele = [];
p.core = [];
p.model = [];
p.network = [];
p.scip = [];
p.all = [];
p.reproduce = [];
p.rng = [];
p.savedfileN = [];
p.savedfileA = [];
p.savedfileR = [];
p.savedpath = [];
p.more = [];

for i = 1:size(info,1)
    if ( ~isempty(info{i,1}) )
        switch info{i,1}
            case 'SNP'
                p.allele = info{i,2};
            case 'CORE'
                p.core = info{i,2};
            case 'GRM'
                p.model = info{i,2};
            case 'NTW'
                p.network = info{i,2};
            case 'MODELPAR'
                p.param = info{i,2};
            case 'REMOVE'
                p.scip = str2double( info{i,2} );
            case 'ALL'
                p.all = str2double( info{i,2} );
            case 'REPRODUCE'
                p.reproduce = str2double( info{i,2} );
            case 'RNG'
                p.rng = str2double( info{i,2} );
            case 'NSAVED'
                p.savedfileN = info{i,2};
            case 'ASAVED'
                p.savedfileA = info{i,2};
            case 'RSAVED'
                p.savedfileR = info{i,2};
            case 'OUTRES'
                p.savedpath = info{i,2};
        end
    end
end

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function [ p ] = fReadNtwPar( filename )

fid = fopen( filename );

tline = '0';
msg = '0';
i = 1;
info = {};
pos = 0;

while ( tline ~= -1 )
    tline = fgetl(fid);
    
    while ( isempty(tline) )
        tline = fgetl(fid);
    end
    
    if ( ~isempty(tline) && tline(1,1) == '%' )
        
        pos = ftell(fid);
        msg = fgetl(fid);
        
        while ( isempty(msg) )
            pos = ftell(fid);
            msg = fgetl(fid);
        end
        if ( msg(1,1) == '%' )
            fseek(fid,pos,'bof');
        else
            if (msg(1,1) ~= -1)
                info{i,2} = msg;
            end
        end
        info{i,1} = tline(2:end);
    end
    i = i + 1;
end

p.arch = [];
p.n1 = [];
p.n2 = [];
p.n3 = [];
p.b1 = [];
p.b2 = [];
p.b3 = [];
p.mrg = [];
p.mode = [];
p.ratio = [];

for i = 1:size(info,1)
    switch info{i,1}
        case 'ARCH'
            p.arch = info{i,2};
        case 'N1'
            p.n1 = str2double( info{i,2} );
        case 'N2'
            p.n2 = str2double( info{i,2} );
        case 'N3'
            p.n3 = str2double( info{i,2} );
        case 'BETA1'
            p.b1 = str2double( info{i,2} );
        case 'BETA2'
            p.b2 = str2double( info{i,2} );
        case 'BETA3'
            p.b3 = str2double( info{i,2} );
        case 'MERGE'
            p.mrg = str2double( info{i,2} );
        case 'MODE'
            p.mode = info{i,2};
        case 'RATIO'
            p.ratio = str2double( info{i,2} );
    end
end

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function [ p, m ] = fReadParam( filename )

fid = fopen( filename );

tline = '0';
msg = '0';
i = 1;
info = {};
pos = 0;

while ( tline ~= -1 )
    tline = fgetl(fid);
    
    while ( isempty(tline) )
        tline = fgetl(fid);
    end
    
    if ( ~isempty(tline) && tline(1,1) == '$' )
        
        pos = ftell(fid);
        msg = fgetl(fid);
        
        while ( isempty(msg) || msg(1,1) == '#' )
            pos = ftell(fid);
            msg = fgetl(fid);
        end
        if ( msg(1,1) == '$' )
            fseek(fid,pos,'bof');
        else
            if (msg(1,1) ~= -1)
                info{i,2} = msg;
            end
        end
        info{i,1} = tline(2:end);
    end
    i = i + 1;
end

p.arch = [];
p.n1 = [];
p.n2 = [];
p.n3 = [];
p.b1 = [];
p.b2 = [];
p.b3 = [];
p.mrg = [];
p.mode = [];
p.ratio = [];

m.raterna = [];
m.ratep = [];
m.degrna = [];
m.degp = [];
m.ebind = [];
m.kbind = [];
m.kact = [];
m.eact = [];
m.krep = [];
m.erep = [];
m.gendiff = [];
m.snpdiff = [];
m.tmax = [];
m.tdil = [];
m.stoch = [];

for i = 1:size(info,1)
    if ( ~isempty(info{i,1}) )
        switch info{i,1}
            case 'ARCH'
                p.arch = info{i,2};
            case 'N1'
                p.n1 = str2double( info{i,2} );
            case 'N2'
                p.n2 = str2double( info{i,2} );
            case 'N3'
                p.n3 = str2double( info{i,2} );
            case 'BETA1'
                p.b1 = str2double( info{i,2} );
            case 'BETA2'
                p.b2 = str2double( info{i,2} );
            case 'BETA3'
                p.b3 = str2double( info{i,2} );
            case 'MERGE'
                p.mrg = str2double( info{i,2} );
            case 'MODE'
                p.mode = info{i,2};
            case 'RATIO'
                p.ratio = str2double( info{i,2} );            
            case 'RATERNA'
                m.raterna = str2double( info{i,2} );
            case 'RATEP'
                m.ratep = str2double( info{i,2} );
            case 'DEGRNA'
                m.degrna = str2double( info{i,2} );
            case 'DEGP'
                m.degp = str2double( info{i,2} );
            case 'EBIND'
                m.ebind = str2double( info{i,2} );
            case 'KBIND'
                m.kbind = str2double( info{i,2} );
            case 'KACT'
                m.kact = str2double( info{i,2} );
            case 'EACT'
                m.eact = str2double( info{i,2} );
            case 'KREP'
                m.krep = str2double( info{i,2} );
            case 'EREP'
                m.erep = str2double( info{i,2} );
            case 'GENDIFF'
                m.gendiff = str2double( info{i,2} );
            case 'SNPDIFF'
                m.snpdiff = str2double( info{i,2} );
            case 'TMAX'
                m.tmax = str2double( info{i,2} );
            case 'TDIL'
                m.tdil = str2double( info{i,2} );
            case 'STOCH'
                m.stoch = str2double( info{i,2} );
        end
    end
end

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function WriteRecodedNetworks( rslt_path, rMap, ntwN, ntwA, ntwR )

ntwSize = size( ntwN );

% write (recoded in snpIDs) N, A and R networks to files

if ( ~rslt_path )
    res_folder = strcat( 'traitsim_results-', date );
else
    res_folder = rslt_path;
end

if ( ~MakeFolder(res_folder) )
    return;
end

fName = strcat( res_folder, '/N_snpID.gntw' );
[row, col, v] = find(ntwN);
for i = 1:size(row,1)
    row( i,1 ) = rMap( row(i,1) );
    col( i,1 ) = rMap( col(i,1) );
end
row( size(row,1)+1, 1 ) = ntwSize(1,1);
col( size(col,1)+1, 1 ) = ntwSize(1,1);
v( size(v,1)+1, 1 ) = 0;
dlmwrite(fName,[row col v], 'delimiter', '\t');

fName = strcat( res_folder, '/A_snpID.gntw' );
[row, col, v] = find(ntwA);
for i = 1:size(row,1)
    row( i,1 ) = rMap( row(i,1) );
    col( i,1 ) = rMap( col(i,1) );
end
row( size(row,1)+1, 1 ) = ntwSize(1,1);
col( size(col,1)+1, 1 ) = ntwSize(1,1);
v( size(v,1)+1, 1 ) = 0;
dlmwrite(fName,[row col v], 'delimiter', '\t');

fName = strcat( res_folder, '/R_snpID.gntw' );
[row, col, v] = find(ntwR);
for i = 1:size(row,1)
    row( i,1 ) = rMap( row(i,1) );
    col( i,1 ) = rMap( col(i,1) );
end
row( size(row,1)+1, 1 ) = ntwSize(1,1);
col( size(col,1)+1, 1 ) = ntwSize(1,1);
v( size(v,1)+1, 1 ) = 0;
dlmwrite(fName,[row col v], 'delimiter', '\t');

end
%---------------------------------------------------------------
