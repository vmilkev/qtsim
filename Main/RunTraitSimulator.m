function RunTraitSimulator()

runtime_log( 1, 'STARTING TraitSim' );

fParam = FindParFile( pwd, '.tsim' );

if ( isempty(fParam) )
    runtime_log( 0, 'ERROR: program cannot find "*.tsim" file!' );
    return;
else
    runtime_log( 0, 'TraitSim parameters file name..........:', fParam );
end

par = ReadTraitSimParam( fParam );

if ( isempty(par) )
    runtime_log( 0, 'ERROR: content of "*.tsim" file is empty, otherwise there are problems while opening file!' );
    return;
end

file_allele = par.allele;
file_core = par.core;
file_model_param = par.model;
file_ntw_param = par.network;
skip_sol = par.scip;
run_all = par.all;
reproduce = par.reproduce;
rngState = par.rng;
fSavedNetwork = par.savedfile;


if (reproduce)
    rng(rngState);
    rng_state = rng;
else
    rng('shuffle');
    rng_state = rng;
end

runtime_log( 0, 'Alleles file name......................:', file_allele );
runtime_log( 0, 'Core alleles file name.................:', file_core );
runtime_log( 0, 'Model parameters file name.............:', file_model_param );
runtime_log( 0, 'Network parameters file name...........:', file_ntw_param );
runtime_log( 0, 'Amount of data to remove from solution.:', num2str(skip_sol) );
runtime_log( 0, 'Simulate ALL genotypes.................:', num2str(run_all) );
runtime_log( 0, 'Reproduce simulation...................:', num2str(reproduce) );
runtime_log( 0, 'State of random number generator.......:', num2str(rngState) );
runtime_log( 0, 'Reused Network file name...............:', fSavedNetwork );

if ( strcmp( fSavedNetwork, 'foo' ) )
    fSavedNetwork = 0;
end

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

[ refGenotype, refParam, rngSeed1 ] = GetRefGenotype( snpMatr, file_model_param, 0.8, true, rng_state, indID );

runtime_log( 0, 'MAKING GENOMIC REGULATORY NETWORK' );

[ N, A, R, N_par, rngSeed2 ] = GetGeneNetwork( file_ntw_param, snp_num, true, rng_state, fSavedNetwork );

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

[ core_map, ncore_map, core_weit ] = GetCoreNonCoreMap( file_core, ntw_genes, ntw_arch );

runtime_log( 0, 'SOLVING FOR REFERENCE GENOTYPE' );

% define structure which contains dynamic solutions
for isol = 1:size(refGenotype.genotypes_id,2)
    
    dynsol(isol).x = [];
    dynsol(isol).y = [];
    
end

Xr = [];
Yr = [];

[ Xr, Yr, rngSeed3, r_weit ] = SolveGRM( refGenotype, refParam, refGenotype.alleles, A, R, core_map, ncore_map, core_weit, true, rng_state );

if ( isempty(Xr) || isempty(Yr) )
    runtime_log( 0, 'STOP: due to GRM solution problem -> solution is empty!' );
    return;
end

[ r_trait, r_vtrait ] = CalcTrait( Xr, Yr, core_map, r_weit, ntw_genes, skip_sol, []);

if ( ~run_all )
    dynsol(1).x = Xr;
    dynsol(1).y = Yr;
    runtime_log( 0, 'WRITING CORE TRAITS AND CORE PROTUCTS TO A FILE' );
    WriteResults( rng_state, r_trait, r_vtrait, core_map, ncore_map, refGenotype.genotypes_id, dynsol, ntw_genes, run_all );
    runtime_log( 0, 'SIMULATION COMPLETED SUCCESSFULLY' );
    return;
end

a_trait = zeros(ind_num,1);

runtime_log( 0, 'SOLVING FOR EACH GENOTYPE' );

parfor ig = 1:ind_num
    
    [ Xi, Yi, rngSeed4, cweit ] = SolveGRM( refGenotype, refParam, snpMatr(ig,:), A, R, core_map, ncore_map, core_weit, true, rng_state );
    dynsol(ig).x = Xi;
    dynsol(ig).y = Yi;
    [ a_trait(ig,1), vtrait(:,ig) ] = CalcTrait( Xi, Yi, core_map, cweit, ntw_genes, skip_sol, r_vtrait);
    
end

runtime_log( 0, 'WRITING CORE TRAITS AND CORE PROTUCTS TO A FILE' );

WriteResults( rng_state, a_trait, vtrait, core_map, ncore_map, refGenotype.genotypes_id, dynsol, ntw_genes );

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
            p.savedfile = char( parameters(i,1) );
    end
end

end
%---------------------------------------------------------------
