function [ RGenotype, ModParam, s ] = GetRefGenotype( SNP, paramDataStr, codingRegionSize, reproduce, rngState, ind_id )

if ~exist('ind_id','var')
    use_indID = false;
else
    use_indID = true;
end

% set the state of random number generator
if (reproduce)
    rng(rngState);
    s = rng;
else
    rng('shuffle');
    s = rng;
end

snpSz = size(SNP);
snp_num = snpSz(1,2);
ind_num = snpSz(1,1);

% (1) get max frequent SNP variants within each snpID
RGenotype.alleles = getRefAlleles( SNP );

% (2) get distribution of polymorphisme location: non-coding|coding
RGenotype.plmph_locat = getPolymorphLocation( snp_num, codingRegionSize );

% (3) Assign the necessary parameters of GR model for each gene in reference genotype
[ RGenotype.grm_param, ModParam ] = setRefGenotModelParam( paramDataStr, snp_num, s );

% (4) create IDs for SNPs and genotypes
if ( use_indID )
    RGenotype.genotypes_id = ind_id;
    RGenotype.snp_id = getDataIDs( snp_num, ind_num );
else
    [ RGenotype.snp_id, RGenotype.genotypes_id ] = getDataIDs( snp_num, ind_num );
end

end
%---------------------------------------------------------------


%---------------------------------------------------------------
% Local Functions
%---------------------------------------------------------------

function [ refGenotypeModelPar, refParam ] = setRefGenotModelParam( parstr, snpNum, rngState )

% Each gene in reference genotype is parametrised by random sampling from
% normal distribution N( base_value, std(base_value) ).

% set the state of random number generator
rng(rngState);

% % read parameters data from *.rgna file
% delimiter = ';';
% 
% formatSpec = '%s%s%*s%[^\n\r]';
% 
% fileID = fopen(fileName,'r');
% 
% dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
% 
% fclose(fileID);
% 
% param = table(dataArray{1:end-1});
% 
% paramSz = size(param);
% for i3 = 2:2:paramSz(1,1)
%     switch i3
%         case 2
%             refGene.kr = str2double( table2array( param(i3,1) ) );
%             refVar.kr = str2double( table2array( param(i3,2) ) );
%         case 4
%             refGene.kp = str2double( table2array( param(i3,1) ) );
%             refVar.kp = str2double( table2array( param(i3,2) ) );
%         case 6
%             refGene.gr = str2double( table2array( param(i3,1) ) );
%             refVar.gr = str2double( table2array( param(i3,2) ) );
%         case 8
%             refGene.gp = str2double( table2array( param(i3,1) ) );
%             refVar.gp = str2double( table2array( param(i3,2) ) );
%         case 10
%             refGene.Gbind = -str2double( table2array( param(i3,1) ) );
%             refVar.Gbind = str2double( table2array( param(i3,2) ) );
%         case 12
%             refGene.kbind = str2double( table2array( param(i3,1) ) );
%             refVar.kbind = str2double( table2array( param(i3,2) ) );
%         case 14
%             refGene.kattr = str2double( table2array( param(i3,1) ) );
%             refVar.kattr = str2double( table2array( param(i3,2) ) );
%         case 16
%             refGene.GattrA = -str2double( table2array( param(i3,1) ) );
%             refVar.GattrA = str2double( table2array( param(i3,2) ) );
%         case 18
%             refGene.krepr = str2double( table2array( param(i3,1) ) );
%             refVar.krepr = str2double( table2array( param(i3,2) ) );
%         case 20
%             refGene.GattrR = -str2double( table2array( param(i3,1) ) );
%             refVar.GattrR = str2double( table2array( param(i3,2) ) );
%         case 22
%             refGene.T = str2double( table2array( param(i3,1) ) );
%             refVar.L = str2double( table2array( param(i3,2) ) );
%             if ( refVar.L == 0.0 )
%                 refVar.L = 0.01;
%             end
%         case 24
%             refGene.stoch = str2double( table2array( param(i3,1) ) );
%             refVar.stoch_var = str2double( table2array( param(i3,2) ) );
%     end
% end
refGene.kr = parstr.raterna;
refVar.kr = parstr.gendiff * parstr.raterna;

refGene.kp = parstr.ratep;
refVar.kp = parstr.gendiff * parstr.ratep;

refGene.gr = parstr.degrna;
refVar.gr = parstr.gendiff * parstr.degrna;

refGene.gp = parstr.degp;
refVar.gp = parstr.gendiff * parstr.degp;

refGene.Gbind = -parstr.ebind;
refVar.Gbind = 0.0;

refGene.kbind = parstr.kbind;
refVar.kbind = 0.0;

refGene.kattr = parstr.kact;
refVar.kattr = 0.0;

refGene.GattrA = -parstr.eact;
refVar.GattrA = 0.0;

refGene.krepr = parstr.krep;
refVar.krepr = 0.0;

refGene.GattrR = -parstr.erep;
refVar.GattrR = 0.0;

refGene.T = parstr.tmax;
refVar.L = parstr.tdil;
if ( refVar.L == 0.0 )
    refVar.L = 0.01;
end

refGene.stoch = ceil(parstr.stoch);
refVar.stoch_var = parstr.stoch;


refGene.GattrA2 = exp(-refGene.GattrA);
refVar.GattrA2 = 0.0;

refParam.mean = refGene;
refParam.std = refVar;

% making report variables to catch incorrect parameters values
parFlag = zeros(15,snpNum);

parChar(1,1) = "kr";
parChar(2,1) = "kp";
parChar(3,1) = "gr";
parChar(4,1) = "gp";
parChar(5,1) = "Gbind";
parChar(6,1) = "kbind";
parChar(7,1) = "kattr";
parChar(8,1) = "krepr";
parChar(9,1) = "GattrA";
parChar(10,1) = "GattrA2";
parChar(11,1) = "GattrR";
parChar(12,1) = "T_max";
parChar(13,1) = "T_d";
parChar(14,1) = "stochastic";
parChar(15,1) = "stochastic_var";
        
% Assign the necessary parameters of GR model for each gene in reference genotype
for i3 = 1:snpNum
    refGenotypeModelPar( i3 ).kr = normrnd( refGene.kr, refVar.kr );
    refGenotypeModelPar( i3 ).kp = normrnd( refGene.kp, refVar.kp );
    refGenotypeModelPar( i3 ).gr = normrnd( refGene.gr, refVar.gr );
    refGenotypeModelPar( i3 ).gp = normrnd( refGene.gp, refVar.gp );
    refGenotypeModelPar( i3 ).Gbind = normrnd( refGene.Gbind, refVar.Gbind );
    refGenotypeModelPar( i3 ).kbind = normrnd( refGene.kbind, refVar.kbind );
    refGenotypeModelPar( i3 ).kattr = normrnd( refGene.kattr, refVar.kattr );
    refGenotypeModelPar( i3 ).krepr = normrnd( refGene.krepr, refVar.krepr );
    refGenotypeModelPar( i3 ).GattrA = normrnd( refGene.GattrA, refVar.GattrA );
    refGenotypeModelPar( i3 ).GattrA2 = normrnd( refGene.GattrA2, refVar.GattrA2 );
    refGenotypeModelPar( i3 ).GattrR = normrnd( refGene.GattrR, refVar.GattrR );
    refGenotypeModelPar( i3 ).T = refGene.T;
    refGenotypeModelPar( i3 ).L = refVar.L;
    refGenotypeModelPar( i3 ).stoch = refGene.stoch;
    refGenotypeModelPar( i3 ).stoch_var = refVar.stoch_var;
    
    % check parameters
    if ( refGenotypeModelPar( i3 ).kr < 0.0 )
        parFlag(1,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).kp < 0.0 )
        parFlag(2,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).gr < 0.0 )
        parFlag(3,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).gp < 0.0 )
        parFlag(4,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).Gbind >= 0.0 )
        parFlag(5,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).kbind <= 0.0 )
        parFlag(6,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).kattr <= 0.0 )
        parFlag(7,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).krepr <= 0.0 )
        parFlag(8,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).GattrA >= 0.0 )
        parFlag(9,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).GattrA2 <= 0.0 )
        parFlag(10,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).GattrR >= 0.0 )
        parFlag(11,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).T <= 0.0 )
        parFlag(12,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).L < 0.0 )
        parFlag(13,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).stoch < 0.0 )
        parFlag(14,i3) = 1;
    end
    if ( refGenotypeModelPar( i3 ).stoch_var < 0.0 )
        parFlag(15,i3) = 1;
    end

end

rowFlag = [];
colFlag = [];

[ rowFlag,colFlag ] = find(parFlag);

if( ~isempty(rowFlag) )
    for ii = 1:size(rowFlag,1)
        runtime_log( 0, 'WARNING: range of parameter............:', parChar(rowFlag(ii,1),1) );
        runtime_log( 0, '                     in SNP............:', num2str( colFlag(ii,1) ) );
    end
    runtime_log( 0, 'WARNING: check the range and|or standard deviation of parameters!' );
end

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function [snpID, indID] = getDataIDs(snpNum, indNum)

snpID = zeros(1, snpNum);
indID = zeros(1, indNum);
for i2 = 1:snpNum
    snpID(1,i2) = i2;
end
for i2 = 1:indNum
    indID(1,i2) = i2;
end

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function [L] = getPolymorphLocation(snpNum, codRate)

% snpNum - is the number of SNPs
% codRate - is the part of all SN polymorphism which fall under coding region

L = zeros(1,snpNum);
for i1 = 1:snpNum
    r = rand;
    if (r <= codRate)
        L(1,i1) = 2; % coding region
    else
        L(1,i1) = 1; % non-coding region
    end
end

end
%---------------------------------------------------------------

%---------------------------------------------------------------
function [P] = getRefAlleles(SNP)

snpSize = size(SNP);

P = zeros( 1,snpSize(1,2) );

for i = 1:snpSize(1,2)
    
    p0 = 0;
    p1 = 0;
    p2 = 0;
    counts = 0;
    for j = 1:snpSize(1,1)
        all = SNP(j,i);
        switch all
            case 0
                p0 = p0 + 1;
                counts = counts + 1;
            case 1
                p1 = p1 + 1;
                counts = counts + 1;
            case 2
                p2 = p2 + 1;
                counts = counts + 1;
        end
    end
    p0 = p0/counts;
    p1 = p1/counts;
    p2 = p2/counts;
    
    pTmp = 0;
    if (p0 >= p1)
        maxAll = 0;
        pTmp = p0;
    else
        maxAll = 1;
        pTmp = p1;
    end
    if (p2 >= pTmp)
        maxAll = 2;
    end
    
    P(1,i) = maxAll;
    
    for j = 1:snpSize(1,1)
        all = SNP(j,i);
        if ( all > 2 )
            SNP(j,i) = maxAll;
        end
    end
    
end

end
%---------------------------------------------------------------
