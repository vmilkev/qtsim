function [ X, Y, s, weit ] = SolveGRM( mut_coef, refGenotype, refParameters, iGenotype, A, R, coreMap, ncoreMap, coreW, reproduce, rngState )

weit = containers.Map(keys(coreW),values(coreW));

if (reproduce)
    rng(rngState);
    s = rng;
else
    rng('shuffle');
    s = rng;
end

k = keys( coreMap );
k_sz = size(k);

%mut_coef = 0.25;

for ik = 1:k_sz(1,2)
    snp = coreMap( k{ik} );
    r_all = refGenotype.alleles( snp );
    i_all = iGenotype( snp );
    mut_type = refGenotype.plmph_locat( snp );
    if (mut_type == 1)
        mut1 = ( mut_coef * (i_all - r_all) + 3 )/3;
        mut2 = 1.0;
    else
        mut2 = ( mut_coef * (i_all - r_all) + 3 )/3;
        mut1 = 1.0;
        if ( mut2 == 1.0 )
            weit( snp ) = coreW( snp ) * mut2;
        else
            % assume that only 10 % of SNP polymorphism in
            % coding siquance leads to changes in gene's product
            if ( rand(1,1) <= 0.1 )
                weit( snp ) = coreW( snp ) * normrnd(1,0.2); % change in core product weight due to mutation
            end
        end
    end
    
    node( k{ik} ).activators = getRegMatInd( A, k{ik} );
    node( k{ik} ).repressors = getRegMatInd( R, k{ik} );
    node( k{ik} ).kr = refGenotype.grm_param( snp ).kr;
    node( k{ik} ).kp = refGenotype.grm_param( snp ).kp;
    node( k{ik} ).gr = refGenotype.grm_param( snp ).gr;
    node( k{ik} ).gp = refGenotype.grm_param( snp ).gp;
    node( k{ik} ).Gbind = mut1 * refGenotype.grm_param( snp ).Gbind;
    node( k{ik} ).kbind = refGenotype.grm_param( snp ).kbind;
    node( k{ik} ).kattr = refGenotype.grm_param( snp ).kattr;
    node( k{ik} ).krepr = refGenotype.grm_param( snp ).krepr;
    node( k{ik} ).GattrA = mut2 * refGenotype.grm_param( snp ).GattrA;
    node( k{ik} ).GattrA2 = refGenotype.grm_param( snp ).GattrA2;
    node( k{ik} ).GattrR = mut2 * refGenotype.grm_param( snp ).GattrR;
end

k = keys( ncoreMap );
k_sz = size(k);

max_snp = max(refGenotype.snp_id);

for ik = 1:k_sz(1,2)
    snp_flag = false;
    snp = ncoreMap( k{ik} );
    if (snp > max_snp)
        snp = max_snp;
        snp_flag = true;
    end
    r_all = refGenotype.alleles( snp );
    i_all = iGenotype( snp );
    mut_type = refGenotype.plmph_locat( snp );
    if (mut_type == 1)
        mut1 = ( mut_coef * (i_all - r_all) + 3 )/3;
        mut2 = 1.0;
    else
        mut2 = ( mut_coef * (i_all - r_all) + 3 )/3;
        mut1 = 1.0;
    end
    
    if (snp_flag) % this is the case if the number of genes in network higher than submitted SNPs
        
        mut1 = 1.0;
        mut2 = 1.0;
        
        % we do not need negative values for these parameters
        sampled_kr = normrnd( refParameters.mean.kr, refParameters.std.kr );
        if ( sampled_kr < 0.0 )
            sampled_kr = 0.0;
        end
        sampled_kp = normrnd( refParameters.mean.kp, refParameters.std.kp );
        if ( sampled_kp < 0.0 )
            sampled_kp = 0.0;
        end
        sampled_gr = normrnd( refParameters.mean.gr, refParameters.std.gr );
        if ( sampled_gr < 0.0 )
            sampled_gr = 0.0;
        end
        sampled_gp = normrnd( refParameters.mean.gp, refParameters.std.gp );
        if ( sampled_gp < 0.0 )
            sampled_gp = 0.0;
        end
        sampled_Gbind = mut1 * normrnd( refParameters.mean.Gbind, refParameters.std.Gbind );
        if ( sampled_Gbind > 0.0 )
            sampled_Gbind = 0.0;
        end

        node( k{ik} ).activators = getRegMatInd( A, k{ik} );
        node( k{ik} ).repressors = getRegMatInd( R, k{ik} );
        node( k{ik} ).kr = sampled_kr;
        node( k{ik} ).kp = sampled_kp;
        node( k{ik} ).gr = sampled_gr;
        node( k{ik} ).gp = sampled_gp;
        node( k{ik} ).Gbind = sampled_Gbind;
        node( k{ik} ).kbind = normrnd( refParameters.mean.kbind, refParameters.std.kbind );
        node( k{ik} ).kattr = normrnd( refParameters.mean.kattr, refParameters.std.kattr );
        node( k{ik} ).krepr = normrnd( refParameters.mean.krepr, refParameters.std.krepr );
        node( k{ik} ).GattrA = mut2 * normrnd( refParameters.mean.GattrA, refParameters.std.GattrA );
        node( k{ik} ).GattrA2 = normrnd( refParameters.mean.GattrA2, refParameters.std.GattrA2 );
        node( k{ik} ).GattrR = mut2 * normrnd( refParameters.mean.GattrR, refParameters.std.GattrR );
    else
        node( k{ik} ).activators = getRegMatInd( A, k{ik} );
        node( k{ik} ).repressors = getRegMatInd( R, k{ik} );
        node( k{ik} ).kr = refGenotype.grm_param( snp ).kr;
        node( k{ik} ).kp = refGenotype.grm_param( snp ).kp;
        node( k{ik} ).gr = refGenotype.grm_param( snp ).gr;
        node( k{ik} ).gp = refGenotype.grm_param( snp ).gp;
        node( k{ik} ).Gbind = mut1 * refGenotype.grm_param( snp ).Gbind;
        node( k{ik} ).kbind = refGenotype.grm_param( snp ).kbind;
        node( k{ik} ).kattr = refGenotype.grm_param( snp ).kattr;
        node( k{ik} ).krepr = refGenotype.grm_param( snp ).krepr;
        node( k{ik} ).GattrA = mut2 * refGenotype.grm_param( snp ).GattrA;
        node( k{ik} ).GattrA2 = refGenotype.grm_param( snp ).GattrA2;
        node( k{ik} ).GattrR = mut2 * refGenotype.grm_param( snp ).GattrR;
    end
end


stochast.perform = refGenotype.grm_param( snp ).stoch;
stochast.var = refGenotype.grm_param( snp ).stoch_var;
stochast.s = s;

t = [0, refGenotype.grm_param( snp ).T];
lags = [refGenotype.grm_param( snp ).L, refGenotype.grm_param( snp ).L];

nodesNum = size(ncoreMap,1)+size(coreMap,1);
y0 = ones(2*nodesNum+1,1);

X = [];
Y = [];

[X,Y] = grm(node, t, y0, lags, stochast);

end


%-------------------------------------------------------
% Local functions
%-------------------------------------------------------

function [reg_ind] = getRegMatInd(reg_matr, which_gene)

j = which_gene;

emptyFlag = true;
ind = 0;
for i5 = 1:numel( reg_matr(:,j) )
    k = reg_matr(i5,j);
    if (k ~= 0)
        emptyFlag = false;
        ind = ind + 1;
        reg_ind(ind,1) = i5;
    end
end

if (emptyFlag)
    reg_ind(1,1) = numel( reg_matr(:,j) )+1;
end

end
%-------------------------------------------------------
