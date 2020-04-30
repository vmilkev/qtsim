function [ X, Y ] = grm ( gene, T, Y0, L, stochastic )

% Function which dynamically solves a Gene Regulation Model (GRM).
%
% Input Parameters:
%
% gene - structure of size of number of genes in a network (i = 1:genesInNetwork); consists of all
%        gene's constants (parameters); it has the following fields:
%        (1) gene(i).activators - vector of genes IDs which are activators for a gene 'i';
%        (2) gene(i).repressors - vector of genes IDs which are repressors for a gene 'i';
%        (3) gene(i).kr - transcription rate for a gene 'i';
%        (4) gene(i).kp - translation rate for a gene 'i';
%        (5) gene(i).gr - mRNA degradation rate for a gene 'i';
%        (6) gene(i).gp - protein degradation rate for a gene 'i';
%        (7) gene(i).Gbind - polymerase_II binding energy for a gene 'i';
%        (8) gene(i).kbind - binding constant for a gene 'i';
%        (9) gene(i).kattr - binding constant for activator TFs in a gene 'i';
%        (10) gene(i).krepr - binding constant for repressor TFs in a gene 'i';
%        (11) gene(i).GattrA - binding energy for activator TFs in a gene 'i';
%        (12) gene(i).GattrA2 - binding energy for activator TFs in a gene 'i';
%        (13) gene(i).GattrR - binding energy for repressor TFs in a gene 'i';
%
% T - time interval for dynamic solution; T[start, end].
%
% Y0 - initial gene's states at T = 0.
%
% L - time delay caused by polymers diffusion;
%     L[1,1] = time delay for mRNA,
%     L[1,2] = time delay for protein.
%
% stochastic - structure served to manage stochastic parameters in GRM; it
%              has the following fields:
%              (1) stochastic.perform = true|false - include|not include
%                  stochastic variation into GRM solution;
%              (2) stochastic.var = variance of normal distribution with a
%                  mean = 0.
%              (3) stochastic.reproduce = true|false - use|do not use a
%                  custom state of a random number generator.
%              (4) stochastic.s = state of a random number generator.
%
% Return Parameters
%
% X - time values of dynymic solution:
%     X = [1,time_steps].
%
% Y - genes solution (molecules concentrations):
%     Y(mRNA) = [1:num_of_genes, 1:time_steps];
%     Y(protein) = [num_of_genes+1:2*num_of_genes, 1:time_steps];
%     Y(virtual) = [2*num_of_genes+1, 1:time_steps] - solution for extra
%     (virtual) variable which should not be taken into account!
%
% s - state of a random number generator; it can be used to reproduce a
%     stochastic solution of GRM; it can be used as following:
%     (1) set stochastic.reproduce = true;
%     (2) set stochastic.s = value of s from previous solution.

rng(stochastic.s);

gSz = size(gene);
n = gSz(1,2);

if (stochastic.perform)
    options = ddeset('RelTol', 1.0e-2, 'NormControl', 'on');
    sol = dde23( @grm_stochastic, L, Y0, T, options );
else
    sol = dde23( @grm, L, Y0, T);
end

X = sol.x;
Y = sol.y;


    function v = grm(t,y,Z)
        
        
        ylag1 = Z(:,1); % lag for mRNA
        ylag2 = Z(:,2); % lag for protein
        
        v = zeros(2*n+1,1);
        
        for i = 1:n
            
            % getting gene's parameters
            indexR = gene(i).repressors+n;
            indexA = gene(i).activators+n;
            kr = gene(i).kr;
            kbind = gene(i).kbind;
            Gbind = gene(i).Gbind;
            krepr = gene(i).krepr;
            GattrR = gene(i).GattrR;
            kattr = gene(i).kattr;
            GattrA = gene(i).GattrA;
            GattrA2 = gene(i).GattrA2;
            gr = gene(i).gr;
            kp = gene(i).kp;
            gp = gene(i).gp;
            
            % equation for mRNA
            v(i) = kr * ...
                   pbind(...
                        kbind, ...
                        Gbind, ...
                        prod( fregr( ylag2(indexR), krepr, GattrR, numel(indexR) ) ).* ...
                        prod( frega( ylag2(indexA), kattr, GattrA, GattrA2, numel(indexA) ) ) ...
                        )- ...
                   gr*y(i);
            
            % equation for protein
            v(i+n) = kp*ylag1(i) - gp*y(i+n);
            
        end
        
        % virtual activator|repressor if there is no real one for a particular gene
        v(2*n+1) = 0;
        
    end


    function v = grm_stochastic(t,y,Z)
        
        
        ylag1 = Z(:,1); % lag for mRNA
        ylag2 = Z(:,2); % lag for protein
        
        v = zeros(2*n+1,1);
        
        for i = 1:n
            
            % getting gene's parameters
            indexR = gene(i).repressors+n;
            indexA = gene(i).activators+n;
            kr = gene(i).kr;
            kbind = gene(i).kbind;
            Gbind = gene(i).Gbind;
            krepr = gene(i).krepr;
            GattrR = gene(i).GattrR;
            kattr = gene(i).kattr;
            GattrA = gene(i).GattrA;
            GattrA2 = gene(i).GattrA2;
            gr = gene(i).gr;
            kp = gene(i).kp;
            gp = gene(i).gp;
            
            % equation for mRNA
            v(i) = kr * ...
                   pbind(...
                        kbind, ...
                        Gbind, ...
                        prod( fregr( ylag2(indexR), krepr, GattrR, numel(indexR) ) ).* ...
                        prod( frega( ylag2(indexA), kattr, GattrA, GattrA2, numel(indexA) ) ) ...
                        )- ...
                   gr*y(i) + normrnd(0,stochastic.var)*y(i);
            
            % equation for protein
            v(i+n) = kp*ylag1(i) - gp*y(i+n) + normrnd(0,stochastic.var)*y(i+n);
            
        end
        
        % virtual activator|repressor if there is no real one for a particular gene
        v(2*n+1) = 0;
        
    end


end