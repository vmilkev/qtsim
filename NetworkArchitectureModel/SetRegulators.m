function [ A, R ] = SetRegulators( N, mode, mode_val, rng_state )

% Function which extract from a directed gene network N(n,n)
% two distinct regulatory sub-networks:
% (1) A(n,n) - gene activators network;
% (2) R(n,n) - gene repressors network.
%
% Input Parameters:
%
% GNA - a directed gene network of size (n,n).
%
% mode - method of the network extraction:
%        mode = 'strict' - a particular gene will have all TFs either
%                          as activators or as repressors;
%        mode = 'mix' - a particular gene will have a mixed amount of
%                       activators and repressors.
%        mode = 'tf' - set a particular genes as activator|repressor
%                      transctiption factor.
% mode_val - parameter indicating a proportion of activators/repressors in GRN:
%            mode_val = 0...1; if mode_val = 0 - all regulators are repressors, if
%            mode_val = 1 - all regulators are activators, if mode_val = 0.5 - amount
%            of activators and repressors are equal.
%            In short, mode_val is a relative amount of A.
%
% reproduce = true|false - indicate whether to use custom state of random number generator
%                          (for reproducibility purposes).
%
% rng_state - saved state of a random number generator.
%
% Output Parameters:
%
% A - activator regulatory network of size (n,n).
%
% R - repressor regulatory network of size (n,n).
%
% s - state of a random number generator.

rng(rng_state);

gSz = size(N);
n = gSz(1,1);

A = zeros(n,n);
R = zeros(n,n);

if ( strcmp(mode, 'strict') )
    for i = 1:n
        which_reg = rand;
        if (which_reg >= mode_val)
            R(:,i) = N(:,i);
        else
            A(:,i) = N(:,i);
        end
    end
elseif ( strcmp(mode, 'mix') )
    for i = 1:n
        for j = 1:i
            if ( N(i,j) ~= 0 )
                which_reg = rand;
                if (which_reg >= mode_val)
                    R(i,j) = N(i,j);
                else
                    A(i,j) = N(i,j);
                end
            end
        end
    end
else
    for i = 1:n
        which_reg = rand;
        if (which_reg >= mode_val)
            R(i,:) = N(i,:);
        else
            A(i,:) = N(i,:);
        end
    end
end

end