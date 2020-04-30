
function [a] = GNd1(...
                    N,...
                    s,...
                    beta,...
                    rngState...
                    )

% The code in this function was borrowed from:
% G. Bianconi and C. Rahmede "Network geometry with flavour: from complexity to quantum geometry".

% s & beta - network parameters; s = -1, 0, 1; beta = [0 10];
% a - adjacency matrix of undirected network (graph).

rng(rngState);

% Initialization
a=sparse(N,N);
a_occ=zeros(1,N);

% Assign energies to the nodes
% If using Poisson and power-law you must define
% the parameters mu, or kappa

% Examples:
 mu=10;
 kappa=1;
for i=1:N
    epsilon(i)=floor(10*rand(1));
    % Alternative energy distributions
    % epsilon(i)=random('Poisson',mu);
    % poisson distribution with average mu
    % epsilon(i)=rand(1)^(1/(kappa+1));
    % power-law distribution with exponent kappa
    
end

% Initial condition: at time t=1 a single link (1,2)
a(1,2)=exp(-beta*(epsilon(1)+epsilon(2)));
a(2,1)=exp(-beta*(epsilon(1)+epsilon(2)));
k(1)=1;
k(2)=1;
a_occ(1)=1;
a_occ(2)=1;

% Addition of new links at time t=in-1 the node in is added to the
% network geometry with flavour

for in=2+1:N
    % Choose the node to which attach a new link
    V=exp(-beta*epsilon).* a_occ;
    norm=sum(V);
    %x=rand(1)*norm;
    x=rand(1)*(norm);
    %j = -1;
    if (norm>0)
        for nj1=1:in-1
            x=x-V(nj1);
            if x<0
                j=nj1;
                break;
            end
        end
    end
    %if(j > -1)
    % Attach the new link between node in and node j
    a(in,j)=exp(-beta*epsilon(in)-beta*epsilon(j));
    a(j,in)=a(in,j);
    a_occ(in)=1;
    a_occ(j)=a_occ(j)+s;
    %a_occ(j)=a_occ(j)*exp(-gamma*a_occ(j));
    %end
end

end