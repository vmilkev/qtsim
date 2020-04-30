function [a] = GNd2(...
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

% Inizialization
a=sparse(N,N);
a_occ=sparse(N,N);
a_occ2=sparse(N,N);

% Assign energies to the nodes
% If using Poisson and power-law you must define
% the parameters mu, or kappa
% Examples:
%mu=10;
% kappa=1
for i=1:N
    epsilon(i)=floor(10*rand(1));
    %epsilon(i)=random('Poisson',mu);
end

% Initial condition at time t=1 including a single triangle between nodes
% 1,2,3
L=0;
for i1=1:3
    for i2=(i1+1):3
        L=L+1;
        a(i1,i2)=exp(-beta*(epsilon(i1)+epsilon(i2)));
        a(i2,i1)=exp(-beta*(epsilon(i1)+epsilon(i2)));
        a_occ(i1,i2)=1;
        a_occ(i2,i1)=1;
        a_occ2(i1,i2)=1;
        a_occ2(i2,i1)=1;
    end
end

% At each time t=in-2 we attach a new triangle
for in=(3+1):N
    % Choose edge (l1,l2) to which we will attach the new triangle
    
    [I,J,V]=find(tril(a.*(a_occ)));
    
    norm=sum(V);
    x=rand(1)*norm;
    if (norm>0)
        for nj1=1:numel(V)
            x=x-V(nj1);
            if x<0
                nj=nj1;
                break;
            end
            
        end
        l1=I(nj);
        l2=J(nj);
        
        
        a_occ(l1,l2)=a_occ(l1,l2)+s;
        a_occ(l2,l1)=a_occ(l2,l1)+s;
        a_occ2(l1,l2)=a_occ2(l1,l2)+1;
        a_occ2(l2,l1)=a_occ2(l2,l1)+1;
        
%         a_occ(l1,l2)=a_occ(l1,l2)*exp(-gamma*a_occ(l1,l2));
%         a_occ(l2,l1)=a_occ(l2,l1)*exp(-gamma*a_occ(l2,l1));
        
        % Attach the new node in to the node l1;
        L=L+1;
        a(in,l1)=exp(-beta*(epsilon(l1)+epsilon(in)));
        a(l1,in)=exp(-beta*(epsilon(l1)+epsilon(in)));
        a_occ(in,l1)=1;
        a_occ(l1,in)=1;
        a_occ2(in,l1)=1;
        a_occ2(l1,in)=1;
        
        % Attach the new node in to the node l2;
        L=L+1;
        a(in,l2)=exp(-beta*(epsilon(l2)+epsilon(in)));
        a(l2,in)=exp(-beta*(epsilon(l2)+epsilon(in)));
        a_occ(in,l2)=1;
        a_occ(l2,in)=1;
        a_occ2(in,l2)=1;
        a_occ2(l2,in)=1;
        
    end
end

end