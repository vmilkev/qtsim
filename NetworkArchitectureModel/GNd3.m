function [a] = GNd3(...
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

nt=0;
% Assign energies to the nodes
% If using Poisson and power-law you must define
% the parameters mu, or kappa
% Examples:
%mu=10;
% kappa=1;
for i=1:N
    epsilon(i)=floor(10*rand(1));
    % Alternative energy distributions
    %epsilon(i)=random('Poisson',mu);
    %poisson distribution with average mu
    %epsilon(i)=rand(1)ˆ(1/(kappa+1));
    %power-law distribution with exponent kappa
    
end

% Initial condition: at time t=1 a single tedrahedron (1,2,3,4)

for i1=1:4
    for i2=(i1+1):4
        a(i1,i2)=1;
        a(i2,i1)=1;
        for i3=(i2+1):4
            nt=nt+1;
            tri(nt,1)=i1;
            tri(nt,2)=i2;
            tri(nt,3)=i3;
            at(nt)=exp(-beta*(epsilon(i1)+epsilon(i2)+epsilon(i3)));
            a_occ(nt)=1;
            a_occ3(nt)=1;
        end
    end
end

% At each time t=in-3 we attach a new tetrahedron

for in=4+1:N
    % Choose triangular face to which to attach the new tetrahedron
    
    [I,J,V]=find(at.* a_occ);
    
    norm=sum(V);
    x=rand(1)*norm;
    for nj1=1:numel(V)
        x=x-V(nj1);
        if x<0
            nj=J(nj1);
            break;
        end
        
    end
    
    l(1)=tri(nj,1);
    l(2)=tri(nj,2);
    l(3)=tri(nj,3);
    
    a_occ(nj)=a_occ(nj)+s;
    a_occ3(nj)=a_occ3(nj)+1;
    
    %a_occ(nj)=a_occ(nj)*exp(-gamma*a_occ(nj));
    
    %Add the tethaedron
    for n=1:3
        a(in,l(n))=1;
        a(l(n),in)=1;
    end
    for n1=1:3
        for n2=n1+1:3
            a(l(n1),l(n2))=a(l(n1),l(n2))+1;
            a(l(n2),l(n1))=a(l(n2),l(n1))+1;
        end
    end
    for n=1:3
        for n2=n+1:3
            nt=nt+1;
            tri(nt,1)=l(n);
            tri(nt,2)=l(n2);
            tri(nt,3)=in;
            at(nt)=exp(-beta*(epsilon(l(n))+epsilon(l(n2))+epsilon(in)));
            %at(nt) = at(nt)*(1/(1+at(nt)/1000));
            a_occ(nt)=1;
            a_occ3(nt)=1;
        end
    end
    
end

end