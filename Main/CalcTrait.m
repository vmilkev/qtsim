function [ trait, v_trait ] = CalcTrait( xi, yi, coreMap, coreWeit, nodes, skip, r_trait )

if (skip ~= 0.0)
    szX = size(xi);
    st = floor(skip*szX(1,2));
    ed = szX(1,2);
end

k = keys( coreMap );
k_sz = size(k);

Yp = zeros( k_sz(1,2), ed-st+1 );
mYp = zeros( k_sz(1,2), 1 );
v_trait = zeros( k_sz(1,2), 1 );
Xp = xi(1,st:ed);

if ( isempty(r_trait) )
    r_trait = ones(k_sz(1,2),1);
end

trait = 0;

for ik = 1:k_sz(1,2)
    ip = k{ik} + nodes;
    Yp(ik,:) = yi(ip,st:ed);
    mYp(ik,1) = mean( yi(ip,st:ed));
    v_trait(ik,1) = ( mYp(ik,1) * coreWeit( coreMap(k{ik}) ) );
    trait = trait + v_trait(ik,1)/r_trait(ik,1);
end

trait = trait/k_sz(1,2);


end
