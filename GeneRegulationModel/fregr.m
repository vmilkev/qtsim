function y = fregr(c_r, Kr, Gr, n)

% y is a regulation factor (F) for repressor
% c_r is a concentration of repressor's molecules
% Gr is a repressor's attracting energy
% Kr is a repressor's attracting energy constant
% n is a total number of different repressors for a gene

    y = ( 1./(1 + (c_r.^2).*exp(-Gr)/Kr) ).^(1/n);

end