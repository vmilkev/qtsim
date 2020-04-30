function y = fregr(R, K, G1, n)

% Regulation factor (F) for repressor

    y = ( 1./(1 + (R.^2).*exp(-G1)/K) ).^(1/n);

end