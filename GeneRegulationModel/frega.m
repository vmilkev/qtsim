function y = frega(A, K, G1, f, n)

% Regulation factor (F) for activator

    y = ( (1 + (A.^2).*exp(-G1)*f/K)./(1 + (A.^2).*exp(-G1)/K) ).^(1/n);

end