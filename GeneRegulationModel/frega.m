function y = frega(c_a, Ka, Ga, f, n)

% y is a regulation factor (F) for activator
% c_a is a concentration of activator's molecules
% K is an activator's constant
% Ga is an acttivator's attracting energy
% f is a between acttivators attracting energy
% n is a total number of different activators for a gene

    y = ( (1 + (c_a.^2).*exp(-Ga)*f/Ka)./(1 + (c_a.^2).*exp(-Ga)/Ka) ).^(1/n);

end