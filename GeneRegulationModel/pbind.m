function y = pbind(K, deltaG, F)

% Function which calculates a binding probability of polymerase II complex.
%
% Parameters:
%
% deltaG -polymerase II binding Gibbs energy, in units [kT].
%
% K - binding constant.
%
% F - combined regulation factor (F = F1*F2*...*Fn).
%
% y - binding probability normalised such as y = [0,2], so
%     rate*y = [0,rate] if there is repressor F and
%     rate*y = [rate,2*rate] if there is activator F.

y = 2/(1+K.*exp(deltaG)./F);

end