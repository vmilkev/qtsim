function [a,b,rsquare,xPred,yPred] = fit_indegree(...
                                                 degree,...
                                                 counts...
                                                 )

% Approximation of in-Connectivity distribution data (degree & counts) by
% exponential fit.
% [a,b,rsquare,xPred,yPred] - parameters of yPred = a*exp(b*xPred);

[xData, yData] = prepareCurveData( degree, counts );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'Bisquare';

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

a = fitresult.a;
b = fitresult.b;
rsquare = gof.rsquare;
xPred = linspace(1,max(degree),100);
yPred = a*exp(b*xPred);

end