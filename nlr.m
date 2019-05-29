% --- Nonlinear regression based on initial estimates from fitting

function [pdfNLR, alphNLR, deltaP] = nlr(alph, C, x, y)

% Function to fit: power law
modelfun = @(c,x) c(1)*x.^(-c(2));
% Initial guesses of coefficients
beta0 = [C, alph];
% Get model by NLR
mdl = fitnlm(x, y, modelfun, beta0);
cfs = mdl.Coefficients{:, 'Estimate'};
alphNLR = cfs(2);
% PDF
pdfNLR = x.^(-cfs(2))/sum(x.^(-cfs(2)));