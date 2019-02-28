function [X,delta,gamma]=TV_alp(A,b,p,q,XX,uu,sita,alpha,gamma,tv_opts,flag)
%% Run TV
% clear opts
% % opts.mu = 2^8;
% % opts.beta = 2^5;
% opts.mu = 2^9;%2^8;
% opts.beta = 2^5;
% 
% opts.tol = 1E-8;
% opts.maxit = 300;
% opts.TVnorm = 1;
% opts.nonneg = false;

t = cputime;
[U, delta,gamma,out] = Modified_ftvcs_alp(A,b,p,q,tv_opts,XX,uu,sita,alpha,gamma,flag);
t = cputime - t;
X=U;