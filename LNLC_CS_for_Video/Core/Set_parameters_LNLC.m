function  par  =  Set_parameters_LNLC(rate,cs_model)  
% mode:1: for SRM+TV and 2 for BCS
par.win        =    6;    % Patch size 6;
par.nblk       =    45;   % number of the blocks;  
par.step       =    5;%min(6,par.win-1);

par.OutIter    =    10;
par.InnerIter  =    10;
par.tvflag     =    1;  %if alpha/delta is update at the outer iteration or in the tv iteration. 

par.sf_alpha   =    0.45; %shrink factor for alpha
par.sf_gamma   =    0.45; %shrink factor fr gamma 
par.sf_lamda   =    0.975;%shrink factor for lamda;

par.mu         =    1.0;
par.sita       =    0.025;  %  
par.beta       =    0.0015;

par.lamda      =   35;
par.max_rank   =11;
    if rate<=0.05
        par.max_rank   = 3;
        par.lamda     =  62;
     elseif rate<=0.1    
         par.lamda     =  42;%42;
         par.max_rank   =  5;
    elseif rate<=0.15                 
         par.lamda     =  28;%ºÃÏóÆ«´ó35
        par.max_rank   =  7;
    elseif rate<=0.2        
        par.lamda     =  15;
        par.max_rank   =  9;
     elseif rate<=0.25  
        par.lamda     =  10;
        par.max_rank   =  11;
    elseif rate<=0.3
        par.lamda      = 5;%2.916;%9;%14%16;%28.097;%28.097;%32;%32; %32;  %28.097;%for 9
        par.max_rank   = 15;%11;%11;%8;%8,30,0.02   
     else
        par.max_rank   = 16;%11;%11;%8;%8,30,0.02
        par.lamda      = 3;%2.916;%9;%14%16;%28.097;%28.097;%32;%32; %32;  %28.097;%for 9
    end    
    
par.o_lamda=par.lamda;
par.o_max_rank=par.max_rank;

% opts.mu = 2^8;
% opts.beta = 2^5;
par.tv_opts.mu =  1.5*2^8;%2^8;
par.tv_opts.beta = 2^5;

par.tv_opts.tol = 1E-8;
par.tv_opts.maxit = 300;
par.tv_opts.TVnorm = 1;
par.tv_opts.nonneg = false;

return;