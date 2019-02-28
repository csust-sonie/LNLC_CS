function   [frames, wei, U_arr,FirstIterRank,AverageRank,ReservedRate]  =  Video_E_LowRank(ims, par, pos_arr, U_arr,iter)

ref_num      =   length(ims);
b            =   par.win;
[h  w ch]    =   size(ims{1});

N            =   h-b+1;
M            =   w-b+1;

r            =   [1:N];
c            =   [1:M];

X            =   Im2Patch( ims{1}, par );
for i=2:ref_num
    X            =   [X Im2Patch( ims{i}, par )];
end;

Ys           =   zeros( size(X) );
W            =   zeros( size(X) );
L            =   size(pos_arr,2);
% fprintf('Average Rank:');

max_rank     =   par.max_rank;
V_arr        =   zeros(par.nblk,par.win^2, L,'single');
Sigma_arr    =    zeros(par.win^2,L, 'single');

ac           =   1.0;

for it=1:2
    CountR          =  0;
    SingularSum0    =  0;
    SingularSum1    =  0;
 

    for  i  =  1 : L
        B          =   X(:, pos_arr(:, i));
        if mod(iter,2) == 1    %mod(iter,5) == 1
            [tmp_y, tmp_w, U_arr(:,:,i),V_arr(:,:,i),Sigma_arr(:,i),R,S0,S1]   =   A_Weighted_SVT( double(B), U_arr(:,:,i),V_arr(:,:,i),Sigma_arr(:,i),par.lamda ,ac,max_rank,it);
        else
            [tmp_y, tmp_w, R,S0,S1]   =   A_Weighted_SVT_fast( double(B),U_arr(:,:,i),par.lamda ,ac,max_rank);
        end;
        CountR                =   CountR+R;
        SingularSum0          =   SingularSum0+S0;
        SingularSum1          =   SingularSum1+S1;
        Ys(:, pos_arr(:,i))   =   Ys(:, pos_arr(:,i)) + tmp_y;
        W(:, pos_arr(:,i))    =   W(:, pos_arr(:,i)) + tmp_w; 
    end
    AverageRank    =  CountR/L;
    if it == 1
        FirstIterRank  =   AverageRank;
    end;
    %      fprintf('Average Rank = %2.2f\n',AverageRank);
    ReservedRate   =   SingularSum1/SingularSum0;
%     fprintf(' rank:%2.2f,singular ratio:%2.3f \n',AverageRank,SingularSum1/SingularSum0);
    ac             =   AverageRank/max_rank;
    %    fprintf('.');
    if ac < 1.1
        break;
    end;
end;
%  fprintf('\nmax_rank:%2.2f\n',max_rank);

Ys_Byframe  =  cell(1,ref_num);
W_Byframe   =  cell(1,ref_num);
fl          =  size(Ys,2)/ref_num;

for ref_no  =  1:ref_num
    Ys_Byframe{ref_no}    =  Ys(:,1+fl*(ref_no-1):fl+fl*(ref_no-1));
    W_Byframe{ref_no}     =  W(:,1+fl*(ref_no-1):fl+fl*(ref_no-1));
end;

frames    =  cell(1,ref_num);
wei       =  cell(1,ref_num);
for ref_no=1:ref_num
    frames{ref_no}     =  zeros(h,w);
    wei{ref_no}        =  zeros(h,w); 
    k    =  0;
    for i  = 1:b
        for j  = 1:b
            k    =  k+1;
            frames{ref_no}(r-1+i,c-1+j)  =  frames{ref_no}(r-1+i,c-1+j) + reshape( Ys_Byframe{ref_no}(k,:)', [N M]);
            wei{ref_no}(r-1+i,c-1+j)  =  wei{ref_no}(r-1+i,c-1+j) + reshape( W_Byframe{ref_no}(k,:)', [N M]);
        end
    end
    
end;

% parallel implement

% function  [Ys, W, U,V,Sigma,CountR,SingularSum0,SingularSum1] ...
%     =    Par_Weighted_SVT(Ys,W,pos_arr,k,L, Y, U,V,Sigma,lamda,ac,rth,it)
% CountR         =  0;
% SingularSum0   =  0;
% SingularSum1=0;
% for j  =  1 : L/4
%     i          =  (k-1)*L/4+j;
%     B          =  Y(:, pos_arr(:, i));
%     [tmp_y, tmp_w, U(:,i),V(:,i),Sigma(:,i),R,S0,S1]   ...
%         =   A_Weighted_SVT( double(B), U(:,i),V(:,i),Sigma(:,i),lamda ,ac,rth,it);
%     CountR                =   CountR+R;
%     SingularSum0          =   SingularSum0+S0;
%     SingularSum1          =   SingularSum1+S1;
%     Ys(:, pos_arr(:,i))   =   Ys(:, pos_arr(:,i)) + tmp_y;
%     W(:, pos_arr(:,i))    =   W(:, pos_arr(:,i)) + tmp_w;
% end;
% return;

% function [Ys, W,U_arr]=Par_A_Weighted_SVT(X, pos_arr,i,par,ac,Ys, W)
% B          =   X(:, pos_arr(:, i));
%
% [tmp_y, tmp_w, U_arr(:,i)]   =   A_Weighted_SVT( double(B), par.lamda ,ac,par.max_rank);
%
% Ys(:, pos_arr(:,i))   =   Ys(:, pos_arr(:,i)) + tmp_y;
% W(:, pos_arr(:,i))    =   W(:, pos_arr(:,i)) + tmp_w;
% return;

%--------------------------------------------------------------------------
%- This function uses the PCA matrixes obtained in the previous iterations
%- to save computational complexity
%--------------------------------------------------------------------------
function  [X,W,r,SingularSum0,SingularSum1]   =   A_Weighted_SVT_fast( Y, U0,lamda, ac,rth)

% n                 =   sqrt(length(U0));
% U0                =   reshape(U0, n, n);
A                 =   U0'*Y;
Sigma0            =   sqrt( sum(A.^2, 2) );
V0                =   (diag(1./Sigma0)*A)';
SingularSum0      =   sum(Sigma0);
S                 =   max( Sigma0.^2/size(Y, 2), 0 );
thr               =   lamda./ ( sqrt(S) + eps );
S                 =   soft(Sigma0, thr);

if ac>1.0
    ac            =  min(2.5,ac);
    r             =  sum(S>0);
    if r > rth
        if r >= 1
            r     =  round(r/ac);
            S(r+1:end)  = 0;
        end;
    end;
end;
SingularSum1      =   sum(S);

r                 =     sum( S>0 );
P                 =     find(S);
X                 =     U0(:,P)*diag(S(P))*V0(:,P)';

W                 =     ones( size(X) );
%X                 =     X*wei;
return;


function  [X, W, U,V,Sigma,r,SingularSum0,SingularSum1]   =   A_Weighted_SVT( Y, U,V,Sigma,lamda,ac,rth,it)
lamda                =   lamda*sqrt(2);
if it == 1
    [U0,Sigma0,V0]    =   svd(full(Y),'econ');
    Sigma0            =   diag(Sigma0);
else
    % SVD isn't needed for second compute
    U0                =   U;
    V0                =   V;
    Sigma0            =   Sigma;
end;

SingularSum0      =   sum(Sigma0);
S                 =   max( Sigma0.^2/size(Y, 2), 0 );
thr               =   lamda./ ( sqrt(S) + eps );
S                 =   soft(Sigma0, thr);

if ac > 1.0
    ac            =   min(2.5,ac);
    r             =   sum(S>0);
    if r > rth
        if  r >= 1
            r           =  round(r/ac);
            S(r+1:end)  =  0;
        end;
    end;
end;
SingularSum1      =   sum(S);
r                 =   sum( S>0 );

U                 =   U0(:,1:r);
V                 =   V0(:,1:r);
X                 =   U*diag(S(1:r))*V';

W                 =   ones( size(X) );
%X                 =   X*wei;%
U                 =   U0 ;
V                 =   V0 ;
Sigma             =   Sigma0(:);
return;
