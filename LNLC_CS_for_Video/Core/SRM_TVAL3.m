function [rec_img, psnr_val, out_y, Phi, Phi_T] = SRM_TVAL3(im, K, trans_mode, rand_type, blk_size, Iteration)
%  ===== Required inputs ============= 
% im : filename of the input image

% K: number of measurements
    % K can be either a vector or a scalar; 
    % If K is a scalar, it will show the reconstructed image;
    % If K is a vector, the function will produce a R-D curve;  
    
% trans_mode: measurement matrix type
    % 'PFFT': Partial FFT without pre-randomization;
    % 'FFT': FFT with pre-randomization;
    % 'BDCT': Block DCT with pre-randomization;
    % 'BWHT': Block Walsh-Hadamard Transform (WHT) with pre-randomization;
    
% rand_type: type of pre-randomization operator. 
    %0=random permutation (global model);
    %1=random flipping the signs of the signal(local model); 

%  ===== Optional Input =============
% blk_size: block size of the measurement operator if the trans_mode is either 
%           'BDCT' or 'BWHT'; 
%           Default: 32 for random permutation (rand_type=0) and row_num of input image
%           for random flipping the signs (rand_type=1);

% **********************************************************************

% ===== Outputs =============
% rec_img: reconstructed image if K is a scalar; []if K is a vector; 
% psnr_val: the PSNR of reconstruced image;  
% out_y: measurement
% Phi/Phit: measurement matrix/ function handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Author: Thong Do, JHU
% modified by: Kan Chang, Guangxi University
% modification: TVAL3 is called here to perform independent reconstruction
% Data: 2014/7/8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Test the input parameters:
if nargin<4,
    help fast_cs2d;
    error('Wrong Number of Input Parameters');
end


% Read the input image;            
x = double(im);
% keep an original copy of the input signal
x0 = x;


[m n] = size(x);
% Total number of pixels;
N = m*n;
x = x(:);

% Substract the mean of the input signal; 
% xmean = mean(x);
% x = x-xmean;
% mean_V = xmean;

% Set the default of blk_size for BDCT or BWHT if it is not specified; 
if (strcmp(trans_mode,'BDCT') || strcmp(trans_mode,'BWHT')) && (nargin<5)
   if rand_type==0,
       blk_size=32;
   else
       blk_size=m;
   end
end

% Initialize the output parameters;
rec_img=[];
psnr_val=[];

% Choose an arbitrary random seed; 
% User can change it
Perm_state = 3587642;
rand('state', Perm_state);
  
% Define the random vector for the input samples: 

% Random permutation;
if rand_type == 0
    if strcmp(trans_mode,'PFFT'),
        rand_vect=1:N; %No random permutation for PFFT;
    else
   % other modes: randomize samples use permuation vector
        rand_vect = randperm(N)';
    end
% Random filipping the sign;    
elseif rand_type == 1
    if strcmp(trans_mode,'PFFT'),
        rand_vect=ones(size(x)); %No sign flipping for PFFT;
    else
        % Other modes: randomize samples using a Bernoulli vector
        rand_vect = 2*round(rand(N,1))-1;
    end
end % of if rand_type

% Main loop to test the reconstruction performance for each K(i);
 for i =1:length(K),
       Ki=round(K(i)); 
       
       % Dense FFT-based Operators;
        if strcmp(trans_mode,'PFFT')|| strcmp(trans_mode,'FFT')
            % Define selected samples
            select_vect = randperm(round(N/2)-1)+1;
            select_vect = select_vect(1:round(Ki/2))';
            % Define Sampling Operator;
            Phi = @(z) fft1d_f(z, select_vect, rand_vect);
            % Define the transpose of the Sampling Operator;
            Phi_T = @(z) fft1d_t(z, N, select_vect, rand_vect);  
            
        % Block-Based Operators:
        elseif strcmp(trans_mode,'BDCT') || strcmp(trans_mode,'BWHT')
                        
            % Define selected samples
            select_vect = randperm(N);
            select_vect = select_vect(1:Ki)'; 
            
            % Define Sampling Operator;    
            Phi = @(z) blk_f1d(z, select_vect, rand_vect, trans_mode, blk_size);
            % Define the transpose of the Sampling Operator;
            Phi_T = @(z) blk_t1d(z, N, select_vect, rand_vect,trans_mode, blk_size); 
            
        end

   % getting measurements
   y = Phi(x);  
   out_y =y;   


   % TVAL3 reconstruction
   clear opts
   opts.mu = 2^12;
   opts.beta = 2^6;
   opts.mu0 = 2^4;       % trigger continuation shceme
   opts.beta0 = 2^-2;    % trigger continuation shceme
   opts.maxcnt = 10;
   opts.tol_inn = 1e-3;
   opts.tol = 1E-6;
   opts.maxit = Iteration;    
 
   A = @(z,mode) dfA(z,mode,Phi,Phi_T);
   [xr, ~] = TVAL3(A,y,m,n,opts);

    psnr_val  = psnr(x0, xr);
 end
 
% Retrun the reconstructed image if K is a scalar; 
if length(K)==1,
    rec_img=uint8(xr);
end

function y = dfA(z,mode,Phi,Phit)
switch mode
    case 1
        y = Phi(z);
    case 2
        y = Phit(z);
    otherwise
        error('Unknown mode passed to f_handleA!');
end