function y= A_DCS_B(X,Phi,h,w)
% [M N] = size(Phi);
% block_size = sqrt(N);
% [im_rows im_cols] = size(im);
% x = im2col(im, [block_size block_size], 'distinct');
% t = Phi * x;
% L=size(t,1)*size(t,2);
% y=reshape(t,L,1);

[M N] = size(Phi);
block_size = sqrt(N);
t = im2col(X, [block_size block_size], 'distinct');
t = Phi * t;
L=size(t,1)*size(t,2);
y=reshape(t,L,1);