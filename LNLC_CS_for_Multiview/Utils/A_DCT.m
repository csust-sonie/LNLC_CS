function y= A_DCT(x,Phi,h,w)
% [M N] = size(Phi);
% block_size = sqrt(N);
% [im_rows im_cols] = size(im);
% x = im2col(im, [block_size block_size], 'distinct');
% t = Phi * x;
% L=size(t,1)*size(t,2);
% y=reshape(t,L,1);

[M N] = size(Phi);
im=reshape(x,h,w);
block_size = sqrt(N);
t = im2col(im, [block_size block_size], 'distinct');
t = Phi * t;
L=size(t,1)*size(t,2);
y=reshape(t,L,1);