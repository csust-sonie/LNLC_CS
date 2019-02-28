function X= At_DCS_B(y,Phi,h,w)
[M N] = size(Phi);
block_size = sqrt(N);
L=length(y);
ty=reshape(y,M,L/M);
tx= Phi'*ty;
L=size(tx,1)*size(tx,2);
x=reshape(tx,L,1);
x=reshape(x,block_size*block_size,(h*w)/(block_size*block_size));
X=col2im(x,[block_size block_size],[h w],'distinct');
X=reshape(X,L,1);

