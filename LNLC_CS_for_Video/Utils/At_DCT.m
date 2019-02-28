function x= At_DCT(y,Phi)
[M N] = size(Phi);
L=length(y);
ty=reshape(y,M,L/M);
tx= Phi'*ty;
L=size(tx,1)*size(tx,2);
x=reshape(tx,L,1);
