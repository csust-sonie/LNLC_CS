function  pos_arr  =  Sequence_Block_matching(ims,par,flag)
%,matched_blk_arr,L,flag=0时会选取本帧内的块，=1时只选GOP相邻帧的块
%%寻找与每个样本块最相似的par.nblk个块
S         =   40;
f         =   par.win;%块大小
f2        =   f^2;
s         =   par.step;%样本点间相隔距离，步长
ref_num=length(ims);
%分块
N         =   size(ims{1},1)-f+1;
M         =   size(ims{1},2)-f+1;
r         =   [1:s:N];
r         =   [r r(end)+1:N];
c         =   [1:s:M];
c         =   [c c(end)+1:M];
L         =   N*M;%块总数
X         =   cell(1,ref_num);
X1         =   zeros(f*f, L, 'single');
% pos_arr_f   =  zeros(par.nblk, N1*M1 );
% Sigma     =   zeros(L,1);
for frameNo=1:ref_num
    k    =  0;
    for i  = 1:f
        for j  = 1:f
            k    =  k+1;
            blk  =  ims{frameNo}(i:end-f+i,j:end-f+j);            
%             det_coefs=conv2(conv2(blk,[-1 1],'valid'),[-1;1],'valid')/(2*0.6745);
%             Sigma(k)=median(abs(det_coefs(:)));            
            X1(k,:) =  blk(:)';         
        end
    end
    X{frameNo}=X1';
   % X1(n,:,:)         =  T';
end;
% Index image
I     =   (1:L);
I     =   reshape(I, N, M);
N1    =   length(r);
M1    =   length(c);
%pos_arr_lst=cell(1,ref_num);

%matched_blk_arr= cell(1,N1*M1);
%X         =  X';
blk_pos=cell(M1,1);
blk_off=cell(M1,1);
pos_arr_t=cell(ref_num,1);
parfor frameNo=1:ref_num
%for frameNo=1:ref_num
     pos_arr_t{frameNo}=Block_Search(X,I,ref_num,frameNo,r,c,S,f2,N,M,N1,M1,par,flag);
%     for  i  =  1 : N1
%         for  j  =  1 : M1
          

%             row     =   r(i);%当前样本块最左上角的点在图像中的行号
%             col     =   c(j);%当前样本块最左上角的点在图像中的列号
% 
%             off     =  (col-1)*N + row;%当前样本块在X中的行号，X中一行就是一块
%             off1    =  (j-1)*N1 + i;%当前样本块在图像中的序号，竖向优先记号
%             
%             %确定对于当前样本块而言，可筛选的相似块的数量及位置
%             rmin    =   max( row-S, 1 );
%             rmax    =   min( row+S, N );
%             cmin    =   max( col-S, 1 );
%             cmax    =   min( col+S, M );
%             
%             idx     =   I(rmin:rmax, cmin:cmax);
%             idx     =   idx(:);%首轮选取的与当前样本块可能相似的所有样本块的序号
%             
%             idx1=idx;
%             B       =   X{n}(idx,:);%选取的与当前样本块可能相似的所有样本块
%             idx=idx+L*(n-1);
%             
%             if(flag==0)
%                 if(ref_num>1)
%                     for t=1:ref_num
%                         if t~=n
%                             T=X{t}(idx1,:);
%                             B =vertcat(B, T);
%                             idx=vertcat(idx,L*(t-1)+idx1);
%                         end;
%                     end;
%                 end;
%             else
%                 if(ref_num>1)
%                     B       = X{n}(off, :);%排除本帧内的其他块，只有当前块了
%                     idx=off+L*(n-1);
%                     for t=1:ref_num
%                         if t~=n
%                             T=X{t}(idx1,:);
%                             B =vertcat(B, T);
%                             idx=vertcat(idx,L*(t-1)+idx1);
%                         end;
%                     end;
%                 end;
%             end
%             v       =   X{n}(off, :);%当前样本块
%             
%             dis     =   (B(:,1) - v(1)).^2;
%             for k = 2:f2
%                 dis   =  dis + (B(:,k) - v(k)).^2;
%             end
%             dis   =  dis./f2;%每一个值对应可能相似块与当前样本块像素值的均方差
%             [val,ind]   =  sort(dis);
%             pos_arr_f(:,off1)  =  idx( ind(1:par.nblk) );%选取与当前样本块最相似的par.nblk-1个块作为最终的匹配块，记录它们在图像中的位置序号（第几个块）

%         end
%     end
end;
for frameNo=1:ref_num
    if frameNo==1
        pos_arr=pos_arr_t{1};
    else
        pos_arr=[pos_arr pos_arr_t{frameNo}];
    end;
end;

