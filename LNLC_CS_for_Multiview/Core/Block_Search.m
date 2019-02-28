function pos_arr_f=Block_Search(X,I,ref_num,frameNo,r,c,S,f2,N,M,N1,M1,par,flag)
L=N*M;
pos_arr_f   =  zeros(par.nblk, N1*M1 );
for  i  =  1 : N1
    for  j  =  1 : M1
        row     =   r(i);%当前样本块最左上角的点在图像中的行号
        col     =   c(j);%当前样本块最左上角的点在图像中的列号
        
        off     =  (col-1)*N + row;%当前样本块在X中的行号，X中一行就是一块
        off1    =  (j-1)*N1 + i;%当前样本块在图像中的序号，竖向优先记号
        
        %确定对于当前样本块而言，可筛选的相似块的数量及位置
        rmin    =   max( row-S, 1 );
        rmax    =   min( row+S, N );
        cmin    =   max( col-S, 1 );
        cmax    =   min( col+S, M );
        
        idx     =   I(rmin:rmax, cmin:cmax);
        idx     =   idx(:);%首轮选取的与当前样本块可能相似的所有样本块的序号
        
        idx1=idx;
        B       =   X{frameNo}(idx,:);%选取的与当前样本块可能相似的所有样本块
        idx=idx+L*(frameNo-1);
        
        if(flag==0)
            if(ref_num>1)
                for t=1:ref_num
                    if t~=frameNo
                        T=X{t}(idx1,:);
                        B =vertcat(B, T);
                        idx=vertcat(idx,L*(t-1)+idx1);
                    end;
                end;
            end;
        else
            if(ref_num>1)
                B       = X{frameNo}(off, :);%排除本帧内的其他块，只有当前块了
                idx=off+L*(frameNo-1); 
                for t=1:ref_num
                    if t~=frameNo
                        T=X{t}(idx1,:);
                        B =vertcat(B, T);
                        idx=vertcat(idx,L*(t-1)+idx1);
                    end;
                end;
            end;
        end
        v       =   X{frameNo}(off, :);%当前样本块
        
        dis     =   (B(:,1) - v(1)).^2;
        for k = 2:f2
            dis   =  dis + (B(:,k) - v(k)).^2;
        end
        %             dis   =  dis./f2;%每一个值对应可能相似块与当前样本块像素值的均方差,没必要吧,2017/12/25
        [val,ind]   =  sort(dis);
         pos_arr_f(:,off1)  =  idx( ind(1:par.nblk) );%选取与当前样本块最相似的par.nblk-1个块作为最终的匹配块，记录它们在图像中的位置序号（第几个块）
    end;
end;
end