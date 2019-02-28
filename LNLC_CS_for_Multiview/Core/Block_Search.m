function pos_arr_f=Block_Search(X,I,ref_num,frameNo,r,c,S,f2,N,M,N1,M1,par,flag)
L=N*M;
pos_arr_f   =  zeros(par.nblk, N1*M1 );
for  i  =  1 : N1
    for  j  =  1 : M1
        row     =   r(i);%��ǰ�����������Ͻǵĵ���ͼ���е��к�
        col     =   c(j);%��ǰ�����������Ͻǵĵ���ͼ���е��к�
        
        off     =  (col-1)*N + row;%��ǰ��������X�е��кţ�X��һ�о���һ��
        off1    =  (j-1)*N1 + i;%��ǰ��������ͼ���е���ţ��������ȼǺ�
        
        %ȷ�����ڵ�ǰ��������ԣ���ɸѡ�����ƿ��������λ��
        rmin    =   max( row-S, 1 );
        rmax    =   min( row+S, N );
        cmin    =   max( col-S, 1 );
        cmax    =   min( col+S, M );
        
        idx     =   I(rmin:rmax, cmin:cmax);
        idx     =   idx(:);%����ѡȡ���뵱ǰ������������Ƶ���������������
        
        idx1=idx;
        B       =   X{frameNo}(idx,:);%ѡȡ���뵱ǰ������������Ƶ�����������
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
                B       = X{frameNo}(off, :);%�ų���֡�ڵ������飬ֻ�е�ǰ����
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
        v       =   X{frameNo}(off, :);%��ǰ������
        
        dis     =   (B(:,1) - v(1)).^2;
        for k = 2:f2
            dis   =  dis + (B(:,k) - v(k)).^2;
        end
        %             dis   =  dis./f2;%ÿһ��ֵ��Ӧ�������ƿ��뵱ǰ����������ֵ�ľ�����,û��Ҫ��,2017/12/25
        [val,ind]   =  sort(dis);
         pos_arr_f(:,off1)  =  idx( ind(1:par.nblk) );%ѡȡ�뵱ǰ�����������Ƶ�par.nblk-1������Ϊ���յ�ƥ��飬��¼������ͼ���е�λ����ţ��ڼ����飩
    end;
end;
end