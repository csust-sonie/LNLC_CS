function  Rec_im0    =   DCT_BCS_CIR( y, par, A, At,blkSize )
ori_im      =    par.ori_im;
[h w]       =    size(ori_im);
x          =    At( y );%����Ҷ��任����������任ǰֵΪ����Ҷ�任ϵ�����ǲ�����ȡ0ֵ
x          =    reshape(x, blkSize*blkSize,(h*w)/(blkSize*blkSize));
im          =   col2im(x,[blkSize blkSize],[h w],'distinct');%����Ҷ��任��Ļָ�ͼ��

lamada      =    1.5;  % 1.8, 1.2-1.7
b           =    par.win*par.win;
D           =    dctmtx(b);%����DCT�任����

for k   =  1:1
    f      =   im;
    for  iter = 1 : 300   
        
        if (mod(iter, 50) == 0)
            if isfield(par,'ori_im')
                PSNR     =   csnr( f, par.ori_im, 0, 0 );                
                fprintf( 'DCT Compressive Image Recovery, Iter %d : PSNR = %f\n', iter, PSNR );
            end
        end
        
        %ѭ���任����任����С�ؽ�ֵ�����ֵ�����
        for ii = 1 : 3
            fb        =   A( f(:) );
            t=At( y-fb );
            t=reshape(t, blkSize*blkSize,(h*w)/(blkSize*blkSize));
            t=col2im(t,[blkSize blkSize],[h w],'distinct');
            f         =   f + lamada.* t;
            %f         =   f + lamada.*reshape(At( y-fb ), h, w);
        end 
        
         %DCT��ֵ������ǿ�ؽ�Ч��
        f          =   DCT_thresholding( f, par, D );
    end
    im     =  f;
end
Rec_im0   =  im;
return;