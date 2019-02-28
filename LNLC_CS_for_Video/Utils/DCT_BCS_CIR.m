function  Rec_im0    =   DCT_BCS_CIR( y, par, A, At,blkSize )
ori_im      =    par.ori_im;
[h w]       =    size(ori_im);
x          =    At( y );%傅里叶逆变换，采样点逆变换前值为傅里叶变换系数，非采样点取0值
x          =    reshape(x, blkSize*blkSize,(h*w)/(blkSize*blkSize));
im          =   col2im(x,[blkSize blkSize],[h w],'distinct');%傅里叶逆变换后的恢复图像

lamada      =    1.5;  % 1.8, 1.2-1.7
b           =    par.win*par.win;
D           =    dctmtx(b);%生成DCT变换矩阵

for k   =  1:1
    f      =   im;
    for  iter = 1 : 300   
        
        if (mod(iter, 50) == 0)
            if isfield(par,'ori_im')
                PSNR     =   csnr( f, par.ori_im, 0, 0 );                
                fprintf( 'DCT Compressive Image Recovery, Iter %d : PSNR = %f\n', iter, PSNR );
            end
        end
        
        %循环变换与逆变换，减小重建值与测量值的误差
        for ii = 1 : 3
            fb        =   A( f(:) );
            t=At( y-fb );
            t=reshape(t, blkSize*blkSize,(h*w)/(blkSize*blkSize));
            t=col2im(t,[blkSize blkSize],[h w],'distinct');
            f         =   f + lamada.* t;
            %f         =   f + lamada.*reshape(At( y-fb ), h, w);
        end 
        
         %DCT阈值法，增强重建效果
        f          =   DCT_thresholding( f, par, D );
    end
    im     =  f;
end
Rec_im0   =  im;
return;