function X=LNLC_CS(filename,subrate,y,ims,A,At,par)

frame_num=length(ims);
[h, w]           =    size( ims{1} );

alpha=cell(frame_num,1);
gamma=cell(frame_num,1);
delta=cell(frame_num,1);

PSNR_LNLC=par.PSNR_LNLC;
SSIM_LNLC=par.SSIM_LNLC;
OutIter=par.OutIter;
InnerIter=par.InnerIter;
beta=par.beta;
sita=par.sita;

mu=par.mu;
RankList=zeros(300,'single');
RankLimitList=zeros(300,'single');

A_At=@(x,mode) getAAt(x,A,At,mode);

for frame_no=1:frame_num
    alpha{frame_no}=zeros(length(y{frame_num}),1);
    gamma{frame_no}=zeros(h,w);
    delta{frame_no}=zeros(h,w);
end;

X=ims;U=ims;Prev_X=ims;

cnt=0;
% AverageErr=zeros(OutIter+1,1);
% ErrList=zeros(OutIter*InnerIter+1,1);
timeBegin=tic;
flag=0;
for  k    =   1 :OutIter
    %     ErrSum=0;
    
    %     if k==4
    %         wkfile=sprintf('.\\temp\\%s_%2.2f_lamda=%2.2f_maxrank=%2.2f_OutIter=%d.mat',...
    %             par.filename,subrate,par.lamda,par.max_rank,k);
    %         save(wkfile);
    %     end;
    
    for it=1:InnerIter
        cnt       =   cnt  +  1;
        
        %TV sub_problem
        %         ErrPerIter=0;
        for frame_no=1:frame_num
            [X{frame_no},~,~]=TV_alp(A_At,y{frame_no},h,w,double(U{frame_no}),double(U{frame_no}),par.beta,alpha{frame_no},gamma{frame_no},par.tv_opts,par.tvflag);
            PSNR_LNLC(cnt,(par.GOPNo-1)*par.GOPSize+frame_no)    =   csnr( X{frame_no}, par.ori_im{frame_no}, 0, 0 );
            SSIM_LNLC(cnt,(par.GOPNo-1)*par.GOPSize+frame_no)    =  cal_ssim(  X{frame_no}, par.ori_im{frame_no}, 0, 0 );
            fprintf( 'LNLC Reconstruction: %s, subrate:%2.2f, frame:%d, Iter:%d, PSNR = %2.4f,SSIM = %2.4f \n',filename,subrate,(par.GOPNo-1)*par.GOPSize+frame_no, cnt,...
                PSNR_LNLC(cnt,(par.GOPNo-1)*par.GOPSize+frame_no) ,SSIM_LNLC(cnt,(par.GOPNo-1)*par.GOPSize+frame_no)  );
            %             Err=norm(X{frame_no}(:) - double(Prev_X{frame_no}(:)))/norm(double(Prev_X{frame_no}(:)));
            %             ErrPerIter=ErrPerIter+Err;
            %             ErrSum=ErrSum+Err;
        end;
        %         ErrList(cnt)=ErrPerIter/frame_num;
        %         Prev_X=X;
        if( it==1 )
            for frame_no=1:frame_num
                alpha{frame_no}=par.sf_alpha*alpha{frame_no}; %zeros(length(y{frame_num}),1);%0.0*alpha{frame_no};
                gamma{frame_no}=par.sf_gamma*gamma{frame_no};%zeros(h,w);
            end;
            pos_arr   = Sequence_Block_matching(X,par,0);
            U_arr        =     zeros(par.win^2, par.win^2,size(pos_arr,2), 'single');
        end;
        
        %  Low rank sub_problem
        [U, wei, U_arr,Rank1,Rank2,ReservedRate]      =   Video_E_LowRank(X, par, pos_arr, U_arr,it);%Extended low rank approximation
        
        if cnt==1        
             par.max_rank=min(Rank1/2.25,par.max_rank);
            if Rank1<10
                par.lamda=par.lamda/4;
                 par.tv_opts.mu=1.5*2^6;
            else
                if Rank1<15
                    par.lamda=par.lamda/2;
                    par.tv_opts.mu=2^7;
                end;
%                 par.max_rank=min(Rank1/2.25,par.max_rank);
            end;
            if Rank1>25
                par.lamda=1.5*par.lamda;
            end;
        end;
        %        RankList(cnt)=Rank2;
        %         RankLimitList(cnt)=Rank2;
        %         fprintf('Rank limit:%2.2f\n',par.max_rank)
        
        if k>4
            par.sf_alpha=1.0;
            par.sf_gamma=1.0;
            par.tv_opts.mu=2^10;
            if  (Rank1<max(par.o_max_rank,par.max_rank))
                par.lamda=par.lamda*par.sf_lamda;
            end;
        end;
        
        for frame_no=1:frame_num
            U{frame_no}     =    (mu*U{frame_no}-gamma{frame_no}+sita*double(X{frame_no}))./(mu*wei{frame_no}+sita);
            PSNR_LNLC_U  =   csnr( U{frame_no}, par.ori_im{frame_no}, 0, 0 );
            SSIM_LNLC_U    =  cal_ssim(  U{frame_no}, par.ori_im{frame_no}, 0, 0 );
            %                                     fprintf( 'fileName= %s, subrate=%2.2f, frame no= %d, LNLC Reconstruction, Iter %d,  U_PSNR = %f,U_SSIM=%2.4f \n',filename,subrate,frame_no, cnt,...
            %                                     PSNR_LNLC_U ,SSIM_LNLC_U  );
        end
%        if(ReservedRate>=0.99)
%            break;
%        end;
        %update Lagarange parameters
        for frame_no=1:frame_num
            alpha{frame_no}=alpha{frame_no}-beta*(y{frame_no}-A(X{frame_no}(:)));
            gamma{frame_no}=gamma{frame_no}-sita*(X{frame_no}-U{frame_no});
            %        delta{frame_no}=delta{frame_no}-sita*(U1{frame_no}-U{frame_no});
        end
    end;
%    if(ReservedRate>=0.99) break;  end;
    %     AverageErr(k+1)=ErrSum/(InnerIter*frame_num);
    %     fprintf( 'ErrSum=%2.6f, Average Err=%2.6f\n',ErrSum,AverageErr(k+1));
    %     save(['.\' filename '_AverageErr.mat'],'AverageErr');
    %     save(['.\result_PSNR\LNLC_' par.filestamp '_results.mat'],'PSNR_LNLC','SSIM_LNLC','PSNR_LNLC_U','SSIM_LNLC_U','par');
end

cnt=cnt+1;
for frame_no=1:frame_num
    [X{frame_no},~,~]=TV_alp(A_At,y{frame_no},h,w,double(U{frame_no}),double(U{frame_no}),par.beta,alpha{frame_no},gamma{frame_no},par.tv_opts,par.tvflag);
    PSNR_LNLC(cnt,(par.GOPNo-1)*par.GOPSize+frame_no)    =   csnr( X{frame_no}, par.ori_im{frame_no}, 0, 0 );
    SSIM_LNLC(cnt,(par.GOPNo-1)*par.GOPSize+frame_no)    =  cal_ssim(  X{frame_no}, par.ori_im{frame_no}, 0, 0 );
    fname=['.\result_img\LNLC_' filename '_' num2str(subrate) '_' num2str((par.GOPNo-1)*par.GOPSize+frame_no) '.bmp'];
    imwrite(uint8(X{frame_no}),fname,'bmp');
    fprintf( 'LNLC Reconstruction, %s:subrate: %2.2f,frame no, %d, Iter %d, PSNR = %f,SSIM = %2.4f \n',filename,subrate,(par.GOPNo-1)*par.GOPSize+frame_no, cnt,...
       PSNR_LNLC(cnt,(par.GOPNo-1)*par.GOPSize+frame_no) ,SSIM_LNLC(cnt,(par.GOPNo-1)*par.GOPSize+frame_no)  );
end;
elapsedTime=toc(timeBegin);
fprintf('\nLNLC Reconstruction total recovery time:%4.2f,recovery time per frame:%4.2f\n,',elapsedTime,elapsedTime/4);
save(['.\result_PSNR\LNLC_' par.filestamp '_results.mat'],'PSNR_LNLC','SSIM_LNLC','PSNR_LNLC_U','SSIM_LNLC_U','par','elapsedTime');
return;

function y=getAAt(x,A,At,mode)
if mode==1
    y=A(x);
else
    y=At(x);
end;
return;
