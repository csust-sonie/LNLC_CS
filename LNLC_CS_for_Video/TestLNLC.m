% Authors: 
% 
%   Yun Song,
%   Associate Professor,
%   School of Computer and Communication Engineering, 
%   Changsha University of Science and Technology
%   Email: sonie@126.com
%
%   Dengyong Zhang,
%   Dr.,
%   School of Computer and Communication Engineering, 
%   Changsha University of Science and Technology
%   Email: zhdy@csust.edu.cn
%
%   Qiang Tang,
%   Dr.,
%   School of Computer and Communication Engineering, 
%   Changsha University of Science and Technology
%   Email: tangqiang@csust.edu.cn
%
%   The code is associated with the following paper: 
%   Yun Song, Dongyong Zhang, and Baoxin Li, "Local and Nonlocal Constraints for Compressed Sensing Recovery," 
%   
%   Please cite this paper if you use this code. 
% 
%   For further information, please contact: sonie@126.com



csmode=1;  %1 for SRM, 2 for BCS

addpath('./yuv');
addpath('./SRM');
addpath('./Utils');
addpath('./Ext_TVAL3');
addpath('./Ext_TVAL3/Solver');
addpath('./Ext_TVAL3/Utilities');
addpath('./Core');

imgHeight=288;
imgWidth=352;

GOPSize=4;
total_num_frames=24;

Sequences={'mother-daughter_cif.yuv','Foreman_cif.yuv','BUS_352x288_420.yuv','mobile_cif.yuv','paris.yuv','football_cif.yuv'};
SequenceName='BUS_352x288_420.yuv';%'Foreman_cif.yuv';%;
subrates={0.05,0.10,0.15,0.20,0.25,0.30};

reconstructedImage=cell(1,GOPSize);
x=cell(1,GOPSize);
y=cell(1,GOPSize);
y_vectored=cell(1,GOPSize);

All_PSNR=zeros(GOPSize);

PSNR_LNLC=zeros(500,total_num_frames);
SSIM_LNLC=zeros(500,total_num_frames);

PSNR_X=zeros(total_num_frames,1);
SSIM_X=zeros(total_num_frames,1);
for s=3:3
    SequenceName=Sequences{s};
    for r=3:3
        subrate=subrates{r};
        switch csmode
            case 1
                N=imgHeight*imgWidth;
                M=round(subrate *N);
            case 2
                blkSize=32;
                load('phi_1024.mat');
                Phi_N=Phi_1024;
                % blkSize=16;
                % load('phi_256.mat');
                % Phi_N=Phi;
                N = blkSize * blkSize;
                M = round(subrate * N);
                Phi = Phi_N(1:M, :);
                A = @(z) A_DCT(z,Phi,imgHeight,imgWidth);
                At=@(z) At_DCS_B(z,Phi,imgHeight,imgWidth);
        end;
        GOPNo=0;
        for k = 1:GOPSize:total_num_frames
            GOPNo=floor(k/4)+1;
            par = Set_parameters_LNLC(subrate,csmode);
            par.GOPNo=GOPNo;
            par.GOPSize=GOPSize;
            originImages=readyuv2(SequenceName,k,GOPSize,imgHeight,imgWidth,'420');
            for i = 1 : GOPSize
                switch csmode
                    case 1  %  %SRM+TVAL3
                        [reconstructedImage{i}, All_PSNR(i), y{i}, A, At] = SRM_TVAL3(originImages{i}, M, 'BWHT', 0, 32, 150);
                        y_vectored{i}=y{i}(:);
                        psnr_name=['.\result_psnr\SRM_TVAL3_' SequenceName '_' num2str(subrate) '_results.mat'];
                        img_name=['.\result_img\SRM_TVAL3_' SequenceName '_' num2str(subrate) '_' num2str((par.GOPNo-1)*par.GOPSize+i) '.bmp'];
                    case 2   %BCS+BCS_SPL
                        x{i}= im2col(originImages{i}, [blkSize blkSize], 'distinct');
                        y{i}=Phi*x{i};%y{i}=A(x);
                        y_vectored{i}=y{i}(:);
                        reconstructedImage{i} = BCS_SPL_DCT_Decoder(y{i}, Phi,imgHeight, imgWidth);
                        psnr_name=['.\result_psnr\BCS_SPL_DCT_' SequenceName '_' num2str(subrate) '_results.mat'];
                        img_name=['.\result_img\BCS_SPL_DCT_' SequenceName '_' num2str(subrate) '_' num2str((par.GOPNo-1)*par.GOPSize+i) '.bmp'];
                end;
                PSNR_X((GOPNo-1)*GOPSize+i)    =   csnr( double(reconstructedImage{i}), originImages{i}, 0, 0 );
                SSIM_X((GOPNo-1)*GOPSize+i)   =  cal_ssim( double(reconstructedImage{i}), originImages{i}, 0, 0 );
                fprintf( 'fileName= %s, subrate=%2.2f, frame no= %d, DCT Reconstruction, X_PSNR = %f,X_SSIM=%2.4f \n',SequenceName,subrate,(par.GOPNo-1)*par.GOPSize+i,PSNR_X((GOPNo-1)*GOPSize+i) ,SSIM_X((GOPNo-1)*GOPSize+i));
               
                imwrite(uint8(reconstructedImage{i}),img_name,'bmp');
            end;
            save(psnr_name,'PSNR_X','SSIM_X');
            
            par.total_num_frames=total_num_frames;
            
            par.PSNR_LNLC=PSNR_LNLC;
            par.SSIM_LNLC=SSIM_LNLC;
            par.filename=SequenceName;
            par.subrate=subrate;
            par.filestamp=sprintf('%s_%2.2f_%d-%d',par.filename,subrate,(GOPNo-1)*GOPSize+1,GOPNo*GOPSize);
            par.ori_im=originImages;
            par.y=y_vectored;
            X=LNLC_CS(SequenceName,subrate,y_vectored,reconstructedImage,A,At,par);
            
        end;
    end;
end;
