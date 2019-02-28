function [Y]=readyuv2(yuvfile,begin,frames,height,width,format)
  fid = fopen(yuvfile,'r');%'akiyo_qcif.yuv'
 % row=height;col=width;    %·´¹ýÀ´
 row=width;col=height;
  Y=cell(1,frames);

 for frame=1:begin-1
  fread(fid,[row,col],'uchar'); 
  switch format
      case '422'
           fread(fid,[row,col],'uchar');
           %U{frame}=fread(fid,[row/2,col/2],'uchar'); 
           %V{frame}=fread(fid,[row/2,col/2],'uchar'); 
      case '420'           
           fread(fid,[row/2,col/2],'uchar'); 
           fread(fid,[row/2,col/2],'uchar');  
     case '444'
           fread(fid,[row,col],'uchar'); 
           fread(fid,[row,col],'uchar'); 
      otherwise
          error('myApp:argChk','invalid yuv format');
  end; 
 end

 for frame=1:frames 
  im= fread(fid,[row,col],'uchar'); 
  Y{frame}=fliplr(imrotate(im,270));
  switch format
      case '422'
           fread(fid,[row,col],'uchar');
           %U{frame}=fread(fid,[row/2,col/2],'uchar'); 
           %V{frame}=fread(fid,[row/2,col/2],'uchar'); 
      case '420'           
           fread(fid,[row/2,col/2],'uchar'); 
           fread(fid,[row/2,col/2],'uchar');  
     case '444'
           fread(fid,[row,col],'uchar'); 
           fread(fid,[row,col],'uchar'); 
      otherwise
          error('myApp:argChk','invalid yuv format');
  end; 
end
fclose(fid);