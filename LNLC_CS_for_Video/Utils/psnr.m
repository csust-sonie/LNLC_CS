function y = psnr(im1,im2)

x1 = double(im1(:));
x2 = double(im2(:));

if x1==x2
    error('error of psnr');
end
err = x1 - x2;

y = 20*log10(255/(sqrt(mean(mean(err.^2)))));