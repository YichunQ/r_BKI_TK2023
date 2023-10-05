clc
clear
close

%% load data
tdata = double(imread('facade.bmp')); 
SNR = 5;
sizes = size(tdata);
dim = numel(sizes);
noisy = 1;
    if noisy
        data = addnoise(tdata,sizes,SNR);
    else 
        data = tdata;
    end
    
%%  denoising
rank = [30,30,3];
e = 0.2;  % sketching over-rate
sksize = min(rank+1/e,sizes);
[X_k,A] = rBKI_TK(data, rank, sksize);

%%  print results

subplot(1,3,1)
imshow(uint8(tdata))
title('clean data')
subplot(1,3,2)
imshow(uint8(data)) 
title('noisy data')
subplot(1,3,3)
imshow(uint8(double(X_k)))
title('denoising data')
