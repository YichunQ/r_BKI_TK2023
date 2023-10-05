function Y = addnoise(X,I,SNR)
%  X is clean data, N is noise, Y is noisy data, I is data-size
%  SNR = 20*log10(norm(tensor(X))/norm(tensor(N)));

norm_X = norm(tensor(X));
noise = randn(I);
sigma = norm_X/(10^(SNR/20)*norm(tensor(noise)));  % noisy
N = sigma.*noise;
Y = X + N;


