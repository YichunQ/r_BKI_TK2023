function [X_k,A] = rBKI_TK(data, rank, sksize)
% A script for randomized Block Krylov Iteration based Tucker decomposition
% data: input data, rank: approximation rank, sksize: sketching size
% X_k: rank_k approximation, A: factor matrices
% If you found this code helpful, please cite the work below:
% title={Towards Efficient and Accurate Approximation: Tensor Decomposition Based on Randomized Block Krylov Iteration}, 
% author={Yichun Qiu and Weijun Sun and Guoxu Zhou and Qibin Zhao}
% Please feel free to contact me at ycqiu@gdut.edu.cn

maxiter = 20;
N = size(data);
dim = numel(N);
A = cell(1,dim);
AAt = cell(1,dim);
mode = 1:dim;

 for i = 1:maxiter    
    for d = 1:dim
        X = double(tenmat(data,d));
        [m,n] = size(X);
        k = rank(d); % lra-parameter
        l = sksize(d); %  sketch-size
        s = min(min(2*k,l),N(d));
        it = 3;    %  num of rBKIterations
        temp = X*X';
        C = X * randn(n, s);
        Krylov = zeros(m, s * it);
        Krylov(:, 1:s) = C;  
            for ii = 2: it
                C = temp * C;   
                [C, ~] = qr(C, 0); 
                Krylov(:, (ii-1)*s+1: ii*s) = C;  
            end
        [Q, ~] = qr(Krylov, 0);
        B = Q'*temp*Q;
        [U, ~, ~] = svd(B,'econ');
        U = Q*U;
        A{d} = U(:,1:k); 
        AAt{d} = A{d}*A{d}'; 
    end                        
    X_k = ttm(tensor(data),AAt,mode);
end

