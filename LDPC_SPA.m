function v = LDPC_SPA(H, LLR_y, iterMax)
% LDPC_SPA - LDPC decoding using the Sum-Product Algorithm (SPA)
% Author: Zhou
%
% Inputs:
%   H       - Parity check matrix
%   LLR_y   - Received log-likelihood ratio vector
%   iterMax - Maximum number of decoding iterations
%
% Output:
%   v       - Decoded binary vector (0/1)

U0i = LLR_y;
Uji = zeros(size(H));
Vij = zeros(size(H'));
[VerificationNodes, VariableNodes] = size(H);
x = zeros(size(LLR_y));

for iter = 1:iterMax
    % Variable node to check node message
    for i = 1:VariableNodes
        idx = find(H(:, i) == 1);
        for k = 1:length(idx)
            Vij(i, idx(k)) = U0i(i) + sum(Uji(idx, i)) - Uji(idx(k), i);
        end
    end

    % Check node to variable node message
    for j = 1:VerificationNodes
        idx = find(H(j, :) == 1);
        for k = 1:length(idx)
            multipleVal = 2 * atanh(prod(tanh(Vij(idx, j) / 2)) / tanh(Vij(idx(k), j) / 2));
            if multipleVal == inf
                Uji(j, idx(k)) = 10;
            elseif multipleVal == -inf
                Uji(j, idx(k)) = -10;
            else
                Uji(j, idx(k)) = multipleVal;
            end
        end
    end

    % Final decision (hard decision from total LLR)
    for i = 1:length(x)
        idx = find(H(:, i) == 1);
        addVal = sum(Uji(idx, i)) + U0i(i);
        x(i) = addVal < 0;
    end

    % Early stopping if syndrome is all zero
    if mod(H * x', 2) == 0
        break;
    end
end

v = x;
end
