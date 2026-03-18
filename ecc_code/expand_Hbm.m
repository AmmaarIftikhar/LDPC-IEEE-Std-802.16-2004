function [H, Hp, Hs] = expand_Hbm(Hbm)
% expand_Hbm - Expand WiMAX base matrix Hbm into full LDPC parity-check matrix
% Author: Zhou
%
% Input:
%   Hbm - Base matrix (12¡Á24) with shift values or -1 for zero blocks
%
% Output:
%   H   - Full parity-check matrix (1152¡Á2304) after expansion (z = 96)
%   Hp  - Left square submatrix (for check bits)
%   Hs  - Right submatrix (for message bits)

z = 96;
[m, n] = size(Hbm);
H = zeros(m*z, n*z);

for i = 1:m
    for j = 1:n
        val = Hbm(i,j);
        row_range = (i-1)*z + (1:z);
        col_range = (j-1)*z + (1:z);

        if val == -1
            block = zeros(z);
        elseif val == 0
            block = eye(z);
        else
            block = circshift(eye(z), [0, mod(val, z)]);
        end

        H(row_range, col_range) = block;
    end
end

Hp = H(1:m*z, 1:m*z);
Hs = H(1:m*z, m*z+1:n*z);
end
