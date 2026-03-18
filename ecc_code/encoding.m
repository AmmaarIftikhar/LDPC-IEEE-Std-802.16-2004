function [x] = encoding(Hs, Hp, s)
% encoding - Encode message s using LDPC generator matrix derived from [Hp | Hs]
% Author: Zhou
%
% Inputs:
%   Hs - Right submatrix of parity-check matrix (used for message bits)
%   Hp - Left submatrix of parity-check matrix (used for check bits)
%   s  - Message vector (1 ¡Á kb*z)
%
% Output:
%   x  - Encoded codeword (1 ¡Á (mb+kb)*z)

mb = 12; kb = 12;
z = 96;

w = zeros(1, mb*z);  % Intermediate result
p = zeros(1, mb*z);  % Parity bits
x = zeros(1, (mb + kb)*z);  % Output codeword

w = s * (Hs');  % Step 1: calculate intermediate vector w

% Step 2: compute parity bits p using Algorithm 2
p(1) = mod(w(1), 2);
idx = 1;
for i = 1:mb*z
    if idx > (mb-1)*z && idx <= mb*z-1
        p(idx - (mb-1)*z + 1) = mod(w(idx - (mb-1)*z + 1) + p(idx), 2);
        idx = idx - (mb-1)*z + 1;
    elseif idx > 0 && idx <= (mb-1)*z
        p(z + idx) = mod(w(z + idx) + p(idx), 2);
        idx = idx + z;
    end
end
p(mb*z) = mod(w(i) + p((mb-1)*z), 2);

% Alternative Algorithm 1 (commented out)
% pb = mod(w * inv(Hp'), 2);

x = [p s];  % Final codeword [parity | message]
end
