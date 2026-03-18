function u = reorder_bits(c, rearranged_cols)
% reorder_bits - Reverse column swaps from Gaussian elimination (Author: Zhou)
% Inputs: 
%   c - codeword vector
%   rearranged_cols - column swap history from H2P
% Output: 
%   u - reordered codeword after reversing swaps

    rows = length(rearranged_cols);
    for i = rows:-1:1
        if rearranged_cols(i) ~= 0
            temp = c(i);
            c(i) = c(rearranged_cols(i));
            c(rearranged_cols(i)) = temp;
        end
    end
    u = c;
end
