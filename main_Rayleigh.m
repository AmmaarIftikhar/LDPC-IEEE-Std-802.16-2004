% LDPC over Rayleigh channel simulation using SPA decoder
% Author: Zhou & Ammaar

clc; clear; close all;

N = 2304;
K = 1152;
R = K / N;

Hbm = Hbm();
[H, Hp, Hs] = expand_Hbm(Hbm);
[P, rearranged_cols] = H2P(H);

EbN0_dB = 0:3:60;
EbN0 = 10 .^ (EbN0_dB / 10);
snrAll = EbN0 * R;
Es = 1;
iterMax = 30;

BER_ray = zeros(1, length(EbN0));
FER_ray = zeros(1, length(EbN0));

for ii = 1:length(EbN0_dB)
    frame = 500 * (EbN0_dB(ii) <= 6) + 1000 * (EbN0_dB(ii) > 6);
    ErrorBits_SP_ray = 0;
    ErrorFrame_SP_ray = 0;
    frame_SP_ray = 0;
    consecutive_no_error = 0; 

    snr = snrAll(ii);
    sigma = sqrt(Es / snr / 2);

    for jj = 1:frame
        fprintf('current snr = %ddB and the frame is %d\n', EbN0_dB(ii), jj);
        s = randi([0, 1], 1, K);
        c1 = mod(P * s', 2);
        u1 = [c1' s];
        x = reorder_bits(u1, rearranged_cols);
        d = 1 - 2 .* x;

        h = (1 / sqrt(2)) * (randn(1, 1) + sqrt(-1) * randn(1, 1));
        n = sigma * (1 / sqrt(2)) * (randn(size(d)) + sqrt(-1) * randn(size(d)));
        y_ray = h .* d + n;

        y_eq = real(conj(h) .* y_ray) ./ (abs(h).^2);
        LLR_ray = 2 * y_eq / (sigma^2);

        v_ray = LDPC_SPA(H, LLR_ray, iterMax);
        v_SP_ray = extract_mesg(v_ray, rearranged_cols);
        errorbits_SP_ray = sum(s ~= v_SP_ray);
        ErrorBits_SP_ray = ErrorBits_SP_ray + errorbits_SP_ray;
        frame_SP_ray = frame_SP_ray + 1;

        if errorbits_SP_ray ~= 0
            ErrorFrame_SP_ray = ErrorFrame_SP_ray + 1;
            consecutive_no_error = 0; 
        else
            consecutive_no_error = consecutive_no_error + 1;
        end

        if consecutive_no_error >= 100
            fprintf('Early stopping at Eb/N0 = %ddB due to 100 consecutive error-free frames.\n', EbN0_dB(ii));
            break;
        end
    end

    BER_ray(ii) = ErrorBits_SP_ray / (K * frame_SP_ray);
    FER_ray(ii) = ErrorFrame_SP_ray / frame_SP_ray;

    if BER_ray(ii) == 0
        BER_ray(ii) = 1e-7;
    end
end

xlswrite('./BERofRayleigh.xlsx', BER_ray);
xlswrite('./FERofRayleigh.xlsx', FER_ray);

figure('numbertitle','off','name','BER of LDPC\_SPA decoding Rayleigh channel');
semilogy(EbN0_dB, BER_ray, 'K-o', 'LineWidth', 1.0, 'MarkerSize', 6); hold on;
BER_BPSK_rayleigh = 0.5 * (1 - sqrt(EbN0 ./ (1 + EbN0)));
semilogy(EbN0_dB, BER_BPSK_rayleigh, 'r-*', 'LineWidth', 1.0, 'MarkerSize', 6);
xlabel('Eb/N0(dB)');
ylabel('BER');
legend('LDPC-Rayleigh-SPA', 'uncoded-Rayleigh-BPSK');
grid on;



