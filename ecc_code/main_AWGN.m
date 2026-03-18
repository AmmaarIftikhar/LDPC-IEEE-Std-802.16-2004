% LDPC Simulation over AWGN Channel with SPA Decoder
% Author: Zhou & Ammaar
clc; clear; close all;

N = 2304;
K = 1152;
R = K / N;
Hbm = Hbm();
[H, Hp, Hs] = expand_Hbm(Hbm);
[P, rearranged_cols] = H2P(H);

EbN0_dB = [-1:0.1:2];
EbN0 = 10 .^ (EbN0_dB / 10);
snrAll = EbN0 * R;
Es = 1;
iterMax = 30;

BER = zeros(1, length(EbN0));
FER = zeros(1, length(EbN0));

target_error_bits = 2000;
max_test_frame = 10000;
consecutive_zero_required = 50;

for ii = 1:length(EbN0_dB)
    ErrorBits_SP = 0;
    ErrorFrame_SP = 0;
    frame_SP = 0;
    consecutive_zero_count = 0;
    snr = snrAll(ii);
    sigma = sqrt(Es / snr / 2);

    while ErrorBits_SP < target_error_bits && frame_SP < max_test_frame
        frame_SP = frame_SP + 1;
        s = randi([0, 1], 1, K);
        c1 = mod(P * s', 2);
        u1 = [c1' s];
        x = reorder_bits(u1, rearranged_cols);
        d = 1 - 2 .* x;
        y = d + sigma * randn(size(d));
        LLR_y = 2 * y / (sigma^2);
        v = LDPC_SPA(H, LLR_y, iterMax);
        v_SP = extract_mesg(v, rearranged_cols);
        errorbits_SP = sum(s ~= v_SP);
        ErrorBits_SP = ErrorBits_SP + errorbits_SP;

        if errorbits_SP ~= 0
            ErrorFrame_SP = ErrorFrame_SP + 1;
            consecutive_zero_count = 0;
        else
            consecutive_zero_count = consecutive_zero_count + 1;
            if consecutive_zero_count >= consecutive_zero_required
                ErrorBits_SP = 0;
                break;
            end
        end
    end

    if ErrorBits_SP == 0
        BER(1, ii) = 0;
    else
        BER(1, ii) = ErrorBits_SP / (K * frame_SP);
    end
    FER(1, ii) = ErrorFrame_SP / frame_SP;

    fprintf('SNR = %.2f dB, Frames = %d, Bit Errors = %d, BER = %.3e\n', ...
        EbN0_dB(ii), frame_SP, ErrorBits_SP, BER(1, ii));
end

xlswrite('./BERofAWGN.xlsx', BER);
xlswrite('./FERofAWGN.xlsx', FER);

BER_plot = BER;
BER_plot(BER_plot == 0) = 1e-8;

figure('numbertitle','off','name','BER of LDPC_SPA decoding AWGN');
semilogy(EbN0_dB, BER_plot, 'K-^', 'LineWidth', 1.0, 'MarkerSize', 6); hold on;
BER_BPSK = qfunc(sqrt(2 * EbN0));
semilogy(EbN0_dB, BER_BPSK, 'r-*', 'LineWidth', 1.0, 'MarkerSize', 6);
xlabel('Eb/N0(dB)');
ylabel('BER');
legend('LDPC-AWGN-SPA', 'uncoded-AWGN-BPSK');
grid on;




