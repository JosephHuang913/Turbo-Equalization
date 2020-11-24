close; 
clear;

load ber_MAPEQ.log
i=1:1:11;
semilogy(ber_MAPEQ(i,1), ber_MAPEQ(i,2), '-.ro');
hold on;
grid on;

load ber_MAPEQ_MAPdec.log;
i=1:1:11;
semilogy(ber_MAPEQ_MAPdec(i,1), ber_MAPEQ_MAPdec(i,2), '-b^');
semilogy(ber_MAPEQ_MAPdec(i,1), ber_MAPEQ_MAPdec(i,3), '-bs');
semilogy(ber_MAPEQ_MAPdec(i,1), ber_MAPEQ_MAPdec(i,4), '-b.');
semilogy(ber_MAPEQ_MAPdec(i,1), ber_MAPEQ_MAPdec(i,6), '-bd');

load ber_BPSK_1_awgn.log;
i=1:1:12;
semilogy(ber_BPSK_1_awgn(i,1), ber_BPSK_1_awgn(i,2), '-.mx');

load ber_CC.log;
i=1:1:11;
semilogy(ber_CC(i,1), ber_CC(i,2), '-.g*');

title('BER Performance of Turbo Equalizer (BCJR-BCJR) in Multipath Fading Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 10 1e-6 1]);

legend('MAP EQ (BCJR)', '0_t_h Iteration', '1_s_t Iteration', '2_n_d Iteration', '4_t_h Iteration', 'BPSK, AWGN', 'SDVA, AWGN', 'Path Weighting Factor:', '\{0.577,0.577,0.577\}', 'MAP EQ + MAP Dec', 'f_d x t = 0.222222', 3);
print -djpeg100 TEQ_MAPEQ_MAP_dec_fading.jpg;
