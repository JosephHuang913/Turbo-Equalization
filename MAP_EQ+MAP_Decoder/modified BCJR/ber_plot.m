close; 
clear;

load ber_MAPEQ.log;
i=1:1:11;
semilogy(ber_MAPEQ(i,1), ber_MAPEQ(i,2), '-.ro');
hold on;
grid on;

load ber.log;
i=1:1:11;
semilogy(ber(i,1), ber(i,2), '-b^');
semilogy(ber(i,1), ber(i,3), '-bs');
semilogy(ber(i,1), ber(i,4), '-b.');
semilogy(ber(i,1), ber(i,12), '-bd');

load ber_BPSK_1_awgn.log;
i=1:1:11;
semilogy(ber_BPSK_1_awgn(i,1), ber_BPSK_1_awgn(i,2), '-.mx');

load ber_CC.log;
i=1:1:11;
semilogy(ber_CC(i,1), ber_CC(i,2), '-.g*');

title('BER Performance of Turbo Equalizer (BCJR-BCJR) in ISI Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 10 1e-6 1]);

legend('MAP EQ (BCJR)', '0_t_h Iteration, MAP EQ + MAP Dec', '1_s_t Iteration', '2_n_d Iteration', '10_t_h Iteration', 'BPSK, AWGN', 'SDVA, AWGN', 'CIR: \{0.407,0.815,0.407\}', 3);
%print -djpeg100 TEQ_MAPEQ_MAP_dec.jpg;
