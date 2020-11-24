close;
clear;

ber1 = fopen('ber0.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:12;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-b^');
fclose(ber1);

clear;
hold on;
ber1 = fopen('ber1.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:13;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-bo');
fclose(ber1);

clear;
ber1 = fopen('ber2.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:13;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-bs');
fclose(ber1);

clear;
ber1 = fopen('ber10.log', 'r');
ber2 = fscanf(ber1, '%e');
i=1:1:9;
Eb_No(i) = ber2(2*i-1);
err_rate(i) = ber2(2*i);
semilogy(Eb_No, err_rate, '-b*');
fclose(ber1);

%grid;
title('BER Performance of Turbo Equalizer in ISI Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');

legend('0_t_h Iteration', '1_t_h Iteration', '2_t_h Iteration', '10_t_h Iteration', 'CIR: \{0.407,0.815,0.407\}');
%print -djpeg100 Turbo_EQ.jpg;
