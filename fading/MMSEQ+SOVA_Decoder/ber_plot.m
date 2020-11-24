close;
clear;

load ber_MMSE_SOVA1.log;

i=1:1:7;
semilogy(ber_MMSE_SOVA1(i,1)-3, ber_MMSE_SOVA1(i,2), '-bs');
hold on;
grid on;
semilogy(ber_MMSE_SOVA1(i,1), ber_MMSE_SOVA1(i,3), '-.rs');

load ber_MMSE_SOVA2.log;

i=1:1:7;
semilogy(ber_MMSE_SOVA2(i,1)-3, ber_MMSE_SOVA2(i,2), '-bo');
hold on;
grid on;
semilogy(ber_MMSE_SOVA2(i,1), ber_MMSE_SOVA2(i,3), '-.ro');

load ber_MMSE_SOVA3.log;

i=1:1:7;
semilogy(ber_MMSE_SOVA3(i,1)-3, ber_MMSE_SOVA3(i,2), '-b^');
hold on;
grid on;
semilogy(ber_MMSE_SOVA3(i,1), ber_MMSE_SOVA3(i,3), '-.r^');

load ber_MMSE_SOVA4.log;

i=1:1:7;
semilogy(ber_MMSE_SOVA4(i,1)-3, ber_MMSE_SOVA4(i,2), '-bv');
hold on;
grid on;
semilogy(ber_MMSE_SOVA4(i,1), ber_MMSE_SOVA4(i,3), '-.rv');

gamma=0:1:22;
gamma_c=10.^(gamma./10)./3;
mju=sqrt(gamma_c./(1+gamma_c));
Pb=((1/2).*(1-mju)).^3.*(1+3.*((1/2).*(1+mju)).^1+6.*((1/2).*(1+mju)).^2);
semilogy(gamma, Pb, '-k.');

r=0:1:12;
Pb=0.5.*erfc(sqrt(10.^(r./10)));
semilogy(r, Pb, '-k');

title('BER Performance of MMSE EQ (RLS) + SOVA Decoder in Multipath Fading Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Probability of Bit Error');
axis([0 30 1e-6 1]);

legend('\{0, 1, 0\}', '\{0, 1, 0\} Coding', '\{0.2096, 0.9551, 0.2096\}', '\{0.2096, 0.9551, 0.2096\} Coding', '\{0.407, 0.815, 0.407\}', '\{0.407, 0.815, 0.407\} Coding', '\{0.577, 0.577, 0.577\}, Fading', '\{0.577, 0.577, 0.577\} Coding', 'Diversity Order: 3', 'BPSK on AWGN', '30 km/h, f_d * t = 0.000056', 1);
%print -djpeg100 MMSEQ+SOVA.jpg;
