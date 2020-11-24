/*************************************************/
/* Author: Chao-wang Huang                       */
/* Date: Wednesday, September 29, 2004           */
/* An (n,k,m) Convolutional code is simulated    */
/* Decoding algorithm: MAP decoding              */
/* Interleaver: S-random interleaver with S=16	 */
/* An MAP Equalizer is implemented               */
/* Channel impulse response: {0.407,0.815,0.407} */
/* TURBO EQUALIZATION									 */
/*************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>
#include <limits.h>
#include <float.h>

void AWGN_noise(float, double, double *);
int Error_count(int, int);
void MAP_EQ_BCJR(double *, double *, double, double *, int, int, int);
void Trellis_diagram(int, int, int, int *);
void CC_Encoder(int *, int *, int, int, int *, int, int);
void MAP_dec_BCJR(double *, double *, double *, int, int, int);

#define Pi 		3.14159265358979
//#define n		2					// (n,k,m) Convolutional Code, # of output symbols
//#define k		1					// (n,k,m) Convolutional Code, # of input symbols
//#define m		2					// (n,k,m) Convolutional Code, # of memories
const int n = 2;
const int k = 1;
const int m = 2;
#define ch_m 	2              // channel memory of ISI channel
#define K		512		  		// packet length of information bit stream
#define N		K*(n/k)	 		// packet length of coded bit stream
#define num_packet 20000		// number of packets simulated
#define Iteration	 10	   		// Iterations of Turbo Equalization ( >= 1 )
#define ch_N_state 4
#define N_state	 4
const int num_state = pow(2,m);		// number of states
const int ch_num_state = pow(2,ch_m);	//number of channel states
const int num_in_sym = pow(2,k);		// number of input symbols
const float CIR[3] = {0.407,0.815,0.407};	// Channel Impulse Response
int gp[2] = {5,7};	// Generator polynomial of CC given in hexadecimal and in reverse order

struct ch_trel_diag			// Data structure of trellis diagram for each branch
		{
      	int from;			// from_state of trellis diagram
         int to;				// to_state of trellis diagram
         int in;				// input data bit of trellis diagram
         float out;			// output codeword symbol of trellis diagram
      };
struct ch_trel_diag ch_Trellis[4][2];			// Trellis[num_state][num_in_sym]

struct trel_diag			// Data structure of trellis diagram for each branch
		{
      	int from;		// from_state of trellis diagram
         int to;			// to_state of trellis diagram
         int in;			// input data bit of trellis diagram
         int out[n];		// output codeword symbol of trellis diagram
      };
struct trel_diag Trellis[4][2];		// Trellis[num_state][num_in_sym]

int main(void)
{
	time_t  t, start, end;
	int i, p, s, input, *data_bit, *coded_bit, *Hk;
   int from_state, to_state, tmp_state, *ch;
   int *S_random, *de_inter, *interleaved, *err_count;
   float ch_out_sym;
   double *bpsk, snr, Eb_No, noise_pwr, noise[2], *Yk, *err_rate;
   double *Dk, *LLR, *intrinsic, *extrinsic;
   FILE *trelis, *ber, *records, *s_inter;

   start = time(NULL);
   printf("BER Performance of TURBO EQUALIZATION in ISI Channel\n");
	printf("Coding Scheme: (2,1,2) Convolutional Code\n");
	printf("Generator polynomials are {5,7} in Octal\n");
   printf("Minimum free distance = 5\n");
   cout << "Decoder: MAP Decoder with modified BCJR algprithm" << endl;
   cout << "Interleaver: S-random interleaver with S = 16" << endl;
   cout << "Channel Impulse Response: {0.407,0.815,0.407}" << endl;
   cout << "Equalizer: MAP Equalizer with modified BCJR Algorithm" << endl;
   printf("Number of bits of simulation = %d\n", K*num_packet);
   printf("Iterations of Turbo Equalization = %d\n\n", Iteration);
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_turbo_EQ.log", "a");
   fprintf(records, "BER Performance of TURBO EQUALIZATION in ISI Channel\n");
   fprintf(records, "Coding Scheme: (2,1,2) Convolutional Code\n");
   fprintf(records, "Generator polynomials are {5,7} in Octal\n");
   fprintf(records, "Minimum free distance = 5\n");
   fprintf(records, "Decoder: MAP Decoder with modified BCJR algprithm\n");
   fprintf(records, "Interleaver: S-random interleaver with S = 16\n");
   fprintf(records, "Channel Impulse Response: {0.407,0.815,0.407}\n");
   fprintf(records, "Equalizer: MAP Equalizer with modified BCJR Algorithm\n");
   fprintf(records, "Number of bits of simulation = %d\n", K*num_packet);
   fprintf(records, "Iterations of Turbo Equalization = %d\n\n", Iteration);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[K];
   coded_bit = new int[N];
   ch = new int[m+1];
   Hk = new int[K];              // Hard Decision
   bpsk = new double[N];
   Yk = new double[N];				// Received signal
   Dk = new double[N];				// Soft Decision
   LLR = new double[N];				// Log-likelihood Ratio
   intrinsic = new double[N];		// Intrinsic information
   extrinsic = new double[N];		// Extrinsic information
   S_random = new int[N];
   de_inter = new int[N];
   interleaved = new int[N];
   err_count = new int[Iteration+1];
   err_rate = new double[Iteration+1];

	srand((unsigned) time(&t));

   Trellis_diagram(n, k, m, &gp[0]);

/****************************************************************/
/* Generate TrellisDiagram for {0.407,0.815,0.407} ISI channel */
/****************************************************************/
	trelis = fopen("ch_Trellis_diagram.log", "w");
   fprintf(trelis, "Trellis diagram of (%.3f,%.3f,%.3f) ISI channel\n", CIR[0], CIR[1], CIR[2]);
   fprintf(trelis, "s(k-1) s(k) input output\n");
   for(from_state=0; from_state<ch_num_state; from_state++) // from_state of trellis diagram
   {
   	for(input=0; input<num_in_sym; input++)		// input symbol (2*input-1) of trellis diagram
      {
      	tmp_state = from_state;
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         ch[0] = 2*(tmp_state & 0x01) - 1;
         ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;
         ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;
         ch_out_sym = ch[0]*CIR[0] + ch[1]*CIR[1] + ch[2]*CIR[2];		// Calculate output symbol
         to_state = tmp_state & (ch_num_state-1); 					// to_state of trellis diagram
         ch_Trellis[from_state][input].from = from_state;
         ch_Trellis[from_state][input].to = to_state;
         ch_Trellis[from_state][input].in = (2*input-1);
         ch_Trellis[from_state][input].out = ch_out_sym;
         fprintf(trelis, "%4d %4d %5d %8.3f\n", from_state, to_state, (2*input-1), ch_out_sym);
      }
   }
   fclose(trelis);

/**************************************************/
/* S-random interleaver and de-interleaver (S=16) */
/**************************************************/
	s_inter = fopen("s_random.txt", "r");
   for(i=0; i<N; i++)		// Interleaver
   	fscanf(s_inter, "%d %d", &de_inter[i], &S_random[i]);
   for(i=0; i<N; i++)		// De-interleaver
   	de_inter[S_random[i]] = i;
   fclose(s_inter);

/************************/
/* main simulation loop */
/************************/
	ber=fopen("ber.log","w");
   for(snr=0; snr<=10; snr++)
   {
      for(s=0; s<=Iteration; s++)
   	  	err_count[s] = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/(((float)k/(float)n)*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);

      from_state = random(4);						// Initial State of Channel
      for(p=0; p<num_packet; p++)
      {
      	for(i=0; i<K-2; i++)
				data_bit[i] = random(2);		// Generate random information bit stream
   		data_bit[K-2] = data_bit[K-1] = 0;

         CC_Encoder(data_bit, coded_bit, n, m, &gp[0], K, num_state);

         // Interleaving
         for(i=0; i<N; i++)
           	interleaved[S_random[i]] = coded_bit[i];

/*****************************************************/
/* BPSK mapping and ISI Channel: {0.407,0.815,0.407} */
/*****************************************************/
   		for(i=0; i<N; i++)
   		{
         	tmp_state = from_state;
            tmp_state = (tmp_state << 1) ^ (interleaved[i] & 0x01);  // read input bit
            ch[0] = 2*(tmp_state & 0x01) - 1;					// input symbol (BPSK)
            ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;        // channel memory (BPSK)
            ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;        // channel memory (BPSK)
            ch_out_sym = ch[0]*CIR[0] + ch[1]*CIR[1] + ch[2]*CIR[2];	// Calculate output symbol
            from_state = tmp_state & (ch_num_state-1); 			// to_state of trellis diagram

            /* AWGN channel */
            AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[i] = ch_out_sym + noise[0];
            // Soft Decision
//            Dk[i] = tanh(Yk[i]/noise_pwr);
            Dk[i] = Yk[i];
         }

/******************************************************************************/
/*                     TURBO EQUALIZATION (BCJR-BCJR)                         */
/******************************************************************************/
			// Initialize a-priori probability
         for(i=0; i<N; i++)
         	intrinsic[i] = 0.0;

         // Iterative Equalization and Decoding
		  	for(s=0; s<=Iteration; s++)
		  	{
         	MAP_EQ_BCJR(Dk, intrinsic, noise_pwr, LLR, N, ch_num_state, num_in_sym);

            // Extrinsic information
				for(i=0; i<N; i++)
            	extrinsic[i] = LLR[i] - intrinsic[i];

            // De-interleaving
				for(i=0; i<N; i++)
               intrinsic[de_inter[i]] = extrinsic[i];

            MAP_dec_BCJR(intrinsic, LLR, extrinsic, K, num_state, num_in_sym);

            // Extrinsic information
            for(i=0; i<N; i++)
            	extrinsic[i] -= intrinsic[i];

            // Interleaving
				for(i=0; i<N; i++)
               intrinsic[S_random[i]] = extrinsic[i];

            for(i=0; i<K; i++)	// Bit error count
            {
            	if(LLR[i]>0.0)
               	Hk[i] = 1;
               else
               	Hk[i] = 0;

             	err_count[s] += Error_count(data_bit[i], Hk[i]);
            }
			}
/******************************************************************************/
		}

      // Statistics and records
      cout << "Error Rate = ";
      for(s=0; s<=Iteration; s++)
      {
      	err_rate[s] = err_count[s] / (double)(K*num_packet);
      	printf("%e, ", err_rate[s]);
      }
      cout << endl;

      fprintf(ber, "%f ", Eb_No);
      fprintf(records, "%f ", Eb_No);
		for(s=0; s<=Iteration; s++)
      {
	      fprintf(ber, "%e ", err_rate[s]);
         fprintf(records, "%e ", err_rate[s]);
      }
      fprintf(ber, "\n");
      fprintf(records, "\n");
      fflush(records);
      fflush(ber);
   }

   delete data_bit;
   delete coded_bit;
   delete ch;
   delete Hk;
   delete bpsk;
   delete Yk;
   delete Dk;
   delete LLR;
   delete S_random;
   delete de_inter;
   delete interleaved;
   delete intrinsic;
   delete extrinsic;
   delete err_count;
   delete err_rate;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
	fclose(ber);
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

void MAP_EQ_BCJR(double *Data_in, double *intrinsic, double noise_pwr,
					  double *LLR, int Packet_length, int Num_state, int Num_in)
{
/******************************************************************/
/* MAP Equalizer (modified BCJR Algorithm) for static ISI channel */
/******************************************************************/
   const float pi = 3.14159265358979;
	int i, j, l;
   double p1, *normal, **alpha, **beta, ***gamma, *delta, min, **a_priori;

   a_priori = new double*[Packet_length];			// a-priori probability, a-priori[time index][input]
   for(i=0; i<Packet_length; i++)
   	a_priori[i] = new double[Num_in];
   alpha = new double*[Packet_length+1];			// alpha[time index][state]
   for(i=0; i<=Packet_length; i++)
      alpha[i] = new double[Num_state];
   beta = new double*[Packet_length+1];       	// beta[time index][state]
   for(i=0; i<=Packet_length; i++)
      beta[i] = new double[Num_state];
   gamma = new double**[Packet_length];			// gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)
      gamma[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         gamma[i][j] = new double[Num_in];
   normal = new double[Packet_length+1];			// renormalization of BCJR
   delta = new double[Num_in];

	// Initialization of alpha and beta
   for(i=1; i<=Packet_length; i++)					// alpha[time index][state]
   	for(j=0; j<Num_state; j++)
      	alpha[i][j] = 0.0;
   for(j=0; j<Num_state; j++)
   	alpha[0][j] = 1/(float)Num_state;

   for(i=0; i<Packet_length; i++)           		// beta[time index][state]
   	for(j=0; j<Num_state; j++)
      	beta[i][j] = 0.0;
   for(j=0; j<Num_state; j++)
   	beta[Packet_length][j] = 1.0;

   // Calculate a-priori probability from intrinsic information
   for(i=0; i<Packet_length; i++)					// time index
   	for(l=0; l<Num_in; l++)							// input symbol
      	a_priori[i][l] = exp(l*intrinsic[i]) / (1 + exp(intrinsic[i]));

   // calculate gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)					// time index
   	for(j=0; j<Num_state; j++)						// state index
      	for(l=0; l<Num_in; l++)						// input symbol
         {
         	p1 = exp(-pow(Data_in[i]-ch_Trellis[j][l].out,2)/(2*noise_pwr))/sqrt(2*pi*noise_pwr);
            gamma[i][j][l] = a_priori[i][l] * p1;		// gamma[time index][state][input]
         }

   // calculate alpha[time index][state]
   for(i=1; i<=Packet_length; i++)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	alpha[i][ch_Trellis[j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

      normal[i] = 0.0;									// for renormalization
      for(j=0; j<Num_state; j++)						// to_state index
      	normal[i] += alpha[i][j];

      for(j=0; j<Num_state; j++)
      	alpha[i][j] = alpha[i][j] / normal[i];
   }
   normal[0] = 1.0;

   // calculate beta[time index][state]
   for(i=Packet_length-1; i>0; i--)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	beta[i][j] += beta[i+1][ch_Trellis[j][l].to] * gamma[i][j][l];

      for(j=0; j<Num_state; j++)
      	beta[i][j] = beta[i][j] / normal[i];
   }

   // calculate conditional LLR
   for(i=0; i<Packet_length; i++)					// time index
   {
   	min = 0.0;		// find the minimum product of alpha*gamma*beta
      for(j=0; j<Num_state; j++)
      	for(l=0; l<Num_in; l++)
         {
         	delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][ch_Trellis[j][l].to];

            if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
            	min = delta[0];
         }

      if(min == 0.0 || min > 1.0)	// if all else fails, make min real small
      	min = 1E-100;

      delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)	// input bit
         	delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][ch_Trellis[j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      LLR[i] = log(delta[1]/delta[0]);
   }

   for(i=0; i<Packet_length; i++)
   	delete a_priori[i];
   delete a_priori;
   for(i=0; i<=Packet_length; i++)
      delete alpha[i];
	delete alpha;
   for(i=0; i<=Packet_length; i++)
       delete beta[i];
   delete beta;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete gamma[i][j];
   for(i=0; i<Packet_length; i++)
      delete gamma[i];
   delete gamma;
   delete normal;
   delete delta;
}
/******************************************************************************/

void Trellis_diagram(int n, int k, int m, int *gp)
{
/**************************************************************/
/* Generate TrellisDiagram for (n=2, k=1, m=2) {5,7} CC code */
/**************************************************************/
	int i, j, input, out_bit, out_sym, from_state, to_state, tmp_state;
   int num_state = pow(2,m), num_in_sym = pow(2,k);

   for(from_state=0; from_state<num_state; from_state++) // from_state of trellis diagram
   {
   	for(input=0; input<num_in_sym; input++)		// input of trellis diagram for (n, k, m)
      {
      	tmp_state = from_state;
         out_sym = 0;                      // output codeword symbol of trellis diagram
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         for(i=0; i<n; i++)
         {
         	out_bit = 0;						// output bit of trellis diagram
            for(j=m; j>=0; j--)
            	out_bit ^= ((tmp_state & gp[i]) >> j) & 1;  	// Calculate output bit

            out_sym = (out_sym << 1) ^ out_bit;				  	// Calculate output symbol
         }
         to_state = tmp_state & (num_state-1); 					// to_state of trellis diagram

         Trellis[from_state][input].from = from_state;
         Trellis[from_state][input].to = to_state;
         Trellis[from_state][input].in = input;
         Trellis[from_state][input].out[0] = ((out_sym>>1)&1);
         Trellis[from_state][input].out[1] = out_sym&1;
      }
   }
}

void CC_Encoder(int *Data_in, int *Code_bit, int n, int m, int *gp, int Packet_length, int Num_state)
{
/**********************************************************************/
/* Convolutional Encoder (n=2, k=1, m=2) Generator polynomial: {5, 7} */
/**********************************************************************/
   int i, j, l, from_state, tmp_state, out_bit;

	from_state = 0;
   for(i=0; i<Packet_length; i++)
   {
   	tmp_state = from_state;
//    out_sym = 0;                      // output codeword symbol of trellis diagram
		tmp_state = (tmp_state << 1) ^ (Data_in[i] & 0x01);  // read input bit
      for(j=0; j<n; j++)
      {
      	out_bit = 0;						// output bit of trellis diagram
         for(l=m; l>=0; l--)
         	out_bit ^= ((tmp_state & gp[j]) >> l) & 1;  	// Calculate output bit

         Code_bit[2*i+j] = out_bit;
//       out_sym = (out_sym << 1) ^ out_bit;				  	// Calculate output symbol
		}
      from_state = tmp_state & (Num_state-1); 					// to_state of trellis diagram
	}
}
/******************************************************************************/

void MAP_dec_BCJR(double *intrinsic, double *LLR, double *extrinsic, int Packet_length, int Num_state, int Num_in)
{
/**********************************************************************/
/* MAP decoder (modified BCJR Algorithm)                              */
/* Convolutional Decoder (n=2, k=1, m=2) Generator polynomial: {5, 7} */
/**********************************************************************/
	int i, j, l;
   double a_priori=0.5, *normal, **alpha, **beta, ***gamma, *delta, min, **code_bit_prob;
   	/* Information Bit */

   code_bit_prob = new double*[2*Packet_length];// a-priori probability, a-priori[time index][input]
   for(i=0; i<2*Packet_length; i++)
   	code_bit_prob[i] = new double[Num_in];
   alpha = new double*[Packet_length+1];			// alpha[time index][state]
   for(i=0; i<=Packet_length; i++)
      alpha[i] = new double[Num_state];
   beta = new double*[Packet_length+1];       	// beta[time index][state]
   for(i=0; i<=Packet_length; i++)
      beta[i] = new double[Num_state];
   gamma = new double**[Packet_length];			// gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)
      gamma[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         gamma[i][j] = new double[Num_in];
   normal = new double[Packet_length+1];			// renormalization of BCJR
   delta = new double[Num_in];

	// Initialization of alpha and beta
   for(i=0; i<=Packet_length; i++)					// alpha[time index][state]
   	for(j=0; j<Num_state; j++)
      	alpha[i][j] = 0.0;
   alpha[0][0] = 1.0;

   for(i=0; i<=Packet_length; i++)           	// beta[time index][state]
   	for(j=0; j<Num_state; j++)
      	beta[i][j] = 0.0;
   beta[Packet_length][0] = 1.0;

   // Calculate a-priori probability of the code bit
   for(i=0; i<2*Packet_length; i++)					// time index
   	for(l=0; l<Num_in; l++)							// code bit
      	code_bit_prob[i][l] = exp(l*intrinsic[i]) / (1 + exp(intrinsic[i]));

   // calculate gamma[time index][state][input]
   for(i=0; i<Packet_length; i++)					// time index
   	for(j=0; j<Num_state; j++)						// state index
      	for(l=0; l<Num_in; l++)						// input symbol
         	gamma[i][j][l] = a_priori * code_bit_prob[2*i][Trellis[j][l].out[0]]
                     				     * code_bit_prob[2*i+1][Trellis[j][l].out[1]];

   // calculate alpha[time index][state]
   for(i=1; i<=Packet_length; i++)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	alpha[i][Trellis[j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

      normal[i] = 0.0;									// for renormalization
      for(j=0; j<Num_state; j++)						// to_state index
      	normal[i] += alpha[i][j];

      for(j=0; j<Num_state; j++)
      	alpha[i][j] = alpha[i][j] / normal[i];
   }
   normal[0] = 1.0;

   // calculate beta[time index][state]
   for(i=Packet_length-1; i>0; i--)					// time index
   {
   	for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	beta[i][j] += beta[i+1][Trellis[j][l].to] * gamma[i][j][l];

      for(j=0; j<Num_state; j++)
      	beta[i][j] = beta[i][j] / normal[i];
   }

   // Calculate conditional LLR of the information bit
   for(i=0; i<Packet_length; i++)					// time index
   {
   	min = 0.0;											// find the minimum product of alpha*gamma*beta
      for(j=0; j<Num_state; j++)
      	for(l=0; l<Num_in; l++)
         {
         	delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

            if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
            	min = delta[0];
         }

      if(min == 0.0 || min > 1.0)					// if all else fails, make min real small
      	min = 1E-100;

      delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)						// input bit
         	delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      LLR[i] = log(delta[1]/delta[0]);

	   // Calculate Extrinsic Information of the 1st code bit
		delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)                // input bit
         	delta[Trellis[j][l].out[0]] += alpha[i][j] * gamma[i][j][l]
                     										 	 * beta[i+1][Trellis[j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      extrinsic[2*i] = log(delta[1]/delta[0]);	// calculate extrinsic information

      // Calculate Extrinsic Information of the 2nd code bit
      delta[0] = delta[1] = 0.0;
      for(j=0; j<Num_state; j++)						// from_state index
      	for(l=0; l<Num_in; l++)                // input bit
         	delta[Trellis[j][l].out[1]] += alpha[i][j] * gamma[i][j][l]
                     										 	 * beta[i+1][Trellis[j][l].to];

      if(delta[1] == 0.0)
      	delta[1] = min;
      if(delta[0] == 0.0)
      	delta[0] = min;

      extrinsic[2*i+1] = log(delta[1]/delta[0]);// calculate extrinsic information
   }

	for(i=0; i<2*Packet_length; i++)
   	delete code_bit_prob[i];
   delete code_bit_prob;
   for(i=0; i<=Packet_length; i++)
      delete alpha[i];
	delete alpha;
   for(i=0; i<=Packet_length; i++)
       delete beta[i];
   delete beta;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete gamma[i][j];
   for(i=0; i<Packet_length; i++)
      delete gamma[i];
   delete gamma;
   delete normal;
   delete delta;
}
/******************************************************************************/

void AWGN_noise(float mu, double variance, double *noise)
{
//	const  float Pi = 3.14159265358979;
   double u1, u2;
   do
   {
   	u1 = (double)rand()/(double)RAND_MAX;
      u2 = (double)rand()/(double)RAND_MAX;
   }
   while(u1 == 0.0 || u2 == 0.0);

   *(noise+0) = (sqrt(-2.0*log(u1))*cos(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
   *(noise+1) = (sqrt(-2.0*log(u1))*sin(2*Pi*u2))*sqrt(variance/2.0)+mu/sqrt(2.0);
}

int Error_count(int x, int y)
{
	if(x == y)
   	return 0;
   else
   	return 1;
}

