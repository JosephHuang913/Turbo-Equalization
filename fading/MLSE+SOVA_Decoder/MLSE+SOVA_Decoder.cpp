/******************************************************/
/* Author: Chao-wang Huang                            */
/* Date: Friday, March 04, 2005                       */
/* An (n,k,m) Convolutional code is simulated         */
/* Decoding algorithm: SOVA decoding                  */
/* Interleaver: S-random interleaver with S=16	      */
/* An MLSE (SOVA) Equalizer is implemented            */
/* 3-path Rayleigh fading channel with equal strength */
/* Multipath weighting factor: {0.577, 0.577, 0.577}  */
/* TURBO EQUALIZATION	         			  	      	*/
/******************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include <iostream.h>
#include <limits.h>
#include <float.h>

void JakesFading(double, float, double, int, double *);
void AWGN_noise(float, double, double *);
int Error_count(int, int);
void MLSE_SOVA(double **, double *, double, int *, double *, int, int, int);
void CH_Trellis_diagram(int, int, float *, double *, double *, double *);
void Trellis_diagram(int, int, int, int *);
void CC_Encoder(int *, int *, int, int, int *, int, int);
void Viterbi_dec_SOVA(double *, int *, double *, int, int, int);

//#define Pi 		3.14159265358979
//#define n		2					// (n,k,m) Convolutional Code, # of output symbols
//#define k		1					// (n,k,m) Convolutional Code, # of input symbols
//#define m		2					// (n,k,m) Convolutional Code, # of memories
const int n = 2;
const int k = 1;
const int m = 2;
#define ch_m 	2              // channel memory of ISI channel
//#define K		512		  		// packet length of information bit stream
//#define N		K*(n/k)	 		// packet length of coded bit stream
const int K = 512;
const int N = 1024;
#define num_packet 20000			// number of packets simulated
#define Iteration	 4	   		// Iterations of Turbo Equalization ( >= 1 )
#define ch_N_state 4
#define N_state	 4
const int num_state = pow(2,m);		// number of states
const int ch_num_state = pow(2,ch_m);	//number of channel states
const int num_in_sym = pow(2,k);		// number of input symbols
int gp[2] = {5,7};	// Generator polynomial of CC given in hexadecimal and in reverse order
const int Num_path = 3;
float CIR[3] = {0.577, 0.577, 0.577};			// Channel Weighting Factor
//float CIR[3] = {0.407, 0.815, 0.407};
const float vc = 120.0; 								/* speed of vehicle in km/hr */
const double C = 3.0E8;  							/* light speed */
const double fc = 2.0e9; 							/* carrier frequency */
const double sym_rate = 1E3; 						/* Symbol rate in symbol/sec */
const double Doppler = (vc*1000.0/3600.0)/(C/fc);  // Maximum Doppler frequency (Hz)

struct ch_trel_diag			// Data structure of trellis diagram for each branch
		{
      	int from;			// from_state of trellis diagram
         int to;				// to_state of trellis diagram
         int in;				// input data symbol of trellis diagram
         float out[2];			// output symbol of trellis diagram
      };
struct ch_trel_diag ch_Trellis[1024][4][2];	// Trellis[time index][num_state][num_in_sym]

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
	int i, p, s, *data_bit, *coded_bit, *Hk;
   int from_state, tmp_state, *ch;
   int *S_random, *de_inter, *interleaved, *err_count;
   double ch_out_sym_I, ch_out_sym_Q, t1, t2, t3, fade1[2], fade2[2], fade3[2];
   double *bpsk, snr, Eb_No, noise_pwr, noise[2], **Yk, *err_rate;
   double *Ak, *LLR, *z, r, *intrinsic, *extrinsic;
   FILE *ber, *records, *s_inter;

   start = time(NULL);
   printf("BER Performance of TURBO EQUALIZATION in ISI Channel\n");
	printf("Coding Scheme: (2,1,2) Convolutional Code\n");
	printf("Generator polynomials are {5,7} in Octal\n");
   printf("Minimum free distance = 5\n");
   cout << "Decoder: Viterbi Decoder with SOVA algprithm" << endl;
   cout << "Interleaver: S-random interleaver with S = 16" << endl;
   cout << "3-path Rayleigh fading channel with equal strength" << endl;
	cout << "Multipath weighting factor: {0.577, 0.577, 0.577}" << endl;
   cout << "Equalizer: MLSE Equalizer with SOVA Algorithm" << endl;
   printf("Speed of the vehicle = %f (km/h)\n", vc);
   printf("Carrier Frequency = %e (Hz)\n", fc);
   printf("Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   printf("Transmission bit Rate = %e (bps)\n", sym_rate*1);
   printf("f_d * t = %f\n", Doppler / sym_rate);
   printf("Number of bits of simulation = %d\n", K*num_packet);
   printf("Iterations of Turbo Equalization = %d\n\n", Iteration);
   printf("This program is running. Don't close, please!\n\n");

   records = fopen("Records_turbo_EQ_SOVA.log", "a");
   fprintf(records, "BER Performance of TURBO EQUALIZATION in ISI Channel\n");
   fprintf(records, "Coding Scheme: (2,1,2) Convolutional Code\n");
   fprintf(records, "Generator polynomials are {5,7} in Octal\n");
   fprintf(records, "Minimum free distance = 5\n");
   fprintf(records, "Decoder: Viterbi Decoder with SOVA algprithm\n");
   fprintf(records, "Interleaver: S-random interleaver with S = 16\n");
   fprintf(records, "3-path Rayleigh fading channel with equal strength\n");
   fprintf(records, "Multipath weighting factor: {0.577, 0.577, 0.577}\n");
   fprintf(records, "Equalizer: MLSE Equalizer with SOVA Algorithm\n");
   fprintf(records, "Speed of the vehicle = %f (km/h)\n", vc);
   fprintf(records, "Carrier Frequency = %e (Hz)\n", fc);
   fprintf(records, "Maximum Doppler Frequency = %f (Hz)\n", Doppler);
   fprintf(records, "Transmission bit Rate = %e (bps)\n", sym_rate*1);
   fprintf(records, "f_d * t = %f\n", Doppler / sym_rate);
   fprintf(records, "Number of bits of simulation = %d\n", K*num_packet);
   fprintf(records, "Iterations of Turbo Equalization = %d\n\n", Iteration);
   fprintf(records, "Eb/No     BER\n");
   fflush(records);

   data_bit = new int[K];
   coded_bit = new int[N];
   ch = new int[m+1];
   Ak = new double[N];				// Soft Decision
   Hk = new int[N];              // Hard output
   bpsk = new double[N];
   Yk = new double*[2];     		// Received signal
   for(i=0; i<2; i++)
   	Yk[i] = new double[N];
   LLR = new double[N];				// Log-likelihood Ratio
   z = new double[N];				// Extrinsic information
   intrinsic = new double[N];		// Intrinsic information
   extrinsic = new double[N];		// Extrinsic information
   S_random = new int[N];
   de_inter = new int[N];
   interleaved = new int[N];
   err_count = new int[Iteration+1];
   err_rate = new double[Iteration+1];

	srand((unsigned) time(&t));

   Trellis_diagram(n, k, m, &gp[0]);

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
	ber=fopen("ber_MLSE_SOVA_dec.log", "w");
   for(snr=0; snr<=10; snr++)
   {
      for(s=0; s<=Iteration; s++)
   	  	err_count[s] = 0;
   	// noise power calculation
      Eb_No = (double)snr;
      noise_pwr = 1.0/(((float)k/(float)n)*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      t1 = 100.0;  	// Initial time of Rayleigh fading pattern
      t2 = 200.0;  	// Initial time of the 2nd path
      t3 = 300.0;		// Innitial time of the 3rd path
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

            JakesFading(fc, vc*1000/3600.0, t1, 2, &fade1[0]);
      		t1 += 1.0/sym_rate;
            JakesFading(fc, vc*1000/3600.0, t2, 2, &fade2[0]);
      		t2 += 1.0/sym_rate;
            JakesFading(fc, vc*1000/3600.0, t3, 2, &fade3[0]);
      		t3 += 1.0/sym_rate;

            // Calculate channel output symbol
            ch_out_sym_I = ch[0]*CIR[0]*fade1[0] + ch[1]*CIR[1]*fade2[0] + ch[2]*CIR[2]*fade3[0];
            ch_out_sym_Q = ch[0]*CIR[0]*fade1[1] + ch[1]*CIR[1]*fade2[1] + ch[2]*CIR[2]*fade3[1];
            from_state = tmp_state & (ch_num_state-1); 			// to_state of trellis diagram

            CH_Trellis_diagram(i, Num_path, CIR, fade1, fade2, fade3);

            /* AWGN channel */
            AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[0][i] = ch_out_sym_I + noise[0];
            Yk[1][i] = ch_out_sym_Q + noise[1];
         }

/******************************************************************************/
/*                     TURBO EQUALIZATION (SOVA-SOVA)                         */
/******************************************************************************/
			// Initialize a-priori probability
         for(i=0; i<N; i++)
         	z[i] = 0.0;
         r = 2 * noise_pwr;

         // Iterative Equalization and Decoding
		  	for(s=0; s<=Iteration; s++)
		  	{
         	MLSE_SOVA(Yk, z, noise_pwr, Hk, LLR, N, ch_N_state, num_in_sym);

            // Extrinsic information, de-interleaving and Soft Decision
				for(i=0; i<N; i++)
            {
            	extrinsic[i] = LLR[i] - r*z[i];
               intrinsic[de_inter[i]] = extrinsic[i];
               Ak[de_inter[i]] = intrinsic[de_inter[i]];
               //Ak[de_inter[i]] = tanh(intrinsic[de_inter[i]]/2.0);
            }

         	Viterbi_dec_SOVA(Ak, Hk, LLR, K, num_state, num_in_sym);

            // Extrinsic information and Interleaving
            for(i=0; i<N; i++)
            {
            	extrinsic[i] = LLR[i] - intrinsic[i];
               z[S_random[i]] = extrinsic[i];
            }

            for(i=0; i<K; i++)	// Bit error count
             	err_count[s] += Error_count(data_bit[i], Hk[i]);
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
   delete Ak;
   delete Hk;
   delete bpsk;
   for(i=0; i<2; i++)
   	delete Yk[i];
   delete Yk;
   delete LLR;
   delete S_random;
   delete de_inter;
   delete interleaved;
   delete z;
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

void MLSE_SOVA(double **Data_in, double *intrinsic, double noise_pwr, int *Data_out,
				   double *Soft_out, int Packet_length, int Num_state, int Num_in)
{
/*********************************************************************/
/* MLSE: Soft Output Viterbi Algorithm (SOVA) for static ISI channel */
/*********************************************************************/
	int i, j, l, q, survivor_state;
   double **mju_f, *mju, *survival_metric, metric, ***branch_metric, **mju_b;
   double r=2*noise_pwr, mju_tmp, survivor_metric;

   struct surv        			// Data structure of survival path
   		{
         	double metric;		// Path metric
            int data_in[N];	// input bit stream, N: packet length, global
            int state[N];		// state transition sequence, N: packet length, global
         };
   struct surv survival_path[N_state], survival_temp[N_state];
   									// Survival_path[num_state], N_state: number of states of channel, globle

   survival_metric = new double[Num_state];
   mju = new double[Num_in];							// minimum path metric
   mju_f = new double*[Packet_length+1];			// forward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_f[i] = new double[Num_state];
   mju_b = new double*[Packet_length+1];       	// backward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_b[i] = new double[Num_state];
   branch_metric = new double**[Packet_length];	// branch[time index][state][input]
   for(i=0; i<Packet_length; i++)
      branch_metric[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         branch_metric[i][j] = new double[Num_in];

   // Initialize survival path
   for(i=0; i<Num_state; i++)
   {
   	survival_path[i].metric = 0.0;
      mju_f[0][i] = 0.0;
   }

/*********************/
/* Forward Recursion */
/*********************/
	for(i=0; i<Packet_length; i++)
   {
   	for(j=0; j<Num_state; j++)					// Initialize the survival path metric
      	survival_metric[j] = DBL_MAX;

      for(j=0; j<Num_state; j++)					// from_state index
      	for(l=0; l<Num_in; l++)					// input bit
         {
         	// branch metric, Euclidean Distance
            branch_metric[i][j][l] = 0.0;
            branch_metric[i][j][l] += (pow(Data_in[0][i]-ch_Trellis[i][j][l].out[0],2) - r*intrinsic[i]*l);
            branch_metric[i][j][l] += (pow(Data_in[1][i]-ch_Trellis[i][j][l].out[1],2) - r*intrinsic[i]*l);
            metric = survival_path[j].metric + branch_metric[i][j][l];

            // find the survival path metric
            if(metric < survival_metric[ch_Trellis[i][j][l].to])
            {
            	survival_metric[ch_Trellis[i][j][l].to] = metric;

               // Record and refresh the survival path
               for(q=0; q<i; q++)
               {
               	survival_temp[ch_Trellis[i][j][l].to].data_in[q] = survival_path[j].data_in[q];
                  survival_temp[ch_Trellis[i][j][l].to].state[q] = survival_path[j].state[q];
               }
               survival_temp[ch_Trellis[i][j][l].to].data_in[i] = l;
               survival_temp[ch_Trellis[i][j][l].to].state[i] = ch_Trellis[i][j][l].to;
            }
         }

      // Record and refresh the survival path
      for(j=0; j<Num_state; j++)					// to_state index
      {
      	survival_path[j].metric = survival_metric[j];
         mju_f[i+1][j] = survival_metric[j];
         for(q=0; q<=i; q++)
         {
         	survival_path[j].data_in[q] = survival_temp[j].data_in[q];
            survival_path[j].state[q] = survival_temp[j].state[q];
         }
      }
   }

   // Find the path with the smallest path metric
   survivor_metric = survival_path[0].metric;
   survivor_state = 0;
   for(j=1; j<Num_state; j++)
   	if(survivor_metric > survival_path[j].metric)
      {
      	survivor_metric = survival_path[j].metric;	// survivor path metric
         survivor_state = j;									// survivor state
      }

/****************************************/
/* Backward Recursion and Soft Decision */
/****************************************/
	// Initialize survival path
   for(j=0; j<Num_state; j++)
      mju_b[Packet_length][j] = 0.0;

   for(i=Packet_length-1; i>=0; i--)
   {
   	for(j=0; j<Num_state; j++)					// Initialize the survival path metric
      	survival_metric[j] = DBL_MAX;

      for(j=0; j<Num_state; j++)					// from_state index
      	for(l=0; l<Num_in; l++)					// input bit
         {
         	metric = mju_b[i+1][ch_Trellis[i][j][l].to] + branch_metric[i][j][l];

            // find the survival path metric
            if(metric < survival_metric[j])
            	survival_metric[j] = metric;
         }

      // Record the survival path metric
      for(j=0; j<Num_state; j++)					// from_state index
      	mju_b[i][j] = survival_metric[j];

      // LLR Calculation
      mju[survival_path[survivor_state].data_in[i]] = survivor_metric;	// mju_f[Packet_length][survivor_state];

      mju[(survival_path[survivor_state].data_in[i]+1)%2] = DBL_MAX;
      for(j=0; j<Num_state; j++)					// from_state index
      {
      	mju_tmp = mju_f[i][j] + branch_metric[i][j][(survival_path[survivor_state].data_in[i]+1)%2]
      				 + mju_b[i+1][ch_Trellis[i][j][(survival_path[survivor_state].data_in[i]+1)%2].to];

         if(mju_tmp < mju[(survival_path[survivor_state].data_in[i]+1)%2])
         	mju[(survival_path[survivor_state].data_in[i]+1)%2] = mju_tmp;
      }

      Soft_out[i] = mju[0] - mju[1];
      Data_out[i] = survival_path[survivor_state].data_in[i];
   }

	delete survival_metric;
   delete mju;
   for(i=0; i<=Packet_length; i++)
      delete mju_f[i];
   delete mju_f;
   for(i=0; i<=Packet_length; i++)
       delete mju_b[i];
   delete mju_b;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete branch_metric[i][j];
   for(i=0; i<Packet_length; i++)
      delete branch_metric[i];
   delete branch_metric;
}
/******************************************************************************/

void CH_Trellis_diagram(int time, int Num_path, float *CIR, double *fade1, double *fade2, double *fade3)
{
/*********************************************************/
/* Generate TrellisDiagram for time varying ISI channel */
/*********************************************************/
	int input, from_state, to_state, tmp_state, ch[3];
   int num_state = pow(2,Num_path-1), num_in_sym = 2;
   double out_sym_I, out_sym_Q;

   for(from_state=0; from_state<num_state; from_state++) // from_state of trellis diagram
   	for(input=0; input<num_in_sym; input++)				// input of trellis diagram
      {
      	tmp_state = from_state;
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         ch[0] = 2*(tmp_state & 0x01) - 1;
         ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;
         ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;
         out_sym_I = ch[0]*CIR[0]*fade1[0] + ch[1]*CIR[1]*fade2[0] + ch[2]*CIR[2]*fade3[0];
         out_sym_Q = ch[0]*CIR[0]*fade1[1] + ch[1]*CIR[1]*fade2[1] + ch[2]*CIR[2]*fade3[1];
         to_state = tmp_state & (num_state-1); 					// to_state of trellis diagram
         ch_Trellis[time][from_state][input].from = from_state;
         ch_Trellis[time][from_state][input].to = to_state;
         ch_Trellis[time][from_state][input].in = (2*input-1);
         ch_Trellis[time][from_state][input].out[0] = out_sym_I;
         ch_Trellis[time][from_state][input].out[1] = out_sym_Q;
      }
}

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

void Viterbi_dec_SOVA(double *Data_in, int *Data_out, double *Soft_out, int Packet_length, int Num_state, int Num_in)
{
/**********************************************************************/
/* Soft Output Viterbi Algorithm decoder (SOVA)                       */
/* Convolutional Decoder (n=2, k=1, m=2) Generator polynomial: {5, 7} */
/**********************************************************************/
	int i, j, l, q;
   double **mju_f, *mju, *survival_metric, metric, ***branch_metric, **mju_b, mju_tmp;

   struct surv        			// Data structure of survival path
   		{
         	double metric;		// Path metric
            int data_in[K];	// input bit stream, K: packet length, global
            int code_bit[N];	// output code bit stream, N: packet length, global
            int state[K];		// state transition sequence, K: packet length, global
         };
   struct surv survival_path[N_state], survival_temp[N_state];
   									// Survival_path[N_state], N_state: number of states of CC code, global

   survival_metric = new double[Num_state];
   mju = new double[Num_in];							// minimum path metric
   mju_f = new double*[Packet_length+1];			// forward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_f[i] = new double[Num_state];
   mju_b = new double*[Packet_length+1];       	// backward path-metric[time index][state]
   for(i=0; i<=Packet_length; i++)
      mju_b[i] = new double[Num_state];
   branch_metric = new double**[Packet_length];	// branch[time index][state][input]
   for(i=0; i<Packet_length; i++)
      branch_metric[i] = new double*[Num_state];
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         branch_metric[i][j] = new double[Num_in];

   // Initialize survival path
   for(i=0; i<Num_state; i++)
   {
// 	survival_path[i].metric = DBL_MAX;	// Initial maximum value for Euclidean distance
		survival_path[i].metric = -DBL_MAX;	// Initial minimum value for cross-correlation
//    mju_f[0][i] = DBL_MAX;					// Initial maximum value for Euclidean distance
//    mju_b[K][i] = DBL_MAX;					// Initial maximum value for Euclidean distance
		mju_f[0][i] = -DBL_MAX;					// Initial minimum value for cross-correlation
      mju_b[K][i] = -DBL_MAX;					// Initial minimum value for cross-correlation
	}
   survival_path[0].metric = 0.0;
   mju_f[0][0] = 0.0;
   mju_b[K][0] = 0.0;

/*********************/
/* Forward Recursion */
/*********************/
	for(i=0; i<Packet_length; i++)
   {
   	for(j=0; j<Num_state; j++)					// Initialize the survival path metric
      	survival_metric[j] = -DBL_MAX;

      for(j=0; j<Num_state; j++)					// from_state index
      	for(l=0; l<Num_in; l++)					// input bit
         {
         	// branch metric, Euclidean Distance
/*          branch_metric[i][j][l] = 0.0;
				branch_metric[i][j][l] += pow(Data_in[2*i]-(2*Trellis[j][l].out[0]-1),2);		// brahch metric
            branch_metric[i][j][l] += pow(Data_in[2*i+1]-(2*Trellis[j][l].out[1]-1),2);	// branch metric
            metric = survival_path[j].metric + branch_metric[i][j][l];
*/
				branch_metric[i][j][l] = 0.0;		// branch metric, Cross-correlation
            branch_metric[i][j][l] += (Data_in[2*i] * (2*Trellis[j][l].out[0]-1));		// brahch metric
            branch_metric[i][j][l] += (Data_in[2*i+1] * (2*Trellis[j][l].out[1]-1));	// branch metric
            metric = survival_path[j].metric + branch_metric[i][j][l];

            // find the survival path metric
//          if(metric < survival_metric[Trellis[j][l].to])	//	Euclidean distance (Minimize)
				if(metric > survival_metric[Trellis[j][l].to])	// Cross-correlation (Maximize)
            {
            	survival_metric[Trellis[j][l].to] = metric;

               // Record and refresh the survival path
               for(q=0; q<i; q++)
               {
               	survival_temp[Trellis[j][l].to].data_in[q] = survival_path[j].data_in[q];
                  survival_temp[Trellis[j][l].to].code_bit[2*q] = survival_path[j].code_bit[2*q];
                  survival_temp[Trellis[j][l].to].code_bit[2*q+1] = survival_path[j].code_bit[2*q+1];
                  survival_temp[Trellis[j][l].to].state[q] = survival_path[j].state[q];
               }
               survival_temp[Trellis[j][l].to].data_in[i] = l;
               survival_temp[Trellis[j][l].to].code_bit[2*i] = Trellis[j][l].out[0];
               survival_temp[Trellis[j][l].to].code_bit[2*i+1] = Trellis[j][l].out[1];
               survival_temp[Trellis[j][l].to].state[i] = Trellis[j][l].to;
            }
			}

		// Record and refresh the survival path
      for(j=0; j<Num_state; j++)		// to_state index
      {
      	survival_path[j].metric = survival_metric[j];
         mju_f[i+1][j] = survival_metric[j];
         for(q=0; q<=i; q++)
         {
         	survival_path[j].data_in[q] = survival_temp[j].data_in[q];
            survival_path[j].code_bit[2*q] = survival_temp[j].code_bit[2*q];
            survival_path[j].code_bit[2*q+1] = survival_temp[j].code_bit[2*q+1];
            survival_path[j].state[q] = survival_temp[j].state[q];
         }
      }
	}

/****************************************/
/* Backward Recursion and Soft Decision */
/****************************************/
	for(i=Packet_length-1; i>=0; i--)
   {
   	for(j=0; j<Num_state; j++)				// Initialize the survival path metric
//    	survival_metric[j] = DBL_MAX;		// Initial maximum value for Euclidean distance
			survival_metric[j] = -DBL_MAX;	// Initial minimum value for cross-correlation

		for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)		// input bit
         {
         	metric = mju_b[i+1][Trellis[j][l].to] + branch_metric[i][j][l];

            // find the survival path metric
//          if(metric < survival_metric[j])	//	Euclidean distance (Minimize)
				if(metric > survival_metric[j])	// Cross-correlation (Maximize)
            	survival_metric[j] = metric;
			}

		// Record the survival path metric
      for(j=0; j<Num_state; j++)		// from_state index
      	mju_b[i][j] = survival_metric[j];
/*
      // LLR Calculation for the information bit
      mju[survival_path[0].data_in[i]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].data_in[i]+1)%2] = DBL_MAX;		//	Euclidean distance (Minimize)
		mju[(survival_path[0].data_in[i]+1)%2] = -DBL_MAX;    // Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      {
      	mju_tmp = mju_f[i][j] + branch_metric[i][j][(survival_path[0].data_in[i]+1)%2]
         			 + mju_b[i+1][Trellis[j][(survival_path[0].data_in[i]+1)%2].to];

//       if(mju_tmp < mju[(survival_path[0].data_in[i]+1)%2])	//	Euclidean distance (Minimize)
			if(mju_tmp > mju[(survival_path[0].data_in[i]+1)%2])	// Cross-correlation (Maximize)
         	mju[(survival_path[0].data_in[i]+1)%2] = mju_tmp;
		}

//    Soft_out[i] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[i] = mju[1] - mju[0];		// Cross-correlation
*/
      // LLR Calculation (for the 1st code bit)
      mju[survival_path[0].code_bit[2*i]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].code_bit[2*i]+1)%2] = DBL_MAX;			//	Euclidean distance (Minimize)
		mju[(survival_path[0].code_bit[2*i]+1)%2] = -DBL_MAX;    	// Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)
         	if(Trellis[j][l].out[0] != survival_path[0].code_bit[2*i])
            {
            	mju_tmp = mju_f[i][j] + branch_metric[i][j][l] + mju_b[i+1][Trellis[j][l].to];
//             if(mju_tmp < mju[(survival_path[0].code_bit[2*i]+1)%2])	//	Euclidean distance (Minimize)
					if(mju_tmp > mju[(survival_path[0].code_bit[2*i]+1)%2])	// Cross-correlation (Maximize)
               	mju[(survival_path[0].code_bit[2*i]+1)%2] = mju_tmp;
				}

//    Soft_out[2*i] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[2*i] = mju[1] - mju[0];		// Cross-correlation

		// LLR Calculation (for the 2nd code bit)
      mju[survival_path[0].code_bit[2*i+1]] = mju_f[Packet_length][0];

//    mju[(survival_path[0].code_bit[2*i+1]+1)%2] = DBL_MAX;		//	Euclidean distance (Minimize)
		mju[(survival_path[0].code_bit[2*i+1]+1)%2] = -DBL_MAX;		// Cross-correlation (Maximize)
      for(j=0; j<Num_state; j++)		// from_state index
      	for(l=0; l<Num_in; l++)
         	if(Trellis[j][l].out[1] != survival_path[0].code_bit[2*i+1])
            {
            	mju_tmp = mju_f[i][j] + branch_metric[i][j][l] + mju_b[i+1][Trellis[j][l].to];
//   	         if(mju_tmp < mju[(survival_path[0].code_bit[2*i+1]+1)%2])	//	Euclidean distance (Minimize)
					if(mju_tmp > mju[(survival_path[0].code_bit[2*i+1]+1)%2])	// Cross-correlation (Maximize)
               	mju[(survival_path[0].code_bit[2*i+1]+1)%2] = mju_tmp;
				}

//    Soft_out[2*i+1] = mju[0] - mju[1];		// Euclidean Distance
		Soft_out[2*i+1] = mju[1] - mju[0];		// Cross-correlation

      Data_out[i] = survival_path[0].data_in[i];
   }

	delete survival_metric;
   delete mju;
   for(i=0; i<=Packet_length; i++)
      delete mju_f[i];
   delete mju_f;
   for(i=0; i<=Packet_length; i++)
       delete mju_b[i];
   delete mju_b;
   for(i=0; i<Packet_length; i++)
   	for(j=0; j<Num_state; j++)
         delete branch_metric[i][j];
   for(i=0; i<Packet_length; i++)
      delete branch_metric[i];
   delete branch_metric;
}
/******************************************************************************/

void JakesFading(double f_c/*Hz*/, float v/*m/s*/, double t/*s*/, int type, double *fade)
{
	const double C = 3.0e8;     // (m/s)
   const float Pi = 3.14159265358979;
   int n, N, N_o = 32;
   double lamda, w_m, beta_n, w_n, alpha, T_c2, T_s2, theta_n;

   lamda = C/f_c;     // wave length (meter)
   w_m = 2.0*Pi*v/lamda;    // maximum Doppler frequency
   N = 2*(2*N_o+1);

   switch(type)
   {
   	case 1:
   		alpha = 0.0;
         T_c2 = (double)N_o;
         T_s2 = (double)N_o + 1.0;
         break;
      case 2:
      	alpha = 0.0;
         T_c2 = (double)N_o + 1.0;
         T_s2 = (double)N_o;
         break;
      case 3:
      	alpha = Pi/4.0;
         T_c2 = (double)N_o + 0.5;
         T_s2 = (double)N_o + 0.5;
         break;
      default:
      	printf("\nInvalid type selection for Jake's fading channel model.\n");
         break;
   }

   if(v == 0.0)
   {
   	*(fade+0) = 1.0;
      *(fade+1) = 0.0;
   }
   else
   {
   	*(fade+0) = sqrt(1.0/T_c2)*cos(alpha)*cos(w_m*t);
      *(fade+1) = sqrt(1.0/T_s2)*sin(alpha)*cos(w_m*t);

      for(n = 1; n <= N_o; n++)
      {
      	switch(type)
         {
         	case 1:
            	beta_n = (double)n*Pi/((double)N_o+1.0);
               break;
            case 2:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
            case 3:
            	beta_n = (double)n*Pi/(double)N_o;
               break;
         	default:
            	break;
         }
         w_n = w_m*cos(2.0*Pi*(double)n/(double)N);
//            theta_n = 2.0*Pi*((double)rand()/(double)RAND_MAX);  // random phase
			theta_n = 0.0;
         *(fade+0) += sqrt(2.0/T_c2)*cos(beta_n)*cos(w_n*t+theta_n);
         *(fade+1) += sqrt(2.0/T_s2)*sin(beta_n)*cos(w_n*t+theta_n);
		}
	}
}

void AWGN_noise(float mu, double variance, double *noise)
{
	const  float Pi = 3.14159265358979;
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

