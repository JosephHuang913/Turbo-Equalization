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

void AWGN_noise(float, double, double *);
int Error_count(int, int);

#define Pi 		3.14159265358979
#define n		2					// (n,k,m) Convolutional Code, # of output symbols
#define k		1					// (n,k,m) Convolutional Code, # of input symbols
#define m		2					// (n,k,m) Convolutional Code, # of memories
#define ch_m 	2              // channel memory of ISI channel
#define K		512		  		// packet length of information bit stream
#define N		K*(n/k)	 		// packet length of coded bit stream
#define num_packet 20000			// number of packets simulated
#define Iteration	 10				// Iterations of Turbo Equalization ( >= 1 )
const int num_state = pow(2,m);		// number of states
const int ch_num_state = pow(2,ch_m);	//number of channel states
const int num_in_sym = pow(2,k);		// number of input symbols
const int gp[2] = {5,7};   // Generator polynomial of CC given in hexadecimal and in reverse order
const float CIR[3] = {0.407,0.815,0.407};	// Channel Impulse Response

int main(void)
{
	time_t  t, start, end;
	int i, j, l, p, q, input, out_sym, out_bit, *data_bit, *coded_bit;
   int from_state, to_state, tmp_state, *Ak, err_count, *ch, state;
   int *S_random, *de_inter, *interleaved;
   double a_priori=0.5, snr, Eb_No, noise_pwr, noise[2], ch_out_sym, err_rate;
   double p1, sum, **alpha, **beta, ***gamma, *delta, *LLR, min, *Yk;
   double *extrinsic[2], *intrinsic[2];
   FILE *trelis, *ber, *records, *s_inter;

   struct trel_diag			// Data structure of trellis diagram for each branch
   		{
         	int from;		// from_state of trellis diagram
            int to;			// to_state of trellis diagram
            int in;			// input data bit of trellis diagram
            int out[n];		// output codeword symbol of trellis diagram
         };
   struct trel_diag Trellis[4][2];			// Trellis[num_state][num_in_sym]

   struct ch_trel_diag		// Data structure of trellis diagram for each branch
   		{
         	int from;		// from_state of trellis diagram
            int to;			// to_state of trellis diagram
            int in;			// input data symbol of trellis diagram
            float out;		// output symbol of trellis diagram
         };
   struct ch_trel_diag ch_Trellis[4][2];	// Trellis[num_state][num_in_sym]

   start = time(NULL);
   printf("BER Performance of TURBO EQUALIZATION in ISI Channel\n");
	printf("Coding Scheme: (2,1,2) Convolutional Code\n");
	printf("Generator polynomials are {5,7} in Octal\n");
   printf("Minimum free distance = 5\n");
   cout << "Decoder: MAP Decoder with modified BCJR algprithm" << endl;
   cout << "Interleaver: S-random interleaver with S = 16" << endl;
   cout << "Channel Impulse Response: {0.407,0.815,0.407}" << endl;
   cout << "Equalizer: MAP Equalizer with modified BCJR Algorithm" << endl;
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
   S_random = new int[N];
   de_inter = new int[N];
   interleaved = new int[N];
   ch = new int[ch_m+1];
   Ak = new int[K];
   Yk = new double[N];
   delta = new double[num_in_sym];
   LLR = new double[N];
   intrinsic[0] = new double[N];
   intrinsic[1] = new double[N];
   extrinsic[0] = new double[N];
   extrinsic[1] = new double[N];
   alpha = new double*[N+1];			  		// alpha[time index][state]
   for(i=0; i<=N; i++)
      alpha[i] = new double[num_state];
   beta = new double*[N+1];       			// beta[time index][state]
   for(i=0; i<=N; i++)
      beta[i] = new double[num_state];
   gamma = new double**[N];					// gamma[time index][state][input]
   for(i=0; i<N; i++)
      gamma[i] = new double*[num_state];
   for(i=0; i<N; i++)
   	for(j=0; j<num_state; j++)
         gamma[i][j] = new double[num_in_sym];

	srand((unsigned) time(&t));

/**************************************************************/
/* Generate TrellisDiagram for (n=2, k=1, m=2) {5,7} CC code */
/**************************************************************/
	trelis = fopen("Trellis_diagram.log", "w");
   fprintf(trelis, "Trellis diagram of (%d,%d,%d) convolutional code\n", n, k, m);
   fprintf(trelis, "Generator polynomials are {%d,%d}\n\n", gp[0], gp[1]);
   fprintf(trelis, "s(k-1) s(k) input out5 out7\n");
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
         fprintf(trelis, "%4d %4d %5d %5d %4d\n", from_state, to_state, input, ((out_sym>>1)&1), out_sym&1);
      }
   }
   fclose(trelis);

/****************************************************************/
/* Generate TrellisDiagram for {0.407,0.815,0.407} ISI channel */
/****************************************************************/
	trelis = fopen("ch_Trellis_diagram.log", "w");
   fprintf(trelis, "Trellis diagram of (%.3f,%.3f,%.3f) ISI channel\n", CIR[0], CIR[1], CIR[2]);
   fprintf(trelis, "s(k-1) s(k) input output\n");
   for(from_state=0; from_state<num_state; from_state++) // from_state of trellis diagram
   {
   	for(input=0; input<num_in_sym; input++)		// input symbol (2*input-1) of trellis diagram
      {
      	tmp_state = from_state;
         tmp_state = (tmp_state << 1) ^ (input & 0x01);  // read input bit
         ch[0] = 2*(tmp_state & 0x01) - 1;
         ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;
         ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;
         ch_out_sym = ch[0]*CIR[0] + ch[1]*CIR[1] + ch[2]*CIR[2];		// Calculate output symbol
         to_state = tmp_state & (num_state-1); 					// to_state of trellis diagram
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

/************************/
/* main simulation loop */
/************************/
   for(snr=-2; snr<=10; snr++)
   {
   	err_count = 0;
      Eb_No = (double)snr;	// noise power calculation
      noise_pwr = 1.0/(((float)k/(float)n)*pow(10.0, Eb_No/10.0));	// BPSK, Nyquist filter assumption
      printf("Eb_No = %f, ", Eb_No);
      state = random(4); // Initial state of channel
      for(p=0; p<num_packet; p++)
      {
   		for(i=0; i<K-2; i++)
				data_bit[i] = random(2);		// Generate random information bit stream
   		data_bit[K-2] = data_bit[K-1] = 0;

/*********************************************************************/
/* Convolutional Encoder (n=2, k=1, m=2) Generator polynomial: {5,7} */
/*********************************************************************/
   		from_state = 0;	// Initial state of convolutional code
   		for(i=0; i<K; i++)
   		{
				tmp_state = from_state;
      		tmp_state = (tmp_state << 1) ^ (data_bit[i] & 0x01);  // read input bit
      		for(j=0; j<n; j++)
      		{
      			out_bit = 0;	// output bit of trellis diagram
         		for(l=m; l>=0; l--)
         			out_bit ^= ((tmp_state & gp[j]) >> l) & 1;	// Calculate output bit

         		coded_bit[2*i+j] = out_bit;
      		}
      		from_state = tmp_state & (num_state-1);	// to_state of trellis diagram
   		}

         // Interleaving
         for(i=0; i<N; i++)
         	interleaved[S_random[i]] = coded_bit[i];

/*****************************************************/
/* BPSK mapping and ISI Channel: {0.407,0.815,0.407} */
/*****************************************************/
   		for(i=0; i<N; i++)
   		{
         	tmp_state = state;
            tmp_state = (tmp_state << 1) ^ (interleaved[i] & 0x01);  // read input bit
            ch[0] = 2*(tmp_state & 0x01) - 1;					// input symbol (BPSK)
            ch[1] = 2*((tmp_state >> 1) & 0x01) - 1;        // channel memory (BPSK)
            ch[2] = 2*((tmp_state >> 2) & 0x01) - 1;        // channel memory (BPSK)
            ch_out_sym = ch[0]*CIR[0] + ch[1]*CIR[1] + ch[2]*CIR[2];	// Calculate output symbol
            state = tmp_state & (ch_num_state-1); 			// to_state of trellis diagram

            /* AWGN channel */
            AWGN_noise(0, noise_pwr, &noise[0]);
         	Yk[i] = ch_out_sym + noise[0];
   		}

/**************************************************/
/*               TURBO EQUALIZATION               */
/**************************************************/
         for(l=0; l<2; l++)	// a-priori probability
         	for(i=0; i<N; i++)
         		intrinsic[l][i] = 0.5;

         q = 0;
         do
         {
/*********************************************/
/* MAP Equalizer (normalized BCJR Algorithm) */
/*********************************************/
				// Initialization of alpha and beta
         	for(i=1; i<=N; i++)				// alpha[time index][state]
         		for(j=0; j<num_state; j++)
         			alpha[i][j] = 0.0;
         	for(j=0; j<num_state; j++)
         		alpha[0][j] = 1/(float)num_state;

         	for(i=0; i<N; i++)           // beta[time index][state]
         		for(j=0; j<num_state; j++)
         			beta[i][j] = 0.0;
         	for(j=0; j<num_state; j++)
         		beta[N][j] = 1.0;

				// calculate gamma[time index][state][input]
         	for(i=0; i<N; i++)	// time index
         		for(j=0; j<num_state; j++)		// state index
            		for(l=0; l<num_in_sym; l++)	// input symbol
               	{
							p1 = exp(-pow(Yk[i]-ch_Trellis[j][l].out,2)/(2*noise_pwr))/sqrt(2*Pi*noise_pwr);
         				gamma[i][j][l] = intrinsic[l][i] * p1;		// gamma[time index][state][input]
               	}

         	// calculate alpha[time index][state]
         	for(i=1; i<=N; i++)		// time index
         	{
         		for(j=0; j<num_state; j++)		// from_state index
            		for(l=0; l<num_in_sym; l++)	// input bit
         				alpha[i][ch_Trellis[j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

            	sum = 0.0;		// for renormalization
	            for(j=0; j<num_state; j++)		// to_state index
   	         	sum += alpha[i][j];

      	      for(j=0; j<num_state; j++)
         	   	alpha[i][j] = alpha[i][j] / sum;
         	}

	         // calculate beta[time index][state]
   	      for(i=N-1; i>0; i--)		// time index
      	   {
         		for(j=0; j<num_state; j++)		// from_state index
            		for(l=0; l<num_in_sym; l++)	// input bit
         				beta[i][j] += beta[i+1][ch_Trellis[j][l].to] * gamma[i][j][l];

	            sum = 0.0;		// for renormalization
   	         for(j=0; j<num_state; j++)		// from_state index
      	      	sum += beta[i][j];

         	   for(j=0; j<num_state; j++)
            		beta[i][j] = beta[i][j] / sum;
	         }

	         // calculate conditional LLR and extrinsic information
   	      for(i=0; i<N; i++)		// time index
      	   {
         		min = 0.0;		// find the minimum product of alpha*gamma*beta
            	for(j=0; j<num_state; j++)
            		for(l=0; l<num_in_sym; l++)
	               {
   	               delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][ch_Trellis[j][l].to];

      	            if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
         	         	min = delta[0];
						}

	            if(min == 0.0 || min > 1.0)	// if all else fails, make min real small
   	         	min = 1E-100;

      	   	delta[0] = delta[1] = 0.0;
         		for(j=0; j<num_state; j++)		// from_state index
            		for(l=0; l<num_in_sym; l++)	// input bit
               		delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][ch_Trellis[j][l].to];

	            if(delta[1] == 0.0)
   	         	delta[1] = min;
//    	        else if(delta[0] == 0.0)
         	   if(delta[0] == 0.0)
            		delta[0] = min;

/***********************************************/
/* log(A/B) = log(A) - log(B)                  */
/* log(d1/d2) - log(p1/p2) = extrinsic         */
/* log((d1/d2)/(p1/p2)) = log((d1/p1)/(d2/p2)) */
/* extrinsic information = (d1/p1) and (d2/p2) */
/***********************************************/

					for(l=0; l<2; l++)	// calculate extrinsic information
               	extrinsic[l][i] = delta[l]/intrinsic[l][i];

//	            LLR[i] = log(delta[1]/delta[0]);
   	      }

            // De-interleave
            for(i=0; i<N; i++)
            	for(l=0; l<2; l++)
	            	intrinsic[l][de_inter[i]] = extrinsic[l][i];

/*******************************************/
/* MAP decoder (normalized BCJR Algorithm) */
/*******************************************/
				// Initialization of alpha and beta
   	      for(i=0; i<=K; i++)				// alpha[time index][state]
      	   	for(j=0; j<num_state; j++)
         			alpha[i][j] = 0.0;
	         alpha[0][0] = 1.0;

   	      for(i=0; i<=K; i++)           // beta[time index][state]
      	   	for(j=0; j<num_state; j++)
         			beta[i][j] = 0.0;
	         beta[K][0] = 1.0;

				// calculate gamma[time index][state][input]
      	   for(i=0; i<K; i++)	// time index
         		for(j=0; j<num_state; j++)		// state index
            		for(l=0; l<num_in_sym; l++)	// input symbol
   	      			gamma[i][j][l] = a_priori * intrinsic[Trellis[j][l].out[0]][2*i]
                     				     * intrinsic[Trellis[j][l].out[1]][2*i+1];

         	// calculate alpha[time index][state]
	         for(i=1; i<=K; i++)		// time index
   	      {
      	   	for(j=0; j<num_state; j++)		// from_state index
         	   	for(l=0; l<num_in_sym; l++)	// input bit
         				alpha[i][Trellis[j][l].to] += alpha[i-1][j] * gamma[i-1][j][l];

	            sum = 0.0;		// for renormalization
   	         for(j=0; j<num_state; j++)		// to_state index
      	      	sum += alpha[i][j];

         	   for(j=0; j<num_state; j++)
            		alpha[i][j] = alpha[i][j] / sum;
	         }

   	      // calculate beta[time index][state]
      	   for(i=K-1; i>0; i--)		// time index
         	{
         		for(j=0; j<num_state; j++)		// from_state index
            		for(l=0; l<num_in_sym; l++)	// input bit
         				beta[i][j] += beta[i+1][Trellis[j][l].to] * gamma[i][j][l];

	            sum = 0.0;		// for renormalization
   	         for(j=0; j<num_state; j++)		// from_state index
      	      	sum += beta[i][j];

         	   for(j=0; j<num_state; j++)
            		beta[i][j] = beta[i][j] / sum;
	         }

            // Calculate Extrinsic Information
      	   for(i=0; i<K; i++)		// time index
         	{
         		min = 0.0;		// find the minimum product of alpha*gamma*beta
	            for(j=0; j<num_state; j++)
   	         	for(l=0; l<num_in_sym; l++)
      	         {
         	         delta[0] = alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

            	      if((delta[0] < min && delta[0] != 0.0) || min == 0.0)
               	   	min = delta[0];
						}

   	         if(min == 0.0 || min > 1.0)	// if all else fails, make min real small
      	      	min = 1E-100;

/***********************************************/
/* log(A/B) = log(A) - log(B)                  */
/* log(d1/d2) - log(p1/p2) = extrinsic         */
/* log((d1/d2)/(p1/p2)) = log((d1/p1)/(d2/p2)) */
/* extrinsic information = (d1/p1) and (d2/p2) */
/***********************************************/

               delta[0] = delta[1] = 0.0;
               for(j=0; j<num_state; j++)		// from_state index
               	for(l=0; l<num_in_sym; l++)
                  	delta[Trellis[j][l].out[0]] += alpha[i][j] * gamma[i][j][l]
                     										 * beta[i+1][Trellis[j][l].to];

               if(delta[1] == 0.0)
   	         	delta[1] = min;
         	   if(delta[0] == 0.0)
            		delta[0] = min;

					for(l=0; l<2; l++)	// calculate extrinsic information
               	extrinsic[l][2*i] = delta[l]/intrinsic[l][2*i];

               delta[0] = delta[1] = 0.0;
               for(j=0; j<num_state; j++)		// from_state index
               	for(l=0; l<num_in_sym; l++)
                  	delta[Trellis[j][l].out[1]] += alpha[i][j] * gamma[i][j][l]
                     										 * beta[i+1][Trellis[j][l].to];

               if(delta[1] == 0.0)
   	         	delta[1] = min;
         	   if(delta[0] == 0.0)
            		delta[0] = min;

               for(l=0; l<2; l++)	// calculate extrinsic information
               	extrinsic[l][2*i+1] = delta[l]/intrinsic[l][2*i+1];

/*         		delta[0] = delta[1] = 0.0;
         		for(j=0; j<num_state; j++)		// from_state index
            		for(l=0; l<num_in_sym; l++)	// input bit
               		delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

	            if(delta[1] == 0.0)
   	         	delta[1] = min;
//    	        else if(delta[0] == 0.0)
         	   if(delta[0] == 0.0)
            		delta[0] = min;

	            LLR[i] = log(delta[1]/delta[0]); */
   	      }

            // Interleave
            for(i=0; i<N; i++)
            	for(l=0; l<2; l++)
	            	intrinsic[l][S_random[i]] = extrinsic[l][i];

            q++;
			}
         while(q<=Iteration);

         // calculate conditional LLR
         for(i=0; i<K; i++)	// time index
         {
         	delta[0] = delta[1] = 0.0;
            for(j=0; j<num_state; j++)		// from_state index
            	for(l=0; l<num_in_sym; l++)	// input bit
               	delta[l] += alpha[i][j] * gamma[i][j][l] * beta[i+1][Trellis[j][l].to];

            if(delta[1] == 0.0)
            	delta[1] = min;
//    	   else if(delta[0] == 0.0)
				if(delta[0] == 0.0)
            	delta[0] = min;

            LLR[i] = log(delta[1]/delta[0]);
         }

/******************************************/
/* Data Decision, Statistics, and Records */
/******************************************/
			for(i=0; i<K; i++)	// data decision
         {
         	if(LLR[i]>=0)
            	Ak[i] = 1;
            else
            	Ak[i] = 0;

            err_count += Error_count(data_bit[i],Ak[i]);
         }
		}

      // Statistics and records
      err_rate = err_count / (double)(K*num_packet);
      printf("Error rate = %e\n", err_rate);
      ber=fopen("ber.log","a");
      fprintf(ber, "%f %e\n", Eb_No, err_rate);
   	fprintf(records, "%f %e\n", Eb_No, err_rate);
      fflush(records);
      fclose(ber);
   }

   delete data_bit;
   delete coded_bit;
   delete S_random;
   delete de_inter;
   delete interleaved;
   delete ch;
   delete Ak;
   delete Yk;
   delete delta;
   delete LLR;
   delete intrinsic[0];
   delete intrinsic[1];
   delete extrinsic[0];
   delete extrinsic[1];
   delete alpha;
   delete beta;
   delete gamma;

   end = time(NULL);
   printf("Total elapsed time: %.0f(sec)\n", difftime(end,start));
   fprintf(records, "\nTotal elapsed time: %.0f(sec)\n\n", difftime(end,start));
   fclose(records);
   printf("This program is ended. Press any key to continue.\n");
   getch();

   return 0;
}

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

