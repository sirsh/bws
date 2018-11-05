//sa comment: different system rules- for binary tree we want powers of 2 + 1, for lattice powers of 2 - 1, keep track

/* GP comment, meeting 16/03
*Branching Wiener Sausage - see readme for details and /scripts for compilation and execution notes*/

//todo document possibly subtle things e.g. we use flags on forloops in writing to decide when to write out things like histograms, lattices, avalanches etc. see for loops
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined linux || defined __APPLE__
#include <getopt.h>
#include <unistd.h>
#endif
#include "lattice_core.h"
#include "bws.h"

int _update_trace(long double key, int L, int D, int ml);
int init_cells(int D);
int get_cell_trace();
int reset();
void lattice_remove(int p);

//write the trace
#define WRITE_TRACE_SIZE(N, L, t) if (__verbose==1 && __trace > 1){printf("#MAX_TRACE(%.3Lf,%d,%d,%f):%d\n", t, L, N,radius_gyr,__trace);}
#define PRINT_VERSION {printf(VERSION); printf("\n");}
#define PRINT_PARAMS printf("#Parameters={ Seed: %i, Hopping rate: %g, Realisations: %i, Chunk size: %i, Dimension: %i, (Max) Lattice size: %i, Graph Type: %d }\n", seed, h, N, CHUNK_SIZE, D, MAX_L,__graph_type__);
#define PRINT_OPTIONS printf("#TODO\n");
//we have predetermined write times at 1,2,3,...9,10,20,....90,100,200,.....900,1000,2000,...
#define WRITE_TIME_FUNCTION(i) (MIN_T * pow(MAX_T/MIN_T, (double)i/(BINS-1)))
//#define INIT_WRITE_TIMES {int i=0; int a =0; int b=0; for(a =-2; a < TIME_WRITE_SCALES-2; a++){for (b = 1; b <= TIME_WRITE_COUNTS; b++){ write_times[i++] = (b)*powl(10,a); }}}
#define INIT_WRITE_TIMES {int i=0; for(i =0;i < BINS; i++){write_times[i] = WRITE_TIME_FUNCTION(i);}}
//moments for just the trace (todo add others) we note that MAX_MOMENTS does not incude the zero'th moment so we add MAX_MOMENTS+1 - this is just a convention
#define INIT_MOMENTS_TRACE {int t; int m; for(m = 0; m < MAX_MOMENTS_ALL; m++){for(t = 0; t < BINS; t++){trace_moments[t][m]=0;  alive_bit[t]=0;}}; INIT_WRITE_TIMES; }
//prepare an array to store a histogram of for trace counts i each sample path - it is assumed to be between 1 and L^2 for practicality - check?
#define INIT_TRACE_HISTOGRAM(L) {int l;  for (l=0;l< HIST_BINS*write_hist; l++ ){__trace_histogram[l]=0;}}
//prepare and set to 0 the avalanche matrix
#define INIT_AVALACNHES {int l, m; 	for (m = 0; m <= MAX_MOMENTS; m++) {for (l = 0; l < 10; l++) {avalanche_moments[m] = 0;}}}
//write each recorded obs(t) - in this case just one - extend with an offset for other moments e.g. 1-8 for trace, 1-8 for population count
#define WRITE_ITEM(x,t,offset)  {int i; for(i = 0; i <= MAX_MOMENTS; i++){trace_moments[t][i+offset] += (long double)pow(x,i);}}
//WRITE_ITEM(p,t,MAX_MOMENTS+1);WRITE_ITEM(cp,t,2*(MAX_MOMENTS+1));
#define WRITE_MOMENTS_TRACE(x,cp,p,t,alive) {alive_bit[t] += alive; WRITE_ITEM(x,t,0);}//
//WRITE_ITEM(p,t,MAX_MOMENTS+1);WRITE_ITEM(cp,t,2*(MAX_MOMENTS+1));} - this was written above but i dont want it
//write the header for the tabular data - make sure it matches what is in COMMIT_MOMENTS_TRACE
//for(i = 0; i <= MAX_MOMENTS; i++){printf("\tP%d",i);} ;for(i = 0; i <= MAX_MOMENTS; i++){printf("\tC%d",i);} 
#define WRITE_HEADER {int i;printf("\nchunk\tL\th\tsig\tD\tBCs\tt"); for(i = 0; i <= MAX_MOMENTS; i++){printf("\tM%d",i);} ;printf("\n");}
//flush out moments to standard out - add some header info that will be used in tabular data handling
//why init movement trace at end ?
#define COMMIT_MOMENTS_TRACE(c,L,h, sig,D,BCs){ int t; int m; for(t = 0; t < BINS; t++){printf("%d\t%d\t%f\t%f\t%d\t%d\t%.3Lf",c,L,h, sig,D,BCs,write_times[t]); for(m = 0; m < MAX_MOMENTS_ALL; m++){printf("\t%.1Lf", trace_moments[t][m]);} printf("\n");} INIT_MOMENTS_TRACE}
//SA:The "Ben Index" is used to map a trace value to a bin - MAXSIZE is currently set to a large number which is something like the expected max value (should maybe change)
#define BIN_TRACE(t,vol, cp,p) {__trace_histogram[(int)floor(HIST_BINS*log(t)/log(vol))]++;}
//flush out the histogram in a comment
#define WRITE_TRACE_HIST(L) {printf("#TRACE_HIST_L=%i: ",L); int l; for(l=0;l<HIST_BINS*write_hist; l++){printf("%d ", __trace_histogram[l]); } printf("\n");}
/* This is such that numbers don't explode. Else, they go with N^n. Rather cut them into pieces here, and then blow them up later */
// At the end of each run in system size L write the final avalanche size into avalanche_moments.
#define COMPUTE_AVALANCHE_STATS(av,N){int i; long double normalised_avalanche = av / ((double)N); for (i = 0; i <= MAX_MOMENTS*write_avalanches; i++) {avalanche_moments[i] += (long double)pow(normalised_avalanche, i);}}
#define WRITE_AVALANCHES(L) {printf("#AVALANCHES_L=%i: ", L); int i; 	for (i = 1; i <= MAX_MOMENTS*write_avalanches; i++) {printf("\t%.1Lf", avalanche_moments[i-1]); } printf("\n");}
//todo: above we should be passng size

//unsigned x= genrand_int32(); while (x >= RNG_MT_MAX - (RNG_MT_MAX % n)) - removed boundary check - needs changing
#define _RANDOM_INT(r, n) {unsigned x = genrand_int32(); x %= n; r = (int)x;}
#define RANDOM_DOUBLE  genrand_real1()
//#define RANDOM_ORIENTATION (genrand_int32() % (2*D))
#define POP_RANDOM_PARTICLE(p)int r,top; _RANDOM_INT(r, stack.top); p = stack.stk[r]; stack.stk[r] = POP(top); p;
#define EXP_WAIT(n) (1.0/(double)n) * (-log(1-RANDOM_DOUBLE))

#define POP(p)  p=stack.stk[--stack.top]; 
#define PUSH(p) stack.stk[stack.top++] = p;
#define PARTICLE_COUNT stack.top
//todo add a stack with a certain capacity and record trace on it unless we exceed threshold in which case?
//#define UPDATE_ROG(p)  radius_gyr += 1;//distance_from_center(p); //test with += 1 to match with trace
//if (__compute_rog == 1) UPDATE_ROG(p);
//when we are an annhilating mode, we only update when there is no exiting trace
#define ADD(p) {if(TRACE_FLAG!=(lattice[p] & TRACE_FLAG)) {__trace++; lattice[p] |= ADD_FLAG; PUSH(p); log_trace(p); }else{ ___annihilations++; if(__is_annhilating__==0){lattice[p] |= ADD_FLAG; PUSH(p);}}}
//log_trace(p)
#define REMOVE(p) {lattice[p] &= ~CURRENT_FLAG;}
//lattice_remove(p);
#define MOVE(from,to) REMOVE(from); ADD(to);
#define STAY(p) ADD(p)
//#define TRY_ADD_IMMOBILE(p)  if(IMMOBILE_FLAG!=(lattice[p] & IMMOBILE_FLAG))__immobileTrace++ ;lattice[p] |= IMMOBILE_FLAG

/*Function declaration */
void print_write_times(void);
/*globals*/
int __sample_n__ = 931;//541; //-1;
int write_hist = 0, write_lattice = 0, write_image = 0, write_hull = 0, write_avalanches = 0;
int __trace = 0, __immobileTrace = 0, __maxParticles = 0, __verbose = 1, __compute_rog = 1, __quit_at_trace_max__ = 1;
int __graph_type__ = 2; //sparse -1
int __is_annhilating__ = 0;

char *lattice;
SSTACK stack;
//temp - idea here is we can configure the resolution of things like write times and historgrams but then we need to malloc
int BIN_N = BINS;
double radius_gyr = 0; //sa change form log

//replace prefact e.g. 1-3 to determine how many obs to capture
#define NUM_OBS 1
int MAX_MOMENTS_ALL = NUM_OBS * (MAX_MOMENTS + 1);
long double avalanche_moments[NUM_OBS*(MAX_MOMENTS + 1)]; // Here, I'm assuming that there will be 10 different system sizes to be checked!
long double trace_moments[BINS][NUM_OBS*(MAX_MOMENTS + 1)];//NB: adding 1 to number of moments to capture 0th moment!! ALL FOREACH shoud terminate at <=MAX_MOMENTS
long long alive_bit[BINS];

long double write_times[BINS]; // bw: changed to double
int __trace_histogram[HIST_BINS];
int CHUNK_SIZE = 0;

//this can be removed - is function of something else
int sizeCount = 0;
int MAX_L = 0;

long double ___events = 0;
long double ___births = 0;
long double ___deaths = 0;
long double ___annihilations = 0;
long double ___deaths_by_dragons = 0;

void test_brs();
void test_brs() {
	init_lattice(0, 127, 2);
	int cen1 = get_center(127, 2);
	int pos1 = 0;
	int next1 = 0;
	ADD(cen1);

	printf("Added %d - Pop is %d", cen1, PARTICLE_COUNT);

	cen1++;
	ADD(cen1);

	printf(" Added %d - Pop is %d", cen1, PARTICLE_COUNT);

	int cold = cen1 - 1;
	printf(" Moving  %d -> %d", cen1, cold);

	POP_RANDOM_PARTICLE(cen1);
	MOVE(cen1, cold);

	printf(" Pop is %d", PARTICLE_COUNT);

}
void spit_out_image(int L, int D, long double t, int n, int new_lines) {
	int a = 0;
	printf("#LATTICE(%.3Lf,%d, %d): ", t, L, n);
	for (a = 0; a < pow(L, D); a++) {
		if (a % L == 0 && new_lines == TRUE)
			printf("\n");
		if (TRACE_FLAG == (lattice[a] & TRACE_FLAG)) {
			if (CURRENT_FLAG == (lattice[a] & CURRENT_FLAG))//if it is there now
				printf("%03d ", -1 * a);
			else
				printf("%03d ", a);
		}
		else {
			//this is for testing and viz
			//printf("    ");
		}
	}
	printf("\n");
}

void spit_out_image_light(int L, int D, long double t, int n, int new_lines);
void spit_out_image_light(int L, int D, long double t, int n, int new_lines) {
	int a = 0;
	FILE *g = NULL;
	fopen_s(&g, "lattice_sp.out", "a");
	for (a = 0; a < L*L; a++) {
		//if (HOLE == (lattice[a] & HOLE)) {
			//fprintf(g, "%.3Lf\t%d\n", t, lattice[a]);
		//}

		//if (TRACE_FLAG == (lattice[a] & TRACE_FLAG)) {
		//	if (CURRENT_FLAG == (lattice[a] & CURRENT_FLAG))//if it is there now
		//		fprintf(g, "%.3Lf\t%d\n", t, -1*a);
		//	else
		//		fprintf(g, "%.3Lf\t%d\n", t, a);
		//}

		//else {
		//	//this is for testing and viz
		//	//printf("    ");
		//}
	}
	fclose(g);
}

int main(int argc, char *argv[])
{
	// Set buffer to zero so that process has to write into stdout immediately
#if defined linux || defined __APPLE__
	setlinebuf(stdout);
#endif

	// generate_binary_tree(600001);
	// generate_binary_tree(200001);
	// generate_binary_tree(300001);
	// generate_binary_tree(400001);
	// generate_binary_tree(800001);

	// return 0;

	double sigma, h;

	int BCs, D, l, Ln, L, C, N, seed, min_l, max_l;
	if (parse_args(&BCs, &D, &L, &C, &Ln, &N, &seed, &min_l, &max_l, &h, argc, argv) == TRUE) {
		//init_genrand(seed);
		CHUNK_SIZE = (int)(N / C);
		sigma = (1 - h) / 2;
		if (__graph_type__ == 2) D = 1; //insist
		PRINT_PARAMS
			printf("#Version: ");
		printf("#");
#if defined linux || defined __APPLE__
		PRINT_VERSION
#endif

			INIT_WRITE_TIMES
			WRITE_HEADER;//SA
			//printf("\nN\tt\th\tsig\ttrace\tevents\tbirths\tdeaths\tann\tdbg\n");


			//printf("n\ttm\tsig\ttrace\tevents\tbirths\tdeaths\tdbd\n");
			// D = 3; //max alloc requirement
			// h = 0.2;
			// L=63;

			// MAX_L = L = 10024;
			//printf("\n#allocating %d", MAX_L);

			//this is important setup - i allocate and generate the maxium graph as an allocation 
			//must be 1D.
			//after, re-alloc will not attempt to enlarge the remory but instead we will write over the memory with increasingly large graphs, confined.
		if (__graph_type__ == 2)
			D = 1;//

		allocate_lattice(&lattice, MAX_L, D, __graph_type__);
		spit_out_image_light(243, 2, 0, 0, 0);

		if (__graph_type__ == 2) generate_graph(MAX_L);

		//run_for_realisations(1, 1001, 1, 0.1, sigma, 0);
		//printf("done");
		//return 0;

		//FOR BINARY TREE I AM MESSING WITH DIMS


				/*
		D = 2;
		__is_annhilating__ = 1;
		allocate_lattice(&lattice, MAX_L, D, __graph_type__);
		if(__graph_type__==2)generate_graph(MAX_L);

		// double sigmas [] = {  0.75, 0.65, 0.75, 0.85, 0.851, 0.852, 853};
		int ds [] = {D};
		int a1 =0, a2=0, a3=1;

		double epsilon = 0.0;
		double start = 0.6;
		double increment = 0.1;
		for (a1 = 0; a1 < 1; a1++){
			for (a2 =0; a2 < 5; a2++){
				//for (a3=3;a3 <= 10; a3++)
				{
					h= 1.0;//0.9 + 0.01 * a3;
					run_for_realisations(N, L, ds[a2], h, start, BCs);
				}
				start+= increment;
			}
		}

		return 0;
		*/
		//**************************************************

	//we either run a single L if specified, otherwise go through the motions
		if (L > 0) { run_for_realisations(N, L, D, h, sigma, BCs); }
		else {
			for (l = min_l; l <= max_l; l++) {
				L = pow(2, l) - 1;
				run_for_realisations(N, L, D, h, sigma, BCs);
			}
		}
	}
}

inline void run_for_realisations(int N, int L, int D, double h, double sigma, int bcs) {
	printf("\n# Running for L = %i\n", L);
	INIT_MOMENTS_TRACE;
	INIT_TRACE_HISTOGRAM(L);
	INIT_AVALACNHES

		double __avalanche = 0; // where should this be? changing SA for per run
	int _n = 0, chunk = 0;
	int ___trace_max = 0;
	int vol = pow(L, D);


	//generate once here and optionally at the end of each chunk
	if (__graph_type__ == 2)
		generate_graph(L);


	//adding test particle
	//
	___events = ___births = ___deaths = ___annihilations = ___deaths_by_dragons = 0;

	//printf("#Generating graph and exiting.... %d", vol);
	//return;

	//test walking the graph - ergodicity test
	int cen1 = get_center(L, D);
	int pos1 = cen1;
	int next1 = 0;

	//For rendering only, we will start at the beginning because easier to test
	//ADD(0);
	ADD(cen1);


	// _update_trace(9021, L, D, 2);
	// _update_trace(8769, L, D, 2);
	//return;
	/*ERGODICITY
	int ch = 0;
	for (ch=0; ch< 10; ch++){
	init_lattice(bcs, L, D);
	stack.top = 0; __trace = 0; __immobileTrace = 0, __maxParticles = 0; __avalanche = 0.0;
	radius_gyr = 0;
	ADD(cen1);
	for (_n=0; _n < 100000; _n++){
		//printf("popping from %d\n", PARTICLE_COUNT);
		POP_RANDOM_PARTICLE(pos1);
		next1 = -1;
		while(next1 == -1){ //reflecting for ergodic tests
			//printf("pos now %d", pos1);
			next1 = diffuse(pos1);
		}
		//printf("#Moving %d -> %d ... trace is %d and %d\n", pos1, next1, __trace, get_cell_trace());
		if (next1 == -1){
			printf("sink\n");
			return;
		}
		//printf(" moving between %d %d \n", pos1, next1);
		MOVE(pos1,next1);
		//ADD(next1);
		//printf("done");
		printf("%d %d %d %d %d\n", ch, _n, __trace, pos1, next1);
		if (__trace == L)break;
	}
	}

	return;
	*/
	//global max time
	long double max_time = 0.;
	for (_n = 0; _n < N; _n++) {
		//printf("%d\n",_n);
		//SA unless we want to total per run but given we store everything we can sum later for now
		___events = ___births = ___deaths = ___annihilations = ___deaths_by_dragons = 0;
		//if (_n == 15034 +1)return;

		long double cumpartcount = 0;

		double rd = 0.0;
		long double time = 0.0, last_time = 0.0;
		int wait_for_for_a_big_one = 1000;
		int write_time_index = 0, pos = 0, ro = 0, next = 0;
		int lattice_flush_index = 0;

		//printf("#INIT...");
		init_lattice(bcs, L, D);
		// printf("\n#Adding test particle");
		// ADD(0);
		// printf("..DONE");

		//return;
		//reset();
		//init_cells();
		stack.top = 0; __trace = 0; __immobileTrace = 0, __maxParticles = 0; __avalanche = 0.0;
		radius_gyr = 0;

		int cen = get_center(L, D);
		//randomize initial condition not on the boundary
		///////////////////////////FOR GRAPH TYPE 2////////////////////////////////
		if (__graph_type__ == 2) {

			_RANDOM_INT(cen, L - 2) + 1;
			//for render may want to start at root
			//cen = 0;//SA TEST - always start at what i think is the root for RT

		}


		ADD(cen);
		__maxParticles = 1;
		// kick out this do while loop

		//printf("\n#SAMPLING....");
		do {
			cumpartcount += PARTICLE_COUNT;
			time += (EXP_WAIT(PARTICLE_COUNT));
			__avalanche += (time - last_time)/*dt*/ * ((double)PARTICLE_COUNT);

			while ((write_times[write_time_index] < time) && (write_time_index <= BINS - 1)) {
				WRITE_MOMENTS_TRACE(__trace, (radius_gyr / __trace), PARTICLE_COUNT, write_time_index, 1); write_time_index++;
				if (__trace > 20) {
					//write the hull IF we have set a particular global option
					if ((write_time_index % 10 == 0) && (__trace > 100)
						&& (
							/*_n == 4122 |
							_n == 8201 |
							_n == 8696 |
							_n == 9323 |
							_n == 11228 |
							_n == 12833 |
							_n == 17576 |
							_n == 19606 |
							_n == 34638 |
							_n == 49745 |*/
							_n == 15034)
						)
					{
						//spit_out_image_light(L, D, time, _n, FALSE);
						//spit_out_image(L, D, time, _n, FALSE);
						//right_turning_walk(L,_n,time);
					}
				}
			}
			if (write_times[BINS - 1] < time) { printf("##Recorded excessive time"); break; }//todo see why this sometimes happens

			POP_RANDOM_PARTICLE(pos);
			rd = RANDOM_DOUBLE;//to choose sub-process...

			___events++;

			// int K = 2, k=0;
			// if (rd <= sigma) {//branch locally but may then move/stay
			// 	for(k=0; k < K; k++){
			// 		//depending on the rate of hopping
			// 		if(RANDOM_DOUBLE <= h){
			// 			next = diffuse(pos);
			// 			while (next == -2) next = diffuse(pos);//reflected from holes - this could *now* be abstracted inside lattice
			// 			if (next == -1) {REMOVE(pos); ___deaths_by_dragons++; continue; }
			// 			else{ ___births++; ADD(next); /*this is a try_move. if annihil: inc counter and dont stack*/ }
			// 		}else{//effectively branch fully locally
			// 			PUSH(pos); //add particle there but dont mess with the lattice
			// 			___births++;//local birth
			// 		}
			// 	}
			// }

			// else {//extinction
			// 	//printf("\n DIE.");
			// 	___deaths++; //- note we infer this type of feath - for us death is by annihilation
			// 	REMOVE(pos);
			// 	continue;
			// }

			if (rd <= sigma) {//branch locally (and stay)
				STAY(pos);
				ADD(pos);
			}

			else if (rd <= sigma + h) { //hop
				next = diffuse(pos);

				while (next == -2)
					next = diffuse(pos);//reflected from holes
				if (next == -1) {//-1 illegal - dead for open boundary - not put back on stack, reset flag on lattice
					REMOVE(pos);
					continue;
				}
				MOVE(pos, next);
			}
			else {//extinction
				REMOVE(pos);
				continue;
			}


			if (PARTICLE_COUNT > (10 * __maxParticles)) {
				__maxParticles = PARTICLE_COUNT;
				//printf("\n#%d - (%d) of (%d) LARGEST PART SO FAR!#############################\n", _n, __maxParticles, vol);
			}

			if (time > (100 * max_time)) {
				max_time = time;
				printf("\n#%d - (%Lf) of (%d) LARGEST TIME SO FAR!#############################\n", _n, time, vol);
			}

			last_time = time;
			if (__quit_at_trace_max__ == 1 && __trace >= vol) {
				break;//this is a safety if trace is the only interesting observable - this prevents particle purgatory 
			}
		} while (PARTICLE_COUNT);

		//if(N==1)return;//special test case

		//custom stat
		//printf("%d\t%Lf\t%f\t%f\t%d\t", _n, last_time,h, sigma,__trace);
		//printf("%Lf\t%Lf\t%Lf\t%Lf\t%Lf\t\n", ___events , ___births, ___deaths, ___annihilations, ___deaths_by_dragons);

		//add on number of births, deaths, deaths by dragon, events

		//TODO: histogram of survival: write ID | LAST_TIME | SIGMA | H 
		//a different thing should is alive at hist: we have that on average but we dont have a N*T life matrix

		//finalise
		radius_gyr /= __trace;

		//printf("%d\n", __trace);

		// if (__trace > ___trace_max){
		// 	___trace_max = __trace;
		// 	printf("\n#%d - (%d,%d) of (%d) LARGEST TRACE SO FAR!#############################\n", _n, __trace,get_cell_trace(), vol);
		// }

		//if N=5732 and L = 8191 go verbose
		// if (_n == 5732 && vol == 8191){
		// 	__verbose = 1;
		// }
		//printf("%f, %d\n",radius_gyr, __trace );
		
			//pad gaps in time for consistent data shape
			// comment that out
		
		while (write_time_index < BINS) { WRITE_MOMENTS_TRACE(__trace, radius_gyr, PARTICLE_COUNT, write_time_index, 0); write_time_index++; }//SA
		BIN_TRACE(__trace, vol, radius_gyr, PARTICLE_COUNT);
		//see ALSO header
		if ((_n + 1) % CHUNK_SIZE == 0 || (_n + 1) == N) { chunk++; printf("#commit trace\n");  COMMIT_MOMENTS_TRACE(chunk, L, h, sigma, D, bcs); if (__graph_type__ == 2) generate_graph(L); } //free graph at chunks -convention - forces rebuild
		//COMPUTE_AVALANCHE_STATS(__avalanche, N) //SA

		printf("#(HIST)#%d\t%d\t%Lf\t%Lf\n", __trace, __maxParticles, cumpartcount, last_time);
		//WRITE_TRACE_SIZE(_n, L, last_time)
	}
	//WRITE_TRACE_HIST(L);//SA
	// WRITE_AVALANCHES(L);//SA
	printf("#okely dokely!");//look for this line int stats out//SA
}


int parse_args(int *bcs, int *D, int *L, int *C, int *Ln, int *N, int *seed, int *min_l, int *max_l, double *h, int argc, char *argv[]) {
	// Define default parameters 5000500
	*seed = 5; *N = 10000; *D = 2; *Ln = -1; *bcs = 0; *L = -1; *C = 100; *h = 0.1;
	*min_l = 16; *max_l =16;

	//test///////////////
	//*L = 243;//binary tree should be odd
	//*h = 0.5;
	//*L = 63730;FB
	//*L = 2559;// 2361;

	//__graph_type__ = 0; //change dim and grap type - if graph type is 2, dim must be 1
	/////////////////////
	// Check if Windows. Then get parameters as "seed, N, Ln (lattice length), D, bcs"
#ifdef _WIN64
	int temp = time(NULL);
	if (argc > 1 && atoi(argv[1]) != -1)
		temp = atoi(argv[1]);
	if (argc == 3) {
		*seed = temp; *N = atoi(argv[2]);
	}

	else if (argc == 6) {
		*seed = temp; *N = atoi(argv[2]); *Ln = atoi(argv[3]); *D = atoi(argv[4]); *bcs = atoi(argv[5]);
	}
	else {
		printf("#Using default arguments - pass either seed and N or all args (seed,N, Ln, D, Bcs)\n");
	}
#endif

	// Check for linux, if so use getopt
#if defined linux || defined __APPLE__
//	int ch;
//	while ((ch = getopt(argc, argv, "S:N:L:D:B:h")) != -1)  // bw change 'R'
//		switch (ch) {
//		case 'S':// seed
//			*seed = atof(optarg);
//			break;
//		case 'N'://realisations
//			*N = atoi(optarg);
//			break;
//		case 'L':
//			*Ln = atoi(optarg);
//			break;
//		case 'D'://dimensions
//			*D = atoi(optarg);
//			break;
//		case 'B':
//			*bcs = atoi(optarg);
//			break;
//		case 'h':
//			printhelp();
//			_exit(0);
//		default:
//			break;
//		}

	int             c;
	int option_index = 0;
	const char    * short_opt = "C:N:L:D:h:";
	struct option   long_opt[] =
	{
	   {"help",          no_argument,       NULL, 0},
	   {"verbose",          no_argument,       NULL, 0 },
	   {"seed",          required_argument, NULL, 0},
	   {"D",          required_argument, NULL, 0},
	   {"L",          required_argument, NULL, 0},
	   {"graph",          required_argument, NULL, 0 },
	   {"specimen",          required_argument, NULL, 0 },
	   {"Ln",          required_argument, NULL, 0},
	   {"BCs",          required_argument, NULL, 0},
	   {"minLn",          required_argument, NULL, 0},
	   {"maxLn",         required_argument, NULL, 0},
	   {"wROG",          no_argument, NULL, 0 },
	   {"wHist",          no_argument, NULL, 0},
	   {"wAval",          no_argument, NULL, 0},
	   {"wImage",          no_argument, NULL, 0},
	   {"wHull",          no_argument, NULL, 0 },
	   {"write_times",          no_argument, NULL, 0},
	   {"version", 			no_argument,NULL, 0},
	   {NULL,    0,                 NULL, 0  }
	};

	while ((c = getopt_long(argc, argv, short_opt, long_opt, &option_index)) != -1)
	{
		switch (c)
		{
		case -1:       /* no more arguments */
		case 0:
			if (long_opt[option_index].flag != 0)
				break;
			if (strcmp(long_opt[option_index].name, "seed") == 0) { *seed = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "wROG") == 0) { __compute_rog = 1; }
			if (strcmp(long_opt[option_index].name, "wHist") == 0) { write_hist = 1; }
			if (strcmp(long_opt[option_index].name, "wAval") == 0) { write_avalanches = 1; }
			if (strcmp(long_opt[option_index].name, "wImage") == 0) { write_image = 1; }
			if (strcmp(long_opt[option_index].name, "wHull") == 0) { write_hull = 1; }
			if (strcmp(long_opt[option_index].name, "maxLn") == 0) { *max_l = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "minLn") == 0) { *min_l = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "L") == 0) { *L = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "graph") == 0) { __graph_type__ = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "specimen") == 0) { __sample_n__ = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "D") == 0) { *D = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "Ln") == 0) { *Ln = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "BCs") == 0) { *bcs = atoi(optarg); }
			if (strcmp(long_opt[option_index].name, "write_times") == 0) { print_write_times(); }
			if (strcmp(long_opt[option_index].name, "verbose") == 0) { __verbose = 1; }
			if (strcmp(long_opt[option_index].name, "version") == 0) { PRINT_VERSION return 0; }
			if (strcmp(long_opt[option_index].name, "help") == 0) {

				printf("__________                             .__    .__                \n");
				printf("\\______   \\____________    ____   ____ |  |__ |__| ____    ____   \n");
				printf(" |    |  _/\\_  __ \\__  \\  /    \\_/ ___\\|  |  \\|  |/    \\  / ___\\    \n");
				printf(" |    |   \\ |  | \\// __ \\|   |  \\  \\___|   Y  \\  |   |  \\/ \\/_\\/>   \n");
				printf(" |______  / |__|  (____  /___|  /\\___  >___|  /__|___|  /\\___  /     \n");
				printf("        \\/             \\/     \\/     \\/     \\/        \\//_____/      \n");
				printf("           __      __.__                                         \n");
				printf("          /  \\    /  \\__| ____   ____   ___________                 \n");
				printf("         \\   \\/\\/   /  |/ __ \\ /    \\_/ __ \\_  __ \\               \n");
				printf("           \\        /|  \\  ___/|   |  \\  ___/|  | \\/             \n");
				printf("            \\__/\\  / |__|\\___  >___|  /\\___  >__|              \n");
				printf("                 \\/          \\/     \\/     \\/          \n");
				printf("\n");

				printf("Simulation of a Branching Wiener Sausage at Critical Point \n");
				printf("Usage: %s [OPTIONS]\n", argv[0]);
				printf("  -h                    set hopping rate\n");
				printf("  -C                    choose the numer of chunks to used. Chunksize will be N/C\n");
				printf("  -N                    choose the numer of realisations to do\n");
				printf("  -L                    choose a specific system size. This takes precedence over other system size parameters\n");
				printf("  -D                    choose the dimension e.g. 1-5. Default is 2\n");

				printf("\n");
				printf("  --seed                choose a sensible seed\n");
				printf("  --Ln                  choose a system size index for L=2^(Ln)-1. This takes precedence over minLn and maxLn\n");
				printf("  --minLn               choose a lower bound system size index for minLn=2^(minLn)-1\n");
				printf("  --maxLn               choose an upper bound system size index for maxLn=2^(maxLn)-1\n");
				printf("  --BCs                 set bitmask to choose specific dimension boundaries. Defaults to 0 i.e. all open boundaries\n");
				//printf("  --chunks              choose chunk size for flush stats (default 1000) \n");

				printf("\n");
				printf("  --wROG                set this flag (without argument) to write radius of gyration data\n");
				printf("  --wHist               set this flag (without argument) to write trace histrogram data\n");
				printf("  --wAval               set this flag (without argument) to write avalanche time integral data\n");
				printf("  --wImage              set this flag (without argument) to write lattice data. WARNING: This can be large!\n");
				printf("  --wHull               set this flag (without argument) to write hull.\n");
				printf("  --verbose             set this flag (without argument) to write additional info e.g (max trace per sample path,..).\n");
				printf("  --specimen            specify sample path N for extended data such as the hull and image. This N value can be determined from sample runs when running with--verbose. Look for #MAX_TRACE(,,index). \n");
				printf("  --graph               choose a particular type of graph. Default to regular lattice or choose sierpinski (1), small-world network(2).\n");

				printf("\n");
				printf("  --write_times         print out the write times used by the program here and now...\n");
				printf("  --version             print out current (git-)version of code\n");
				printf("  --help                print help and exit\n");
				printf("\n");
				return 0;
			}
			break;

		case 'h':
			*h = atof(optarg);
			break;

		case 'L':
			*L = atoi(optarg);
			break;

		case 'N':
			*N = atoi(optarg);
			break;

		case 'C':
			*C = atoi(optarg);
			break;

		case 'D':
			*D = atoi(optarg);
			break;

		case ':':
		case '?':
			fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
			return(-2);

		default:
			fprintf(stderr, "%s: invalid option -- %c\n", argv[0], c);
			fprintf(stderr, "Try `%s --help' for more information.\n", argv[0]);
			return(-2);
		};
	}
#endif

	//if Ln is specified, use that as the log range - note, if L also specified, L will take precedence in main function
	if (*Ln > 0 && *L < 0) {
		*min_l = *Ln;
		*max_l = *Ln;
	}
	if (*L > 0)MAX_L = *L;
	else MAX_L = pow(2, *max_l) - 1;

	return TRUE;
}

void printhelp() {
	printf("Simulation of a Branching Wiener Sausage at Critical Point \n\
	Correct syntax (UNIX-only): ./bws -S (Seed) -N (#Realisations) -L (L <= 2^input) -D (dim) -B (Boundary-conditions as bitmask) -h [help]\n ");

	//printf("#AVALCOMMENT below are avalanche moment <s^n>. However, read <s^n> = <s^0>^n*<s^n>! Moments have to be inflated by N^n.\n");
}

void print_write_times() {
	INIT_WRITE_TIMES
		int i = 0;
	for (i = 0; i < BINS; i++) {
		printf("%0.1Lf\n", write_times[i]);
	}
}
