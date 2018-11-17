//#basic bosonic L^D lattice
#include "lattice_core.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NODE_TYPE 0
#define EDGE_TYPE 1

#define RANDU  rand() / (RAND_MAX + 1.) 

int maxl = 0;
int _D = 0, _L = 0, _volume = 0, _BCS = 0, _center = 0;
int prev_n = -1;
int lattice_type = 0;
int sparse_mode = 0;

int cache[CACHE_CAPACITY];
int capacity_exceeded = 0; // to disable sttucture
int cache_size = 0;
char* lattice;
//graph
static unsigned long last_malloc_sizes[2];
static unsigned long graph_allocated_bytes = 0;
static unsigned int * degree_dist;
static unsigned int * offsets;
static unsigned int * adjacency;
static long long non_one_degree_nodes[10000000];//for test fixed size - i will allocate this as L-ks later
static unsigned int degree_sum;

int lattice_actions[20];
int dimension_increment_maps[20];
int wrap_increment_maps[20];

/*CELL BITS
*/
#define cell_stack_capacity 10
#define tilestack_size 2000000
cell* cells[cell_stack_capacity];
trace128 tiles[tilestack_size];

static int ___trace = 0;
static int __flushes = 0;
static int __max_ll_length = 1;
static double __avg_ll_length = 1;
static int __complete_cells = 0;
static int stack_top =0;
static int counter = 0;
//static trace128 done;
int iknowd = -1;
///////

#define TRACE_FLAG_LAT (1<<2) //redefined from bws.c

#define FALSE 0
#define TRUE 1
#define BOUNDARY_FLAG(d) (1<<d)
#define IS_BOUNDARY_CLOSED(d,flag) (BOUNDARY_FLAG(d) == (flag & BOUNDARY_FLAG(d)))
#define IS_ALLOWED_TRANSITION(i,inext,dim_scale) (PARTITION_OF(i, dim_scale*_L) == PARTITION_OF(inext, dim_scale*_L))
#define PARTITION_OF(i,j) (int)( floor(i/((double)j)))
#define IS_ACCESSIBLE(p) (HOLE != (lattice[p] & HOLE))
//#define RANDOM_EDGE(loc) (genrand_int32() % (2*_D))

inline static void print_adj(int L){
     //print the non one degree helper
    int _of=0;
    int pairings = 0;
    int alt =0;
    printf("\n");
    for(alt =0; alt < L; alt++){
        int alt_degree = degree_dist[alt];
        int alt_offset = offsets[alt];
        printf("%d(%d, %d):", alt, alt_offset, alt_degree);
        for(_of=alt_offset; _of<alt_offset+alt_degree;_of++ ){
            printf(" %d ", adjacency[_of]);
        }
        printf("\n");
    }
}

inline static int _power_sample(double xmax, double degree) {
    double xmin = 1.; //weight a little away from one and below L in this distributiocle
    double pmin = pow(xmin, degree);
    double pmax = pow(xmax, degree);
    double v = pmin + (pmax - pmin) * RANDU;
    return (int)pow(v, 1.0/(degree));
}

inline static void _allocate(unsigned int ** arr, int size, int arr_offset){
    unsigned long req = size*sizeof(int);
    int last_size = last_malloc_sizes[arr_offset];
    if (last_size == 0){
		//printf("#First time allocation...");
        *arr = (unsigned int*)malloc(req);
		if(*arr!= NULL){ graph_allocated_bytes += req; }//otherwise complain like hell
		else{printf("#FATAL: could not allocate structure...");}
    }else if (size > last_size){
		printf("#Realloc %d->%d", last_size, size);
        *arr  = (unsigned int*)realloc(*arr , req);
		if(*arr != NULL){graph_allocated_bytes += (req-last_size*sizeof(unsigned int));}//otherwise complain
		else{printf("#FATAL: could not allocate structure...");}
    }
}

inline static int power_law_dist(int L, int stats[1]){
    _allocate(&degree_dist, L, NODE_TYPE);
    _allocate(&offsets, L, NODE_TYPE);
    last_malloc_sizes[0] = L;
    int i = 0, val =0, sum =0, count_ones=0;
    int last_non_one = 0;
    for(i =0; i < L; i ++){ 
        //SCALE TESTING -------------------- THIS NEEDS TO BE CONDITIONED
        val =  _power_sample(L, -1.35);
        sum += val;
        if( val== 1)count_ones++;
        else{last_non_one=i;}
        degree_dist[i] = val;
    } 
    //ensure sum is of even degree by random correction
    if (sum % 2 != 0){  degree_dist[last_non_one]++; sum++; }
    //printf("distribution generated\n");
    stats[0] = count_ones;
    return sum;
}

inline static void init_graph(unsigned int L,int type){
	printf("\n#Allocating PL..");
    int stats[1];
    int i=0,num_ks = stats[0], rest=0;
    degree_sum = power_law_dist(L,stats);
    num_ks = stats[0];
    rest = L - num_ks;

	int chain = 2*(L-num_ks-1);
	printf("\n#Allocating ADJ..");
	//printf("\nChain will cost %d leaving %d for arbitrary pairings inluding %d for 1 degrees", chain, degree_sum-chain, 2*num_ks);
	_allocate(&adjacency, degree_sum, EDGE_TYPE);
	last_malloc_sizes[1] = degree_sum;
	for (i = 0; i < degree_sum; i++) { adjacency[i] = -1; }
	printf("\n#Allocated adjacency");
	configure_edges(L, num_ks);
}

int build_edge_distribution(int n) {
	//actually these could be allocated to a fixed size - and maybe we could get away with allocating 2n for the edges
	if (prev_n == -1 && !degree_dist) {
		printf("# allocating degree first time");
		degree_dist = malloc(maxl * sizeof(double));
		offsets = malloc(maxl * sizeof(int));
	}
	// else if (prev_n < n){
	// 	printf("# re-allocating degree");
	// 	degree_dist = realloc(degree_dist, n*sizeof(double));
	// 	adj_offset = realloc(adj_offset, n*sizeof(double));
	// }

	if (!degree_dist) {
		printf("\nFatal error: an attempt to malloc degree failed!");
	}

	int i = 0;
	if (lattice_type == 2) {//powerlaw
		for (i = 0; i < n; i++) { //MAKE POWER LAW ORDER N FOR TEST
			degree_dist[i] = _power_sample(n/*scale parameter describes say a tail*/, -2);
		}
	}

	int sum = 0;

	for (i = 0; i < n; i++) {
		sum += degree_dist[i];
	}

	if (prev_n == -1 && !adjacency) {
		printf("# allocating adj first time"); //i cannot see why the re-alloc is not working sometimes.
		//#this graph edge factor can start at 2 and we can re-attempt to re-alloc on an incremented factor if we learnsomethig
		adjacency = malloc(/*sum*/(2 * maxl) * sizeof(int));
	}
	// else if (prev_n < n){
	// 	printf("# re-allocating adj for adjacency of length %d", sum);
	// 	adjacency = realloc(adjacency,sum*sizeof(int));
	// }

	if (!adjacency) {
		printf("\nFatal error: an attempt to malloc adj failed!");
	}


	for (i = 0; i < sum; i++) { adjacency[i] = -1; }

	//track allocated 
	//prev_n = n;

	return sum;
}

//trees normally call for recursion to implement statefullness on a stack
//but not best for optimized code because you are not in control of the stack 
//if you can use problem specific information to manage the stack (some state will be required to expand the tree postfix, prefix or infix) then do so
//here we stack "branched parents" so that on visiting downstream nodes, they can pop the parents from the stack 
//-we do not care what the order is but we must consistently pair - note there can only every be a small number of these as any child must choose an available parent
void generate_binary_tree(long long L,int allocate) {
	//preamble
	long long n = 0;//the root
	long long i = 0;
	long long _n_alive_stack = 0;
	unsigned long long degree_sum = 2 * (L - 1);
	printf("\n#generating binary tree L= %d",L);


	if (allocate == 1){
		printf("re-allocating");
		//printf("generating node elements....\n");
		_allocate(&degree_dist, L, NODE_TYPE);
		_allocate(&offsets, L, NODE_TYPE);
		if (L > last_malloc_sizes[0]) last_malloc_sizes[0] = L;
		printf("\n#adjacency size %llu for lattice of size %d\n", degree_sum, L);
		_allocate(&adjacency, degree_sum, EDGE_TYPE);
		if (degree_sum > last_malloc_sizes[1]) last_malloc_sizes[1] = degree_sum;
		////////////////////////////////////////
	}
	//temp - allocate once then come in here and just return somewhere
	//else {return;}

	for (i = 0; i < L; i++) { offsets[i] = -1; }
	for (i = 0; i < degree_sum; i++) { adjacency[i] = -1; }
	
	
	int e_vals[11];

	//algo
	long long branchings = 0;
	int parent_id=0, parity_param=0, par_adj=0;
	int offset = -1;//start at -1 to deal with root exception
	int E = 1;//can never return to 0
	int terminate = ((L - 1) / 2);
	int qFront = 0;
	int qRear = -2;//subtle: we want it to be -1 except for the root, we do not want to increment at end of loop (exception)
	i = 0;//visited node id init
	
	while (i<L) {
		//e_vals[i] = E;
		
		offsets[i] = offset;
		parity_param = ((i+1) % 2) + 1;
		if (i == 0) offsets[i] = 0;//root exception
		if (qRear >= -1) {//key step in terms of pairing logic here
			//peak the last parent on the FIFO queue (queue works for breadth-first visitation for example)
			parent_id = non_one_degree_nodes[qRear+1];
			//i will be pinned to my parent from their offset at position 1 or 2, (parity depending) -> [[0], [1], [2]]
			par_adj = offsets[parent_id] + parity_param;
			if (parent_id == 0)par_adj--;//root exception
		}
		//branching process
		if (genrand_real1() > 0.5 || i ==0/*force parent branch*/) {
			if (branchings >= terminate)continue;// we continue the other path on clean up - we have branched as much as we 
			//printf("(%d, %d)", 1, E);
			non_one_degree_nodes[qFront++] = i;//ENQUEUE next visited child will use this
			non_one_degree_nodes[qFront++] = i;//ENQUEUE next child visited after will use this
			branchings++;
			E++;
			degree_dist[i] = 3; 
			if (i == 0) {degree_dist[i] = 2; }
			else { adjacency[par_adj] = i; adjacency[offset] = parent_id; }//all important stub pairing
			//printf("1");
			offset += 3;
		}
		else {
			/*assuming you are not the very last, second check*/
			if ( ((E -1) == 0) && (i+1 != L) ) continue; //violation, something needs to branch
			//printf("(%d)", 0);
			E--;
			degree_dist[i] = 1; //I (terminal) only point to my parent, contributing 1 half-edge
			adjacency[par_adj] = i; adjacency[offset] = parent_id; //all important stub pairing
			//printf("0");
			offset += 1;
		}
		//printf("%d..", E);	
		i++;	
		qRear++;//DEQUEUE parent-pointer which we only peaked at before
	}
	
	// FILE *f = NULL;
	// fopen_s(&f, "rt_graph_adj", "w");
	// if (f) {
	// 	for (i = 0; i < offset; i++) {
	// 		fprintf(f, "%d\n", adjacency[i]);
	// 	}
	// 	fclose(f);
	// }

	// FILE *g = NULL;
	// fopen_s(&g, "rt_graph_nodes", "w");
	// if (g) {
	// 	for (i = 0; i < L; i++) {
	// 		fprintf(g, "%d\n", degree_dist[i]);
	// 	}
	// 	fclose(g);
	// }

	printf("\n#*************random tree graph generated********* \n");

	/*
	printf("\nvalence: ");
	for (i = 0; i < L; i++) {
		printf(" %d ", degree_dist[i]);
	}
	
	printf("\noffsets: ");
	for (i = 0; i < L; i++) {
		printf(" %d ", offsets[i]);
	}

	printf("\nadjacency: ");
	for (i = 0; i < offset; i++) {//
		printf(" %d ", adjacency[i]);
	}
	*/
}


void generate_graph_PA(int n){
	//pref attachment

	printf("\n#generating graph %d",n);

	int uedge_max = 5*n;
	int i =0, _j =0;

	//node_t * adj = malloc(n*sizeof(node_t));
	if (!degree_dist){degree_dist = malloc(n*sizeof(int));}
	else{degree_dist = realloc(degree_dist, n*sizeof(int));}

	if (!offsets){offsets = malloc(n*sizeof(int));}
	else{offsets = realloc(offsets, n*sizeof(int));}
	
	int* adjacency_temp = NULL;
	if (!adjacency_temp){adjacency_temp = malloc(uedge_max*sizeof(int));}
	else{adjacency_temp = realloc(adjacency_temp, uedge_max*sizeof(int));}

	for (i=0; i < uedge_max;i++){ adjacency_temp[i] = -1;}
	for (i=0; i < n;i++){  degree_dist[i] = 0; }
	
	int sum = 1;
	int max = 1;
	int uedges = 0;
	int offset = 0;

	//printf("\nlength %d", n);

	for (i =1; i < n; i++){
		//printf("\nadding %d", i);
		int added = 0;
		
		//double rd = genrand_real1(); //<-faster to use just one RV
		for(_j =i; _j >=0  ; _j--){
		//while(added == 0){ //not sure if we want to assume we can only connect to one or we can connect to many
			int j = _j;// % i; //only look up to added and wrap
			double rd = genrand_real1();
			if (i != j){
				double weight = degree_dist[j]/(double)sum; 
				if (sum ==1){weight = 1; sum =0;}//i spoof the start
				//printf(" comp: %d... %f <= %f ", j, rd, weight);
				if (weight == 1 || rd <= weight){
					sum += 2;

					adjacency_temp[uedges] = j;
					adjacency_temp[(uedge_max-1) - uedges] = i;
					
					degree_dist[i]++;
					degree_dist[j]++;
					uedges++;
	
					added = 1; 
				}
			}
			if(_j ==0 && added == 0){_j = i;} //go again
		}
	}
	
	//printf("\n....");
	//initialise the offset of the nodes
	for (i =0; i < n; i++){ offsets[i] = offset; offset += degree_dist[i]; }
	int j =0;
	//create proper adj

	printf("\n#degree sum is %d" , sum);

	if (!adjacency){adjacency = malloc(sum*sizeof(int));}
	else{adjacency = realloc(adjacency, sum*sizeof(int));}

	for (i=0; i < sum;i++){ adjacency[i] = -1;}

	//printf("\nchecking temp edges");
	for (i=0; i < uedges; i++){//double check size
		j = (uedge_max - 1) - i;
		int a = adjacency_temp[i];
		int b = adjacency_temp[j];
		// printf("\n--(%d, %d)", a, b);
		// printf("\nadding update to (%d, %d) ",adj_offset[a],adj_offset[b] );
		//the majic is in our ability to know where to place the numbers
		adjacency[offsets[a]] = b;
		adjacency[offsets[b]] = a;
		//printf("\n linking %d <-> %d @ (%d,%d)",a,b, adj_offset[a],adj_offset[b]);
		//trust this never gets messed up if above is correct
		offsets[a]++;
		offsets[b]++;
	}

	//printf("\n");
	//printf("\n freeing adj temp");
	free(adjacency_temp);
	
	//reset offset
	offset = 0;
	for (i =0; i < n; i++){ offsets[i] = offset; offset += degree_dist[i]; }


	FILE *f = NULL;
	fopen_s(&f, "pa_graph_adj", "w");
	if (f) {
		for (i = 0; i < sum; i++) {
			fprintf(f, "%d\n", adjacency[i]);
		}
		fclose(f);
	}

	FILE *g = NULL;
	fopen_s(&g, "pa_graph_nodes", "w");
	if (g) {
		for (i = 0; i < n; i++) {
			fprintf(g, "%d\n", degree_dist[i]);
		}
		fclose(g);
	}

	printf("\n#PA graph generated \n");

	// printf("\n");
	// for (i=0; i < n;i++){ printf("%d ", degree_dist[i]);}
	// printf("\n");
	// printf("\n");
	// for (i=0; i < sum;i++){ printf("%d ", adjacency[i]);}

	//commented out on win compiler
}


void read_in_network()
{
	FILE* stream = NULL;
	
	int a,b;
	char line[1024];
	int counter = 0;
	fopen_s(&stream, "degrees_n1.dat", "r");
	while (fgets(line, 1024, stream))
	{
		//sscanf_s(&line, "%d,%d", &a, &b);
		counter++;
	}
	fclose(stream);

	_allocate(&degree_dist, counter, NODE_TYPE);
	_allocate(&offsets, counter, NODE_TYPE);
	//malloc degree and adj for counter
	fopen_s(&stream, "degrees_n1.dat", "r");
	counter = 0;
	while (fgets(line, 1024, stream))
	{
		sscanf_s(&line, "%d,%d", &a, &b);

		degree_dist[counter] = a;
		offsets[counter] = b;

		counter++;
	}
	fclose(stream);

	counter = 0;

	fopen_s(&stream, "adj_bi_n1.dat", "r");
	while (fgets(line, 1024, stream))
	{
		//sscanf_s(&line, "%d,%d", &a, &b);
		counter++;
	}
	fclose(stream);
	_allocate(&adjacency, counter, EDGE_TYPE);
	counter = 0;
	//malloc degree and adj for counter
	fopen_s(&stream, "adj_bi_n1.dat", "r");
	while (fgets(line, 1024, stream))
	{
		sscanf_s(&line, "%d,%d", &a, &b);
		adjacency[counter] = b;
		counter++;
	}
	fclose(stream);

}

void generate_graph(int n,int allocate){
	//init_graph(n,0);

	generate_binary_tree(n,allocate);
	//generate_graph_PA(n,allocate);
	//read_in_network();
}

int inline Xc(int i, int d) {
	int _d = 0, x = 0, fact = 0;
	for (_d = _D - 1; _d >= 0;_d--) {
		fact = pow(_L, _d);
		x = (int)floor(i / fact);
		if (_d == d)return x;
		i = i % fact;
	}
}

//this code is very literal and crude at the moment and should be optimized and done more elegantly
int inline get_site_default(int i) {
	if (lattice_type != 1)return 0;
	if (lattice_type == 1) {//sierpinski b=3, m=1
		int b = 3;//move out
		int c = 0, p = 0, p_dash, _p = 0;
		int partitions = log(_L) / log(b);//move out

		for (_p = 0; _p < partitions; _p++)
		{
			int invalid = 1;
			p = pow(b, _p + 1);
			p_dash = pow(b, _p);
			for (c = 0; c < _D; c++) {
				if (((Xc(i, c) % p) / p_dash) != 1) {
					invalid = 0;
					break;
				}
			}
			if (invalid == 1) {//we could not prove it is NOT a hole in D
				//printf("%d,%d\n", Xc(i, 0), Xc(i, 1));
				return 0 | HOLE;
			}
		}
		return 0;
	}
	return 0;
}


int count_full_resets = 0;
int count_cache_resets = 0;

void __init(int L, int D) {
	_L = L; _D = D;
	_volume = (int)pow(_L, D);

	int bisectAt = floor(L / 2);
	int center = 0, i = 0;
	for (i = 0; i < _D; i++) {
		center += (bisectAt * (int)(pow(L, i)));
	}
	_center = center;
	//printf("#CENTER OF LATTICE %d", _center);
}

int allocate_lattice(char **buffer, int L, int D, int type) {
	lattice_type = type;
	maxl = L;
	__init(L, D);
	int choice;
	for (choice = 0; choice < 2 * D; choice++) {
		int dim = (int)floor(choice / 2);
		int sign = 1; if (choice % 2 == 0)sign = -1;
		//these could be computed once and added to a list
		int j = (int)pow(L, dim); //the jump for moving around the lattice structure e.g. 1, L, L^2 etc. see hier. pend.
		dimension_increment_maps[choice] = j;
		lattice_actions[choice] = sign * j;
		wrap_increment_maps[choice] = (-1 * sign) * (j*L - j);
	}

	char* myBuf = (char*)*buffer;
	int size = _volume * sizeof(char), i = 0;

	if (type == -1){
		printf("#running sparse mode...\n");
		sparse_mode = 1;
		//init_cells(D);
	}else{ 
		//printf("\nMallocing %d", size);
	//TDOD(EX1)! FOR EXCESSIVE ALLOCATION - WE DONT DO THIS
		lattice = myBuf = (char*)malloc(size);
		for (i = 0;i < _volume;i++) { lattice[i] = get_site_default(i); }//first time init, we may not always run this
	}
	//int test1 = Xc(12, 0);
	//int test2 = Xc(12, 1);
	//int nacc = IS_ACCESSIBLE(30);
	//int yacc = IS_ACCESSIBLE(29);

	if (*buffer == NULL) { printf("\nERROR: Memory allocation did not complete successfully! Required size %i bytes", size); }
	return 0;
}

void init_lattice(int boundary_conditions, int L, int D)
{
	_BCS = boundary_conditions;
	__init(L, D);
	int i = 0;

	//if in excessive memory mode - we should call super_cells.init() - and this should use the cache
	if (sparse_mode == 1){
		//we do need to clean up the cells - for now call cell init though costly
		printf("#costly call to cells init....\n");
		init_cells(D);
		return;
	}

	//reset logic follows
	capacity_exceeded == -1;
	if (capacity_exceeded != 1) {//turn off this behaviour using cap_exceeded == -1
		for (i = 0;i < cache_size;i++) {
			lattice[cache[i]] = 0;//get_site_default(i); <- if it was visited, it should be zero on reset
		}
		count_cache_resets++;
	}
	else
	{
		for (i = 0;i < _volume;i++) { lattice[i] = get_site_default(i); }//;
		count_full_resets++;
	}

	cache_size = capacity_exceeded = 0;
	//mark bit for boundary - todo
	//check BC for closed boundaries

}

int diffuse(int pos) {
	//if the graph is of a network type for now we assume no boundary - later we could create some infinity nodes
	//so this step is more straightforward - we assume the an adjacency has been allocated and that the lattice stores metadata
	
	//if(pos< 0){ printf("\n###-1 encountered in the structure - this should never happen"); }
	//if(pos>= _L){ 	printf("\n###Excess encountered pos(%d) > L (%d) in structure - this should never happen", pos, _L); }

	if (lattice_type >= 2){ //scale free etc.
		int new_pos= adjacency[  offsets[pos] + (genrand_int32() % degree_dist[pos]) ];
		// if (new_pos == 0 || new_pos == _volume-1){
		// 	return -1;//poor man's boundary conditions
		// }
		
		return new_pos;
	}

	//regular
	int choice = (genrand_int32() % (2*_D));

	int pos_next = pos + lattice_actions[choice];
	int d_scale = dimension_increment_maps[choice];
	if (IS_ALLOWED_TRANSITION(pos, pos_next, d_scale)==0)
		return -1;
	if (IS_BOUNDARY_CLOSED(choice, _BCS) == TRUE) {
		//return pos;//reflecting
		pos_next= pos + wrap_increment_maps[choice];//wrap?
	}
 if (IS_ACCESSIBLE(pos_next)==0) {  return -2;  } //reflect and surface the issue, process can go again
	//if we are wrapping?
	return pos_next;
}

int get_center(int L, int D) {
	__init(L, D);

	if (lattice_type == 1) {
		int random_dir = abs(genrand_int32() % (2 * _D));
		int random_dir_offset = lattice_actions[random_dir];
		int sierpinski_tile_offset = random_dir_offset * (int)ceil(L / (double)(3 * 2));//see notes S3,8 always has central tile of side sqrt(L). So we place ourparticle beside it
		
		//printf("center at %d\n", sierpinski_tile_offset);
		return _center + sierpinski_tile_offset; //must be a random offset from center - not a single one
	}

	return _center;
}

int probe(int p, int orientation, int offsets[]) {
	//consider lattice boundary cond.
	int L = offsets[4];
	if (orientation != -1) {
		//first check that the offset is legal for this position - if not return 0
		int temp = p;
		int offset = offsets[orientation];
		temp += offset;
		if (IS_ALLOWED_TRANSITION(p, temp, L) == FALSE) {//2d only
			return 0;
		}
		p = temp;
	}
	if (TRACE_FLAG_LAT == (lattice[p] & TRACE_FLAG_LAT))
		return 1;
	return 0;
}


void right_turning_walk(int L, int n, long double t) {
	int pos = L / 2, orientation = 0;
	int offsets[8] = { -L,-L + 1,1,L + 1,+L,+L - 1,-1,-L - 1 };
	while (probe(pos, -1, offsets) == FALSE) {
		pos += L;
	} //move down
	int ring_start = pos;
	printf("#HULL(%.3Lf,%d, %d): ", t, L, n);
	do {
		//turn until we find something
		while (probe(pos, orientation, offsets) == FALSE) {
			orientation = (orientation + 1) % 8;
		}
		printf("%d ", pos);
		//move the cursor to this something
		pos += offsets[orientation];
		//re-orient - by looking "backwards" to where you came from + 1 to the right. There should be nothing there
		orientation = (orientation + 5) % 8;
		//terminate when we come back around
	} while (ring_start != pos);
	printf("\n");
}

double distance_from_center(int i) {
	int _d = 0;
	long double sdist = 0;
	for (_d = 0; _d < _D; _d++) {
		sdist += pow(Xc(_center, _d) - Xc(i, _d), 2);
	}
	return sdist;//sqrt(dist);
}

void get_done_bits(int d, trace128 t[1]);

inline int update_trace_bit(int tile_id, int tile_bit);

int init_cells(int D){
	iknowd = D;
    stack_top = counter = 0;
	int i =0;
	
    for (i=0; i < cell_stack_capacity;i++){
        cells[i] = NULL;
    }
    for (i=0; i < tilestack_size;i++){
        tiles[i].elements[0] = tiles[i].elements[1] = 0;
    }
    printf("\ninitialized");
	
	// for(i =0; i < 8; i++){
	// 	update_trace_bit(0, i);
	// }
	// printf(" %d, %llu ", tiles[0].elements[0], tiles[0].elements[0]);

	return 42;
}


void get_done_bits(int d, trace128 t [1]){
	int tile_bit =0;
	for (tile_bit=0; tile_bit < pow(2, d);tile_bit++){
		int temp =0;
		if (tile_bit> 63)
			t[0].elements[1] = 1ULL << (tile_bit % 64); 
		else
			t[0].elements[0] = 1ULL << (tile_bit); 
		}
}

inline int update_trace_bit(int tile_id, int tile_bit){
    //or greater again ..... 
    int temp =0;
    if (tile_bit> 63){
        temp = 1ULL << (tile_bit % 64); 
        if ((tiles[tile_id].elements[1] & temp) == temp ) return 0;
        tiles[tile_id].elements[1] |= temp;    
    }else{
        temp = 1ULL << tile_bit; 
        if ((tiles[tile_id].elements[0] & temp) == temp ) return 0;
        tiles[tile_id].elements[0] |= temp;   
    }
    return 1;
}

int get_cell_trace(){
	return ___trace;
}

int reset(){
	___trace = 0;
}


static int flush(){
	int i =0;
	trace128 done[1];
	get_done_bits(iknowd, done);
	long long a = done[0].elements[0];
	long long b = done[0].elements[1];

	for (i=0; i < tilestack_size;i++){
		//cap the max_times at 63 - but fill from the left
		//printf("In %d, [%d, %d] comp [%d,%d]\n", iknowd, tiles[i].elements[0] ,tiles[i].elements[1],a,b);
		if((tiles[i].elements[0] == a) 
		&&  (tiles[i].elements[1] == b)) {
				printf(" ...........> %d removed\n", i);
				tiles[i] = tiles[stack_top-1];
				stack_top--;
		}
	}
	
	__flushes++;
}


static int update_cell_trace(int index, int offset, int tile_index){
    cell* c = cells[index];
    int trace_change =0;
    int mode = 0;
    int stack_location = 0;
    if (c != NULL){//we found something
       int _n =0;
        for (_n=c->n; _n >0; _n--){
            if (c->div == offset){//seek
                stack_location = c->trace_ref;
                mode = 1;
                //printf("\nfound existing at %d", offset);
                break;
            }
            c = c->next;
        }
        if (mode == 0)c = cells[index]; //unfortunate - if we introduce macros we can clean this buisness up
    }

    //if no exact match, 
    if (mode != 1){   
        cell * temp = (cell*)malloc(sizeof(cell));
        temp->n = 1;
        temp->div=offset;
        if (c != NULL){ 
            temp->n = c->n + 1;
            temp->next = c;
            mode = 2;
            //printf("\nappending new at pos %d", temp->n);
            if (temp->n > __max_ll_length){
                __max_ll_length = temp->n;
                //printf("updating to %d",__max_ll_length);
            }
        }
        else{
            //printf("\nnew obj");
        }
        
        if ((stack_top + 1) < tilestack_size){ 
            stack_top++;
            stack_location = stack_top;
        }
        else{
            printf("#flushing....\n");
			flush();
			return 0;

        }

        c = temp;
        c->trace_ref = stack_location;
        cells[index] = c;
        counter++;
    }

    trace_change = update_trace_bit(stack_location, tile_index); 
    return trace_change;
}

static int i_decomp(long double i, int M, int D, int L, int rating[1]){//
    long double temp = i;
    int _d = 0;
    long sum =0;
    int _rating = 0;
    for (_d = D-1; _d >=0; _d--){
        int rem = (int)floorl((temp/pow(L,_d)));//fllot to integer - it is originla basis so should be integer
        temp -= rem*pow(L,_d);
        sum += (long)(rem/M)*pow(L/M,_d); //base factor is L/M (renorm)
        _rating += ((rem % M)*pow(M,_d));
    }

    rating[0] = _rating;
    return sum;
}

int _update_trace(long double key, int L, int D, int ml){
    int tile_index[1];
    double cell_id =  i_decomp(key,ml, D, L, tile_index);
    long long cell_offset = (long long)(cell_id/cell_stack_capacity);
    long double ck = fmodl(cell_id, (long double)cell_stack_capacity);
    int index = (int)ck;
	//printf(" - %.0Lf, %f, %d, %d (%d)- \n", key, cell_id, index, tile_index[0], stack_top);
	printf("%d\t%d\n", ___trace, stack_top);
    ___trace += update_cell_trace(index, cell_offset, tile_index[0]);
    return 0;//trace_change;
}
//static int _update_trace(long double key, int L, int D, int ml);

void lattice_remove(int p){
	lattice[p] &= ~CURRENT_FLAG;
}

void log_trace(int i) {

	//hint_stack.push(i)
	//if (sparse_mode == 1){  _update_trace(i, _L, _D, 2); }//in sparse mode use tiles, no lattice is allocated
	//else 
	// {
	// 	if(TRACE_FLAG!=(lattice[i] & TRACE_FLAG)) ___trace++; 	lattice[i] |= ADD_FLAG;
	// } 

	if (capacity_exceeded == -1)return;//disabled
	if ((cache_size) == CACHE_CAPACITY) {
		capacity_exceeded = 1;
		return;
	}
	cache[cache_size] = i;
	cache_size += 1;
}


inline static void configure_edges(unsigned int L, int ks){
    //chain loosely and then do the check sum
    //keep a stack of 1s and 2s which we can use to manipulate the graph after

    int tracker = degree_sum;
    int tree_seed =-1;
    int threshold = 0;
    int i =0, j=0, k=0;
    int tree_seed_adj_tail;
    int twos = 0;
    int others = 0;
    int total_tree_stubs = 0;
    int offset = 0;

    int test_array_top=0;

    for(i=0; i < L; i++){
        int my_degree = degree_dist[i];
        offsets[i] = offset;

        if (my_degree==2)twos++;

        //seek tree seed
        if(tree_seed == -1){
            if(my_degree>2){
                tree_seed = i;threshold = my_degree;
                tree_seed_adj_tail = my_degree + offset - 1;
                total_tree_stubs += (my_degree);
                others++;
                //test_array[test_array_top++] = i;
            }
        }
        else if (my_degree > 2){
            //printf("\n(%d) im now the king pin - im going to join the tree and consume edges",i);
            adjacency[tree_seed_adj_tail] = i;
            adjacency[offset] = tree_seed;//proof of being in the tree is the first element for quick lookup
            tracker-=2;
            tree_seed = i;threshold = my_degree;
            tree_seed_adj_tail = my_degree + offset - 1;
            //go back - if there is a free edge take it, if there is not, fold in - make sure the node is not the last seed
            total_tree_stubs += (my_degree-2);
            others++;
            //test_array[test_array_top++] = i;
        }
        offset += my_degree;
    }

    //in this second pass i need to make say, 640 pairings and then i can add thes 1s and 2s anywhere i want - but twos first

    threshold = ks;//this is the threshold for stubs on trees, not total stubs!!
    //if(degree_dist[0]==1)threshold--; //this is because we allowed the tree_seed to be 1 degree so there is actualy 1 less 1-degree node

    //printf("\n\nStage 1: %d ones %d twos %d others %d tree stubs", ks, twos, others, total_tree_stubs);
    //print_stub_count(3);

    int j_watermark =0;

    //print_adj();

    //non_one_degree_nodes <- stack these assuming they have residual degree and correct threshold
    int poly_top =0;
    int alt, alt_degree, alt_offset, alt_tail;
    for(i=0; i < L;i++){
        if(total_tree_stubs <= threshold)break;//threshold condition
        //now we are going to build a queue of things for inter-stubs until the threshold is met
        int my_degree = degree_dist[i];
        int my_offset = offsets[i];
        int my_edge =my_offset;
        int my_capacity = my_degree -2;
        //if (poly_top==5)break;
        if(my_degree > 2){
            if(poly_top==0){
                non_one_degree_nodes[poly_top++] = i;
            }
            else if (poly_top > 0){//fill from behind but build up stack when nothing to do
                int cursor = 0;//this is outside because we only scan once per edges
                int accepted_edges = 0; 
                //todo check how big the stack gets - nice to know          
                while(my_edge < my_offset+my_degree){
                    int accept = 0;
                    while(cursor++ < poly_top){
                        alt = non_one_degree_nodes[cursor-1];
                        alt_degree = degree_dist[alt];
                        alt_offset = offsets[alt];  
                        alt_tail = alt_offset+alt_degree-1;
                        //break out on 'full' condition but also remove from stack
                        if(adjacency[alt_tail-1]!=-1){ //printf("\n**********remove %d  @ %d************\n", non_one_degree_nodes[cursor-1],i);
                                                 non_one_degree_nodes[cursor-1] =non_one_degree_nodes[--poly_top]; continue;}
                        //avoid self loop and multiple edges (chain condition)
                        if(adjacency[alt_tail]== i || alt==i)continue; //we would never accept this for any of our edges
                        //printf("\n>My i is %d and %d is the value of the pen-last of alt %d (adj index =%d)",i, adjacency[alt_tail], alt,alt_tail );
                        accept = 1; 
                        break;//passes all tests but we need to scan it below
                    }

                    if(accept){
                        //printf("\nFound possible stub matching %d-%d..cursor @ %d out of %d", alt, i, cursor, poly_top);
                        //here we should have a pointer to a candidate but we need to check we are not connected while scanning for a stub
                        int alt_edge =alt_offset;//start here
                        while(alt_edge < alt_offset+alt_degree){
                            //printf("\n@starting search at adj %d,%d", alt_edge,my_edge);
                            if(adjacency[alt_edge]==i)break;//already connected to one in our past
                            //printf(".");
                            if(total_tree_stubs <= threshold)break;//threshold condition - could detect this in other places
                            //printf(".");
                            if((adjacency[my_edge]==-1) && (adjacency[alt_edge]==-1)){//stub conditions
                                //printf("\n+ adj@ %d,%d linking %d-%d, current vals %d, %d",(my_edge),alt_edge, i, alt, adjacency[my_edge],adjacency[alt_edge]);
                                //if(alt==-1 || i ==-1){  printf("\n**********************FATAL ERROR - we have a negative index*************************"); }
                                adjacency[my_edge] = alt; 
                                adjacency[alt_edge] = i;
                                total_tree_stubs-=2;
                                //print_stub_count(3);
                                accepted_edges++;
                                break;
                            }
                            alt_edge++;
                        }
                    }
                    my_edge++;
                }
                //unless im done, add to the stack
                if(accepted_edges != my_capacity){non_one_degree_nodes[poly_top++] = i; }
            }
        }
    }

    //printf("\n\nStage 2: %d ones %d twos %d others %d tree stubs", ks, twos, others, total_tree_stubs);
    //print_stub_count(3);

    //print_adj();
    //return;

    /*MODE: IF WE HAVE AN EXCESS, WE CAN COBURN SOME BASED ON OUR TWO CAPACITY BUT FOR LARGE N MAY NOT BE NECESSARY*/
    //A TEST HERE SHOULD SHOW THAT *ALL* 2s are consumed - check last
    j_watermark =0; //- if we can only attach to > 2, watermark is safe
    for(i=0;i<L;i++){
        int my_degree = degree_dist[i];
        int my_offset = offsets[i];  
        if (my_degree==2){
            for(j=j_watermark/*IS THIS SAFE WATERMARK - SUPPOSED TO SAY WHERE WE HAVE FILLED TO?*/;j<L;j++ ){
                if(adjacency[my_offset]!=-1)break;
                int alt_degree = degree_dist[j];
                int alt_offset = offsets[j];
                int alt_edge = my_offset;
                int edge_found = 0;
                if (alt_degree > 2){//SHOULD THIS BE GREATER THAN 1 OR 2? if i make it two i need 
                     for(alt_edge=alt_offset; alt_edge<alt_offset+alt_degree;alt_edge++){
                         if(adjacency[alt_edge]==-1){
                             adjacency[my_offset] = j; 
                             adjacency[alt_edge] = i;
                             //total_tree_stubs-=2; IM REMOVING ONE BUT ADDING ONE (in another mode i could remove TWO)
                             edge_found =1;
                             break;
                         }
                     }
                }
                if(edge_found==1)break;
                else{j_watermark =j;} //- have not tested this but i think its fine - we are actually functionally the same as 1 unless we are coburning
            }
        }
    }

    //printf("\n\nStage 3: %d ones %d twos %d others %d tree stubs", ks, twos, others, total_tree_stubs);
    //print_stub_count(2);

    //now just add the ones but include newly added twos from the start with a new watermark
    j_watermark =0;
    for(i=0;i<L;i++){
        int my_degree = degree_dist[i];
        int my_offset = offsets[i];
        if (my_degree==1){
            for(j=j_watermark/*check safety*/;j<L;j++ ){
                if(adjacency[my_offset]!=-1)break;
                int alt_degree = degree_dist[j];
                int alt_offset = offsets[j];
                int alt_edge = my_offset;
                int edge_found = 0;
                if (alt_degree > 1){
                     for(alt_edge=alt_offset; alt_edge<alt_offset+alt_degree;alt_edge++){
                         if(adjacency[alt_edge]==-1){
                            adjacency[my_offset] = j; 
                            adjacency[alt_edge] = i;
                            total_tree_stubs--;
                            edge_found = 1;
                            break;
                         }
                     }
                }
                if(edge_found==1)break;    
                else{j_watermark =j;}//alt-j is alwats in degree excess of 1, and if we cant join now, noone can. so increment
            }
        }
    }

    if(total_tree_stubs != 0){
        printf("\nFATAL ERROR - COULD NOT CONSUME ALL STUBS!\n");
    }else{
        printf("\n#GRAPH GENERATED: 0 STUBS\n");
    }
    //print_adj();
    //printf("\n\nStage 4: %d ones %d twos %d others %d tree stubs", ks, twos, others, total_tree_stubs);
    //print_stub_count(0);
}
