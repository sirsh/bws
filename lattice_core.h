#pragma once

#include <math.h>

#define HOLE (1<<6)
#define CACHE_CAPACITY 100000
#define TRACE_FLAG (1<<2)
#define CURRENT_FLAG (1<<3)
#define IMMOBILE_FLAG (1<<4)
#define ADD_FLAG ( TRACE_FLAG | CURRENT_FLAG )
//int inline get_site_default(int i);
void log_trace(int i);
double distance_from_center(int i);
int diffuse(int pos);
void init_lattice(int boundary_conditions, int L, int D);
int allocate_lattice(char** buffer, int L, int D, int type);
int get_center(int L, int D);
int probe(int loc, int dir, int offsets[]);
void right_turning_walk(int L, int n, long double time);
#define RNG_TYPE unsigned long
RNG_TYPE genrand_int32(void);
double genrand_real1(void);//[0,1]
double genrand_real2(void);//[0,1)

void generate_graph(int N,int allocate);
int build_edge_distribution(int n);
int init_cells(int D);
void lattice_remove(int p);

inline static int power_law_dist(int L, int stats[1]);
inline static int _power_sample(double xmax, double degree);
inline static void _allocate(unsigned int ** arr, int size, int arr_offset);
inline static void init_graph(unsigned int L, int type);
inline static void configure_edges(unsigned int L, int ks);
void generate_binary_tree(long long L,int allocate);
//int inline Xc(int i, int d);
//void regenerate_graph(int n);

typedef struct trace128{
unsigned long long elements[2];
} trace128;

typedef struct cell {
    int n;
    int div;
    int invalidated:1;
    int persisted:1;
    int trace_ref;
    struct cell* next;
} cell;
