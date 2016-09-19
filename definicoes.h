#ifndef __definicoes_h_
#define __definicoes_h_

#include <assert.h>

#define MAX_POP_SIZE	8100
#define MAX_X			100
#define MAX_Y			100
#define MAX_VARIABLES   100

/////////////////////////////////////////////////////// STRUCTURES
typedef float PARTICLE[MAX_VARIABLES];

typedef struct {
	int x;
	int y;
	int new_x;
	int new_y;
	int old_mark_x;
	int old_mark_y;
	long double fitness;
	int mutated;
	int neighbors;
	int mutatedLastIteration;
	long double best_fitness_so_far;
	long double informants_best_fitness_so_far;
	PARTICLE informants_best_position_so_far;
	PARTICLE best_position_so_far;
	PARTICLE position;
	PARTICLE velocity;
} INDIVIDUAL;   

typedef INDIVIDUAL SWARM[MAX_POP_SIZE];

typedef struct {
	int occupied;      // if 1, the cell is occupied by at least one ant
	int mark;          // if > 0, the cell is marked by the ant with the id
	long double fitnessMark;
	PARTICLE positionMark;
    int ant;		   // if ocupied, ant is the id, else, -1
    long double fitness;     // fitness of the ant in the current cell		
	
} HABITAT;

typedef HABITAT CELL[MAX_X][MAX_Y];

typedef struct {
    PARTICLE best_position;
	PARTICLE best_position_so_far;    
	SWARM ant;
	CELL cell;
	float k;
	float distance;
	float minDistance;
	float maxDistance;
	long double minFitness;
	int movingAnts;
	int neighborsDistribution[10];
	unsigned int evaluations;
	long double average_fitness;
	long double best_fitness;
    long double best_so_far;
	long double worst_fitness;
	long double worst_so_far;
	int best_so_far_id;
    int worst_id;
	int worst_so_far_id;
} MODEL;

#endif
