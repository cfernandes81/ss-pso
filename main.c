///////////////////////////////////////////////////////
// Graphic ANT COLONY
// Carlos Fernandes
// 

#include <malloc.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <winbgim.h>

#include "definicoes.h"
#include "aloca.h"
#include "functions.h"

// GRID CONSTANTS
#define WIDTH 1000
#define HEIGHT 700
#define GRID_MAX_SIZE 600

// *************** Global Variables *****************

// System's parameters
int max_x, max_y; 
int max_t;                        // habitat x*y
int max_evaluations;			                      // Maximum number of iterations
int n_ants;
int n_runs;			                      // Number of particles in the swarm at time t=0	
int position[8][2];                       // Array used to check the occupation of each ant' surroundig sites
int mutation;
int bakSneppen;                           // Run Bak-Sneppen model on the structure
int pso;                                    // run PSO on the structure
int randomStructure;  
int neighborhood; 
int gbest;   
int memory;                        
int problem;
float Xmax;
float Vmax;
float chi;
float omega;
float c;
int numberVariables;
int assyInitialization;
float initialXmin;
float initialXmax;
float crit;
int iWeightStrategy,cStrategy;
int parametros, semilla, pausa;


// ****************** Function Headings ************************
extern float **aloc_matrizf(int linhas , int colunas);
extern int **aloc_matrizi(int linhas , int colunas);
extern int *aloc_vetori(int linhas);
extern long *aloc_vetorl(int linhas);
extern float *aloc_vetorf(int linhas);
extern long double *aloc_vetorld(int linhas);
extern void desaloc_matrizf(float **Matriz , int linhas);
extern void desaloc_matrizi(int **Matriz , int linhas);
extern long double evaluate (PARTICLE position);
void DrawGrid(int fils, int cols);
void DrawParticles(unsigned int z, MODEL *antColony);

////////////////////////////// FUNCTIONS //////////////////////////////
void DrawGrid(int fils, int cols) {
int i,j;
int minx, miny, posx, posy;
int tamc;
int posax,posay;
    tamc = GRID_MAX_SIZE/fils; // CELLS SIZE
    minx = 20;
    miny = 40;
    setfillstyle(SOLID_FILL,WHITE); 
   // bar(minx,miny,minx+GRID_MAX_SIZE,miny+GRID_MAX_SIZE);
    setcolor(WHITE);
    for(i = 0;      i < fils;   i++) { // Drawing the Grid
       for(j = 0; j < cols; j++) {
          posx = (minx+(j*tamc));
          posy = (miny+(i*tamc));
          rectangle(posx,posy,(posx+tamc),(posy+tamc));
       }
    }
}
///////////////////////////////////////////////////////
void DrawParticles(unsigned int t, MODEL *antColony) {
int i,j;
int **auxcell;
int count;
int minx, miny, posx, posy;
int tamc;
int posax,posay;

    tamc= GRID_MAX_SIZE/max_y; // CELLS SIZE
    minx = 20;
    miny = 40;
    DrawGrid(max_y,max_x); // Drawing the empty grid
    
    for (i = 0;  i < max_x;  i++) {
		for (j = 0;  j < max_y;   j++) { // Drawing the particles
            
            if (antColony->cell[i][j].occupied == 1) {
               if (antColony->ant[antColony->cell[i][j].ant-1].mutated == 1)
                  setfillstyle(SOLID_FILL, RED);
               else
                   setfillstyle(SOLID_FILL, BLACK);
            }
            else 
               setfillstyle(SOLID_FILL, WHITE);
            
            posx= (minx+(j*tamc))+2;
            posy= (miny+(i*tamc))+2;
            bar(posx,posy,(posx+tamc-3),(posy+tamc-3));
		}
	}
}
/////RANDOM NUMBERS AND MATH FUNCTIONS
////////////////////////////////////////////////////////////
long double	real_al_entre_a_b(long double a,long double b) {
	return  a+(rand()*(b-a)) / (pow(2.0,15.0)-1.0);
}
////////////////////////////////////////////////////////////
long	aleatorio_entre_a_b(long double a,long double b) {
	return	(long) floor(a+(rand()*(b+1-a)) / (pow(2.0,15.0)));
}
////////////////////////////////////////////////////////////
double square(float x) {
	return (x*x);
}
/////////////////////////////////////////////////////////////


///////////////// INITIALIZE THE SWARM AND THE HABITAT
void initialize (MODEL *antColony) {
int i, j, z;
float aux;
    for (i = 0;   i < max_x;  ++i) {                 //INITIALIZE EMPTY GRID
        for (j = 0; j < max_y;     ++j) {
            antColony->cell[i][j].occupied = 0;
            antColony->cell[i][j].mark = 0;
            antColony->cell[i][j].ant = 0;
            antColony->cell[i][j].fitnessMark = 0;
            antColony->cell[i][j].fitness = 0;
            for (z = 0;   z < numberVariables;   ++z)
                 antColony->cell[i][j].positionMark[z] = 0;
        }
    } 
    for (i = 0;   i < n_ants;  ++i)  // Assign a random fitness value to each particle
        antColony->ant[i].fitness = real_al_entre_a_b(0.00001,1.0);
    for (i = 0;   i < n_ants-1;  ++i) { // rank by increasing fitness
        for (j = i+1; j <= n_ants-1;   ++j) {
            if (antColony->ant[i].fitness > antColony->ant[j].fitness) {
               aux = antColony->ant[j].fitness;
               antColony->ant[j].fitness = antColony->ant[i].fitness;
               antColony->ant[i].fitness = aux;
            }
        }
    }     
    for (i = 0;   i < n_ants;  ++i) { // Distribute the particles randomly through the grid
        do {
            antColony->ant[i].x = aleatorio_entre_a_b (0, max_x-1);
            antColony->ant[i].y = aleatorio_entre_a_b (0, max_y-1);
        } while (antColony->cell[antColony->ant[i].x][antColony->ant[i].y].occupied == 1);
        antColony->cell[antColony->ant[i].x][antColony->ant[i].y].occupied = 1;
        antColony->cell[antColony->ant[i].x][antColony->ant[i].y].ant = i+1;  // ant id
        antColony->cell[antColony->ant[i].x][antColony->ant[i].y].fitness = antColony->ant[i].fitness;  // keep fitness value in the cell
        antColony->ant[i].old_mark_x = antColony->ant[i].x;
        antColony->ant[i].old_mark_y = antColony->ant[i].y;   
        antColony->ant[i].new_x = antColony->ant[i].x;
        antColony->ant[i].new_y = antColony->ant[i].y;   
        antColony->ant[i].mutated = 0;
        antColony->ant[i].mutatedLastIteration = 0;
    }
    antColony->evaluations = 0;
    antColony->k = 0;
    antColony->distance = 0;
    antColony->minFitness = antColony->ant[0].fitness;
}           
///////////////////// BAK-SNEPPEN MUTATION
void bs(unsigned int t, MODEL *antColony) {
int i, j, ii, jj, index;
int minx, maxx, miny, maxy;
int mutation;
int id;
FILE *out1;
float minFitness;
int avalanche;
float aux;
int aux2; 
int mutations;
    avalanche = 0;
    mutations = 0;
    minFitness = 2;
    for (i = 0;  i < n_ants;  ++i)  // reinitialize the swarm: no mutation in this iteration
        antColony->ant[i].mutated = 0;
    for (i = 0;  i < n_ants;  ++i) { // Finds the particle with lower fitness
        if (antColony->ant[i].fitness < minFitness) {
           minFitness = antColony->ant[i].fitness;
           index = i;
        }
    }
    out1 = fopen("minFitness.txt","a");
    fprintf(out1,"%f\n", antColony->minFitness);
    fclose(out1);
    minx = antColony->ant[index].x-1;
    miny = antColony->ant[index].y-1;
    maxx = antColony->ant[index].x+1;
    maxy = antColony->ant[index].y+1;
    for (i = minx;	i <= maxx;	++i) { 
        for (j = miny;	j <= maxy;	++j) {
	       ii = i;
	       jj = j;
	       if (i < 0)      ii = max_x-1;
           if (i >= max_x) ii = 0;
           if (j < 0)      jj = max_y-1;
           if (j >= max_y) jj = 0;
           id = antColony->cell[ii][jj].ant;
           if (id > 0) { // there is a particle in the site, mutate it
              mutations = mutations+1;
              antColony->ant[id-1].fitness = real_al_entre_a_b(0.00001,1.0);
              antColony->cell[ii][jj].fitness = antColony->ant[id-1].fitness;
              antColony->ant[id-1].mutated = 1;
              antColony->ant[id-1].mutatedLastIteration = 1;
           }
        }
     }
     for (i = 0;  i < n_ants;  ++i) { // Finds the particle with lower fitness
        if (antColony->ant[i].fitness < minFitness)
           minFitness = antColony->ant[i].fitness;
    }  
    if (antColony->minFitness < minFitness) antColony->minFitness = minFitness;
    out1 = fopen("BAK-SNEPPEN.txt","a");
    fprintf(out1,"%f\t%d\t%d\n", antColony->minFitness, avalanche, mutations);
    fclose(out1);
}

////////////////////////////////////////////////////////////////////////////////
void rank(unsigned int t, MODEL *antColony) {
int i, j, z;
long double aux;
int aux2; 
float aux3;
   for (i = 0;   i < n_ants-1;  ++i) { //rank by fitness
        for (j = i+1; j <= n_ants-1;   ++j) {
            if (antColony->ant[i].fitness >= antColony->ant[j].fitness) {              
               //swap fitness
               aux = antColony->ant[j].fitness;
               antColony->ant[j].fitness = antColony->ant[i].fitness;
               antColony->ant[i].fitness = aux;
               //swap position
               aux2 = antColony->ant[j].x;
               antColony->ant[j].x = antColony->ant[i].x;
               antColony->ant[i].x = aux2;
               aux2 = antColony->ant[j].y;
               antColony->ant[j].y = antColony->ant[i].y;
               antColony->ant[i].y = aux2;
               //swap mutation condition
               aux2 = antColony->ant[j].mutated;
               antColony->ant[j].mutated = antColony->ant[i].mutated;
               antColony->ant[i].mutated = aux2;
               //swap ids in the cells
               aux2 = antColony->cell[antColony->ant[j].x][antColony->ant[j].y].ant;
               antColony->cell[antColony->ant[j].x][antColony->ant[j].y].ant = antColony->cell[antColony->ant[i].x][antColony->ant[i].y].ant;
               antColony->cell[antColony->ant[i].x][antColony->ant[i].y].ant = aux2;
               //change positions
               antColony->ant[i].new_x = antColony->ant[i].x;
               antColony->ant[i].new_y = antColony->ant[i].y;
               antColony->ant[j].new_x = antColony->ant[j].x;
               antColony->ant[j].new_y = antColony->ant[j].y;  
               // UPDATE PARTICLE SWARM
               aux = antColony->ant[j].best_fitness_so_far;
               antColony->ant[j].best_fitness_so_far = antColony->ant[i].best_fitness_so_far;
               antColony->ant[i].best_fitness_so_far = aux;
               aux = antColony->ant[j].informants_best_fitness_so_far;
               antColony->ant[j].informants_best_fitness_so_far = antColony->ant[i].informants_best_fitness_so_far;
               antColony->ant[i].informants_best_fitness_so_far = aux;
               for (z = 0;    z < numberVariables;   ++z) {
                   // swap position
                   aux3 = antColony->ant[j].position[z];
                   antColony->ant[j].position[z] = antColony->ant[i].position[z];
                   antColony->ant[i].position[z] = aux3;
                   // swap velocity
                   aux3 = antColony->ant[j].velocity[z];
                   antColony->ant[j].velocity[z] = antColony->ant[i].velocity[z];
                   antColony->ant[i].velocity[z] = aux3;
                   // swap particle's best position so far
                   aux3 = antColony->ant[j].best_position_so_far[z];
                   antColony->ant[j].best_position_so_far[z] = antColony->ant[i].best_position_so_far[z];
                   antColony->ant[i].best_position_so_far[z] = aux3;
                   // swap particle's informant best position so far
                   aux3 = antColony->ant[j].informants_best_position_so_far[z];
                   antColony->ant[j].informants_best_position_so_far[z] = antColony->ant[i].informants_best_position_so_far[z];
                   antColony->ant[i].informants_best_position_so_far[z] = aux3;
               }                                 
            }
        }
    }   
}

///////////////////// MUTATE 
void mutate(unsigned int t, MODEL *antColony) {
int i, j;
float aux;
int aux2;     
     for (i =  0; i < mutation;    ++i)
         antColony->ant[i].fitness = real_al_entre_a_b(0.00001,1.0);  
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////// PARTICLE SWARM OPTIMIZATION    ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
void initializeParticles (MODEL *antColony) {
int i,j;
float xmin, xmax;
    //Initialize Position and Velocity
    for (i = 0;  i < n_ants;    ++i) {
        if (assyInitialization == 1) { // Assymetric initialization of the population 
           xmin = initialXmin;
           xmax = initialXmax;
        }
        else {                         // Normal initialization
             xmin = -Xmax;
             xmax = Xmax;
        }
        for (j = 0;  j < numberVariables;   ++j) { 
            antColony->ant[i].position[j] = real_al_entre_a_b(xmin, xmax);
            antColony->ant[i].best_position_so_far[j] = antColony->ant[i].position[j];  
            antColony->ant[i].informants_best_position_so_far[j] = antColony->ant[i].position[j];          
            antColony->ant[i].velocity[j] = real_al_entre_a_b(-Xmax, Xmax)*(0.5-real_al_entre_a_b(0,1.0));
        }
        antColony->ant[i].fitness = evaluate (antColony->ant[i].position);
        antColony->ant[i].best_fitness_so_far = antColony->ant[i].fitness;
        antColony->ant[i].informants_best_fitness_so_far = antColony->ant[i].fitness; 
    }
    antColony->best_so_far = antColony->ant[0].fitness;
    antColony->best_so_far_id = 0;
}

///////////////////////////////////
void updatePopulationData (MODEL *antColony) {
int i,j;
FILE *out1;
    antColony->best_fitness = antColony->ant[0].fitness;
    antColony->worst_fitness = 0;
    antColony->average_fitness = 0;
    for (i = 0;  i < n_ants;     ++i) {
        // Updates worst in population
        if (antColony->ant[i].fitness > antColony->worst_fitness) {
           antColony->worst_fitness = antColony->ant[i].fitness;
           antColony->worst_id = i;
        }
        // Updates best_so_far in population
        if (antColony->ant[i].fitness < antColony->best_so_far) {
           antColony->best_so_far = antColony->ant[i].fitness;
           for (j = 0;   j < numberVariables;   ++j)
               antColony->best_position_so_far[j] = antColony->ant[i].position[j];
        }
        // Updates best in current population
        if (antColony->ant[i].fitness < antColony->best_fitness) {
           antColony->best_fitness = antColony->ant[i].fitness;
           for (j = 0;   j < numberVariables;   ++j)
               antColony->best_position[j] = antColony->ant[i].position[j];
        }
        // Updates particle's best position 
        if (antColony->ant[i].fitness < antColony->ant[i].best_fitness_so_far) {
           antColony->ant[i].best_fitness_so_far = antColony->ant[i].fitness;
           for (j = 0;  j < numberVariables;  ++j)
               antColony->ant[i].best_position_so_far[j] = antColony->ant[i].position[j];
        }
        // Updates best informant
        if (antColony->ant[i].fitness < antColony->ant[i].informants_best_fitness_so_far) {
           antColony->ant[i].informants_best_fitness_so_far = antColony->ant[i].fitness;
           for (j = 0;  j < numberVariables;   ++j)
               antColony->ant[i].informants_best_position_so_far[j] = antColony->ant[i].position[j];
        }
        antColony->average_fitness = antColony->average_fitness+antColony->ant[i].fitness;
    }
    // Updates informants best position and fitness
    antColony->average_fitness = antColony->average_fitness/n_ants;
}

/////////////////////////////////////    
void updateParticles (unsigned int t, MODEL *antColony, int i){
int j, k;
float v, x;
float pi, pg;
long double phi1, phi2;
float c1, c2;
float maxIW, minIW;
FILE *out1;
c1 = c;
    c2 = c;
    maxIW = 0.9;
    minIW = 0.4;
  // if (antColony->ant[i].neighbors > 1) {
          if (iWeightStrategy == 1)     // TVIW-PSO
             omega = minIW+(maxIW-minIW)*(((float)max_t-(float)t)/(float)max_t);     
          if (cStrategy == 1) {         //TVAC-PSO
             c1 = (0.5-2.5)*((float)t/(float)max_t)+2.5;
             c2 = (2.5-0.5)*((float)t/(float)max_t)+0.5; 
          }
          for (j = 0;  j < numberVariables;    ++j) {  
              if (gbest == 0) pg = antColony->ant[i].informants_best_position_so_far[j];
              if (gbest == 1) pg = antColony->best_position_so_far[j]; 
              pi = antColony->ant[i].best_position_so_far[j];
              v = antColony->ant[i].velocity[j];
              x = antColony->ant[i].position[j];
              phi1 = real_al_entre_a_b (0.0, c1);
              phi2 = real_al_entre_a_b (0.0, c2);
              // Update Velocity
              v = omega*v+(float)phi1*(pi-x)+(float)phi2*(pg-x);
              if (v > Vmax) v = Vmax;
              if (v < -Vmax) v = -Vmax;
              // Update Position
              x = x+v;
              if (x > Xmax) {
                 x = Xmax;
                 v = 0;
              }
              if (x < -Xmax) {
                 x = -Xmax;   
                 v = 0;
              }
              antColony->ant[i].position[j] = x;  
              antColony->ant[i].velocity[j] = v; 
          }
          
          if (antColony->ant[i].neighbors > 0) {
            antColony->ant[i].fitness = evaluate (antColony->ant[i].position);   
            antColony->evaluations = antColony->evaluations+1;
            if (antColony->evaluations == 49000 || antColony->evaluations == 147000 || antColony->evaluations == 294000 || antColony->evaluations == 490000) {
             out1=fopen("INTERMEDIARY.DAT","a");
             fprintf(out1,"%.50f\t", (float)antColony->best_so_far);	
             fclose (out1);
          }
    }
}
// END PARTICLE SWARM OPTIMIZATION

//////////////////////MOVE THE PARTICLES
void move(unsigned int t, MODEL *antColony) {
int i, a, j, w, ii, jj, z;
int minx, maxx, miny, maxy;
int **position;
int freeSlots;
int destinationCell_x, destinationCell_y;
int currentMark;
long double currentFitnessMark;
long double diference;
long double currentMaxDistance;
int neighborAnt;
float distance;
long rol;
int move;
int nei;
long double distanceToMark;
int update;
     antColony->k = 0;
     antColony->distance = 0;
     position = aloc_matrizi(10,2);
     antColony->minDistance = 1.0;
     antColony->maxDistance = 0.0;
     move = 0;
     for (i = 0;   i < 9;   ++i)
         antColony->neighborsDistribution[i] = 0;
     for (a = 0;   a < n_ants;       ++a) {
         update = 0;
         freeSlots = 0;
         currentFitnessMark = 0;
         currentMaxDistance = 10000000;
         antColony->ant[a].neighbors = 0;
         distance = 0;
         if (memory == 0) {
             antColony->ant[a].informants_best_fitness_so_far = antColony->ant[a].best_fitness_so_far;
             for (i = 0;    i < numberVariables;   ++i) 
                   antColony->ant[a].informants_best_position_so_far[i] = antColony->ant[a].best_position_so_far[i];        
         }                        
         minx = antColony->ant[a].x-1;
         miny = antColony->ant[a].y-1;
         maxx = antColony->ant[a].x+1;
         maxy = antColony->ant[a].y+1;
         destinationCell_x = antColony->ant[a].x;
         destinationCell_y = antColony->ant[a].y;
         for (i = minx;	i <= maxx;	++i) { 
	         for (j = miny;	j <= maxy;	++j) {
		         ii = i;
		         jj = j;
		         if (i < 0)      ii = max_x-1;    // ISTO ESTA ERRADO!!!!!!
                 if (i >= max_x) ii = 0;
                 if (j < 0)      jj = max_y-1;    // ISTO ESTA ERRADO!!!!!!
                 if (j >= max_y) jj = 0; 
                 if (antColony->cell[ii][jj].occupied == 0) { // cell is not ocupied
                    position[freeSlots][0] = ii;
                    position[freeSlots][1] = jj;
                    freeSlots = freeSlots+1;
                    if (antColony->cell[ii][jj]. mark > 0) { // cell is ocupied and has a mark
                       // compute euclidean distance between ant and mark
                       diference = 0;
                       for (z = 0; z < numberVariables;  ++z)
                           diference = diference+(sqrt(square(antColony->ant[a].position[z]-antColony->cell[ii][jj].positionMark[z])));
                       diference = diference/numberVariables;
                       //printf("\nDiference = %f", (float)diference);
                       if (diference >= currentFitnessMark) {   
                          currentFitnessMark = diference;
                          destinationCell_x = ii;
                          destinationCell_y = jj;
                       }
                    }
                 }
                 if (antColony->cell[ii][jj].occupied == 1) {  // theres is a neighbor
                                    
                    // Updates best neighbor 
                    if (neighborhood == 0 || (i == minx+1 && j == miny+1) || (i == minx+1 && j == miny) || (i == minx && j == miny+1) || (i == minx+1 && j == maxy) || (i == maxx && j == miny+1)) {
                       if (antColony->cell[ii][jj].ant == n_ants) // the worst ant is a neighbor 
                             update = 1;    // mark particle for updating; used for new PSO based on Bak-Sneppen strategy 
                       antColony->ant[a].neighbors = antColony->ant[a].neighbors+1;
                       neighborAnt = antColony->cell[ii][jj].ant-1;
                       if (antColony->ant[neighborAnt].best_fitness_so_far < antColony->ant[a].informants_best_fitness_so_far) {
                          antColony->ant[a].informants_best_fitness_so_far = antColony->ant[neighborAnt].best_fitness_so_far;
                          for (z = 0;   z < numberVariables; ++z)
                              antColony->ant[a].informants_best_position_so_far[z] = antColony->ant[neighborAnt].best_position_so_far[z];
                       }
                    }
                    /////////////////////////// measures distance between neighboring ants                                         
                    distance = distance+sqrt(square(antColony->cell[ii][jj].fitness-antColony->ant[a].fitness));
                    if (sqrt(square(antColony->cell[ii][jj].fitness-antColony->ant[a].fitness)) > antColony->maxDistance)
                       antColony->maxDistance = sqrt(square(antColony->cell[ii][jj].fitness-antColony->ant[a].fitness));
                    if (sqrt(square(antColony->cell[ii][jj].fitness-antColony->ant[a].fitness)) < antColony->minDistance)
                       antColony->minDistance = sqrt(square(antColony->cell[ii][jj].fitness-antColony->ant[a].fitness));  
                 }       
           }
         }
         nei = antColony->ant[a].neighbors;
         antColony->neighborsDistribution[nei] = antColony->neighborsDistribution[nei]+1;
         antColony->k = antColony->k+(float)nei;
         if (antColony->ant[a].neighbors > 1)
            antColony->distance = antColony->distance+(distance/(float)(nei));          
         if (freeSlots == 0) {   // the ants is surrounded and cannot move
            antColony->ant[a].new_x = antColony->ant[a].x;
	        antColony->ant[a].new_y = antColony->ant[a].y;
         }
         else {   
              move = move+1;
              if (currentFitnessMark < 10000000 && randomStructure == 0) {
                 antColony->ant[a].new_x = destinationCell_x;
	             antColony->ant[a].new_y = destinationCell_y;
              }
              else {
                   rol = aleatorio_entre_a_b (0.0, (long double) (freeSlots-1));
		           antColony->ant[a].new_x = position[rol][0];
	               antColony->ant[a].new_y = position[rol][1];
              }
	          antColony->cell[antColony->ant[a].x][antColony->ant[a].y].occupied = 0;                  // Frees the cell
	          antColony->cell[antColony->ant[a].x][antColony->ant[a].y].fitness = 0;
	          antColony->cell[antColony->ant[a].x][antColony->ant[a].y].ant = 0;                      
	          antColony->cell[antColony->ant[a].x][antColony->ant[a].y].mark = a+1;                    // Mark old cell
	          antColony->cell[antColony->ant[a].x][antColony->ant[a].y].fitnessMark = antColony->ant[a].fitness;
	          for (z = 0; z < numberVariables;  ++z)
	              antColony->cell[antColony->ant[a].x][antColony->ant[a].y].positionMark[z] = antColony->ant[a].position[z];
	          
	          antColony->cell[antColony->ant[a].old_mark_x][antColony->ant[a].old_mark_y].mark = 0;   // remove old mark
	          antColony->cell[antColony->ant[a].old_mark_x][antColony->ant[a].old_mark_y].fitnessMark = 0;   // remove old fitness mark
              antColony->ant[a].old_mark_x = antColony->ant[a].x;                                     // saves latest marked position
              antColony->ant[a].old_mark_y = antColony->ant[a].y;
              
              antColony->cell[antColony->ant[a].new_x][antColony->ant[a].new_y].occupied = 1;	 // Occupy new cell
              antColony->cell[antColony->ant[a].new_x][antColony->ant[a].new_y].fitness = antColony->ant[a].fitness;
              antColony->cell[antColony->ant[a].new_x][antColony->ant[a].new_y].ant = a+1;       // Occupy new cell - ant id
              antColony->cell[antColony->ant[a].new_x][antColony->ant[a].new_y].mark = 0;        // clears new cell 
              antColony->cell[antColony->ant[a].new_x][antColony->ant[a].new_y].fitnessMark = 0; // no mark
              for (z = 0; z < numberVariables;  ++z)
	              antColony->cell[antColony->ant[a].new_x][antColony->ant[a].new_y].positionMark[z] = 0;
              
              antColony->ant[a].x = antColony->ant[a].new_x;                                      // updates position
              antColony->ant[a].y = antColony->ant[a].new_y;
              // Update PSO fitness
              
        }
        if (pso == 1)
            updateParticles (t, antColony, a);
        if (pso == 2) { // new PSO; only the worst and its neigbobors are updated
           if (update == 1)
              updateParticles (t, antColony, a);
        }   
     }
     //printf ("\n kkkkkkkkkk");
     antColony->k = antColony->k/((float)n_ants);
     antColony->distance = antColony->distance/n_ants;
     antColony->movingAnts = move;
     desaloc_matrizi (position, 10);
}

//################################################################## MAIN ALGORITHM
int main(int argc, char* argv[]) {
FILE *out1;
int i,j,n,a,w;
MODEL *antColony;
unsigned int z, counter;
char cadena[50];
float minPosDistance;
FILE *in1;
long double *averageBestsofar;
int flag;
    // Read parameters ////////////
    in1=fopen("INPUT.txt","r");
	fscanf(in1,"%d", &max_x);
	fscanf(in1,"%d", &max_y);
	fscanf(in1,"%d", &n_runs);
    fscanf(in1,"%d", &max_t);
	fscanf(in1,"%d", &max_evaluations);
	fscanf(in1,"%d", &n_ants);
	fscanf(in1,"%d", &mutation);
	fscanf(in1,"%d", &bakSneppen);
	fscanf(in1,"%d", &pso);
	fscanf(in1,"%d", &randomStructure);
	fscanf(in1,"%d", &neighborhood);
	fscanf(in1,"%d", &gbest);
	fscanf(in1,"%d", &memory);
	fscanf(in1,"%d", &problem);
	fscanf(in1,"%f", &Xmax);
	fscanf(in1,"%f", &Vmax);
	fscanf(in1,"%f", &chi);
	fscanf(in1,"%f", &omega);
	fscanf(in1,"%f", &c);
	fscanf(in1,"%u", &numberVariables);
	fscanf(in1,"%d", &iWeightStrategy);
	fscanf(in1,"%d", &cStrategy);
	fscanf(in1,"%d", &assyInitialization);
	fscanf(in1,"%f", &initialXmin);
	fscanf(in1,"%f", &initialXmax);
	fscanf(in1,"%f", &crit);
	fclose (in1);
	//////////////////////////////
	averageBestsofar = aloc_vetorld(50000);
    for (i = 0;      i < 50000;           ++i)
        averageBestsofar[i] = 0.0;
    initwindow(WIDTH,HEIGHT); 
    settextstyle(0,0,2);
    setcolor(WHITE);
    sprintf(cadena,"%d PARTICLES ON A %dx%d GRID", n_ants, max_x, max_y);
    outtextxy(60,15,cadena);
    DrawGrid(max_y,max_x);
    setfillstyle(SOLID_FILL,BLACK);
    bar(690,130,990,630);
	    
        
    for (w = 0;   w < n_runs;   ++w) {
	    minPosDistance = (2*(1/(float)n_ants)+2*2*(1/(float)n_ants)+2*3*(1/(float)n_ants)+2*4*(1/(float)n_ants))/8;
        if((antColony = (MODEL *) calloc(1,sizeof(MODEL))) == NULL) {
            printf("\nERROR: Out of Memory");
            exit(0);
        }
        flag = 0;
        z = 0;
        counter = 1;
        initialize(antColony);
        DrawParticles(z,antColony);
        setfillstyle(SOLID_FILL,BLACK);
        bar(690,130,990,630);
        if (pso >= 1) 
           initializeParticles (antColony);
        // *************** MAIN CYCLE *****************
        do {
            z = z+1;
            if (bakSneppen == 1) bs (z, antColony);
            if (mutation > 0) mutate (z, antColony);
            if (bakSneppen == 1 || mutation > 0) rank (z, antColony);
            if (pso >= 1) {
               updatePopulationData (antColony);
               rank (z, antColony);
            }
            
            if (pso == 0) DrawParticles(z, antColony);
          //  DrawParticles(z, antColony);
            move(z, antColony);
            if (antColony->evaluations > counter*100) {
               averageBestsofar[counter] = averageBestsofar[counter]+antColony->best_so_far/(long double)n_runs;
               counter = counter+1;
            }
            if (antColony->best_so_far < crit && flag == 0) {
                out1=fopen("AES.DAT","a");
                fprintf(out1,"\n%d", antColony->evaluations);	
                fclose (out1);
                //z = max_t;
                flag = 1;
            }  
            setcolor(WHITE);
            settextstyle(0,0,2);
            sprintf(cadena,"t = %d / %d", z, max_t);
            outtextxy(650,50,cadena);  
            sprintf(cadena,"k = %f", antColony->k);
            outtextxy(650,100,cadena);
            sprintf(cadena,"dist ants = %f", minPosDistance);
            outtextxy(650,120,cadena); 
            sprintf(cadena,"aver Dist = %f", antColony->distance);
            outtextxy(650,140,cadena);  
            sprintf(cadena,"min Dist = %f", antColony->minDistance);
            outtextxy(650,160,cadena);
            sprintf(cadena,"max Dist = %f", antColony->maxDistance);
            outtextxy(650,180,cadena);
            sprintf(cadena,"moving = %d", antColony->movingAnts);
            outtextxy(650,220,cadena);
            sprintf(cadena,"bestfit = %.5f", (float)antColony->best_so_far);
            outtextxy(650,240,cadena);
         /*   out1 = fopen("K.txt","a");
            fprintf(out1,"%f\t%d\n", antColony->k, antColony->movingAnts);
            fclose(out1);*/
            if (n_runs == 1) {
               /* if (max_x > 1) {
                  out1 = fopen("antsDistribution.DAT","w");
                       for (j = 0;  j < max_x;  ++j) {
                           for (i = 0;  i < max_y;  ++i) {
                               if (antColony->cell[j][i].occupied == 0)
                                  fprintf(out1,"%d\n", 255);
                               else
                                   fprintf(out1,"%d\n", 0);
                           }
                       }
                  fclose(out1);
               }*/
               if (max_x == 1) {
                  out1 = fopen("antsDistribution.DAT","a");
                  for (i = 0;  i < max_y;  ++i) {
                      if (antColony->cell[0][i].occupied == 0)
                         fprintf(out1,"%d\n", 0);
                      else
                          fprintf(out1,"%d\n", 255);
                  }
                  for (i = max_y; i <1000;  ++i)
                      fprintf(out1,"%d\n", 0);
                  fclose(out1);
               }
               out1 = fopen("FITNESS_MAPS.DAT","w");
               for (j = 0;   j < max_x;     ++j) {
                   for (i = 0;  i < max_y;  ++i) {
                       if (antColony->cell[j][i].occupied == 0)
                          fprintf(out1,"%d\n", 255);
                       else
                           fprintf(out1,"%d\n", 255-(int)(antColony->ant[antColony->cell[j][i].ant-1].fitness*255));
                   }
               }
               fclose(out1);
               out1 = fopen("Kdistribution.txt","a");
               for (i = 0;  i <= 9;  ++i)
                   fprintf(out1,"%d\t", antColony->neighborsDistribution[i]);
               fprintf(out1,"\n");
               fclose(out1);
               out1 = fopen("DISTANCE.txt","a");
               fprintf(out1,"%f\t%f\t%f\n", antColony->distance, antColony->minDistance, antColony->maxDistance);
               fclose(out1);
            }
            /*if (antColony->evaluations == 1000 || z == 3000 || z == 6000 || z == 10000) {
               out1=fopen("INTERMEDIARY.DAT","a");
               fprintf(out1,"%.41f\t", (float)antColony->best_so_far);	
               fclose (out1);
            }*/
        } while (antColony->evaluations < max_evaluations);
        out1=fopen("INTERMEDIARY.DAT","a");
        fprintf(out1,"\n");	
        fclose (out1);
        out1=fopen("FINAL.DAT","a");
        fprintf(out1,"%.45f\n", (float)antColony->best_so_far);	
        fclose (out1);
        //getch();
        free (antColony);
    } 
    out1=fopen("AVE_BESTSOFAR.DAT","a");
	for (i = 1;	i < counter+1;	++i) 
	   fprintf(out1,"%.40f\n", (float)averageBestsofar[i]);	
	fclose (out1);
    free(averageBestsofar);
    return 0;
} // main
/////////////////////////////////////////////////////// END ANT ALGORITHM
