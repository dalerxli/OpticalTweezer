#include <iostream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;
int ID=1;


const gsl_rng_type *T = gsl_rng_taus;
gsl_rng *r = gsl_rng_alloc(T);

double mean = 5;
double sigma = 0.5;

class Particle
{
	public:
	 double x,y,z;
	 double vx,vy,vz;
	 Particle();
};

Particle::Particle()
{
				x = 0;
				y = 0;
				z = 0;
				vx = 0;
				vy = 0;
				vz = 0;
}

struct Node
{
				Particle *part;
				int id;
				Node *next;
};

typedef struct Node s_Node;

Particle *testparticle = new Particle[20];

//initialization 
void initNode(struct Node *head)
{
				head -> part = new Particle();
				head -> id = ID;
				ID+=1;
				head -> part -> x = gsl_rng_uniform(r);
				head -> part -> y = gsl_rng_uniform(r);
				head -> part -> z = gsl_rng_uniform(r);
				head -> next = NULL;
}

//appending
void addNode(struct Node *head)
{
				Node *newNode = new Node;
				newNode -> part = new Particle();
				newNode -> id = ID;
				ID+=1;
				newNode -> next = NULL;

				Node *current = head;
				while(current)
				{
								if(current -> next == NULL)
								{
												current -> next = newNode;
												return;
								}
								current = current -> next;
				}
}

bool deleteNode(struct Node  **head, Node *toDelete)
{
				Node *current = *head;
				if(toDelete == *head)
				{
								cout << "Cannot delete Head" << endl;
								return false;
				}
				//if(toDele6te == *head)
				//{
								//*head = current -> next;
								//delete toDelete;
								//return true;
				//}

				while(current)
				{
								if(current -> next == toDelete)
								{
												current -> next = toDelete -> next;
												delete toDelete;
												return true;
								}
								current = current -> next;
				}
				return false;
}

void printList(struct Node *head)
{
				Node *current = head;
				while(current)
				{
								cout << "Particle " << current -> id << 
												"("<< current -> part -> x << "," << current -> part -> y << ","<< current -> part -> z << ")"<<endl;
								current = current -> next;
				}
				cout << "\n\n";
}

struct Node *searchNode(struct Node *head, int id)
{
				Node *current = head;
				while(current)
				{
								if(current -> id == id)
												return current;
								current = current -> next;
				}
				cout << "ID nicht gefunden!\n";
}

struct Node *LastNode(struct Node *head)
{
				Node *current = head;
				while(current)
				{
								if(current -> next == NULL)
												return current;
								current = current -> next;
				}
}

void InitializePositions(struct Node *head, struct Node *start)
{
				Node *current = head;
				while(current)
				{
								if(current == start)
								{
												break;
								}
								current = current -> next;
				}
				cout << "Starting Position: Particle " << current -> id << endl;
				while(current)
				{
								current -> part -> x = gsl_rng_uniform(r);
								current -> part -> y = gsl_rng_uniform(r);
								current -> part -> z = gsl_rng_uniform(r);
								current = current -> next;
				}
}

int InitializeVelocities(struct Node *head, struct Node *start,char dir)
{
				Node *current = head;
				while(current)
				{
								if(current == start)
								{
												break;
								}
								current = current -> next;
				}
				cout << "Starting Position: Particle " << current -> id << endl;
				while(current) // FEHLT: -x, -y, -z
				{
								if(dir == 'x')
								{
												current -> part -> x = gsl_ran_rayleigh(r,sigma);
												current -> part -> y = gsl_ran_gaussian(r,sigma);
												current -> part -> z = gsl_ran_gaussian(r,sigma);
								}
								else if(dir == 'y')
								{
												current -> part -> y = gsl_ran_rayleigh(r,sigma);
												current -> part -> x = gsl_ran_gaussian(r,sigma);
												current -> part -> z = gsl_ran_gaussian(r,sigma);
								}
								else if(dir =='z')
								{
												current -> part -> z = gsl_ran_rayleigh(r,sigma);
												current -> part -> y = gsl_ran_gaussian(r,sigma);
												current -> part -> x = gsl_ran_gaussian(r,sigma);
								}
								else
								{
												cout << "ERROR! Direction not correctly specified!" << endl;
												return 217;
								}

								current = current -> next;
				}
				return 0;
}
void newGasParticles(struct Node *head)
{
				struct Node *last = new Node;
				last = LastNode(head);
				unsigned int Nfac = gsl_ran_poisson(r,mean);
				for(int i=0;i<Nfac;i++)
				{
								addNode(head);
				}
				if(last->next != NULL)
								last = last -> next;
				InitializePositions(head,last);
}

void initParticle()
{
				for(int i=0;i<20;i++)
				{
								testparticle[i].x = gsl_rng_uniform(r);
								testparticle[i].y = gsl_rng_uniform(r);
								testparticle[i].z = gsl_rng_uniform(r);
								cout << "Particle " << i << ": " << testparticle[i].x << testparticle[i].y << testparticle[i].z << endl;
				}
				cout << endl;
}

void interact(struct Node *head)
{
				int i;
				double *r = new double[3];
				struct Node *current = head;
				while(current)
				{
								for(i=0;i<20;i++)
								{
										r[0] = current->part->x - testparticle[i].x;
										r[1] = current->part->y - testparticle[i].y;
										r[2] = current->part->z - testparticle[i].z;
										cout << "Distance Particle " << current->id << " to Testparticle " << i << ": " << sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])<<endl;
								}
								current = current -> next;
				}
}

void clearParticles(struct Node *head)
{
				
				struct Node *current = head;
				while(current)
				{
								if(abs(current -> part -> x) > 0.5 || abs(current -> part -> y) > 0.5 || abs(current -> part -> z) > 0.5)
								{
												cout << "Deleting Node " << current -> id << endl;
												deleteNode(&head,current);
												if(current -> next == NULL)
																break;
								}
								current = current -> next;
				}
}


int main()
{
				gsl_rng_set(r,132);
				initParticle();

				struct Node *particle_list = new Node;
				initNode(particle_list);

				newGasParticles(particle_list);
				printList(particle_list);
				interact(particle_list);
				newGasParticles(particle_list);
				printList(particle_list);
				newGasParticles(particle_list);
				clearParticles(particle_list);
				printList(particle_list);
				newGasParticles(particle_list);
				printList(particle_list);

				struct Node *last = new Node;
				last = LastNode(particle_list);
				cout << "Last Node: Particle " << last -> id << endl;


/*
 *        addNode(particle_list);
 *        for(int i=0;i<10;i++)
 *        {
 *                addNode(particle_list);
 *        }
 *        printList(particle_list);
 *        Node *delnode = searchNode(particle_list,5);
 *        if(deleteNode(&particle_list,delnode))
 *                cout << "Particle 5 deleted! \n";
 *        printList(particle_list);
 *
 *        struct Node *iterator = new Node;
 *        iterator = particle_list;
 *        while(iterator)
 *        {
 *                iterator->part->x=1;
 *                cout << "Particle " << iterator->id << " -> x = " << iterator->part->x << endl;
 *                iterator = iterator -> next;
 *        }
 *
 *        Node *part1 = searchNode(particle_list,1);
 *        part1 -> part -> y = 5;
 */


				return 0;
}
