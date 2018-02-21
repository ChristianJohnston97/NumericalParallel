// Christian Johnston, mdsw22
// Numerical Algorithms and Parallel Computing Summative: Spaceboddies
// Submitted Wednesday 21st Feburary 2018 


// Translate this file with
// g++ -O3 --std=c++11 spaceboddies.c -o spaceboddies

// Run it with
// ./spaceboddies time px py pz vx vy vz m


// e.g. 

// two boddie collision
// ./spaceboddies 10.0 0 0 0  0 0 0  2     0.001 0 0  -10 0 0  2

// 3 body in a row collision 
// ./spaceboddies 10.0 0 0 0 1 0 0 10 1 0 0 -1 0 0 20 2 0 0 -0.1 0 0 5


// Open MP compile command 
//g++-7 -O3 -fopenmp --std=c++11 spaceboddies.c -o spaceboddies 
// run in the same way

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <random>

//open MP
#include <omp.h>


// Global variables
double t = 0;
double tFinal = 0;

int NumberOfBodies = 0;

double** x;
double** v;
double*  mass;

double** force;

void setUp(int argc, char** argv) {

// if only one input argument, take this to be the number of bodies 
// how it runs for 10, 100, 1000 and 10000 particles
  if(argc == 2)
  {
    // convert first argument to an integer
    NumberOfBodies = atoi(argv[1]);
    x = new double*[NumberOfBodies];
    v = new double*[NumberOfBodies];
    mass = new double [NumberOfBodies];
    force = new double*[NumberOfBodies];

    // tFinal 
    tFinal = 0.1;

    for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];
    force[i] = new double[3];

    // random positions and velocities 
    std::random_device rd;
    std::default_random_engine e2(rd());
    // position range -10 --> 10
    std::uniform_real_distribution<> pos(-10,10);
    // velocity range -50 --> 50
    std::uniform_real_distribution<> vel(-50,50);


    // Create random numbers for position and velocity of all particles
    x[i][0] = pos(e2);
    x[i][1] = pos(e2);
    x[i][2] = pos(e2);

    v[i][0] = vel(e2);
    v[i][1] = vel(e2);
    v[i][2] = vel(e2);

    // initialise the force to 0 
    force[i][0] = 0.0;
    force[i][1] = 0.0;
    force[i][2] = 0.0;

    // constant mass of 10 for all particles
    mass[i] = 10;
    }
  }
  else
  {
    // else do normally
    // run-time command line arguments 
    // tFinal px py pz vx vy vz m for one body
    NumberOfBodies = (argc-2) / 7;

    x = new double*[NumberOfBodies];
    v = new double*[NumberOfBodies];
    force = new double*[NumberOfBodies];
    mass = new double [NumberOfBodies];


    int readArgument = 1;

    tFinal = std::stof(argv[readArgument]); readArgument++;


    for (int i=0; i<NumberOfBodies; i++) {
      x[i] = new double[3];
      v[i] = new double[3];
      force[i] = new double[3];

      x[i][0] = std::stof(argv[readArgument]); readArgument++;
      x[i][1] = std::stof(argv[readArgument]); readArgument++;
      x[i][2] = std::stof(argv[readArgument]); readArgument++;

      v[i][0] = std::stof(argv[readArgument]); readArgument++;
      v[i][1] = std::stof(argv[readArgument]); readArgument++;
      v[i][2] = std::stof(argv[readArgument]); readArgument++;

      mass[i] = std::stof(argv[readArgument]); readArgument++;

      if (mass[i]<=0.0 ) {
        std::cerr << "invalid mass for body " << i << std::endl;
        exit(-2);
        }
      force[i][0] = 0.0;
      force[i][1] = 0.0;
      force[i][2] = 0.0;

      }
    std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
    }
}


std::ofstream videoFile;


void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}



void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


void printParaviewSnapshot(int counter) {
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";


  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}

// function to print important details about particles
void printparticles () {
  printf("Printing particles for time %f.\n",t);
  for (int i=0; i<NumberOfBodies; i++) {
    printf(" NumberOfBodies: %d\n", NumberOfBodies);
    printf(" Particle %d:\n",i);
    printf(" Time: %f.\n",t);
    printf(" Pos: %f, %f, %f\n",x[i][0],x[i][1],x[i][2]);
    printf(" Vel: %f, %f, %f\n",v[i][0],v[i][1],v[i][2]);
    printf(" Mass: %f\n\n",mass[i]);
  }
}


void updateBody() {

// initial time step size, can change.
double timeStepSize = 1e-5;


// machine precision constant.
const double eps = 0.0001;
 
 // begin parallelisation region

#pragma omp parallel
{
  #pragma omp for 
    for (int i=0; i<NumberOfBodies; i++) {
    force[i][0] = 0.0;
    force[i][1] = 0.0;
    force[i][2] = 0.0;
  }
   
int i,j;
#pragma omp for
  for (i=0; i<NumberOfBodies; i++) {
    for (j=0; j<NumberOfBodies; j++) {

      // ignore if they are the same particle
      if(i != j)
      {
        // distance between i and j 
        double distance = sqrt(
        (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
        (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
        (x[i][2]-x[j][2]) * (x[i][2]-x[j][2])
        );
        //printf("Distance: %f \n ",distance);

        if(distance < 1e-8)
          // collision occures
          {
            // Combine all of particle i's and particle j's data into i
            // set j's information to NumberOfBodies-1 's
            // effectively 'removing' other collided particle

            printf("***COLLISION*** between particle:  %d and %d \n ",i, j);

            // Momentum conservation- inelastic collision
            v[i][0] = (mass[i]*v[i][0] + mass[j]*v[j][0])/(mass[i] + mass[j]);
            v[i][1] = (mass[i]*v[i][1] + mass[j]*v[j][1])/(mass[i] + mass[j]);
            v[i][2] = (mass[i]*v[i][2] + mass[j]*v[j][2])/(mass[i] + mass[j]);

            // Combined mass 
            mass[i] += mass[j];

            // Averaged position  
            x[i][0] = (x[i][0] + x[j][0]) / 2;
            x[i][1] = (x[i][1] + x[j][1]) / 2;
            x[i][2] = (x[i][2] + x[j][2]) / 2;
 
            // remove the other collided particle
            // set j's information to NumberOfBodies-1 's
                
              mass[j] = mass[NumberOfBodies-1];
              
               for(int direction=0; direction<3; direction++)
              {
                force[j][direction] = force[NumberOfBodies-1][direction];
                x[j][direction] = x[NumberOfBodies-1][direction];
                v[j][direction] = v[NumberOfBodies-1][direction];
              } 

            // recalculate the distance between i and the 'new' j
            double distance = sqrt(
            (x[i][0]-x[j][0]) * (x[i][0]-x[j][0]) +
            (x[i][1]-x[j][1]) * (x[i][1]-x[j][1]) +
            (x[i][2]-x[j][2]) * (x[i][2]-x[j][2])
            );
            // keeps track of the number of bodies over time
            // critical region to ensure only one thread enters this region and executes 
            // the commands
            #pragma omp critical(myCritical)
            {
              NumberOfBodies--;
              printparticles();
               //j--;
            }    
        }

        // recalculate distance and force
        double dist_cubed = distance * distance * distance;

        
        force[i][0] += (x[j][0]-x[i][0]) * mass[i]*mass[j] / dist_cubed;
        force[i][1] += (x[j][1]-x[i][1]) * mass[i]*mass[j] / dist_cubed;
        force[i][2] += (x[j][2]-x[i][2]) * mass[i]*mass[j] / dist_cubed;
      }
    }
    // Adaptive time-stepping (can turn off when required)

      while ((v[i][0]>0) && (std::abs(timeStepSize * force[i][0] / mass[i]) / v[i][0] > eps)
        ||  (v[i][1]>0) && (std::abs(timeStepSize * force[i][1] / mass[i]) / v[i][1] > eps)
        ||  (v[i][2]>0) && (std::abs(timeStepSize * force[i][2] / mass[i]) / v[i][2] > eps))
      {
        timeStepSize = timeStepSize / 2.0;
        // code below is for large number of particles, stops its getting stuck in with a tiny time step size 
        
        if(timeStepSize < 1e-6)
        {
          timeStepSize = 1e-6;
          break;
        }
      }
  }
  
    // Update particles position and velocity depending on time-step size 
#pragma omp for 
  for (int i=0; i<NumberOfBodies; i++) {
      x[i][0] = x[i][0] + timeStepSize * v[i][0];
      x[i][1] = x[i][1] + timeStepSize * v[i][1];
      x[i][2] = x[i][2] + timeStepSize * v[i][2];

      v[i][0] = v[i][0] + timeStepSize * force[i][0] / mass[i];
      v[i][1] = v[i][1] + timeStepSize * force[i][1] / mass[i];
      v[i][2] = v[i][2] + timeStepSize * force[i][2] / mass[i];
  }
}
  // end of parallel region

  // increment the time by the time step size
t += timeStepSize;
}



// Main method
int main(int argc, char** argv) {

    // We can call the clock function at the beginning and end of the code 
    // for which we measure time, subtract the values
    // Open MP function to calculate time when multiple threads used 

  double start_time = omp_get_wtime();

  if (argc==1) {
    std::cerr << "please add the final time plus a list of object configurations as tuples px py pz vx vy vz m" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
              << "100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-2)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }

  setUp(argc,argv);

  openParaviewVideoFile();

  printParaviewSnapshot(0);

// SET NUMBER OF THREADS
  omp_set_num_threads(2);

  int timeStepsSinceLastPlot = 0;
  const int plotEveryKthStep = 100;
    int i = 0;

 
    // this is the loop that keeps on going round, calling update body
    while (t<tFinal) {
      updateBody();
      // only print particles every 10,000 timesteps
      if (i++%10000 == 0) {
      printparticles();
      }
    
    timeStepsSinceLastPlot++;
    if (timeStepsSinceLastPlot%plotEveryKthStep==0) {
      printParaviewSnapshot(timeStepsSinceLastPlot/plotEveryKthStep);
    }
  }

  // Open MP function to calculate time when multiple threads used 
  double time = omp_get_wtime() - start_time;
  printf("Time taken: %f\n",time);

  closeParaviewVideoFile();

  return 0;
}

