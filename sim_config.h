#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using olb::util::Randomizer;
using olb::util::Timer;
using namespace particles;
using namespace particles::subgrid;
using namespace particles::communication;
using namespace particles::dynamics;
using namespace particles::creators;
using namespace particles::io;

using T = FLOATING_POINT_TYPE;
typedef D3Q19<> DESCRIPTOR;
typedef SubgridParticle3D PARTICLETYPE;
typedef BGKdynamics<T,DESCRIPTOR> BulkDynamics;

#define BOUZIDI

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Physical Settings
const T tau = 0.51;
const T rhoP = 998.2; // Physical Density (kg/m3)
const T muP = 0.02;
const T nuP = muP / rhoP; // Physical Nu
const T nuL = (tau - nuP) / 3; // Lattice Nu (LU)


// Sim Resolution Determination
const int N = 25;                         
const T charMinL =0.005; // m       
const T deltaX = charMinL / N; // LU 
const T deltaT = (nuL * deltaX * deltaX) / nuP; // s
const T maxAllowableLatticeVel = 0.4;

// Sim Time Settings
const T fluidMaxPhysT = T( 1 );     // max. fluid simulation time in s, SI unit
const T particleMaxPhysT = T( 10 ); // max. particle simulation time in s, SI unit  

// Average Velocity Determination
const T flowRate = 2. * 1.0e-6; // m3/s
const T radInlet = 0.005; // m
const T avgVel = flowRate / (M_PI * radInlet* radInlet); // m/s
const T avgLVel = (avgVel*deltaT) / deltaX; // LU

// Particle Settings
std::size_t noOfParticles = 2000;   // total number of inserted particles   
const T radius = 7.5e-5;            // particles radius
const T partRho = 998.2; 
const T particleInjectionX = 0.001;
const T injectionRadius = 3.0; // LU away from the wall
const T injectionP = 0.1; // LU into the geometry

//Set capture method:
// materialCapture: based on material number
// wallCapture:     based on more accurate stl description
typedef enum {materialCapture, wallCapture} ParticleDynamicsSetup;
const ParticleDynamicsSetup particleDynamicsSetup = wallCapture;

// center of inflow and outflow regions [m]
Vector<T, 3> inletCenter( 0.000, 0.000019, 0.000025 );
Vector<T, 3> outletCenter( 0.020, -0.000386, 0.000322 );


// normals of inflow and outflow regions
Vector<T, 3> inletNormal(T(-1), T(0), T(0) );
Vector<T, 3> outletNormal( T(1), T(0), T(0) );

//writing data constants
char vtkFileName[] = "rtdVal";
char vtkParticleFileName[] = "particle";
const T physVTKiter = 0.1;
const T physStatiter = 0.002;