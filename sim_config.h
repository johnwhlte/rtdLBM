#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;
using namespace particles;
using namespace particles::subgrid;
using namespace particles::communication;
using namespace particles::dynamics;
using namespace particles::creators;
using namespace particles::io;

using T = FLOATING_POINT_TYPE;
typedef D3Q19<> DESCRIPTOR;
typedef SubgridParticle3D PARTICLETYPE;

#define BOUZIDI

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Physical Settings
const T tau = 0.51;
const T nuP = 0.002; // Physical Nu
const T nuL = (tau - nuP) / 3; // Lattice Nu (LU)

// Sim Resolution Determination
int N = 40;                         
const T charMinL =0.001; // m       
const T deltaX = charMinL / N; // LU 
const T deltaT = (nuL * pow(deltaX,2)) / nuP; // s

// Sim Time Settings
const T fluidMaxPhysT = T( 5 );     // max. fluid simulation time in s, SI unit
const T particleMaxPhysT = T( 20 ); // max. particle simulation time in s, SI unit  

// Average Velocity Determination
const T flowRate = 2.; // mL/s
const T radInlet = 0.005; // m
const T avgVel = flowRate / (M_PI * pow(radInlet,2)); // m/s
const T avgLVel = (avgVel*deltaT) / deltaX; // LU

// Particle Settings
std::size_t noOfParticles = 10000;   // total number of inserted particles   
const T radius = 1.5e-4;            // particles radius
const T partRho = 998.2; 

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