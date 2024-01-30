/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Louis Kronberg, Christian Vorwerk, Bastian Schäffauer
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* rtdVal.cpp
Validation case for discrete particle RTD in OpenLB validated against
a previously run RTD in Ansys Fluent.
*/
#include "sim_config.h"



void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter, 
                      IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, 
                      SuperGeometry<T,3>& superGeometry)
{
    OstreamManager clout( std::cout,"prepareGeometry" );
    clout << "Prepare Geometry ..." << std::endl;

    //take all values in the 3d volume of the cylinder.stl and assign them 
    //material number 2
    superGeometry.rename(0, 2, indicator);
    //Now rename all values on the surface of the stl to material 1
    superGeometry.rename(2, 1, stlReader);
    //cleans auxillary voxels
    superGeometry.clean();

    //setup the inlet and outlet faces

    IndicatorCircle3D<T> inflow(inletCenter, inletNormal, radInlet);
    IndicatorCylinder3D<T> layerInflow(inflow, 2. * converter.getConversionFactorLength() );
    superGeometry.rename(2,3,1,layerInflow);

    IndicatorCircle3D<T> outflow(outletCenter, outletNormal, radInlet);
    IndicatorCylinder3D<T> layerOutflow(outflow, 2. * converter.getConversionFactorLength() );
    superGeometry.rename(2, 4, 1, layerOutflow);

    superGeometry.clean();
    superGeometry.innerClean( 3 );
    superGeometry.checkForErrors();

    superGeometry.print();
    clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(SuperLattice<T, DESCRIPTOR>& lattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry )


{
    OstreamManager clout( std::cout,"prepareLattice" );
    clout << "Prepare Lattice ..." << std::endl;

    const T omega = converter.getLatticeRelaxationFrequency();

    //set material 1 to bulk dynamics
    lattice.defineDynamics<BulkDynamics>(superGeometry, 1);

    setBounceBackBoundary(lattice, superGeometry, 2);
    setInterpolatedVelocityBoundary<T, DESCRIPTOR>(lattice, omega, superGeometry, 3);
    lattice.defineDynamics<BulkDynamics>(superGeometry.getMaterialIndicator(4));
    setInterpolatedPressureBoundary<T,DESCRIPTOR>(lattice, omega, superGeometry.getMaterialIndicator(4));

    // Initial conditions
    AnalyticalConst3D<T, T> rhoF( 1 );
    std::vector<T> velocity(3, T());
    AnalyticalConst3D<T, T> uF( velocity );

    lattice.defineRhoU( superGeometry.getMaterialIndicator({1, 3, 4}),rhoF,uF );
    lattice.iniEquilibrium( superGeometry.getMaterialIndicator({1, 3, 4}),rhoF,uF );
    lattice.setParameter<descriptors::OMEGA>(omega);

    lattice.initialize();

    clout << "Prepare Lattice ... OK" << std::endl;


}

void prepareParticles( UnitConverter<T,DESCRIPTOR> const& converter,
  SuperParticleSystem<T,PARTICLETYPE>& superParticleSystem,
  SolidBoundary<T,3>& wall,
  SuperIndicatorMaterial<T,3>& materialIndicator,
  ParticleDynamicsSetup particleDynamicsSetup,
  Randomizer<T>& randomizer)
{
  OstreamManager clout( std::cout, "prepareParticles" );
  clout << "Prepare Particles ..." << std::endl;

  //Add selected particle dynamics
  if (particleDynamicsSetup==wallCapture){
    //Create verlet dynamics with material aware wall capture
    superParticleSystem.defineDynamics<
      VerletParticleDynamicsMaterialAwareWallCapture<T,PARTICLETYPE>>(
        wall, materialIndicator);
  } else {
    //Create verlet dynamics with material capture
    superParticleSystem.defineDynamics<
      VerletParticleDynamicsMaterialCapture<T,PARTICLETYPE>>(materialIndicator);
  }

  // particles generation at inlet3
  Vector<T, 3> c( inletCenter );
  c[2] = 0.074;
  IndicatorCircle3D<T> inflowCircle( c, inletNormal, radInlet -
                                     converter.getConversionFactorLength() * 2.5 );
  IndicatorCylinder3D<T> inletCylinder( inflowCircle, 0.01 *
                                        converter.getConversionFactorLength() );

  //Add particles
  addParticles( superParticleSystem, inletCylinder, partRho, radius, noOfParticles, randomizer );

  //Print super particle system summary
  superParticleSystem.print();
  clout << "Prepare Particles ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T,DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,3>& superGeometry )
{
    OstreamManager clout( std::cout,"setBoundaryValues" );

    // No of time steps for smooth start-up
    const int iTmaxStart = converter.getLatticeTime( fluidMaxPhysT*0.2 );
    const int iTupdate = converter.getLatticeTime( 0.01 );

    if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
        sLattice.setProcessingContext(ProcessingContext::Evaluation);

        // Smooth start curve, sinus
        // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

        // Smooth start curve, polynomial
        PolynomialStartScale<T,int> startScale( iTmaxStart, T( 1 ) );

        // Creates and sets the Poiseuille inflow profile using functors
        int iTvec[1]= {iT};
        T frac[1]= {};
        startScale( frac,iTvec );
        std::vector<T> maxVelocity( 3,0 );
        maxVelocity[0] = frac[0]*converter.getLatticeVelocity(avgVel);

        T distance2Wall = converter.getConversionFactorLength()/2.;
        CirclePoiseuille3D<T> poiseuilleU( superGeometry, 3, maxVelocity[0], distance2Wall );
        sLattice.defineU( superGeometry, 3, poiseuilleU );
        sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
        ProcessingContext::Simulation);

        if ( iT % (10*iTupdate) == 0 && iT <= iTmaxStart ) {
        clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;
        }

        sLattice.setProcessingContext(ProcessingContext::Simulation);
    }

}

bool getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR>& converter, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer, STLreader<T>& stlReader )
{
  OstreamManager clout( std::cout, "getResults" );
  const int vtkIter  = converter.getLatticeTime( physVTKiter );
  const int statIter = converter.getLatticeTime( physStatiter );

  if ( iT==0 ) {
    SuperVTMwriter3D<T> vtmWriter(vtkFileName);

    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }
  // Writes output to vtk files for viewing in Paraview
  if ( iT%vtkIter==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);
    sLattice.scheduleBackgroundOutputVTK([&,iT](auto task) {
      SuperVTMwriter3D<T> vtmWriter(vtkFileName);
      SuperLatticePhysVelocity3D velocity(sLattice, converter);
      SuperLatticePhysPressure3D pressure(sLattice, converter);
      vtmWriter.addFunctor(velocity);
      vtmWriter.addFunctor(pressure);
      task(vtmWriter, iT);
    });
  }
  // Writes output on the console
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    clout << sLattice.getStatistics().getMaxU() << std::endl;
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Flux at the inflow and outflow region
    std::vector<int> materials = { 1, 3, 4};

    IndicatorCircle3D<T> inflow(inletCenter[0], inletCenter[1]-2.*converter.getConversionFactorLength(), inletCenter[2], inletNormal[0], inletNormal[1], inletNormal[2], radInlet+2.*converter.getConversionFactorLength() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( sLattice, converter, superGeometry, inflow, materials, BlockDataReductionMode::Discrete );
    vFluxInflow.print( "inflow","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( sLattice, converter, superGeometry, inflow, materials, BlockDataReductionMode::Discrete );
    pFluxInflow.print( "inflow","N","mmHg" );

    IndicatorCircle3D<T> outflow(outletCenter[0], outletCenter[1]+2.*converter.getConversionFactorLength(), outletCenter[2],  outletNormal[0], outletNormal[1], outletNormal[2], radInlet+2.*converter.getConversionFactorLength() );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow( sLattice, converter, superGeometry, outflow, materials, BlockDataReductionMode::Discrete );
    vFluxOutflow.print( "outflow","ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow( sLattice, converter, superGeometry, outflow, materials, BlockDataReductionMode::Discrete );
    pFluxOutflow.print( "outflow","N","mmHg" );
  }


  if ( sLattice.getStatistics().getMaxU() > 0.7 ) {
    clout << "PROBLEM uMax=" << sLattice.getStatistics().getMaxU() << std::endl;
    std::exit(0);
  }
    return true;
}

int main( int argc, char* argv[] )
{
    // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    (T) N,            // Number of lattice cells across the char length defined below
    (T) tau,          // relaxation time tau
    (T) charMinL,     // charPhysLength: reference length of simulation geometry
    (T) avgVel,    // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T) nuP,          // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T) rhoP          // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write(vtkFileName);

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  STLreader<T> stlReader( "myCylinder.stl", converter.getConversionFactorLength(), 0.001, 0, true );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getConversionFactorLength() );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = util::min(16*N, 8*singleton::mpi().getSize());
#else
  const int noOfCuboids = 2;
#endif
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(), noOfCuboids, "volume" );
  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  util::Timer<T> timer1( converter.getLatticeTime( fluidMaxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer1.start();

  prepareLattice( sLattice, converter, stlReader, superGeometry );

  timer1.stop();
  timer1.printSummary();

  // // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( fluidMaxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT <= converter.getLatticeTime( fluidMaxPhysT ); iT++ ) {
  //   // === 5th Step: Definition of Initial and Boundary Conditions ===
     setBoundaryValues( sLattice, converter, iT, superGeometry );

  //   // === 6th Step: Collide and Stream Execution ===
     sLattice.collideAndStream();

  //   // === 7th Step: Computation and Output of the Results ===
     getResults( sLattice, converter, iT, superGeometry, timer, stlReader );
   }

   timer.stop();
   timer.printSummary();
    return 0;
}