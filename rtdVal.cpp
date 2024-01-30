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



void prepareGeometry(IndicatorF3D<T>& indicator,
                        STLreader<T>& stlReader,
                        SuperGeometry<T,3>& superGeometry,
                        UnitConverter<T, DESCRIPTOR> const& converter)
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

void setBoundaryValues()
{

}

bool getResults()
{
    return true;
}

int main( int argc, char* argv[] )
{
    return 0;
}