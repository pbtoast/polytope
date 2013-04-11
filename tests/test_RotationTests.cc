// RotationTests
//
// A collection of rotation fields on the unit square
//
// Flow fields:
// 1. Constant Vorticity / Rigid Body Rotation 
// 2. Single "Vortex in a Box" (Bell, Colella, Glas, JCP, 1989)
// 3. Taylor-Green Vortex
// 4. 16-Vortex Deformation

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>
#include <sstream>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

#if HAVE_MPI
#include "mpi.h"
#endif

using namespace std;
using namespace polytope;


// -----------------------------------------------------------------------
// computeConstantVorticityFlow
// -----------------------------------------------------------------------
void computeConstantVorticityFlow(const vector<double>& points,
                                  vector<double>& velocities) {
   const unsigned numGenerators = points.size()/2;
   for (unsigned i = 0; i != numGenerators; ++i) {
      velocities[2*i  ] = (0.5 - points[2*i+1]) * M_PI / 314;
      velocities[2*i+1] = (points[2*i  ] - 0.5) * M_PI / 314;
   }
}

// -----------------------------------------------------------------------
// computeSingleVortexFlow
// -----------------------------------------------------------------------
void computeSingleVortexFlow(const vector<double>& points,
                             vector<double>& velocities) {
   const unsigned numGenerators = points.size()/2;
   for (unsigned i = 0; i != numGenerators; ++i) {
      velocities[2*i  ] =  sin(M_PI*points[2*i]) * cos(M_PI*points[2*i+1]);
      velocities[2*i+1] = -cos(M_PI*points[2*i]) * sin(M_PI*points[2*i+1]);
   }
}

// -----------------------------------------------------------------------
// computeTaylorGreenVortexFlow
// -----------------------------------------------------------------------
void computeTaylorGreenVortexFlow(const vector<double>& points,
                                  vector<double>& velocities) {
   const unsigned numGenerators = points.size()/2;
   for (unsigned i = 0; i != numGenerators; ++i) {
      velocities[2*i  ] =  0.5 * (cos(2*M_PI*(points[2*i  ]-0.25)) * 
                                  sin(2*M_PI*(points[2*i+1]-0.25)) );
      velocities[2*i+1] = -0.5 * (sin(2*M_PI*(points[2*i  ]-0.25)) * 
                                  cos(2*M_PI*(points[2*i+1]-0.25)) );
   }
}

// -----------------------------------------------------------------------
// computeDeformationFlow
// -----------------------------------------------------------------------
void computeDeformationFlow(const vector<double>& points,
                             vector<double>& velocities) {
   const unsigned numGenerators = points.size()/2;
   double x,y;
   for (unsigned i = 0; i != numGenerators; ++i) {
      x = points[2*i  ];
      y = points[2*i+1];
      velocities[2*i  ] = -sin(4*M_PI*(x + 0.5)) * sin(4*M_PI*(y + 0.5));
      velocities[2*i+1] = -cos(4*M_PI*(x + 0.5)) * cos(4*M_PI*(y + 0.5));
   }
}

// -----------------------------------------------------------------------
// getVelocities
// -----------------------------------------------------------------------
void getVelocities(const vector<double>& points,
                   const unsigned flowType,
                   vector<double>& velocities) {
   POLY_ASSERT(points.size() == velocities.size());
   switch(flowType){
   case 1:
      computeConstantVorticityFlow(points, velocities);
      break;
   case 2:
      computeSingleVortexFlow(points, velocities);
      break;
   case 3:
      computeTaylorGreenVortexFlow(points, velocities);
      break;
   case 4:
      computeDeformationFlow(points, velocities);
      break;
   }
}


// -----------------------------------------------------------------------
// outputMesh
// -----------------------------------------------------------------------
void outputMesh(Tessellation<2,double>& mesh, 
                const vector<double>& points,
                const unsigned flowType,
                int nstep,
                double time) {
#if HAVE_SILO
   vector<double> index(mesh.cells.size());
   vector<double> genx (mesh.cells.size());
   vector<double> geny (mesh.cells.size());
   for (int i = 0; i < mesh.cells.size(); ++i){
      index[i] = double(i);
      genx[i] = points[2*i  ];
      geny[i] = points[2*i+1];
   }
   map<string,double*> nodeFields, edgeFields, faceFields, cellFields;
   cellFields["cell_index"   ] = &index[0];
   cellFields["cell_center_x"] = &genx[0];
   cellFields["cell_center_y"] = &geny[0];
   ostringstream os;
   os << "test_RotationTests_test_" << flowType;
   polytope::SiloWriter<2, double>::write(mesh, nodeFields, edgeFields, 
                                          faceFields, cellFields, os.str(),
                                          nstep, time);
#endif
}


// -----------------------------------------------------------------------
// runTest
// -----------------------------------------------------------------------
void runTest(const unsigned flowType) {
   POLY_ASSERT(flowType >= 1 and flowType <= 4);

   // Boundary size parameters
   const unsigned nx = 50;
   const double xmin = 0.0, xmax = 1.0;
   const double ymin = 0.0, ymax = 1.0;
   const double dx = (xmax-xmin)/nx,  dy = (ymax-ymin)/nx;
   
   // Time stepping and point-resizing stuff
   double dt, Tmax, scaleFactor;
   switch(flowType){
   case 1:
      dt = 2.0;
      Tmax = 628.0/2;
      scaleFactor = sqrt(2);
      cout << "\nTest 1: Solid Rotation Flow\n" << endl;
      break;
   case 2:
      dt = sqrt(2)*dx;
      Tmax = 4.0;
      scaleFactor = 1.0;
      cout << "\nTest 2: Single Vortex Flow\n" << endl;
      break;
   case 3:
      dt = sqrt(2)*dx;
      Tmax = 4.0;
      scaleFactor = 2.0;
      cout << "\nTest 3: Taylor-Green (4-Vortex) Flow\n" << endl;
      break;
   case 4:
      dt = sqrt(2)*dx;
      Tmax = 2.0;
      scaleFactor = sqrt(2);
      cout << "\nTest 4: Deformation (16-Vortex) Flow\n" << endl;
      break;
   }

   // The PLC points for a unit square
   vector<double> PLCpoints;
   PLCpoints.push_back(xmin);  PLCpoints.push_back(ymin);
   PLCpoints.push_back(xmax);  PLCpoints.push_back(ymin);
   PLCpoints.push_back(xmax);  PLCpoints.push_back(ymax);
   PLCpoints.push_back(xmin);  PLCpoints.push_back(ymax);

   // The unit square PLC facets
   PLC<2,double> boundary;
   boundary.facets.resize(4, vector<int>(2));
   for (int i = 0; i != 4; ++i) {
      boundary.facets[i][0] = i;
      boundary.facets[i][1] = (i+1)%4;
   }
   
   // The generator set
   vector<double> points;
   unsigned ix, iy;
   double xi, yi;
   for (iy = 0; iy != nx; ++iy) {
      yi = ymin + (iy + 0.5)*dy;
      for (ix = 0; ix != nx; ++ix) {
         xi = xmin + (ix + 0.5)*dx;
         points.push_back(xi);
         points.push_back(yi);
      }
   }

   // Resize the generator so we dont' fling them out of the boundary
   for (unsigned i = 0; i != points.size(); ++i) {
      points[i] = 0.5 + (points[i]-0.5)/scaleFactor;  
   }

   // The tessellator
   TriangleTessellator<double> tessellator;

   // The velocity field
   vector<double> velocityField(points.size());


   // The initial tessellation
   unsigned step = 0;
   double time = 0.0;
   Tessellation<2,double> mesh;
   MeshEditor<2,double> meshEditor(mesh);
   tessellator.tessellate(points, PLCpoints, boundary, mesh);
   outputMesh(mesh, points, flowType, step, time);
   
   // Update the point positions and generate the mesh
   vector<double> halfTimePositions(points.size());
   while (time < Tmax) {
      if (step % 5 == 0) cout << (time/Tmax)*100 << "%" << endl;
      mesh.clear();
      getVelocities(points, flowType, velocityField);
      for (unsigned i = 0; i != points.size(); ++i) {
         halfTimePositions[i] = points[i] + dt*velocityField[i];
      }
      getVelocities(halfTimePositions, flowType, velocityField);
      for (unsigned i = 0; i != points.size(); ++i) {
         points[i] += dt*velocityField[i];
	 POLY_ASSERT(xmin <= points[i] and points[i] <= xmax);
      }
      time += dt;
      ++step;
      tessellator.tessellate(points, PLCpoints, boundary, mesh);
      meshEditor.cleanEdges(0.001);
      outputMesh(mesh, points, flowType, step, time);
   }
}


// -----------------------------------------------------------------------
// main
// -----------------------------------------------------------------------
int main(int argc, char** argv)
{
#if HAVE_MPI
   MPI_Init(&argc, &argv);
#endif

   for (unsigned flowTest = 1; flowTest < 5; ++flowTest) runTest(flowTest); 

   cout << "PASS" << endl;

#if HAVE_MPI
   MPI_Finalize();
#endif
   return 0;
}
