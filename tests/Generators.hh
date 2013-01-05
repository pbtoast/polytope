#ifndef POLYTOPE_GENERATORS_HH
#define POLYTOPE_GENERATORS_HH

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>

#include "polytope.hh"
#include "Boundary2D.hh"

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------
template<int Dimension, typename RealType>
class Generators
{
public:
   // -------------- Public member variables and routines ---------------- //

   // Number of generators
   unsigned nPoints;
   vector<RealType> mGenerators;
   Boundary2D<RealType>& mBoundary;
   
   
   //------------------------------------------------------------------------
   // Constructor, destructor
   //------------------------------------------------------------------------
   Generators(Boundary2D<RealType>& boundary):
      nPoints(0),
      mGenerators(0),
      mBoundary(boundary) {};

   ~Generators() {};


   //------------------------------------------------------------------------
   // Place random generators into spatial domain
   //------------------------------------------------------------------------
   void randomPoints(unsigned nGenerators)
   {
      mGenerators.clear();
      nPoints = nGenerators;

      mBoundary.getBoundingBox();
      POLY_ASSERT( mBoundary.mLow  != 0 );
      POLY_ASSERT( mBoundary.mHigh != 0 );
      
      for (unsigned iter = 0; iter < nGenerators; ++iter ){
         RealType pos[Dimension];
         bool inside = false;
         while( !inside )
         {
            for (unsigned n = 0; n < Dimension; ++n){
               pos[n] = (mBoundary.mHigh[n]-mBoundary.mLow[n]) * 
                  RealType(::random())/RAND_MAX + mBoundary.mLow[n];
            }
            inside = mBoundary.testInside( pos );
         }
         mGenerators.insert( mGenerators.end(), pos, pos + Dimension );
      }
      POLY_ASSERT( mGenerators.size()/Dimension == nGenerators );
   }

   
   //------------------------------------------------------------------------
   // Place Cartesian points of constant mesh spacing
   //------------------------------------------------------------------------
   void cartesianPoints(vector<unsigned> nCellsPerDimension)
   {
      mGenerators.clear();
      POLY_ASSERT( nCellsPerDimension.size() == Dimension );

      mBoundary.getBoundingBox();
      POLY_ASSERT( mBoundary.mLow  != 0 );
      POLY_ASSERT( mBoundary.mHigh != 0 );

      if( Dimension == 2 ){
         cartesian2D( nCellsPerDimension[0], nCellsPerDimension[1] );
      }else{
         cartesian3D( nCellsPerDimension[0], nCellsPerDimension[1], nCellsPerDimension[2] );
      }
      nPoints = mGenerators.size()/Dimension;
   }
    

   //------------------------------------------------------------------------
   // Place radial generators about the center specified in Boundary2D
   //------------------------------------------------------------------------
   void radialPoints(const unsigned nr)
   {
      mGenerators.clear();
      POLY_ASSERT( Dimension == 2 );
      RealType maxDistance;
      mBoundary.getBoundingRadius( maxDistance );
      POLY_ASSERT( maxDistance > 0 );
      
      RealType dRadius = maxDistance/nr;
      for( unsigned i = 0; i != nr; ++i ){
         RealType rad = (i+0.5)*dRadius;

         // This is supposed to befloor(2*pi*i), however 2*floor(pi)*i=6*i 
         // is found to work better
         unsigned nArcs = 6*i;
         RealType* pos;
         for( unsigned j = 0; j != nArcs; ++j ){
            RealType theta = 2*M_PI*j/nArcs;
            pos[0] = mBoundary.mCenter[0] + rad*cos(theta);
            pos[1] = mBoundary.mCenter[1] + rad*sin(theta);
            if( mBoundary.testInside( pos ) ){
               mGenerators.push_back( pos[0] );
               mGenerators.push_back( pos[1] );
            }
         }
      }
      nPoints = mGenerators.size()/Dimension;
   }
   
   
   //------------------------------------------------------------------------
   // add a point to the generator set
   //------------------------------------------------------------------------
   void addGenerator(RealType& pos)
   {
      POLY_ASSERT( mBoundary.testInside( pos ) );
      for( unsigned n=0; n < Dimension; ++n ){
         mGenerators.push_back( pos[n] );
      }
   }

   //------------------------------------------------------------------------
   // for 2D problems
   //------------------------------------------------------------------------
   void cartesian2D(const unsigned nx, const unsigned ny)
   {
      RealType x, y;
      RealType dx = (mBoundary.mHigh[0] - mBoundary.mLow[0]) / nx;
      RealType dy = (mBoundary.mHigh[1] - mBoundary.mLow[1]) / ny;
      RealType pos[2];
      for (unsigned iy = 0; iy != ny; ++iy){
         y = mBoundary.mLow[1] + (iy + 0.5)*dy;
         for (unsigned ix = 0; ix != nx; ++ix){
            x = mBoundary.mLow[0] + (ix + 0.5)*dx;
            pos[0] = x;   pos[1] = y;
            if( mBoundary.testInside( pos ) ){
               mGenerators.push_back( x );
               mGenerators.push_back( y );
            }
         }
      }
   }

   //------------------------------------------------------------------------
   // for 3D problems
   //------------------------------------------------------------------------
   void cartesian3D(const unsigned nx, const unsigned ny, const unsigned nz)
   {
      RealType x, y, z;
      RealType dx = (mBoundary.mHigh[0] - mBoundary.mLow[0]) / nx;
      RealType dy = (mBoundary.mHigh[1] - mBoundary.mLow[1]) / ny;
      RealType dz = (mBoundary.mHigh[2] - mBoundary.mLow[2]) / nz;
      for (unsigned iz = 0; iz != nz; ++iz){
         z = mBoundary.mLow[2] + (iz + 0.5)*dz;
         for (unsigned iy = 0; iy != ny; ++iy){
            y = mBoundary.mLow[1] + (iy + 0.5)*dy;
            for (unsigned ix = 0; ix != nx; ++ix){
               x = mBoundary.mLow[0] + (ix + 0.5)*dx;
               mGenerators.push_back( x );
               mGenerators.push_back( y );
               mGenerators.push_back( z );
            }
         }
      }
   }

   //------------------------------------------------------------------------
   // perturb the locations of the generators by +/- epsilon/2 in each direction
   //------------------------------------------------------------------------
   void perturb(RealType epsilon)
   {
      for (unsigned i = 0; i < mGenerators.size()/Dimension; ++i){
         for (unsigned n = 0; n < Dimension; ++n){
            mGenerators[n*i] += epsilon*( RealType(::random())/RAND_MAX - 0.5 );
         }
      }
   }
};


#endif
