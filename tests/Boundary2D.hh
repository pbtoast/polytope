#ifndef POLYTOPE_BOUNDARY2D_HH
#define POLYTOPE_BOUNDARY2D_HH

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>

#include "polytope.hh"

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------
template<typename RealType>
class Boundary2D
{
public:
   // -------------- Public member variables and routines ---------------- //
   
   // Number of dimensions
   unsigned nDims;
   
   // Piecewise linear construct to define the boundary facets + holes
   PLC<2, RealType> mPLC;
   // Vector of generators to define the boundary
   vector<RealType> mGens;
   // Prescribed center of the computational domain
   RealType mCenter[2];
   // Ranges of bounding box
   RealType mLow[2], mHigh[2];

   
   // Define enum to keep track fo the type of boundary called for
   enum BoundaryType{
      box                = 0,
      circle             = 1,
      torus              = 2,
      mwithholes         = 3,
      funkystar          = 4,
      circlewithstarhole = 5,
      cardioid           = 6,
      cardioidwithhole   = 7,
   };
   
   // Boundary type
   mutable BoundaryType mType;

   
   //------------------------------------------------------------------------
   // Constructor, destructor
   //------------------------------------------------------------------------
   Boundary2D():
      nDims(2),
      mLow(),
      mHigh(),
      mCenter(),
      mType(box)
   {
   };

   ~Boundary2D() {};

   void clear()
   {
      mPLC.clear();
      mGens.clear();
      mCenter(0);
      mLow(0);
      mHigh(0);
   }
   
   //------------------------------------------------------------------------
   // computeBoundary
   //------------------------------------------------------------------------
   void computeDefaultBoundary(int bType)
   {
      switch(bType){
      case box:
         this->unitSquare();
         break;
      case circle:
         this->unitCircle();
         break;
      case torus:
         this->donut();
         break;
      case mwithholes:
         this->MWithHoles();
         break;
      case funkystar:
         this->funkyStar();
         break;
      case circlewithstarhole:
         this->circleWithStarHole();
         break;
      case cardioid:
         this->cardioidBoundary();
         break;
      case cardioidwithhole:
         this->cardioidWithHole();
         break;
      }
      getBoundingBox();
   }


   //------------------------------------------------------------------------
   // unitSquare
   //------------------------------------------------------------------------
   void unitSquare()
   {
      const RealType x1 = mCenter[0] - 0.5;
      const RealType x2 = mCenter[0] + 0.5;
      const RealType y1 = mCenter[1] - 0.5;
      const RealType y2 = mCenter[1] + 0.5;
      
      mGens.push_back( x1 );   mGens.push_back( y1 );
      mGens.push_back( x2 );   mGens.push_back( y1 );
      mGens.push_back( x2 );   mGens.push_back( y2 );
      mGens.push_back( x1 );   mGens.push_back( y2 );
      
      mPLC.facets.resize(4);
      for (unsigned f = 0; f < 4; ++f)
      {
         mPLC.facets[f].resize(2);
         mPLC.facets[f][0] = f;
         mPLC.facets[f][1] = (f+1) % 4;
      }
      mType = box;
   }
   
   
   //------------------------------------------------------------------------
   // unitCircle
   //------------------------------------------------------------------------
   void unitCircle()
   {
      // Boundary generators.
      unsigned Nb = 90; // 4-degree resolution.
      for (unsigned b = 0; b < Nb; ++b)
      {
         RealType theta = 2.0*M_PI*b/(Nb+1);
         RealType x = mCenter[0] + cos(theta);
         RealType y = mCenter[1] + sin(theta);
         mGens.push_back(x);
         mGens.push_back(y);
      }

      // Facets.
      mPLC.facets.resize(Nb); 
      for (unsigned f = 0; f < Nb; ++f)
      {
         unsigned fBegin =  mGens.size()/2 - Nb + f;
         unsigned fEnd   = (mGens.size()/2 - Nb + f + 1) % Nb;
         mPLC.facets[f].resize(2);
         mPLC.facets[f][0] = fBegin;
         mPLC.facets[f][1] = fEnd;
      }
      mType = circle;
   }

   
   //------------------------------------------------------------------------
   // donut
   // Unit circle with a hole of prescribed radius
   //------------------------------------------------------------------------
   void donut( RealType innerRadius=0.25 )
   {
      POLY_ASSERT2( innerRadius > 0, "Must provide a positive inner radius" );
      POLY_ASSERT2( innerRadius < 1, "Inner radius may not exceed outer (unit) radius" );

      mPLC.holes = vector< vector< vector<int> > >(1);
      // The outer circle
      unitCircle();
      
      // Inner circle.
      unsigned Nb = mGens.size()/2;
      unsigned Nh = Nb;
      for (unsigned b = 0; b < Nh; ++b)
      {
         RealType theta = 2.0*M_PI*(1.0 - RealType(b)/RealType(Nb+1));
         RealType x = mCenter[0] + innerRadius*cos(theta);
         RealType y = mCenter[1] + innerRadius*sin(theta);
         mGens.push_back(x);
         mGens.push_back(y);
      }

      // Facets on the inner circle
      mPLC.holes[0].resize(Nh);
      for (unsigned f = 0; f < Nh; ++f)
      {
         unsigned fBegin = mGens.size()/2 - Nh + f;
         unsigned fEnd   = mGens.size()/2 - Nh + (f + 1) % Nh;
         mPLC.holes[0][f].resize(2);
         mPLC.holes[0][f][0] = fBegin;
         mPLC.holes[0][f][1] = fEnd;
      }
      mType = torus;
   }

   
   //------------------------------------------------------------------------
   // MWithHoles
   // M-shaped domain with two square holes
   //------------------------------------------------------------------------
   void MWithHoles()
   {
      // Outer boundary of the M-shape
      mGens.push_back(0.0); mGens.push_back(0.0);
      mGens.push_back(2.0); mGens.push_back(0.0);
      mGens.push_back(2.0); mGens.push_back(2.0);
      mGens.push_back(1.0); mGens.push_back(1.0);
      mGens.push_back(0.0); mGens.push_back(2.0);

      int nSides = mGens.size()/2;
      mPLC.facets.resize( nSides, vector<int>(2) );
      for (unsigned i = 0; i != nSides; ++i )
      {
         mPLC.facets[i][0] = i;
         mPLC.facets[i][1] = (i+1) % nSides;
      }

      int nHoles = 2;
      
      // Square hole #1
      mGens.push_back(0.25); mGens.push_back(0.25);
      mGens.push_back(0.25); mGens.push_back(0.75);
      mGens.push_back(0.75); mGens.push_back(0.75);
      mGens.push_back(0.75); mGens.push_back(0.25);
      
      // Square hole #2
      mGens.push_back(1.25); mGens.push_back(0.25);
      mGens.push_back(1.25); mGens.push_back(0.75);
      mGens.push_back(1.75); mGens.push_back(0.75);
      mGens.push_back(1.75); mGens.push_back(0.25);

      mPLC.holes.resize(nHoles, vector<vector<int> >(4, vector<int>(nHoles)));
      for (unsigned i = 0; i != 4; ++i)
      {
         mPLC.holes[0][i][0] = nSides + i;
         mPLC.holes[0][i][1] = nSides + ((i+1) % 4);
         mPLC.holes[1][i][0] = nSides + 4 + i;
         mPLC.holes[1][i][1] = nSides + 4 + ((i+1) % 4);
      }
      mType = mwithholes;  
   }

   
   //------------------------------------------------------------------------
   // funkyStar
   // Star-shaped(-ish) region from Misha Shashkov's Voronoi test suite
   //------------------------------------------------------------------------
   void funkyStar()
   {
      // Get the boundary points, organized counterclockwise
      int Nsides = 10;
      for (int i = Nsides-1; i >= 0; --i )
      {
         RealType rad = 2.0 + 0.5*sin(12.0* M_PI*RealType(i)/(RealType(Nsides)-1.0));
         RealType x = rad * sin(M_PI*RealType(i)/(RealType(Nsides)-1.0));
         RealType y = rad * cos(M_PI*RealType(i)/(RealType(Nsides)-1.0));
         mGens.push_back( x );
         mGens.push_back( y );
      }

      // Connect the boundary facets
      mPLC.facets.resize(Nsides);
      for (unsigned f = 0; f < Nsides; ++f )
      {
         mPLC.facets[f].resize(2);
         mPLC.facets[f][0] = f;
         mPLC.facets[f][1] = (f+1) % Nsides;
      }
      mType = funkystar;
   }

   
   //------------------------------------------------------------------------
   // circleWithStarHole
   // Unit circle with a hole shaped like a regular n-pointed star
   //------------------------------------------------------------------------
   void circleWithStarHole( int nPoints = 5 )
   {
      // The outer boundary
      unitCircle();
      int Nb = mGens.size()/2;
      
      RealType theta0 = 2*M_PI/nPoints;
      RealType outerRadius = 0.75;
      RealType innerRadius = outerRadius*( sin(theta0/4.0) / sin(3*theta0/4.0) );
         
      RealType theta;
      for (unsigned p = 0; p < nPoints; ++p )
      {
         // For the pointy bits of the star
         theta = M_PI/2 - p*theta0;
         mGens.push_back( mCenter[0] + outerRadius*cos(theta) );
         mGens.push_back( mCenter[1] + outerRadius*sin(theta) );
         
         // For the concave bits of the star
         theta = M_PI/2 - p*theta0 - theta0/2.0;
         mGens.push_back( mCenter[0] + innerRadius*cos(theta) );
         mGens.push_back( mCenter[1] + innerRadius*sin(theta) );
      }
      
      // Facets on the inner circle
      mPLC.holes = vector< vector< vector<int> > >(1);      
      mPLC.holes[0].resize(2*nPoints);
      for (unsigned f = 0; f < 2*nPoints; ++f)
      {
         unsigned fBegin = mGens.size()/2 - 2*nPoints + f;
         unsigned fEnd   = mGens.size()/2 - 2*nPoints + ((f + 1) % (2*nPoints));
         mPLC.holes[0][f].resize(2);
         mPLC.holes[0][f][0] = fBegin;
         mPLC.holes[0][f][1] = fEnd;
      }

      mType = circlewithstarhole;
   }
   
   
   //------------------------------------------------------------------------
   // cardioid
   // Cardioid with parametric equations
   //   x(t) = z*cos(t) + cos(2*t)
   //   y(t) = z*sin(t) + sin(2*t)
   // for t in [0,2*pi] and coefficient z
   //
   // NOTE:   z=2 : defines the unit cardioid with cusp pointing inward
   //       0<z<2 : cusp pushes inward and creates a hole
   //         z>2 : cusp less pronounced, cardioid more circular, non-unit
   //------------------------------------------------------------------------
   void cardioidBoundary( RealType z = 2 )
   {
      POLY_ASSERT2( z > 0, "Must provide a positive coefficient for the cardioid" );

      // Boundary generators.
      unsigned Nb = 90; // 4-degree resolution.
      for (unsigned b = 0; b < Nb; ++b)
      {
         RealType theta = 2.0*M_PI*b/(Nb+1);
         RealType x = z*cos(theta) - cos(2*theta);
         RealType y = z*sin(theta) - sin(2*theta);
         mGens.push_back(x);
         mGens.push_back(y);
      }

      // Facets.
      mPLC.facets.resize(Nb); 
      for (unsigned f = 0; f < Nb; ++f)
      {
         unsigned fBegin =  mGens.size()/2 - Nb + f;
         unsigned fEnd   = (mGens.size()/2 - Nb + f + 1) % Nb;
         mPLC.facets[f].resize(2);
         mPLC.facets[f][0] = fBegin;
         mPLC.facets[f][1] = fEnd;
      }
      mType = cardioid;
   }

   
   //------------------------------------------------------------------------
   // cardioidWithHole
   // The z=1 cardioid above
   //------------------------------------------------------------------------
   void cardioidWithHole()
   {
      mType = cardioidwithhole;
      cardioidBoundary( 1.0 );
   }

   
   // //------------------------------------------------------------------------
   // // testInside
   // // Tests if a given point (x,y) lies in the interior of the domain
   // //------------------------------------------------------------------------
   // bool testInside( RealType* pos )
   // {
   //    RealType x = pos[0];
   //    RealType y = pos[1];
   //    POLY_ASSERT( mGens.size() > 0 );
   //    //const unsigned nSides = mGens.size()/2;
   //    const unsigned nSides = mPLC.facets.size();
   //    unsigned j = nSides - 1;
   //    bool isInside = false;
   //    for (unsigned i = 0; i < nSides; ++i )
   //    {
   //       unsigned ix = 2*i  ,  jx = 2*j;
   //       unsigned iy = 2*i+1,  jy = 2*j+1;  
   //       if( (mGens[2*i+1] <  y && mGens[2*j+1] >= y ||
   //            mGens[2*j+1] <  y && mGens[2*i+1] >= y) &&
   //           (mGens[2*i  ] <= x || mGens[2*j  ] <= x) )
   //       {
   //          isInside ^= ( mGens[2*i] + ( y            - mGens[2*i+1] ) /
   //                                     ( mGens[2*j+1] - mGens[2*i+1] ) *
   //                                     ( mGens[2*j  ] - mGens[2*i  ] ) < x );
   //       }
   //       j = i;
   //    }
   //    return isInside;
   // }

   



   //------------------------------------------------------------------------
   // testInside
   // Tests if a given point (x,y) lies in the interior of the domain.
   // If holes are defined in mPLC, then each hole is tested separately.
   // Returns true if the point is inside the facets but outside the holes.
   //------------------------------------------------------------------------
   bool testInside( RealType* pos )
   {
      RealType x = pos[0];
      RealType y = pos[1];
      unsigned offset = 0;
      POLY_ASSERT( mGens.size() > 0 );
      const unsigned nSides = mPLC.facets.size();
      bool isInside = inside(x,y,nSides,offset);
      
      for (unsigned hIt = 0; hIt != mPLC.holes.size(); ++hIt )
      {
         isInside ^= inside(x,y,mPLC.holes[hIt].size(),offset);
      }
      return isInside;
   }


   //------------------------------------------------------------------------
   // inside
   // Tests if (x,y) is inside the nSide-sided polygon defined by the ordered
   // set of points in mGen starting at index 'offset'
   //------------------------------------------------------------------------
   bool inside( const RealType x, const RealType y, 
                const unsigned nSides ,unsigned& offset )
   {
      unsigned j = nSides - 1;
      bool isInside = false;
      for (unsigned i = 0; i < nSides; ++i )
      {
	 unsigned ix = 2*(i+offset),   iy = 2*(i+offset)+1;
	 unsigned jx = 2*(j+offset),   jy = 2*(j+offset)+1;
	 if( (mGens[iy] <  y  &&  mGens[jy] >= y  ||
	      mGens[jy] <  y  &&  mGens[iy] >= y) &&
	     (mGens[ix] <= x  ||  mGens[jx] <= x) )
	 {
	    isInside ^= ( mGens[ix] + ( y         - mGens[iy] ) /
                                      ( mGens[jy] - mGens[iy] ) *
                                      ( mGens[jx] - mGens[ix] ) < x );
	 }
	 j = i;
      }
      offset += nSides;
      return isInside;
   }




   //------------------------------------------------------------------------
   // getBoundingBox
   // Get the maximal L-infinity norm of the boundary generator set
   // NOTE: method is general to 2D or 3D. nDims is set to 2 here
   //------------------------------------------------------------------------
   void getBoundingBox()
   {
      for (unsigned i = 0; i < mGens.size()/nDims; ++i ){
         for (unsigned n = 0; n < nDims; ++n ){
            mLow [n] = min( mLow [n], mGens[nDims*i+n] );
            mHigh[n] = max( mHigh[n], mGens[nDims*i+n] );
         }
      }
   }


   //------------------------------------------------------------------------
   // getBoundingCircle
   // Get the maximal L-2 norm of the boundary generator set about center pt
   // NOTE: method is general to 2D or 3D. nDims is set to 2 here
   //------------------------------------------------------------------------
   void getBoundingRadius(RealType& radius)
   {
      POLY_ASSERT( mCenter != 0 );
      radius = 0;
      for (unsigned i = 0; i < mGens.size()/nDims; ++i ){
         RealType distance = 0;
         for (unsigned n = 0; n < nDims; ++n ){
            distance += (mGens[nDims*i+n] - mCenter[n]) *
               (mGens[nDims*i+n] - mCenter[n]);
         }
         radius = max( radius, sqrt( distance ) );
      }
   } 
};

#endif
