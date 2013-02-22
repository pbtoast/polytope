#ifndef POLYTOPE_BOUNDARY2D_HH
#define POLYTOPE_BOUNDARY2D_HH

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <cstdlib>

#include "polytope.hh"
#include "polytope_test_utilities.hh"

// We use the Boost.Geometry library to handle polygon intersections and such.
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>

using namespace std;
using namespace polytope;

//------------------------------------------------------------------------
template<typename RealType>
class Boundary2D
{
public:
   // ------------- Some handy typedefs for Boost.Geometry --------------- //
   typedef boost::geometry::model::point<RealType, 2, boost::geometry::cs::cartesian>
      BGpoint;
   typedef boost::geometry::model::polygon<BGpoint,false> 
      BGpolygon;
   typedef boost::geometry::model::ring<BGpoint,false>
      BGring;

   // -------------- Public member variables and routines ---------------- //
   
   // Number of dimensions
   unsigned Dimension;
   
   // Piecewise linear construct to define the boundary facets + holes
   PLC<2, RealType> mPLC;
   // Vector of generators to define the boundary
   std::vector<RealType> mPLCpoints;
   // Prescribed center of the computational domain
   std::vector<RealType> mCenter;
   // Ranges of bounding box
   RealType mLow[2], mHigh[2], mArea;
   
   BGpolygon mBGboundary;
   
   
   // Define enum to keep track fo the type of boundary called for
   enum BoundaryType{
      box                = 0,
      circle             = 1,
      torus              = 2,
      mwithholes         = 3,
      funkystar          = 4,
      circlewithstarhole = 5,
      cardioid           = 6,
      burninator         = 7
   };
   
   // Boundary type
   mutable BoundaryType mType;

   
   //------------------------------------------------------------------------
   // Constructor, destructor
   //------------------------------------------------------------------------
   Boundary2D():
      Dimension(2),
      mType(box)
   {
      mCenter = std::vector<RealType>(Dimension,0);
   };

   ~Boundary2D() {};

   void clear()
   {
      mPLC.clear();
      mPLCpoints.clear();
      mCenter(0);
      mLow(0);
      mHigh(0);
      boost::geometry::clear(mBGboundary);
   }

   void finalize()
   {
      getBoundingBox();
      boostMyBoundary();
      mArea = boost::geometry::area( mBGboundary );
   }
   
   //------------------------------------------------------------------------
   // computeBoundary
   //------------------------------------------------------------------------
   void computeDefaultBoundary(int bType)
   {
      //mCenter = std::vector<RealType>(Dimension,0);
      
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
      case burninator:
         this->trogdor();
         break;
      }

      this->finalize();
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
      
      mPLCpoints.push_back( x1 );   mPLCpoints.push_back( y1 );
      mPLCpoints.push_back( x2 );   mPLCpoints.push_back( y1 );
      mPLCpoints.push_back( x2 );   mPLCpoints.push_back( y2 );
      mPLCpoints.push_back( x1 );   mPLCpoints.push_back( y2 );
      
      mPLC.facets.resize(4);
      for (unsigned f = 0; f < 4; ++f)
      {
         mPLC.facets[f].resize(2);
         mPLC.facets[f][0] = f;
         mPLC.facets[f][1] = (f+1) % 4;
      }
      mType = box;
      this->finalize();
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
         mPLCpoints.push_back(x);
         mPLCpoints.push_back(y);
      }

      // Facets.
      mPLC.facets.resize(Nb); 
      for (unsigned f = 0; f < Nb; ++f)
      {
         unsigned fBegin =  mPLCpoints.size()/2 - Nb + f;
         unsigned fEnd   = (mPLCpoints.size()/2 - Nb + f + 1) % Nb;
         mPLC.facets[f].resize(2);
         mPLC.facets[f][0] = fBegin;
         mPLC.facets[f][1] = fEnd;
      }
      mType = circle;
      this->finalize();
   }

   
   //------------------------------------------------------------------------
   // donut
   // Unit circle with a hole of prescribed radius
   //------------------------------------------------------------------------
   void donut( RealType innerRadius=0.25 )
   {
      POLY_ASSERT2( innerRadius > 0, "Must provide a positive inner radius" );
      POLY_ASSERT2( innerRadius < 1, "Inner radius may not exceed outer (unit) radius" );

      // The outer circle
      unitCircle();
            
      // Inner circle.
      unsigned Nb = mPLCpoints.size()/2;
      unsigned Nh = Nb;
      for (unsigned b = 0; b < Nh; ++b)
      {
         RealType theta = 2.0*M_PI*(1.0 - RealType(b)/RealType(Nb+1));
         RealType x = mCenter[0] + innerRadius*cos(theta);
         RealType y = mCenter[1] + innerRadius*sin(theta);
         mPLCpoints.push_back(x);
         mPLCpoints.push_back(y);
      }

      // Facets on the inner circle
      mPLC.holes = std::vector< std::vector< std::vector<int> > >(1);
      mPLC.holes[0].resize(Nh);
      for (unsigned f = 0; f < Nh; ++f)
      {
         unsigned fBegin = mPLCpoints.size()/2 - Nh + f;
         unsigned fEnd   = mPLCpoints.size()/2 - Nh + (f + 1) % Nh;
         mPLC.holes[0][f].resize(2);
         mPLC.holes[0][f][0] = fBegin;
         mPLC.holes[0][f][1] = fEnd;
      }
      mType = torus;
      this->finalize();
   }

   
   //------------------------------------------------------------------------
   // MWithHoles
   // M-shaped domain with two square holes
   //------------------------------------------------------------------------
   void MWithHoles()
   {
      // Outer boundary of the M-shape
      mPLCpoints.push_back(0.0); mPLCpoints.push_back(0.0);
      mPLCpoints.push_back(2.0); mPLCpoints.push_back(0.0);
      mPLCpoints.push_back(2.0); mPLCpoints.push_back(2.0);
      mPLCpoints.push_back(1.0); mPLCpoints.push_back(1.0);
      mPLCpoints.push_back(0.0); mPLCpoints.push_back(2.0);

      int nSides = mPLCpoints.size()/2;
      mPLC.facets.resize( nSides, std::vector<int>(2) );
      for (unsigned i = 0; i != nSides; ++i )
      {
         mPLC.facets[i][0] = i;
         mPLC.facets[i][1] = (i+1) % nSides;
      }

      int nHoles = 2;
      
      // Square hole #1
      mPLCpoints.push_back(0.25); mPLCpoints.push_back(0.25);
      mPLCpoints.push_back(0.25); mPLCpoints.push_back(0.75);
      mPLCpoints.push_back(0.75); mPLCpoints.push_back(0.75);
      mPLCpoints.push_back(0.75); mPLCpoints.push_back(0.25);
      
      // Square hole #2
      mPLCpoints.push_back(1.25); mPLCpoints.push_back(0.25);
      mPLCpoints.push_back(1.25); mPLCpoints.push_back(0.75);
      mPLCpoints.push_back(1.75); mPLCpoints.push_back(0.75);
      mPLCpoints.push_back(1.75); mPLCpoints.push_back(0.25);

      mPLC.holes.resize(nHoles, std::vector<std::vector<int> >(4, std::vector<int>(nHoles)));
      for (unsigned i = 0; i != 4; ++i)
      {
         mPLC.holes[0][i][0] = nSides + i;
         mPLC.holes[0][i][1] = nSides + ((i+1) % 4);
         mPLC.holes[1][i][0] = nSides + 4 + i;
         mPLC.holes[1][i][1] = nSides + 4 + ((i+1) % 4);
      }
      mType = mwithholes;
      this->finalize();
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
         mPLCpoints.push_back( x );
         mPLCpoints.push_back( y );
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
      this->finalize();
   }

   
   //------------------------------------------------------------------------
   // circleWithStarHole
   // Unit circle with a hole shaped like a regular n-pointed star
   //------------------------------------------------------------------------
   void circleWithStarHole( int nPoints = 5 )
   {
      // The outer boundary
      unitCircle();
      
      RealType theta0 = 2*M_PI/nPoints;
      RealType outerRadius = 0.75;
      RealType innerRadius = outerRadius*( sin(theta0/4.0) / sin(3*theta0/4.0) );
         
      RealType theta;
      for (unsigned p = 0; p < nPoints; ++p )
      {
         // For the pointy bits of the star
         theta = M_PI/2 - p*theta0;
         mPLCpoints.push_back( mCenter[0] + outerRadius*cos(theta) );
         mPLCpoints.push_back( mCenter[1] + outerRadius*sin(theta) );
         
         // For the concave bits of the star
         theta = M_PI/2 - p*theta0 - theta0/2.0;
         mPLCpoints.push_back( mCenter[0] + innerRadius*cos(theta) );
         mPLCpoints.push_back( mCenter[1] + innerRadius*sin(theta) );
      }
      
      // Facets on the inner circle
      mPLC.holes = std::vector< std::vector< std::vector<int> > >(1);      
      mPLC.holes[0].resize(2*nPoints);
      for (unsigned f = 0; f < 2*nPoints; ++f)
      {
         unsigned fBegin = mPLCpoints.size()/2 - 2*nPoints + f;
         unsigned fEnd   = mPLCpoints.size()/2 - 2*nPoints + ((f + 1) % (2*nPoints));
         mPLC.holes[0][f].resize(2);
         mPLC.holes[0][f][0] = fBegin;
         mPLC.holes[0][f][1] = fEnd;
      }

      mType = circlewithstarhole;
      this->finalize();
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
         mPLCpoints.push_back(x);
         mPLCpoints.push_back(y);
      }

      // Facets.
      mPLC.facets.resize(Nb); 
      for (unsigned f = 0; f < Nb; ++f)
      {
         unsigned fBegin =  mPLCpoints.size()/2 - Nb + f;
         unsigned fEnd   = (mPLCpoints.size()/2 - Nb + f + 1) % Nb;
         mPLC.facets[f].resize(2);
         mPLC.facets[f][0] = fBegin;
         mPLC.facets[f][1] = fEnd;
      }
      mType = cardioid;
      this->finalize();
   }

   
   //------------------------------------------------------------------------
   // trogdor
   // No explanation necessary
   //------------------------------------------------------------------------
   void trogdor()
   {
     const unsigned nSides = 30;
     const double points[60] = {2, 9,
                                4, 8.9000000000000004,
                                5, 9.1999999999999993,
                                6.5, 8.8000000000000007,
                                7, 8,
                                6.5, 7,
                                5, 6.2999999999999998,
                                4, 5.5,
                                3.7000000000000002, 4.5999999999999996,
                                4, 3.5,
                                5.5, 2.7999999999999998,
                                6.7999999999999998, 3.3999999999999999,
                                5.5, 2.5,
                                4, 2.6000000000000001,
                                3, 3,
                                2.5, 4,
                                2.3999999999999999, 4.5,
                                2.5, 5.2999999999999998,
                                3, 5.9000000000000004,
                                4.5, 7,
                                4.9000000000000004, 7.2999999999999998,
                                5.0999999999999996, 7.7000000000000002,
                                5, 8,
                                4.5, 8,
                                3.5, 7.5,
                                2, 7,
                                2, 7.4000000000000004,
                                3.2999999999999998, 8,
                                1.75, 8,
                                1.75, 9.1999999999999993};
     
     mPLC.facets.resize( nSides, std::vector<int>(2) );
     for (unsigned i = 0; i != nSides; ++i){
        mPLCpoints.push_back(points[2*i  ]);
        mPLCpoints.push_back(points[2*i+1]);
        mPLC.facets[i][0] = i;
        mPLC.facets[i][1] = (i+1) % nSides;
     }
           
     mType = burninator;
     this->finalize();
   }

   
   // //------------------------------------------------------------------------
   // // testInside
   // // Tests if a given point (x,y) lies in the interior of the domain
   // //------------------------------------------------------------------------
   // bool testInside( RealType* pos )
   // {
   //    RealType x = pos[0];
   //    RealType y = pos[1];
   //    POLY_ASSERT( mPLCpoints.size() > 0 );
   //    //const unsigned nSides = mPLCpoints.size()/2;
   //    const unsigned nSides = mPLC.facets.size();
   //    unsigned j = nSides - 1;
   //    bool isInside = false;
   //    for (unsigned i = 0; i < nSides; ++i )
   //    {
   //       unsigned ix = 2*i  ,  jx = 2*j;
   //       unsigned iy = 2*i+1,  jy = 2*j+1;  
   //       if( (mPLCpoints[2*i+1] <  y && mPLCpoints[2*j+1] >= y ||
   //            mPLCpoints[2*j+1] <  y && mPLCpoints[2*i+1] >= y) &&
   //           (mPLCpoints[2*i  ] <= x || mPLCpoints[2*j  ] <= x) )
   //       {
   //          isInside ^= ( mPLCpoints[2*i] + ( y            - mPLCpoints[2*i+1] ) /
   //                                     ( mPLCpoints[2*j+1] - mPLCpoints[2*i+1] ) *
   //                                     ( mPLCpoints[2*j  ] - mPLCpoints[2*i  ] ) < x );
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
      POLY_ASSERT( mPLCpoints.size() > 0 );
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
   // set of points in mPLCpoints starting at index 'offset'
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
	 if( ((mPLCpoints[iy] <  y  &&  mPLCpoints[jy] >= y)  ||
	      (mPLCpoints[jy] <  y  &&  mPLCpoints[iy] >= y)) &&
	     (mPLCpoints[ix] <= x  ||  mPLCpoints[jx] <= x) )
	 {
	    isInside ^= ( mPLCpoints[ix] + ( y         - mPLCpoints[iy] ) /
                                      ( mPLCpoints[jy] - mPLCpoints[iy] ) *
                                      ( mPLCpoints[jx] - mPLCpoints[ix] ) < x );
	 }
	 j = i;
      }
      offset += nSides;
      return isInside;
   }




   //------------------------------------------------------------------------
   // getBoundingBox
   // Get the maximal L-infinity norm of the boundary generator set
   // NOTE: method is general to 2D or 3D. Dimension is set to 2 here
   //------------------------------------------------------------------------
   void getBoundingBox()
   {
      for (unsigned n = 0; n < Dimension; ++n){
         mLow[n]  =   numeric_limits<RealType>::max();
         mHigh[n] =  (numeric_limits<RealType>::is_signed ? -mLow[n] : 
                      numeric_limits<RealType>::min());
      }
      for (unsigned i = 0; i < mPLCpoints.size()/Dimension; ++i ){
         for (unsigned n = 0; n < Dimension; ++n ){
            mLow [n] = min( mLow [n], mPLCpoints[Dimension*i+n] );
            mHigh[n] = max( mHigh[n], mPLCpoints[Dimension*i+n] );
         }
      }
   }


   //------------------------------------------------------------------------
   // getBoundingCircle
   // Get the maximal L-2 norm of the boundary generator set about center pt
   // NOTE: method is general to 2D or 3D. Dimension is set to 2 here
   //------------------------------------------------------------------------
   void getBoundingRadius(RealType& radius)
   {
      POLY_ASSERT( mCenter != 0 );
      radius = 0;
      for (unsigned i = 0; i < mPLCpoints.size()/Dimension; ++i ){
         RealType distance = 0;
         for (unsigned n = 0; n < Dimension; ++n ){
            distance += (mPLCpoints[Dimension*i+n] - mCenter[n]) *
               (mPLCpoints[Dimension*i+n] - mCenter[n]);
         }
         radius = max( radius, sqrt( distance ) );
      }
   }


   //------------------------------------------------------------------------
   // BGmyBoundary
   // Store the boundary info as a Boost.Geometry polygon
   //------------------------------------------------------------------------
   void boostMyBoundary()
   {
      boost::geometry::clear(mBGboundary);
      mBGboundary = makePolygon<RealType>( mPLC, mPLCpoints );
   }

};

#endif
