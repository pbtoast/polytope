// Try tessellating a simple lattice of generators in a box.

#include <iostream>
#include <vector>
#include "polytope.hh"

// We really want cassert for the test!
#include <cassert>

using namespace std;

int main() {
  const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
  const double x2 = 100.0, y2 = 100.0, z2 = 100.0;

  // Try tessellating increasing numbers of generators.
  for (unsigned nx = 2; nx != 11; ++nx) {
    cout << "Testing nx=" << nx << endl;

    // Create the generators.
    vector<double> generators;
    const double dx = (x2 - x1)/nx, dy = (y2 - y1)/nx, dz = (z2 - z1)/nx;
    unsigned ix, iy, iz;
    double xi, yi, zi;
    for (iz = 0; iz != nx; ++iz) {
      zi = std::max(z1, std::min(z2, z1 + (iz + 0.5)*dz));
      for (iy = 0; iy != nx; ++iy) {
        yi = std::max(y1, std::min(y2, y1 + (iy + 0.5)*dy));
        for (ix = 0; ix != nx; ++ix) {
          xi = std::max(x1, std::min(x2, x1 + (ix + 0.5)*dx));
          generators.push_back(xi);
          generators.push_back(yi);
          generators.push_back(zi);
        }
      }
    }

    // Create the tessellation.
    double xmin[3] = { x1, y1, z1 };
    double xmax[3] = { x2, y2, z2 };
    polytope::Tessellation<3, double> mesh;
    polytope::VoroPP_3d<double> voro;
    voro.tessellate(generators, xmin, xmax, mesh);

    // Spew the mesh statistics.
    cout << "   num mesh nodes : " << mesh.nodes.size()/3 << endl;
    cout << "   num mesh cells : " << mesh.cells.size() << endl;
    cout << "   num mesh faces : " << mesh.faces.size() << endl;
    cout << "Node positions: " << endl;
    for (unsigned i = 0; i != mesh.nodes.size()/3; ++i) {
      cout << "   Node " << i << " @ (" << mesh.nodes[3*i] << " " << mesh.nodes[3*i + 1] << " " << mesh.nodes[3*i + 2] << ")" << endl;
    }
    cout << "Face node sets: " << endl;
    for (unsigned i = 0; i != nx*nx*nx; ++i) {
      cout << "   FACES for mesh cell " << i << " :";
      for (unsigned j = 0; j != mesh.cells[i].size(); ++j) cout << " " << mesh.cells[i][j];
      cout << endl;
    }
    for (unsigned i = 0; i != mesh.faces.size(); ++i) {
      double xf = 0.0, yf = 0.0, zf = 0.0;
      cout << "   NODES for mesh face " << i << " :";
      for (unsigned j = 0; j != mesh.faces[i].size(); ++j) {
        unsigned k = mesh.faces[i][j];
        cout << " " << k;
        xf += mesh.nodes[3*k];
        yf += mesh.nodes[3*k + 1];
        zf += mesh.nodes[3*k + 2];
      }
      xf /= mesh.faces[i].size();
      yf /= mesh.faces[i].size();
      zf /= mesh.faces[i].size();
      cout << " @ (" << xf << " " << yf << " " << zf << ")"  << endl;
    }

    // Now do the checks.
    assert(mesh.nodes.size() == 3*(nx + 1)*(nx + 1));
    assert(mesh.cells.size() == nx*nx*nx);
    for (unsigned i = 0; i != nx*nx*nx; ++i) assert(mesh.cells[i].size() == 6);
    // assert(mesh.faces.size() == 2*nx*nx*(nx + 1));
  }

  cout << "PASS" << endl;
  return 0;
}
