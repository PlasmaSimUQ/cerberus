#include "MFP_interp.H"
#include "MFP_global.H"
#include "MFP_source.H"

#include <fstream>
#include <sstream>

#include "delaunator.hpp"


//#define DBL_EPSILON 2.2204460492503131e-16
#define EPS 100 * DBL_EPSILON

Interp2D::Interp2D() {}

void Interp2D::init(const std::vector<double> xs, const std::vector<double> ys,
                    const std::vector<double> zs) {
  std::vector<double> coords;

  vals = zs;

  for (size_t i = 0; i < xs.size(); i++) {
    coords.push_back(xs[i]);
    coords.push_back(ys[i]);
  }

  delaunator::Delaunator d(coords);

  d_triangles = d.triangles;
  d_coords = d.coords;
  d_halfedges = d.halfedges;

  is_valid = true;
}

Interp2D::Interp2D(const std::vector<double> xs, const std::vector<double> ys,
                   const std::vector<double> zs) {
  init(xs, ys, zs);
}

Interp2D::Interp2D(const std::string ifile) {
  std::vector<double> coords;
  std::vector<double> xs;
  std::vector<double> ys;
  std::vector<double> zs;

  std::ifstream input_file;
  input_file.open(ifile);

  std::string line;
  std::string x_st, y_st, z_st;

  while (std::getline(input_file, line, '\n')) {
    std::istringstream line_stream(line);
    getline(line_stream, x_st, ',');
    xs.push_back(stod(x_st));

    getline(line_stream, y_st, ',');
    ys.push_back(stod(y_st));

    getline(line_stream, z_st, ',');
    zs.push_back(stod(z_st));
  }

  x_min = xs[0];
  y_min = ys[0];
  x_max = xs[xs.size() - 1];
  y_max = ys[ys.size() - 1];

  input_file.close();

  init(xs, ys, zs);
}

std::array<size_t, 3> Interp2D::edges_of_triangle(size_t t) const {
  return {3 * t, 3 * t + 1, 3 * t + 2};
}

size_t Interp2D::triangle_of_edge(size_t e) const { return floor(e / 3); }

std::vector<size_t> Interp2D::triangles_adjacent_to_triangle(size_t t) const {
  std::vector<size_t> adjacent_triangles;
  for (auto e : edges_of_triangle(t)) {
    std::size_t opposite = d_halfedges[e];
    if (opposite >= 0) {
      adjacent_triangles.push_back(triangle_of_edge(e));
    }
  }
  return adjacent_triangles;
}

// double Interp2D::distance_along_z(double x1, double y1, double z1, double x2,
//                                   double y2, double z2, double x3, double y3,
//                                   double z3, double x, double y, bool broad)
//                                   const {
//   double detT = ((y2 - y3) * (x1 - x3)) + ((x3 - x2) * (y1 - y3));
//   double z = std::numeric_limits<double>::infinity();
//   double eps_min = 100*DBL_EPSILON;
//   double eps_max = 100*DBL_EPSILON;
//   if (broad) eps_min = epsilon2;

//   if (std::abs(detT) >= 1e-16) {
//     double lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT;
//     double lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT;

//     if (lambda1 >= 0 - eps_min && lambda1 <= 1 + eps_max && lambda2 >= 0 -
//     eps_min &&
//         lambda2 <= 1 + eps_max) {
//       // printf("2\n");
//       double sum = lambda1 + lambda2;
//       if (sum <= 1) {
//         // printf("3\n");
//         double lambda3 = 1 - sum;
//         z = (lambda1 * z1) + (lambda2 * z2) + (lambda3 * z3);
//       }
//     }
//   }
//   return z;
// }

double Interp2D::distance_along_z(double x1, double y1, double z1, double x2,
                                  double y2, double z2, double x3, double y3,
                                  double z3, double x, double y) const {
  double detT = ((y2 - y3) * (x1 - x3)) + ((x3 - x2) * (y1 - y3));
  double z = std::numeric_limits<double>::infinity();

  if (std::abs(detT) >= 1e-16) {
    double lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / detT;
    double lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / detT;
    double lambda3 = 1 - lambda2 - lambda1;

    if (lambda1 >= 0 - EPS && lambda1 <= 1 + EPS && lambda2 >= 0 - EPS &&
        lambda2 <= 1 + EPS && lambda3 >= 0 - EPS && lambda3 <= 1 + EPS) {
      // printf("3\n");
      z = (lambda1 * z1) + (lambda2 * z2) + (lambda3 * z3);
    }
  }
  return z;
}

double Interp2D::point_in_triel(size_t triel, double x, double y) const {
  double d;
  int i, j, ztest = 0;
  double xi, yi;
  double xj, yj;

  for (int i = 2, j = i - 1; i > -1; i--, j--) {
    if (j == -1) j = 3;

    auto ijpoint = d_triangles[triel + j];
    auto iipoint = d_triangles[triel + i];

    xj = d_coords[2 * ijpoint];
    yj = d_coords[2 * ijpoint + 1];

    xi = d_coords[2 * iipoint];
    yi = d_coords[2 * iipoint + 1];

    d = (xj - xi) * (y - yi) - (yj - yi) * (x - xi);
    if (d < 0.0) return -1;
    if (d == 0.0) ztest = 1;
  }
  return (ztest ? 0 : 1);
}

double Interp2D::interpolate(double x, double y) const {
  size_t triel;
  int triel_found = 0;
  double sum, ans = 0.0;
  std::vector<std::size_t> neighbours;

  // use nearest co-ordinate for out of bounds
  if (x < x_min) {
    x = x_min;
  } else if (x > x_max) {
    x = x_max;
  }

  if (y < y_min) {
    y = y_min;
  } else if (y > y_max) {
    y = y_max;
  }

  for (triel = 0; triel < d_triangles.size(); triel += 3) {
    auto ai = d_triangles[triel];
    auto bi = d_triangles[triel + 1];
    auto ci = d_triangles[triel + 2];

    double ax = d_coords[2 * ai];
    double ay = d_coords[2 * ai + 1];
    double az = vals[ai];

    double bx = d_coords[2 * bi];
    double by = d_coords[2 * bi + 1];
    double bz = vals[bi];

    double cx = d_coords[2 * ci];
    double cy = d_coords[2 * ci + 1];
    double cz = vals[ci];

    double z1 = distance_along_z(ax, ay, az, bx, by, bz, cx, cy, cz, x, y);

    if (z1 != std::numeric_limits<double>::infinity()) {
      return z1;
    }
  }
  std::cout << triel / 3 << "\n";
  return std::numeric_limits<double>::infinity();
}
