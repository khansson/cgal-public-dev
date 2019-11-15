#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Barycentric_coordinates_2/Harmonic_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/analytic_coordinates_2.h>

// Typedefs.
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = Kernel::FT;
using Point_2 = Kernel::Point_2;

using Points_2 = std::vector<Point_2>;

/*
using Harmonic = CGAL::Barycentric_coordinates::Harmonic_coordinates_2<Points_2, Kernel>;
using Policy   = CGAL::Barycentric_coordinates::Computation_policy; */

int main() {  

}