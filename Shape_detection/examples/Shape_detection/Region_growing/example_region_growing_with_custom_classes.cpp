// STL includes.
#include <vector>
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>

namespace SD = CGAL::Shape_detection;
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT      = typename Kernel::FT;
using Point_3 = typename Kernel::Point_3;

// Custom Connectivity and Conditions classes for region growing.
namespace custom {

  // A 3D sphere of the fixed radius 
  // that stores its center point and 
  // indices of all adjacent to it other spheres.
  struct Sphere {
    const FT radius = FT(1);

    Point_3 center;
    std::vector<std::size_t> neighbors;
  };

  // A range of spheres.
  using Spheres = std::vector<Sphere>;

  // Connectivity class, where the function get_neighbors()
  // simply returns indices of all neighbors stored in
  // the sphere struct above. This is the only function that 
  // we have to define.
  class Connectivity {
    const Spheres &m_spheres;

  public:
    Connectivity(Spheres& spheres) : 
    m_spheres(spheres) 
    { }

    template<typename OutputIterator>
    void get_neighbors(
      const std::size_t query_index, 
      OutputIterator neighbors) const {

      for (std::size_t i = 0; i < m_spheres[query_index].neighbors.size(); ++i)
        *(neighbors++) = m_spheres[query_index].neighbors[i];
    }
  };

  // Conditions class, where the function belongs_to_region() checks
  // a very specific condition that the first and second spheres in the
  // range are in fact neighbors; is_valid_region() function always 
  // returns true after the first call to the update() function.
  // These are the only functions that we have to define.
  class Conditions {
    bool m_is_valid = false;

  public:
    Conditions() { }

    bool belongs_to_region(
      const std::size_t query_index,
      const std::vector<std::size_t>& region) const {

      if (region.size() == 0) 
        return false;

      const std::size_t index = region[0]; 
      return (index == 0 && query_index == 1);
    }

    inline bool is_valid_region(const std::vector<std::size_t>&) const {
      return m_is_valid;
    }

    void update(const std::vector<std::size_t>&) {
      m_is_valid = true;
    }
  };
}

// Type declarations.
using Spheres      = custom::Spheres;
using Connectivity = custom::Connectivity;
using Conditions   = custom::Conditions;

using Region_growing = SD::Region_growing<Spheres, Connectivity, Conditions>;

int main(int argc, char *argv[]) { 
  
  std::cout << std::endl << 
    "region_growing_with_custom_classes example started" 
  << std::endl << std::endl;

  // Define a range of spheres, where the first two spheres form
  // the first region, while the thrid sphere forms the second region.
  Spheres spheres(3);

  // Region 1.
  spheres[0].center = Point_3(FT(0), FT(0), FT(0));
  spheres[0].neighbors.resize(1, 1);

  spheres[1].center = Point_3(FT(1), FT(1), FT(1));
  spheres[1].neighbors.resize(1, 0);

  // Region 2.
  spheres[2].center = Point_3(FT(4), FT(4), FT(4));

  // Create instances of the classes Connectivity and Conditions.
  Connectivity connectivity = Connectivity(spheres);
  Conditions   conditions   = Conditions();

  // Create an instance of the region growing class.
  Region_growing region_growing(spheres, connectivity, conditions);

  // Run the algorithm.
  region_growing.detect();

  // Print the number of found regions. It must be two regions.
  std::cout << "* " << 
  region_growing.number_of_regions() << 
    " regions have been found among 3 spheres" 
  << std::endl;
  
  std::cout << std::endl << 
    "region_growing_with_custom_classes example finished" 
  << std::endl << std::endl;

  return EXIT_SUCCESS;
}