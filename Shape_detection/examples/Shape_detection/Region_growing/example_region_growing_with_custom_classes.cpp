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

namespace custom {

  struct Sphere {
    const FT radius = FT(1);

    Point_3 center;
    std::vector<std::size_t> neighbors;
  };

  using Spheres = std::vector<Sphere>;

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

  class Conditions {
    bool m_is_valid = false;

  public:
    Conditions() { }

    template<typename ItemRange>
    bool belongs_to_region(
      const std::size_t query_index,
      const ItemRange& region) const {

      if (region.size() == 0) 
        return false;

      const std::size_t index = region[0]; 
      return (index == 0 && query_index == 1);
    }

    template<typename ItemRange>
    inline bool is_valid_region(const ItemRange&) const {
      return m_is_valid;
    }

    template<typename ItemRange>
    void update(const ItemRange&) {
      m_is_valid = true;
    }
  };
}

using Spheres      = custom::Spheres;
using Connectivity = custom::Connectivity;
using Conditions   = custom::Conditions;

using Region_growing = SD::Region_growing<Spheres, Connectivity, Conditions>;

int main(int argc, char *argv[]) { 
  
  std::cout << std::endl << 
    "region_growing_with_custom_classes example started" 
  << std::endl << std::endl;

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

  // Print the number of found regions.
  std::cout << "* " << 
  region_growing.number_of_regions() << 
    " regions have been found among 3 spheres" 
  << std::endl;
  
  std::cout << std::endl << 
    "region_growing_with_custom_classes example finished" 
  << std::endl << std::endl;

  return EXIT_SUCCESS;
}