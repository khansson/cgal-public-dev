// STL includes.
#include <map>
#include <vector>
#include <string>
#include <iostream>

// CGAL includes.
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>

namespace SD = CGAL::Shape_detection;

// Custom Connectivity, Conditions, and Seed_map classes for region growing.
namespace custom {

  // An object that stores indices of all adjacent to it other objects.
  struct Object {
    std::vector<std::size_t> neighbors;
  };

  // A range of objects.
  using Objects = std::vector<Object>;

  // Connectivity class, where the function neighbors()
  // simply returns indices of all neighbors stored in
  // the object struct above. This is the only function that 
  // we have to define.
  class Connectivity {
    const Objects& m_objects;

  public:
    Connectivity(Objects& objects) : 
    m_objects(objects) 
    { }

    void neighbors(
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) const {
      
      neighbors = m_objects[query_index].neighbors;
    }
  };

  // Conditions class, where the function belongs_to_region() checks
  // a very specific condition that the first and second objects in the
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

      if (index == 0 && query_index == 1) return true;
      if (query_index == 0 && index == 1) return true;
      
      return false;
    }

    inline bool is_valid_region(const std::vector<std::size_t>&) const {
      return m_is_valid;
    }

    void update(const std::vector<std::size_t>&) {
      m_is_valid = true;
    }
  };

  // A seed map with the minimum requirements that is using the m_objects_map 
  // to define the seeding order of objects.
  class Seed_map {
                        
  public:
    using key_type = std::size_t;
    using value_type = std::size_t;
    using category = boost::lvalue_property_map_tag;

    Seed_map(const std::map<std::size_t, std::size_t>& objects_map) : 
    m_objects_map(objects_map) 
    { }

    value_type operator[](const key_type key) const { 
      return m_objects_map.at(key);
    }

    friend value_type get(
      const Seed_map& seed_map, 
      const key_type key) { 
      
      return seed_map[key];
    }

  private:
    const std::map<std::size_t, std::size_t>& m_objects_map;
  };
}

// Type declarations.
using Objects      = custom::Objects;
using Connectivity = custom::Connectivity;
using Conditions   = custom::Conditions;
using Seed_map     = custom::Seed_map;

using Region_growing = 
SD::Region_growing<Objects, Connectivity, Conditions, Seed_map>;

int main(int argc, char *argv[]) { 
  
  std::cout << std::endl << 
    "region_growing_with_custom_classes example started" 
  << std::endl << std::endl;

  // Define a range of objects, where the first two objects form
  // the first region, while the third object forms the second region.
  Objects objects(4);

  // Region 1.
  objects[0].neighbors.resize(1, 1);
  objects[1].neighbors.resize(1, 0);

  // Region 2.
  objects[2].neighbors.resize(0);

  // Extra object to skip.
  objects[3].neighbors.resize(0);

  // Create instances of the classes Connectivity and Conditions.
  Connectivity connectivity = Connectivity(objects);
  Conditions   conditions   = Conditions();

  // Create a seed map.
  std::map<std::size_t, std::size_t> objects_map;
  objects_map[0] = 1; // the order is swapped with the next object
  objects_map[1] = 0;
  objects_map[2] = 2; // the default order
  objects_map[3] = std::size_t(-1); // skip this object
  const Seed_map seed_map(objects_map);

  // Create an instance of the region growing class.
  Region_growing region_growing(objects, connectivity, conditions, seed_map);

  // Run the algorithm.
  region_growing.detect();

  // Print the number of found regions. It must be two regions.
  std::cout << "* " << 
  region_growing.number_of_regions() << 
    " regions have been found among " << objects.size() <<  " objects" 
  << std::endl;
  
  // Release all internal memory that is used to store regions and
  // other related data.
  region_growing.release_memory();

  std::cout << std::endl << 
    "region_growing_with_custom_classes example finished" 
  << std::endl << std::endl;

  return EXIT_SUCCESS;
}
