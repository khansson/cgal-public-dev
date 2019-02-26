// STL includes.
#include <map>
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <iterator>

// CGAL includes.
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>

// Custom Neighbor_query, Region_type, and Seed_map classes for region growing.
namespace Custom {

  // An object that stores indices of all adjacent to it other objects.
  struct Object {
    std::vector<std::size_t> neighbors;
  };

  // A range of objects.
  using Objects = std::list<Object>;

  // Neighbor_query class that
  // simply returns indices of all neighbors stored in
  // the object struct above.
  class Neighbor_query {
    const Objects& m_objects;

  public:
    Neighbor_query(Objects& objects) : 
    m_objects(objects) 
    { }

    void operator()(
      const std::size_t query_index, 
      std::vector<std::size_t>& neighbors) const {
      
      std::size_t i = 0;
      for (const auto& object : m_objects) {
        if (i == query_index) {
          
          neighbors = object.neighbors;
          return;
        } ++i;
      }
    }
  };

  // Region_type class, where the function is_part_of_region() checks
  // a very specific condition that the first and second objects in the
  // range are in fact neighbors; is_valid_region() function always 
  // returns true after the first call to the update() function.
  // These are the only functions that we have to define.
  class Region_type {
    bool m_is_valid = false;

  public:
    Region_type() { }

    bool is_part_of_region(
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
      return m_objects_map.find(key)->second;
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
using Object         = Custom::Object;
using Objects        = Custom::Objects;
using Neighbor_query = Custom::Neighbor_query;
using Region_type    = Custom::Region_type;
using Seed_map       = Custom::Seed_map;

using Region_growing = 
CGAL::Shape_detection::Region_growing<Objects, Neighbor_query, Region_type, Seed_map>;

int main(int argc, char *argv[]) { 
  
  std::cout << std::endl << 
    "region_growing_with_custom_classes example started" 
  << std::endl << std::endl;

  // Define a range of objects, where the first two objects form
  // the first region, while the third object forms the second region.
  Objects objects;

  // Region 1.
  Object object1;
  object1.neighbors.resize(1, 1); 
  objects.push_back(object1);
  
  Object object2;
  object2.neighbors.resize(1, 0); 
  objects.push_back(object2);

  // Region 2.
  Object object3;
  objects.push_back(object3);

  // Extra object to skip.
  Object object4;
  objects.push_back(object4);

  // Create instances of the classes Neighbor_query and Region_type.
  Neighbor_query neighbor_query = Neighbor_query(objects);
  Region_type    region_type    = Region_type();

  // Create a seed map.
  std::map<std::size_t, std::size_t> objects_map;
  objects_map[0] = 1; // the order is swapped with the next object
  objects_map[1] = 0;
  objects_map[2] = 2; // the default order
  objects_map[3] = std::size_t(-1); // skip this object
  const Seed_map seed_map(objects_map);

  // Create an instance of the region growing class.
  Region_growing region_growing(objects, neighbor_query, region_type, seed_map);

  // Run the algorithm.
  using Region = std::vector<std::size_t>;
  std::list<Region> regions;
  region_growing.detect(std::back_inserter(regions));

  // Print the number of found regions. It must be two regions.
  std::cout << "* " << regions.size() << 
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
