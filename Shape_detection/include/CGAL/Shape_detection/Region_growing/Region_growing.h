#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_H

// STL includes.
#include <queue>
#include <vector>

// CGAL includes.
#include <CGAL/Iterator_range.h>

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief Region growing algorithm
            \tparam InputRange An arbitrary range with user-defined items.
            \tparam Connectivity Model of `RegionGrowingConnectivity`
            \tparam Conditions Model of `RegionGrowingPropagationConditions`
        */
        template<class InputRange, class Connectivity, class Conditions>
        class Region_growing {

        public:

            using Input_range             = InputRange;
            ///< An arbitrary range with user-defined items. The range must implement the function size().
            
            using Item_index              = std::size_t;
            ///< Index of a given item.

            ///< \cond SKIP_IN_MANUAL
            using Visited_items           = std::vector<bool>;

            using Running_queue           = std::queue<Item_index>;

            using Items                   = std::vector<Item_index>;
            ///< \endcond

            using Neighbors               = Items;
            ///< Indices of the neighbors of the given item that should be returned by the call to the function `RegionGrowingConnectivity::get_neighbors()`.

            using Region                  = Items;
            ///< A random access container with indices of all items that belong to a region.

            using Regions                 = std::vector<Region>;
            ///< All found regions.

            using Region_range            = CGAL::Iterator_range<typename Regions::const_iterator>;
            ///< An `Iterator_range` of the iterators in `CGAL::Shape_detection::Region_growing::Regions`.

            using Unclassified_items      = Items;
            ///< A random access container with indices of all unclassified items.

            using Unclassified_range      = CGAL::Iterator_range<typename Unclassified_items::const_iterator>;
            ///< An `Iterator_range` of the iterators in `CGAL::Shape_detection::Region_growing::Unclassified_items`. Here, we store indices of all unclassified items.

            /*!
                The constructor requires an input range and instances of the Connectivity class and Conditions class.
            */
            Region_growing(const Input_range &input_range, Connectivity &connectivity, Conditions &conditions) :
                m_input_range(input_range),
                m_connectivity(connectivity),
                m_conditions(conditions)
                { }

            /*!
                Perform the region growing algorithm over the input range, using the Connectivity class to find neighbors and the Conditions class to validate regions.
            */
            void detect() {

                clear();
                Region region;

                for (Item_index seed_index = 0; seed_index < m_input_range.size(); ++seed_index) {
                    
                    // Try to grow a new region from the seed item.
                    if (!m_visited[seed_index]) {
                        propagate(seed_index, region);

                        // Check global conditions.
                        if (!m_conditions.are_valid(region)) 
                            revert(region);
                        else 
                            m_regions.push_back(region);
                    }
                }
                m_output_regions = Region_range(m_regions.begin(), m_regions.end());

                // Return unclassified items.
                for (Item_index item_index = 0; item_index < m_input_range.size(); ++item_index)
                    if (!m_visited[item_index]) m_unclassified.push_back(item_index);

                m_output_unclassified = Unclassified_range(m_unclassified.begin(), m_unclassified.end());
            }

            /*!
                Return a pair of begin iterator and pass-the-end iterator of the container with found regions. If the function `CGAL::Shape_detection::Region_growing::detect()` has not been called,
                the first and second of the pair will be the same, which implies an empty container.
            */
            const Region_range &regions() {
                return m_output_regions;
            }

            const Unclassified_range &unclassified_items() {
                return m_output_unclassified;
            }

            /*!
                Return the number of found regions.
            */
            size_t number_of_regions() {
                return m_regions.size();
            }

            /*!
                Return the number of unclassified items.
            */
            size_t number_of_unclassified_items() {
                return m_unclassified.size();
            }

            /*!
                Clear all internal data structures.
            */
            void clear() {
                
                m_visited.clear();
                m_regions.clear();

                m_unclassified.clear();
                m_visited.resize(m_input_range.size(), false);
            }

        private:

            void propagate(const Item_index seed_index, Region &region) {
                region.clear();

                // Use two queues, while running on this queue, push to the other queue;
                // When the queue is done, update the shape of the current region and swap to the other queue;
                // depth_index is the index of the queue we are using.
                Running_queue running_queue[2];
                bool depth_index = 0;

                // Once an item is pushed to the queue, it is pushed to the region too.
                m_visited[seed_index] = true;
                running_queue[depth_index].push(seed_index);
                region.push_back(seed_index);

                // Update internal properties of the propagating region.
                m_conditions.update(region);

                while (!running_queue[depth_index].empty() || !running_queue[!depth_index].empty()) {

                    // Call the next item of the queue and remove it from the queue.
                    const Item_index item_index = running_queue[depth_index].front();
                    running_queue[depth_index].pop();

                    // Get neighbors of the current item.
                    Neighbors neighbors;
                    m_connectivity.get_neighbors(item_index, neighbors);

                    // Visit the neighbors.
                    for (Item_index i = 0; i < neighbors.size(); ++i) {
                        const Item_index neighbor_index = neighbors[i];
                        
                        if (!m_visited[neighbor_index] && m_conditions.belongs_to_region(neighbor_index, region)) {

                            // Add this neighbor to the other queue so that we can visit it later.
                            m_visited[neighbor_index] = true;
                            running_queue[!depth_index].push(neighbor_index);
                            region.push_back(neighbor_index);
                        }
                    }

                    // Update internal properties of the propagating region.
                    if (running_queue[depth_index].empty()) {

                        m_conditions.update(region);
                        depth_index = !depth_index;
                    }
                }
            }

            void revert(const Region &region) {

                for (Item_index i = 0; i < region.size(); ++i)
                    m_visited[region[i]] = false;
            }

            // Fields.
            const Input_range       &m_input_range;
            Connectivity            &m_connectivity;
            Conditions              &m_conditions;

            Visited_items           m_visited;
            Regions                 m_regions;
            Unclassified_items      m_unclassified;

            Region_range            m_output_regions      = Region_range(m_regions.begin(), m_regions.end());
            Unclassified_range      m_output_unclassified = Unclassified_range(m_unclassified.begin(), m_unclassified.end());
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_H
