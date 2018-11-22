#ifndef CGAL_SHAPE_DETECTION_REGION_GROWING_H
#define CGAL_SHAPE_DETECTION_REGION_GROWING_H

// STL includes.
#include <map>
#include <list>
#include <queue>

// CGAL includes.
#include <CGAL/Iterator_range.h>

// Todo:
// Reduce to one queue? - no, because we have to know when we need to update the shape.
// Change Element to Element_with_properties in m_visited? - this may be fixed with unordered_map, because it does not have dependency on < operator.
// Element_with_properties neighbor = *it; and Element_with_properties ewp = ewp_queue[depth_index].front(); - add const and reference.

namespace CGAL {

    namespace Shape_detection {

        /*!
            \ingroup PkgShapeDetection
            \brief Region growing algorithm
            \tparam RegionGrowingTraits `CGAL::Shape_detection::Region_growing_traits`
            \tparam Connectivity Model of `RegionGrowingConnectivity`
            \tparam Conditions Model of `RegionGrowingPropagationConditions`
        */

        template<class RegionGrowingTraits, class Connectivity, class Conditions>
        class Region_growing {

        public:

            using Region_growing_traits   = RegionGrowingTraits;
            ///< Type storing region growing traits.

            using Input_range             = typename Region_growing_traits::Input_range;
            ///< Type storing user-defined elements.

            using Element                 = typename Region_growing_traits::Element;
            ///< The primary geometric element on which the region growing algorithm will run
            ///< The operator '<' must be defined for Element so that `std::map<Element>` can be used.
            
            using Element_map             = typename Region_growing_traits::Element_map;
            ///< An `LvaluePropertyMap` that maps to type `CGAL::Shape_detection::Region_growing::Element`.
            
            using Element_with_properties = typename Element_map::key_type;
            ///< Value type of the iterator in the input range.

            using Neighbors               = std::list<Element_with_properties>;
            ///< Type storing items, used in `RegionGrowingConnectivity::get_neighbors()`.

            using Region                  = std::list<Element_with_properties>;
            ///< Type storing items that represent a region.

            using Regions                 = std::list<Region>;
            ///< Type storing regions.

            using Region_range            = CGAL::Iterator_range<typename Regions::const_iterator>;
            ///< An `Iterator_range` of iterators in `CGAL::Shape_detection::Region_growing::Regions`.

            /*!
                The constructor requires an input range and instances of the Connectivity class and Conditions class.
            */
            Region_growing(const Input_range &input_range, Connectivity &connectivity, Conditions &conditions) :
                m_input_range(input_range),
                m_connectivity(connectivity),
                m_conditions(conditions) 
                { }

            /*!
                Perform the region growing algorithm over the input range, using the Connectivity class to find neighbors and the PropagationConditions class to validate elements and regions.
            */
            void find_regions() {

                m_regions.clear();
                Region region;

                for (auto it = m_input_range.begin(); it != m_input_range.end(); ++it) {
                    if (!m_visited[get(m_elem_map, *it)]) { // element is available

                        // Grow a region from that element.
                        region.clear();
                        grow_region(*it, region);

                        // Check global conditions.
                        if (!m_conditions.is_valid(region)) revert(region);
                        else m_regions.push_back(region);
                    }
                }
                m_output = Region_range(m_regions.begin(), m_regions.end());
            }

            /*!
                Return a pair of begin iterator and pass-the-end iterator of the list of regions found. If `CGAL::Shape_detection::Region_growing::find_regions()` has not been called, 
                the first and second of the pair will be the same, which implies an empty list.
            */
            const Region_range &regions() {
                return m_output;
            }

            /*!
                Return the number of regions found.
            */
            size_t number_of_regions() {
                return m_regions.size();
            }

        private:
            void grow_region(const Element_with_properties &seed, Region &region) {
                region.clear();

                // Use two queues, while running on this queue, push to the other queue;
                // When the queue is done, update the shape of the current region and swap to the other queue;
                // depth_index is the index of the queue we're using.
                std::queue<Element_with_properties> ewp_queue[2];
                bool depth_index = 0;

                // Once an element is pushed to the queue, it is pushed to the region too.
                ewp_queue[depth_index].push(seed);
                m_visited[get(m_elem_map, seed)] = true;
                region.push_back(seed);

                // Update internal properties of the growing shape if needed.
                m_conditions.update_shape(region);

                while (!ewp_queue[depth_index].empty() || !ewp_queue[!depth_index].empty()) {

                    // Call the next element of the queue and remove it from the queue
                    // but the element is not removed from the region.
                    Element_with_properties ewp = ewp_queue[depth_index].front();
                    ewp_queue[depth_index].pop();

                    // Get neighbors.
                    Neighbors neighbors;
                    m_connectivity.get_neighbors(ewp, neighbors);

                    // Visit the neighbors.
                    for (auto it = neighbors.begin(); it != neighbors.end(); ++it) {
                        Element_with_properties neighbor = *it;

                        if (!m_visited[get(m_elem_map, neighbor)] && m_conditions.is_in_same_region(ewp, neighbor, region)) {

                            // Add to the other queue the neighbor which doesn't belong to any regions
                            // so that we can visit them later.
                            ewp_queue[!depth_index].push(neighbor);
                            region.push_back(neighbor);
                            m_visited[get(m_elem_map, neighbor)] = true;
                        }
                    }

                    if (ewp_queue[depth_index].empty()) {

                        m_conditions.update_shape(region);
                        depth_index = !depth_index;
                    }
                }
            }

            void revert(const Region &region) {
                for (auto it = region.begin(); it != region.end(); ++it)
                    m_visited[get(m_elem_map, *it)] = false;
            }

        private:
            const Input_range              &m_input_range;
            Connectivity                   &m_connectivity;
            Conditions                     &m_conditions;
            
            Element_map                     m_elem_map;
            Regions                         m_regions;
            std::map<Element, bool>         m_visited;
            Region_range                    m_output = Region_range(m_regions.begin(), m_regions.end());
        };

    } // namespace Shape_detection

} // namespace CGAL

#endif // CGAL_SHAPE_DETECTION_REGION_GROWING_H
