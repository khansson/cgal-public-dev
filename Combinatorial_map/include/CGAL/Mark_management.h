#ifndef CGAL_MARK_MANAGEMENT_H
#define CGAL_MARK_MANAGEMENT_H 1

#include <CGAL/config.h>
#include <CGAL/disable_warnings.h>
#include <array>
#include <bitset>

namespace CGAL
{

    template <int NB_MARKS, typename size_type = unsigned int>
    struct mmark_stack
    {

        using self_type = mmark_stack<NB_MARKS,size_type>;


        std::array<size_type, NB_MARKS> free_marks;
        std::array<size_type, NB_MARKS> used_marks;
        std::array<size_type, NB_MARKS> mark_indecies;
        size_type nb_used_marks;
        size_type min_mark;
        size_type max_mark;

        mmark_stack() : min_mark(0), max_mark(0), nb_used_marks(0) 
        {

            for( auto& elem :  mark_indecies)
            {
                elem =  NB_MARKS;
            }
        };

        void clear()
        {
        for( auto& elem :  mark_indecies)
            {
                elem = NB_MARKS;
            }
        };
        bool is_free(size_type min)
        {

            for(size_type i = min; i < max_mark; i++)
            {
                
                if (mark_indecies.at(i) != NB_MARKS)
                {
            
                std::cerr << "Not enough Boolean marks: "
                             "increase NB_MARKS in item class."
                          << std::endl;
                std::cerr << "  (exception launched)" << std::endl;
                return false;
                }
            
            }

            //ensure that none of the given marks are considered free
            for(size_type i = 0; i < min; i++)
            {
                
                if (free_marks.at(i) >= min)
                {
                    std::cerr << "Not enough Boolean marks: "
                             "increase NB_MARKS in item class."
                          << std::endl;
                std::cerr << "  (exception launched)" << std::endl;

                return false;
                }
            
                
            
            }


            return true;
        };

        void initialize_mark_range(size_type min, size_type number)
        {
            CGAL_assertion(nb_used_marks == 0);

            for (int i = 0; i < number; i++)
            {
                free_marks.at(i) = min + i;
            }

            for (auto &elem : used_marks)
            {
                elem = 0;
            }
            min_mark = min;
            max_mark = min + number;
        };


        void give_marks(self_type &other, size_type count)
        {

            CGAL_assertion(other.min_mark == other.max_mark);
            


            size_type cutoff = max_mark - count;
            CGAL_assertion(cutoff > min_mark);

            CGAL_assertion(is_free(cutoff));


            


            max_mark = cutoff;

            other.initialize_mark_range(cutoff, count);



        };

        void absorb_marks(self_type &other)
        {

            CGAL_assertion(other.is_free(other.min_mark));
            CGAL_assertion(max_mark == other.min_mark);

            max_mark = other.max_mark;

            size_type mark_count = max_mark-min_mark;
            
            for(size_type i = 0; i < other.max_mark - other.min_mark; i++)
            {

                free_marks.at(mark_count+i) = other.free_marks.at(i);

            }



        };

        size_type create_mark()
        {

            size_type created_mark = NB_MARKS;
            if (nb_used_marks == max_mark - min_mark)
            {
                std::cerr << "Not enough Boolean marks: "
                             "increase NB_MARKS in item class."
                          << std::endl;
                std::cerr << "  (exception launched)" << std::endl;
                return created_mark;
            }

            created_mark = free_marks.at(nb_used_marks);
            mark_indecies.at(created_mark) = nb_used_marks;
            nb_used_marks++;
            return created_mark;
        };

        void free_mark(size_type amark)
        {
            CGAL_assertion(min_mark <= amark && amark < max_mark);

            nb_used_marks--;

            used_marks.at(mark_indecies.at(amark)) = used_marks.at(nb_used_marks);

            free_marks.at(nb_used_marks) = amark;
            mark_indecies.at(amark) = NB_MARKS;
        };
    };

    template <int PARTITION_SIZE, int NUM_PARTITIONS>
    struct PartitionedBitset
    {

        std::array<std::bitset<PARTITION_SIZE>, NUM_PARTITIONS> data;

        auto operator[](std::size_t i) -> decltype(data.at(0)[0])
        {
            return data.at(i / PARTITION_SIZE)[i % PARTITION_SIZE];
        }

        void flip(int i) { data.at(i / PARTITION_SIZE).flip(i % PARTITION_SIZE); }

        void reset()
        {
            for (auto &bit : data)
            {
                bit.reset();
            }
        }

        void reset(std::size_t i) { data.at(i / PARTITION_SIZE).reset(i % PARTITION_SIZE); }

        void set(std::size_t i, bool val = true) { data.at(i / PARTITION_SIZE).set(i % PARTITION_SIZE, val); }
    };

    template <int PARTITION_SIZE>
    struct PartitionedBitset<PARTITION_SIZE, 1> : public std::bitset<PARTITION_SIZE>
    {
    };

    template <int PARTITION_SIZE, int NUM_PARTITIONS>
    PartitionedBitset<PARTITION_SIZE, NUM_PARTITIONS>
    operator^(const PartitionedBitset<PARTITION_SIZE, NUM_PARTITIONS> one,
              const PartitionedBitset<PARTITION_SIZE, NUM_PARTITIONS> two)
    {
        PartitionedBitset<PARTITION_SIZE, NUM_PARTITIONS> to_return;
        for (std::size_t i = 0; i < NUM_PARTITIONS; i++)
        {

            to_return.data.at(i) = one.data.at(i) ^ two.data.at(i);
        }

        return to_return;
    }

} // namespace CGAL
#endif