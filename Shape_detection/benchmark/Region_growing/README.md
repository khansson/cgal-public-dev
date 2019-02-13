# This code generates the benchmark table in the user manual.

First, compile the file `benchmark_region_growing_on_points_2.cpp` using the provided cmake file:

```bash
$ cd /path/to/build/directory
$ cmake -DCGAL_DIR=/path/to/cgal/release/build -DCMAKE_BUILD_TYPE=Release /path/to/benchmark/Region_growing
$ make benchmark_region_growing_on_points_2
```

The program uses the data file `data/points_2.xyz`, run the benchmark program as follows:

```bash
$ ./benchmark_region_growing_on_points_2 /path/to/data/points_2.xyz
```

Results:

```
Test #1
  search_radius = 1
  min_region_size = 5
  max_distance_to_line = 4.5
  normal_threshold = 0.7
  -----
  Time elapsed: 0.138201
  Number of detected regions: 796
  Number of assigned points: 4491
  Number of unassigned points: 63277


Test #2
  search_radius = 3
  min_region_size = 5
  max_distance_to_line = 4.5
  normal_threshold = 0.7
  -----
  Time elapsed: 0.072401
  Number of detected regions: 3054
  Number of assigned points: 63154
  Number of unassigned points: 4614


Test #3
  search_radius = 6
  min_region_size = 5
  max_distance_to_line = 4.5
  normal_threshold = 0.7
  -----
  Time elapsed: 0.080463
  Number of detected regions: 2483
  Number of assigned points: 64977
  Number of unassigned points: 2791


Test #4
  search_radius = 9
  min_region_size = 5
  max_distance_to_line = 4.5
  normal_threshold = 0.7
  -----
  Time elapsed: 0.093556
  Number of detected regions: 2282
  Number of assigned points: 65353
  Number of unassigned points: 2415
```
