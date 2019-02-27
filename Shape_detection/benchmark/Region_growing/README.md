# This code generates the benchmark table in the user manual.

First, compile the file `benchmark_region_growing_on_point_set_2.cpp` using the provided cmake file:

```bash
$ cd /path/to/build/directory
$ cmake -DCGAL_DIR=/path/to/cgal/release/build -DCMAKE_BUILD_TYPE=Release /path/to/benchmark/Region_growing
$ make benchmark_region_growing_on_point_set_2
```

The program uses the data file `data/point_set_2.xyz`.
You can run the benchmark program as follows:

```bash
$ ./benchmark_region_growing_on_point_set_2 /path/to/data/point_set_2.xyz
```

Results:

```
Test #1
  sphere_radius = 1
  min_region_size = 5
  distance_threshold = 4.5
  angle_threshold = 45
  -----
  Time elapsed: 0.138831
  Number of detected regions: 794
  Number of assigned points: 4483
  Number of unassigned points: 63285


Test #2
  sphere_radius = 3
  min_region_size = 5
  distance_threshold = 4.5
  angle_threshold = 45
  -----
  Time elapsed: 0.069098
  Number of detected regions: 3063
  Number of assigned points: 63038
  Number of unassigned points: 4730


Test #3
  sphere_radius = 6
  min_region_size = 5
  distance_threshold = 4.5
  angle_threshold = 45
  -----
  Time elapsed: 0.077703
  Number of detected regions: 2508
  Number of assigned points: 64906
  Number of unassigned points: 2862


Test #4
  sphere_radius = 9
  min_region_size = 5
  distance_threshold = 4.5
  angle_threshold = 45
  -----
  Time elapsed: 0.093415
  Number of detected regions: 2302
  Number of assigned points: 65334
  Number of unassigned points: 2434
```
