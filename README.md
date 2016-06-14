
## Implementations of different k-means algorithms in C.


## Using guide:

### Installation

With having gcc compiler and Make installing the program should be as easy as writing make. An executable "clbin" will be created.

### Flags

```
Help menu. Argument and its value must be separated by a space:
-h  - help menu
-f  - input data points file name
-o  - output file name
-ci - centers input file name
-a  - algorithm choice (lloyd, elkan, hamerly, macqueen, hartigan, closest)
-i  - initialization method choice (kpp, forgy, partition, furthest, firstn)
-m  - metric choice (euclidean, manhattan)
-k  - clusters count, -ci flag has a higher priority
-s  - random seed nr, otherwise uses current time as the seed. Used to confirm clustering results
-n  - iteration count (default 100)
```

## Contains the following:

### Choosing initial cluster centers:

1. k-means++
2. Forgy
3. Partition (assigns points to random cluster and then finds means of these assignments as centers)
4. Furthest first
5. firstn (chooses k first points from input file as initial cluster centers).

### Clustering algorithms:

1. Lloyd
2. Elkan
3. Hamerly
4. MacQueen
5. Hartigan-Wong
6. Closest (just assigns points to closest centers and stops)

### Supports following metrics

1. Euclidean
2. Manhattan

It should be rather easy to add more of them.

### Input data format

1. First line consists of two integers: Point count and data dimensionality. Follows n rows with d doubles on each of them.

### Output data formats

Outputs 2 files

1. File consisting of numbers to which cluster a point belongs to
2. File consisting of means' vectors. First row has an integer which shows how many means follow.