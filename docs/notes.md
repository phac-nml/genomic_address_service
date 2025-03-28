# Single Linkage Clustering Assignment Example

There are two parts to assigning a clustering label: generating a hierarchical linkage and flattening that linkage into clusters (with labels).

## Input Distance Matrix

```
dists  A  B  C  D  E  F  G  H  I  J
A      0  1  2  2  5  5  6  6  9  9
B      1  0  3  3  6  6  7  7  9  9
C      2  3  0  0  3  3  6  6  9  9
D      2  3  0  0  3  3  6  6  9  9
E      5  6  3  3  0  0  6  6  9  9
F      5  6  3  3  0  0  6  6  9  9
G      6  7  6  6  6  6  0  0  3  3
H      6  7  6  6  6  6  0  0  3  3
I      9  9  9  9  9  9  3  3  0  1
J      9  9  9  9  9  9  3  3  1  0
```

## 1-D Condensed Distance Matrix

`[1. 2. 2. 5. 5. 6. 6. 9. 9. 3. 3. 6. 6. 7. 7. 9. 9. 0. 3. 3. 6. 6. 9. 9. 3. 3. 6. 6. 9. 9. 0. 6. 6. 9. 9. 6. 6. 9. 9. 0. 3. 3. 3. 3. 1.]`

Notice that this is functionally the top right part of the distance matrix (above the main diagonal) all in one row.

## Linkage

`scipy.cluster.hierarchy.linkage` will generate this linkage matrix Z (my rendition):

```
10: [ 2.  3.  0.  2.] -> (C,D)
11: [ 4.  5.  0.  2.] -> (E,F)
12: [ 6.  7.  0.  2.] -> (G,H)
13: [ 0.  1.  1.  2.] -> (A,B)
14: [ 8.  9.  1.  2.] -> (I,J)
15: [10. 13.  2.  4.] -> ((C,D), (A,B))
16: [11. 15.  3.  6.] -> ((E,F), ((C,D), (A,B)))
17: [12. 14.  3.  4.] -> ((G,H), (I,J))
18: [16. 17.  6. 10.] -> (((E,F), ((C,D), (A,B))), ((G,H), (I,J)))
```

With the following extended implicit linkage matrix:

```
0: [ 0.  0.  0.  1.] -> (A)
1: [ 1.  1.  0.  1.] -> (B)
2: [ 2.  2.  0.  1.] -> (C)
3: [ 3.  3.  0.  1.] -> (D)
4: [ 4.  4.  0.  1.] -> (E)
5: [ 5.  5.  0.  1.] -> (F)
6: [ 6.  6.  0.  1.] -> (G)
7: [ 7.  7.  0.  1.] -> (H)
8: [ 8.  8.  0.  1.] -> (I)
9: [ 9.  9.  0.  1.] -> (J)
10: [ 2.  3.  0.  2.] -> (C,D)
11: [ 4.  5.  0.  2.] -> (E,F)
12: [ 6.  7.  0.  2.] -> (G,H)
13: [ 0.  1.  1.  2.] -> (A,B)
14: [ 8.  9.  1.  2.] -> (I,J)
15: [10. 13.  2.  4.] -> ((C,D), (A,B))
16: [11. 15.  3.  6.] -> ((E,F), ((C,D), (A,B)))
17: [12. 14.  3.  4.] -> ((G,H), (I,J))
18: [16. 17.  6. 10.] -> (((E,F), ((C,D), (A,B))), ((G,H), (I,J)))
```

The form I'm showing it is like this:

`cluster_index: [ cluster_index_1.  cluster_index_2.  distance.  num_elements.] -> (representation)`

Where `cluster_index` is composed of two other clusters (`cluster_index_1.` and `cluster_index_2.`) with a distance between them of `distance.` and has `num_elements.` elements. `cluster_index` 0-(n-1) are implicit and are clusters with one element each (themselves). So this:

`10: [ 2.  3.  0.  2.] -> (C,D)`

means cluster index 10 is made from cluster index 2 and 3, has a distance of 0 (between children), and contains 2 original elements (C, D). And this:

`15: [10. 13.  2.  4.] -> ((C,D), (A,B))`

means cluster index 15 is composed of 10 and 13, has a distance of 2, and contains 4 elements ((C,D), (A,B)). Etc.

## Generating Flat Clusters from Linkage

Flattening the hierarchical linkage into clusters just descends down the linkage and stops once the two children of the node have a distance that is less than or equal to the threshold value. The same linkage is used for various thresholds, since the thresholds only change how the cluster is flattened, not how it was originally hierarchically organized.

### Example

A threshold of `2` would mean:

```
18 (root) distance is 6 > 2, descend to children (16, 17)

16 distance is 3 > 2, descend to children (11, 15)
17 distance is 3 > 2, descend to children (12, 14)

11 distance is 0 <= 2, label together (E,F)
15 distance is 2 <= 2, label together ((C,D), (A,B))
12 distance is 0 <= 2, label together (G,H)
14 distance is 1 <= 2, label together (I,J)
```

## Labeling Clusters

Since labels don't mean anything other than showing the group, this matches the program output:

```
id      address level_1
A       2       2
B       2       2
C       2       2
D       2       2
E       1       1
F       1       1
G       3       3
H       3       3
I       4       4
J       4       4
```

# Complete Linkage Example

## Distance Matrix

```
dists    A    B    C    D    E
A        0    1    2    3    4
B        1    0    5    6    7
C        2    5    0    8    9
D        3    6    8    0    10
E        4    7    9    10   0
```

## Manual Complete Linkage

### Round 1

The smallest distance is (A,B). Cluster and recalculate:

```
dists    (A,B)    C    D    E
(A,B)    0        5    6    7
C        5        0    8    9
D        6        8    0    10
E        7        9    10   0
```

### Round 2

The smallest distance is ((A,B),C). Cluster and recalculate:

```
dists        ((A,B),C)    D    E
((A,B),C)    0            8    9
D            8            0    10
E            9            10   0
```

### Round 3

The smallest distance is (((A,B),C),D). Cluster and recalculate:

```
dists            (((A,B),C),D)    E
(((A,B),C),D)    0                10
E                10               0
```

### Round 4

The smallest distance and only remaining cluster operation is ((((A,B),C),D),E).

## Generated Linkage (Programmatically)

```
5: [0. 1. 1. 2.] -> (A,B)
6: [2. 5. 5. 3.] -> (C,(A,B))
7: [3. 6. 8. 4.] -> (D,(C,(A,B)))
8: [4. 7. 10. 5.] -> (E,(D,(C,(A,B))))
```

Which matches the manually calculated linkage.

## Generating Flat Clusters from Linkage

We can see from above (third column of the linkage) that `{1, 5, 8, 10}` are distances between clusters in the linkage. So we should choose thresholds that are aware of these (for demonstration and testing purposes).

### threshold=10

Root: `8: [4. 7. 10. 5.] -> (E,(D,(C,(A,B))))`

The whole linkage's distance is `10<=10`, so everything will cluster together as the same label.

Result: `{A,B,C,D,E}`

### threshold=8

Root: `8: [4. 7. 10. 5.] -> (E,(D,(C,(A,B))))`

The distance is too great (`10>8`). We need to descend to the children (4 and 7):

`4: [4. 4. 0. 1.] -> (E) (implicit)`  
`7: [3. 6. 8. 4.] -> (D,(C,(A,B)))`

4 is implicit and has a distance of `0<=8`, so it forms one cluster: E
7 has a distance of `8<=8`, so all of it forms one cluster: A,B,C,D

Result: {A,B,C,D}, {E}

### threshold=5

Root: `8: [4. 7. 10. 5.] -> (E,(D,(C,(A,B))))`

The distance is too great (`10>5`). We need to descend to the children (4 and 7):

`4: [4. 4. 0. 1.] -> (E) (implicit)`  
`7: [3. 6. 8. 4.] -> (D,(C,(A,B)))`

4 is implicit and has a distance of `0<=5`, so it forms one cluster: E
7 has a distance of `8>5`, so we need to descend to the children (3 and 6).

`3: [3. 3. 0. 1.] -> (D) (implicit)`  
`6: [2. 5. 5. 3.] -> (C,(A,B))`

3 is implicit and has a distance of `0<=5`, so it forms one cluster: D
6 has a distance of `5<=5`, so all of it forms one cluster: A,B,C

Result: `{A,B,C}, {D}, {E}`

## GAS mcluster / fcluster Output

```
id      address level_1 level_2 level_3
A       1.1.1   1       1       1
B       1.1.1   1       1       1
C       1.1.1   1       1       1
D       1.1.2   1       1       2
E       1.2.3   1       2       3
```

Which matches the manual expectation calculated above.

# Average Linkage Example

## Distance Matrix

```
dists    A    B    C
A        0    4    8
B        4    0    12
C        8    12    0
```

## Manual Complete Linkage

### Round 1

The smallest distance is (A,B). Cluster and recalculate:

```
dists    (A,B)    C
(A,B)    0        10
C        10       0
```

```
dist((A,B), C) = (dist(A, C) + dist(B, C)) / (size((A,B)) * size(C))
dist((A,B), C) = (8 + 12) / (2 * 1)
dist((A,B), C) = 10
```

## Round 2

The only remaining option is clustering (A,B) with C: ((A,B), C).

## Generated Linkage (Programmatically)

3: [ 0.  1.  4.  2.] -> (A,B)
4: [ 2.  3. 10.  3.] -> (C, (A,B))

Which matches our manual calculation above.

# Links

https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html

https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html
