# k-vertex-critical 
This repository contains code for exhaustively generating k-vertex-critical (P_t,F)-free graphs (see e.g. the paper "Critical (P_5,Dart)-free graphs" at https://arxiv.org/pdf/2308.03414.pdf).

The program can be compiled by using
```bash
make
```

Some examples of how the program can be used:
```bash
./generateKVertexCriticalGraphs -l5 -k4 25 < dartGraph.g6
```
Generates all 5-vertex-critical (P5,dart)-free graphs on at most 25 vertices

```bash
./generateKVertexCriticalGraphs -l5 -k5 25 < dartGraph.g6
```
Generates all 6-vertex-critical (P5,dart)-free graphs on at most 25 vertices
