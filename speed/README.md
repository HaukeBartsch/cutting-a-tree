# Speed up processing of graphs

Prune is a reimplementation of the javascript code in the parent directory. 

Here is an example on how to call prune:

```bash
./prune --verbose -s 10 \
    -e data/CD1-41-bulbs-edges.csv \
    -n data/CD1-41-bulbs-nodes.csv \
    -o data/output.csv \
    -r data/remove_vertices.json
```

We switch on verbose mode and remove 10 edges from the graph. Edges and nodes are read from csv files and the assignment of classes is stored in output.csv (one row per node). The 10 edges are stored in the remove_vertices.json file. If that file already exists the included edge indices are removed before processing starts. 

> [!NOTE] 
> Edge indices are incremental, e.g. you need to remove the first index before you can use the index of the second entry.

![example run](https://github.com/HaukeBartsch/cutting-a-tree/raw/main/speed/images/example_run.gif)

## How it works

Assumption: All edges are the same irrespective of parameters like length and thickness.

Objective function: Remove edges iteratively and ensure that a) resulting graph is at least equally or more 'balanced' compared to before and b) the graph is still fully connected after removal of an edge.

### 'balanced'

A diffusion algorithm is used to distribute arterial and venous input through edges in the graph based on an initial assignment of two nodes to artery and vein status (boundary contition). All other nodes receive an artery or a vein assignment based on the maximum value of either incoming arterial or incoming venous input. 

A graph is balanced if across all levels an equal number of arteries and veines are present. 

Levels are computed recursively using an octree-based space separation. The graphs (3D) bouding box is split into eight equally large bounding boxes that are in-turn split again. The space separation ends if a bounding box does not contain any nodes or if a maximum number of level splits are performed (4).


### Build

The project is build using cmake and tbb (threads building blocks).

```bash
cmake -DCMAKE_BUILD_TYPE=Debug .
make
```

Test the build by calling prune:

```bash
❯ ./prune 
Allowed options:
  -h [ --help ]                prune: Run a diffusion quad-tree pruning on a 
                               given graph.
  -V [ --version ]             Print the version number.
  -v [ --verbose ]             Print more verbose output during processing.
  -e [ --edges ] arg           The edges csv file.
  -n [ --nodes ] arg           The nodes csv file.
  -o [ --output ] arg          Path to output csv file.
  -c [ --stopping ] arg        When to assume diffusion solution has converged 
                               based on summed overall change [default 0.3].
  -r [ --remove_vertices ] arg Path to a json that contains vertices that 
                               should be removed initially. This file will be 
                               updated continuously.
  -s [ --steps ] arg           Number of connections to remove.

```

