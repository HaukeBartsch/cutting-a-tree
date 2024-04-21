function readDataVert(callback) {
    var deferred = new jQuery.Deferred();
    jQuery.get( 'data/DLA1-edges.csv', function(csv) {
        
        var header = [];
        var data = {};
        var vertices = [];
        csv.split("\n").forEach(function(line) {
            if (header.length == 0) {
                header = line.split(",");
                header = header.map(function(a) { return a.trim(); });
                return true;
            }
            // assume we have a header now
            dat = line.split(",");
            dat = dat.map(function(a) { return a.trim(); });
            for (var i = 0; i < dat.length; i++) {
                var key = "new";
                if (i < header.length) {
                    key = header[i];
                    if (typeof data[key] == 'undefined')
                    data[key] = [];
                    data[key].push(dat[i]);
                }
            }
        });
        // get the nodes from avx_1, avx_2, and avx_3
        data["EndNodes_1"].map(function(a, idx) {
            var p = [
                parseInt(data["EndNodes_1"][idx]), 
                parseInt(data["EndNodes_2"][idx])
            ];
            if (isFinite(p[0]) && isFinite(p[1]))
                vertices.push(p);
        });
        deferred.resolve(vertices);
    });
    return deferred.promise();
}


function readDataPos() {
    var deferred = new jQuery.Deferred();
    
    jQuery.get( 'data/DLA1-nodes.csv', function(csv) {
        // convert the csv data to a structure
        var header = [];
        var data = {};
        var pos = [];
        csv.split("\n").forEach(function(line) {
            if (header.length == 0) {
                header = line.split(",");
                header = header.map(function(a) { return a.trim(); });
                return true;
            }
            // assume we have a header now
            dat = line.split(",");
            dat = dat.map(function(a) { return a.trim(); });
            for (var i = 0; i < dat.length; i++) {
                var key = "new";
                if (i < header.length) {
                    key = header[i];
                    if (typeof data[key] == 'undefined')
                    data[key] = [];
                    data[key].push(dat[i]);
                }
            }
        });
        // get the nodes from avx_1, avx_2, and avx_3
        data["avx_1"].map(function(a, idx) {
            var p = [parseFloat(data["avx_2"][idx]), 
            parseFloat(data["avx_3"][idx]), 
            parseFloat(data["avx_1"][idx])];
            pos.push(p);
        });
        var mean = [0,0,0];
        for(var i = 0; i < pos.length; i++) {
            mean[0] += pos[i][0];
            mean[1] += pos[i][1];
            mean[2] += pos[i][2];
        }
        mean[0] /= pos.length;
        mean[1] /= pos.length;
        mean[2] /= pos.length;
        
        for (var i = 0; i < pos.length; i++) {
            pos[i][0] -= mean[0];
            pos[i][1] -= mean[1];
            pos[i][2] -= mean[2];
        }
        // get the belongroot array (node_type)
        var node_type = [];
        pos.map(function(a,idx) { 
            node_type[idx] = "";
        });
        data["avx_1"].map(function(a, idx) {
            node_type[idx] = data["belongroot"][idx];
            if (node_type[idx] == "none")
                node_type[idx] = "";
        });

        deferred.resolve([pos, node_type]);
    });
    return deferred.promise();
}

function computeBounds(points) {
    var increase_by = 1.0/10.0; // percent
    // compute bounds
    var bounds = [points[0][0],points[0][1],points[0][2],points[0][0],points[0][1],points[0][2]];
    for (var i = 0; i < points.length; i++) {
        if (points[i][0] < bounds[0])
            bounds[0] = points[i][0];
        if (points[i][1] < bounds[1])
            bounds[1] = points[i][1];
        if (points[i][2] < bounds[2])
            bounds[2] = points[i][2];
        if (points[i][0] > bounds[3])
            bounds[3] = points[i][0];
        if (points[i][1] > bounds[4])
            bounds[4] = points[i][1];
        if (points[i][2] > bounds[5])
            bounds[5] = points[i][2];
    }
    var size = Math.min(bounds[3]-bounds[0],bounds[4]-bounds[1],bounds[5]-bounds[2]);

    bounds[0] -= size*increase_by; 
    bounds[1] -= size*increase_by;
    bounds[2] -= size*increase_by;
    bounds[3] += size*increase_by; 
    bounds[4] += size*increase_by;
    bounds[5] += size*increase_by;
     
    return bounds;
}

function inside(bounds, p) {
    // is point p inside bounds?
    if (p[0] >= bounds[0] && p[0] <= bounds[3] && p[1] >= bounds[1] && p[1] <= bounds[4] && p[2] >= bounds[2] && p[2] <= bounds[5])
        return true;
    return false;
}

// create a quad-tree for the points
// based on bounds points can be inside more than one branch
function getQuadTree( points, tree, startLevel ) {
    // the quad tree, key is path into the tree, length is depth
    // all nodes on level 0
    if (Object.keys(tree).length == 0) { // init the tree
        tree = { "0": { children: Array.from(Array(points.length).keys()), bounds: computeBounds(points) } };
        startLevel = "0";
    }
    // split the startLevel now and create 8 more entries based on bounds
    var indices = tree[startLevel].children;
    var bounds = tree[startLevel].bounds;
    var widths = [bounds[3]-bounds[0], bounds[4]-bounds[1], bounds[5]-bounds[2]];

    var newK = startLevel + "0";
    var newChildren = []; // array with indices
    var bb = [bounds[0], bounds[1], bounds[2], 0, 0, 0]; // start index
    bb = [bb[0], bb[1], bb[2], 
        bb[0] + widths[0]/2.0, 
        bb[1] + widths[1]/2.0, 
        bb[2] + widths[2]/2.0];
    for (var i = 0; i < indices.length; i++) {
        var p = points[indices[i]];
        if (inside(bb, p)) {
            newChildren.push(indices[i]);
        }
    }
    if (newChildren.length > 0)
        tree[newK] = { children: newChildren, bounds: bb };

    var newK = startLevel + "1";
    var newChildren = []; // array with indices
    var bb = [bounds[0] + widths[0]/2, bounds[1], bounds[2], 0, 0, 0]; // start index
    bb = [bb[0], bb[1], bb[2], 
        bb[0] + widths[0]/2.0, 
        bb[1] + widths[1]/2.0, 
        bb[2] + widths[2]/2.0];
    for (var i = 0; i < indices.length; i++) {
        var p = points[indices[i]];
        if (inside(bb, p)) {
            newChildren.push(indices[i]);
        }
    }
    if (newChildren.length > 0)
        tree[newK] = { children: newChildren, bounds: bb };


    var newK = startLevel + "2";
    var newChildren = []; // array with indices
    var bb = [bounds[0], bounds[1] + widths[1]/2, bounds[2], 0, 0, 0]; // start index
    bb = [bb[0], bb[1], bb[2], 
        bb[0] + widths[0]/2.0, 
        bb[1] + widths[1]/2.0, 
        bb[2] + widths[2]/2.0];
    for (var i = 0; i < indices.length; i++) {
        var p = points[indices[i]];
        if (inside(bb, p)) {
            newChildren.push(indices[i]);
        }
    }
    if (newChildren.length > 0)
        tree[newK] = { children: newChildren, bounds: bb };

    var newK = startLevel + "3";
    var newChildren = []; // array with indices
    var bb = [bounds[0] + widths[0]/2, bounds[1] + widths[1]/2, bounds[2], 0, 0, 0]; // start index
    bb = [bb[0], bb[1], bb[2], 
        bb[0] + widths[0]/2.0, 
        bb[1] + widths[1]/2.0, 
        bb[2] + widths[2]/2.0];
    for (var i = 0; i < indices.length; i++) {
        var p = points[indices[i]];
        if (inside(bb, p)) {
            newChildren.push(indices[i]);
        }
    }
    if (newChildren.length > 0)
        tree[newK] = { children: newChildren, bounds: bb };

    //
    // next level counting again, same order
    //
    var newK = startLevel + "4";
    var newChildren = []; // array with indices
    var bb = [bounds[0], bounds[1], bounds[2] + widths[2]/2, 0, 0, 0]; // start index
    bb = [bb[0], bb[1], bb[2], 
        bb[0] + widths[0]/2.0, 
        bb[1] + widths[1]/2.0, 
        bb[2] + widths[2]/2.0];
    for (var i = 0; i < indices.length; i++) {
        var p = points[indices[i]];
        if (inside(bb, p)) {
            newChildren.push(indices[i]);
        }
    }
    if (newChildren.length > 0)
        tree[newK] = { children: newChildren, bounds: bb };

    var newK = startLevel + "5";
    var newChildren = []; // array with indices
    var bb = [bounds[0] + widths[0]/2, bounds[1], bounds[2] + widths[2]/2, 0, 0, 0]; // start index
    bb = [bb[0], bb[1], bb[2], 
        bb[0] + widths[0]/2.0, 
        bb[1] + widths[1]/2.0, 
        bb[2] + widths[2]/2.0];
    for (var i = 0; i < indices.length; i++) {
        var p = points[indices[i]];
        if (inside(bb, p)) {
            newChildren.push(indices[i]);
        }
    }
    if (newChildren.length > 0)
        tree[newK] = { children: newChildren, bounds: bb };


    var newK = startLevel + "6";
    var newChildren = []; // array with indices
    var bb = [bounds[0], bounds[1] + widths[1]/2, bounds[2] + widths[2]/2, 0, 0, 0]; // start index
    bb = [bb[0], bb[1], bb[2], 
        bb[0] + widths[0]/2.0, 
        bb[1] + widths[1]/2.0, 
        bb[2] + widths[2]/2.0];
    for (var i = 0; i < indices.length; i++) {
        var p = points[indices[i]];
        if (inside(bb, p)) {
            newChildren.push(indices[i]);
        }
    }
    if (newChildren.length > 0)
        tree[newK] = { children: newChildren, bounds: bb };

    var newK = startLevel + "7";
    var newChildren = []; // array with indices
    var bb = [bounds[0] + widths[0]/2, bounds[1] + widths[1]/2, bounds[2] + widths[2]/2, 0, 0, 0]; // start index
    bb = [bb[0], bb[1], bb[2], 
        bb[0] + widths[0]/2.0, 
        bb[1] + widths[1]/2.0, 
        bb[2] + widths[2]/2.0];
    for (var i = 0; i < indices.length; i++) {
        var p = points[indices[i]];
        if (inside(bb, p)) {
            newChildren.push(indices[i]);
        }
    }
    if (newChildren.length > 0)
        tree[newK] = { children: newChildren, bounds: bb };


    return tree;
}

// based on the level we can make the box larger or smaller (in pixel)
function bounds2lines( b, offset ) {
    var tiny = 0.168/800.0;
    // for the first bounding box the offset will be 0
    var c = b.map(function(x, idx) { 
        if (idx < 3)
            return x + (offset * tiny);
        return x - (offset * tiny);
    });

    var ar = [];
    ar.push( [c[0], c[1], c[2]], [c[0], c[1], c[5]]);
    ar.push( [c[0], c[1], c[2]], [c[3], c[1], c[2]]);
    ar.push( [c[0], c[1], c[2]], [c[0], c[4], c[2]]);

    ar.push( [c[3], c[4], c[5]], [c[0], c[4], c[5]]);
    ar.push( [c[3], c[4], c[5]], [c[3], c[1], c[5]]);
    ar.push( [c[3], c[4], c[5]], [c[3], c[4], c[2]]);

    ar.push( [c[0], c[1], c[5]], [c[3], c[1], c[5]]);
    ar.push( [c[0], c[1], c[5]], [c[0], c[4], c[5]]);

    ar.push( [c[0], c[4], c[5]], [c[0], c[4], c[2]]);
    ar.push( [c[3], c[1], c[5]], [c[3], c[1], c[2]]);

    ar.push( [c[0], c[4], c[2]], [c[3], c[4], c[2]]);
    ar.push( [c[3], c[4], c[2]], [c[3], c[1], c[2]]);

    return ar;
}

function occupancy( key, tree, types) {
    // compute and returns the occupancy of this octree node
    if (typeof tree[key] == 'undefined')
        return -1;

    // types should be as long as positions

    var occupancy = { "none": 0, "arterial": 0, "venous": 0 };
    for (var i = 0; i < tree[key].children.length; i++) {
        // compute the number of none, arterious and venous elements
        var node = tree[key].children[i];
        if (types[node] == "") {
            occupancy["none"]++;
        } else if (types[node] == "venous") {
            occupancy["venous"]++;
        } else {
            occupancy['arterial']++;
        }
    }

    return occupancy;
}

function shuffle(array) {
    let currentIndex = array.length;
  
    // While there remain elements to shuffle...
    while (currentIndex != 0) {
  
      // Pick a remaining element...
      let randomIndex = Math.floor(Math.random() * currentIndex);
      currentIndex--;
  
      // And swap it with the current element.
      [array[currentIndex], array[randomIndex]] = [
        array[randomIndex], array[currentIndex]];
    }
}
  

// Compute the probable assignment of nodes based on neighbors with a fixed type
// in array 'types'.
// Output is a probability of being arterial or venous
function diffuse(types, vertices, numConnectionByPosition, diffusionSteps=1) {
    // in probs we have [arterial prob, venous prob, 1==is fixed]
    var probs = Array.from(Array(types.length)).map(function(x,idx) { 
        if (types[idx] == "")
            return [0.0,0.0,0];
        else if (types[idx] == "venous") {
            return [0.0, 1.0, 1];
        } else {
            return [1.0, 0.0, 1];// arterial
        }
    });
    var lr = 0.01;
    var probs_tmp = probs.map(function(a) {
        return [a[0],a[1]];
    });

    // we need to iterate here a while to converge
    // we will not change any fixed once
    for (var iter = 0; iter < types.length*diffusionSteps; iter++) {
        // diffuse along the connections based on edges
        var indices = Array.from(Array(vertices.length)).map(function(a,idx) { return idx; });
        shuffle(indices);
        for (var i_idx = 0; i_idx < vertices.length; i_idx++) { // add random access update
            var i = indices[i_idx];
            // if we have a probs distribute to both sides (only if not fixed)
            var p1 = vertices[i][0]-1; // index starts with 1
            var p2 = vertices[i][1]-1;

            // weigth arterial and venous amounts by number of outgoing connections
            var sA = lr * ((probs_tmp[p1][0] / numConnectionByPosition[p1]) + (probs_tmp[p2][0] / numConnectionByPosition[p2]));
            var sV = lr * ((probs_tmp[p1][1] / numConnectionByPosition[p1]) + (probs_tmp[p2][1] / numConnectionByPosition[p2]));
            if (probs[p1][2] != 1 && sA > 0) { // not fixed
                probs[p1][0] += sA/2.0;
            }
            if (probs[p2][2] != 1 && sA > 0) { // not fixed
                probs[p2][0] += sA/2.0;
            }
            if (probs[p1][2] != 1 && sV > 0) { // not fixed
                probs[p1][1] += sV/2.0; // distribute to both sides
            }
            if (probs[p2][2] != 1 && sV > 0) { // not fixed
                probs[p2][1] += sV/2.0;
            }
        }
        // normalize probs
        for (var i = 0; i < probs.length; i++) {
            // make sure that probabilities sum up to 1 for each node
            var s = (probs[i][0] + probs[i][1])/10.0;
            if (s > 0 && probs[i][2] != 1) { // not fixed 
                probs[i][0] /= s;
                probs[i][1] /= s;
            }
        }


        // switch arrays for next iteration
        probs.map(function(a,idx) {
            probs_tmp[idx] = [a[0],a[1]];
        });
    }
    // normalize probs
    for (var i = 0; i < probs.length; i++) {
        // make sure that probabilities sum up to 1 for each node
        var s = probs[i][0] + probs[i][1];
        if (s > 0 && probs[i][2] != 1) { // not fixed 
            probs[i][0] /= s;
            probs[i][1] /= s;
        }
    }

    return probs;
}

// Compute the number of outgoing connections for every point based on
// vertices. A node with many connections needs to diffuse less material.
function getNumConnectionByPosition(positions, vertices) {
    var numConnectionByPosition = Array.from(Array(positions.length)).map(function(a) { return 0; });
    vertices.map(function(v) {
        numConnectionByPosition[v[0]]++;
        numConnectionByPosition[v[1]]++;
    });

    return numConnectionByPosition;
}

// Return new types given the existing types (fixed values)
// and the probabilites based on vertices (from diffuse).
function sampleProbs(probs, types) {
    var eps = 1e-4;
    return types.map(function(a,idx) {
        // take higher probability entry
        var diff = probs[idx][0] - probs[idx][1];
        if (diff > eps)
            return "arterial";
        else if (Math.abs(diff) <= eps) // if they are both the same take a random entry
            return ""; // (Math.random() > 0.5)?"arterial":"venous";
        else
            return "venous";
    }); // make a copy
}

// add sample probability for each quad tree box
function computeOccupancy(tree, newTypes) {
    var keys = Object.keys(tree);
    for (var i = 0; i < keys.length; i++) {
        var node = tree[keys[i]];
        var probs = [0, 0]; // [ arterious, venous ]
        for (var j = 0; j < node.children.length; j++) {
            var typ = newTypes[node.children[j]];
            if (typ == "arterial")
                probs[0]++;
            else if (typ == "venous")
                probs[1]++;
            else { // an undecided node, count for both, count for nothing?
                console.log("Should not happen, unknown node type");
            }
        }
        var s = probs[0] + probs[1];
        if (s > 0) {
            probs[0] /= s;
            probs[1] /= s;
        }
        tree[keys[i]].probs = probs;
    }
}

function renderGraph(positions, types, vertices, diffusionSteps) {

        // now in probs we have [probability for arterious, probability for venous, is_fixed_prob]
        probs = diffuse(types, vertices, numConnectionByPosition, diffusionSteps);
        var newTypes = sampleProbs(probs, types);

        //
        // DISPLAY
        //

        jQuery('#graph01').children().remove();

        const element = document.querySelector("#graph01");
    
        const options = {
            element: element,
            plugins: ["core", "controls", "cursor"],
            controls: {
                klass: THREE.OrbitControls
            },
        };
        if (root) {
            delete root;
        }
        root = MathBox.mathBox(options);
        if (root.fallback) 
            throw "WebGL not supported";
        
        root.set({ scale: 300, focus: 1 });
        
        var three = root.three;
        three.renderer.setClearColor(new THREE.Color(0xcccccc), 1.0);
        
        camera = root.camera({
            proxy: true,
            position: [0,2,0],
        });
        
        view = root
        .cartesian({
            range: [
                [-0.1, 0.1],
                [-0.1, 0.1],
                [-0.1, 0.1],
            ],
            scale: [1, 1, 1],
        });
    
        view.array({
            items: 1,
            width: newTypes.length,
            channels: 4,
            live: false,
            expr: function(emit, idx) {
                if (newTypes[idx] == "") {
                    emit(.9,.9,.9,1);
                } else if (newTypes[idx] == "venous") {
                    emit(0.1,0.1,1,1);
                } else { // types[idx] == "arterial"
                    emit(1,0,0,1);
                }
            },
            id: 'type'
        });

        view.array({
            items: 1,
            width: positions.length,
            channels: 3,
            live: false,
            data: positions,
            id: "position"
        }).swizzle({
            order: "xyz", // swizzle is done in the data reader
        })
        .transform({
            scale: [1, 1, 1],
            position: [0, 0, 0],
        })
        .point({
            colors: "#type",
            color: 0xffffff, // all colors are pre-multiplied with this color (0..255)
            size: 5
        });
        
        var lines = [];
        for (var i = 0; i < vertices.length; i++) {
            lines.push(positions[vertices[i][0]-1], positions[vertices[i][1]-1]);
        }
        
        view.array({
            width: lines.length / 2,
            items: 2,
            channels: 3,
            live: false,
            data: lines
        })
        .vector({
            color: 0x909091,
            join: "bevel",
            width: 1
        });

        // a sequential colormap
        // var cols = [ [255,255,229,1], [255,247,188,1], [254,227,145,1], [254,196,79,1], [254,153,41,1], [236,112,20,1], [204,76,2,1], [153,52,4,1], [102,37,6,1] ];
        // a diverging colormap (red to blue)
        var cols = [ [165,0,38], [215,48,39],[244,109,67],[253,174,97],[254,224,144],[255,255,191],[224,243,248],[171,217,233],[116,173,209],[69,117,180],[49,54,149] ];
        cols = cols.map(function(a) {
            return [a[0]/255.0, a[1]/255, a[2]/255, 0.5];
        });

        // draw boxes with a given color, based on bounds in the tree
        var boxes = []; // a box is a list of lines (wireframe box)
        var box_colors = [];
        for (var i = 0; i < Object.keys(tree).length; i++) {
            var key = Object.keys(tree)[i];
            var lineset = bounds2lines(tree[key].bounds, key.length-1);
            boxes.push(lineset); // make deeper boxes smaller by length of key
            var occ = occupancy(key, tree, newTypes);
            if (occ["none"] + occ["venous"] + occ["arterial"] == 0) {
                box_colors.push( Array.from(Array(lineset.length)).map(function(a,idx) { return [0.1, 0.1, 0.1, 0.5]; }) );
            } else {
                var part = occ["venous"]/(/*occ["none"] +*/ occ["venous"] + occ["arterial"]);
                if (part > 1)
                    part = 1;
                box_colors.push( Array.from(Array(lineset.length)).map(function(a,idx) { return cols[parseInt(part*(cols.length-1))]; }) ); // colors from 0 to 1
            }
        }
        for (var i = 0; i < box_colors.length; i++) {
            view.array({
                id: 'box_colors' + i,
                channels: 4,
                data: box_colors[i]
            });
        }

        for (var i = 0; i < boxes.length; i++) {
            view.array({
                width: boxes[i].length / 2,
                items: 2,
                channels: 3,
                live: false,
                data: boxes[i]
            }).vector({
                color: 0xffffff,
                colors: '#box_colors' + i,
                join: "bevel",
                width: 0.5
            });
        }
}

function isFullyConnected(vertices, edge2remove) {
    var connectionByPoint = {};
    for (var i = 0; i < vertices.length; i++) {
        if (edge2remove == i)
            continue; // ignore this edge
        var p1 = vertices[i][0];
        var p2 = vertices[i][1];
        if (typeof connectionByPoint[p1] == 'undefined')
            connectionByPoint[p1] = [];
        if (connectionByPoint[p1].indexOf(p2) == -1)    
            connectionByPoint[p1].push(p2);
        if (typeof connectionByPoint[p2] == 'undefined')
            connectionByPoint[p2] = [];
        if (connectionByPoint[p2].indexOf(p1) == -1)    
            connectionByPoint[p2].push(p1);
    }

    // we have all points
    var numPoints = positions.length;
    // do a tree traversal
    var queue = [ vertices[0][0] ]; // start with one point
    var visitedPoints = [ vertices[0][0] ];

    while(queue.length > 0) {
        var node = queue.shift();
        for (var i = 0; i < connectionByPoint[node].length; i++) {
            var np = connectionByPoint[node][i];
            if (visitedPoints.indexOf(np) == -1) {
                queue.push(connectionByPoint[node][i]);
                // we visited this point now
                visitedPoints.push(np);
            }
        }
    }
    if (visitedPoints.length == numPoints)
        return true;
    return false;
}


// we could remove each one and pick the option with the lower summed_occupancy_score
// we should not remove a vertex if the graph is afterwards not connected anymore
function step( vertices, tree) {
    // our tree has some occupancy score before we remove a vertex
    numConnectionByPositionBefore = getNumConnectionByPosition(positions, vertices);

    probs = diffuse(types, vertices, numConnectionByPositionBefore, 1);
    var newTypes = sampleProbs(probs, types);
    computeOccupancy(tree, newTypes);
    var summed_occupancy_score_before = Object.keys(tree).reduce(function(s, a) {
        if (a.length < 2) // ignore the tiny boxes, equal prob tends to 0
            return s + (Math.abs(tree[a].probs[0] - tree[a].probs[1]));
        return s;
    }, 0.0);

    // change a vertex now and compute the overall balance
    var idx_2_remove = Math.floor(Math.random() * (vertices.length-1));

    // is the graph still fully connected without that vertex?
    // if its disconnected now we should have vertices with "" as type (no input)
    if (!isFullyConnected(vertices, idx_2_remove)) { // lets pick another one ++ But we want to have two graphs that are not connected with each other...
        console.log("graph would no longer be fully connected");
        numConnectionByPosition = numConnectionByPositionBefore;
        return vertices; 
    }

    verticesNew = vertices.map(function(a) { return a; }); // make a copy
    verticesNew.splice(idx_2_remove, 1);

    // Use vertices to compute arterial and venous probs
    numConnectionByPositionAfter = getNumConnectionByPosition(positions, vertices);
    probs = diffuse(types, verticesNew, numConnectionByPositionAfter, 1);
    // use probs to update node types (winner takes all)
    var newTypes = sampleProbs(probs, types);
    // use the new types to update quad tree box probs
    computeOccupancy(tree, newTypes);
    // each but the smallest box should have the same probs (0.5, 0.5)
    var box_ids = Object.keys(tree);
    var summed_occupancy_score = box_ids.reduce(function(s, a) {
        if (a.length < 2) // ignore the tiny boxes, equal prob tends to 0
            return s + (Math.abs(tree[a].probs[0] - tree[a].probs[1]));
        return s;
    }, 0.0);
    // minimize occupancy_score
    jQuery('#occupancy_score').text(summed_occupancy_score.toFixed(3) + " (" + summed_occupancy_score_before.toFixed(3) + ")");

    // this is all done in render()
    // use the original types array here (const global)
    // probs = diffuse(types, verticesNew);
    // var newTypes = sampleProbs(probs, types);
    if (summed_occupancy_score <= summed_occupancy_score_before) {
        numConnectionByPosition = numConnectionByPositionAfter;
        return verticesNew;
    }
    numConnectionByPosition = numConnectionByPositionBefore;
    console.log("Removing vertex id: " + idx_2_remove + " would not decrease occupancy score");
    return vertices;
}

var root = null;
var view;
var camera = null;
var tree = {};

var types;
var positions;
var vertices;
var numConnectionByPosition; // array needs to change if a connection has been removed
jQuery(document).ready(function() {

    jQuery('#step').on('click', function() {
        // lets change the graph and draw again
        vertices = step(vertices, tree); // overwrites our current vertices with a new version
        jQuery('#num_vertices').text(vertices.length);
        renderGraph(positions, types, vertices, 1);
    });

    jQuery('#step10').on('click', function() {
        // lets change the graph and draw again
        var numDeleted = 0;
        var max_attempts = 1000;
        var attempts = 0;
        while (numDeleted < 10 && attempts < max_attempts) {
            var numVerticesBefore = vertices.length;
            vertices = step(vertices, tree); // overwrites our current vertices with a new version
            var numVerticesAfter = vertices.length;
            numDeleted += Math.abs(numVerticesBefore - numVerticesAfter);
            attempts++;
        }
        jQuery('#num_vertices').text(vertices.length);
        renderGraph(positions, types, vertices, 1);
    });


    const element = document.querySelector("#graph01");
    
    const options = {
        element: element,
        plugins: ["core", "controls", "cursor"],
        controls: {
            klass: THREE.OrbitControls
        },
    };
    root = MathBox.mathBox(options);
    if (root.fallback) 
        throw "WebGL not supported";
    
    root.set({ scale: 300, focus: 1 });
    
    var three = root.three;
    three.renderer.setClearColor(new THREE.Color(0xcccccc), 1.0);
    
    camera = root.camera({
        proxy: true,
        position: [0, 2, 0],
    });
    
    view = root
    .cartesian({
        range: [
            [-0.1, 0.1],
            [-0.1, 0.1],
            [-0.1, 0.1],
        ],
        scale: [1, 1, 1],
    });
    /*.axis({
        axis: 1,
    })
    .axis({
        axis: 3,
    })
    .axis({
        axis: 2,
    })*/
    //.grid({axes: [1, 3]});
    
    // root.select('axis').set('color', 'black');
    
    var colors = {
        x: 0xff4136, // red
        y: 0xffdc00, // yellow
        z: 0x0074d9, // blue
        xy: 0xff851b, // orange
        xz: 0xb10dc9, // purple
        yz: 0x2ecc40, // green
        xyz: 0x654321, // brown
    };
    
    Promise.all([
        readDataPos(),
        readDataVert()  
    ]).then(function(data) {
        // data should be an array now:
        positions = data[0][0];
        types = data[0][1]; // fixed assignment we start with only one node arterious only one node is venous
        vertices = data[1];
        numConnectionByPosition = getNumConnectionByPosition(positions, vertices);

        tree = getQuadTree(positions, tree, ""); // start with an empty tree
        // expand all level 2 elements, afterwards expand all level 3 elements
        [2, 3].forEach(function(a) {
            for (var i = 0; i < Object.keys(tree).length; i++) {
                var k = Object.keys(tree)[i];
                if (k.length == a)
                    tree = getQuadTree(positions, tree, k); // add more levels
            }    
        });

        // for this initial rendering of the graph we don't want to diffuse anything (0)
        renderGraph(positions, types, vertices, 0);
    });
    
})
