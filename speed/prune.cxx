
#include <boost/filesystem.hpp>
#include <boost/date_time.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <map>
#include <codecvt>
#include <locale> // wstring_convert
#include "json.hpp"
#include <iostream>
#include <ctime>
#include <random>
#include <tbb/parallel_for.h>

#include <filesystem>
namespace fs = std::filesystem;
using namespace boost::filesystem;

namespace po = boost::program_options;
using json = nlohmann::json;

typedef struct node {
    float avx_1;
    float avx_2;
    float avx_3;
    std::string type; // "none", "venous", "arterial"
} node_t;

// a QuadTree is a map (key is path in tree)
typedef struct {
    std::vector<int> children;
    std::array<float, 6> bounds;
    std::array<float, 2> probs;
} tree_node_t;


typedef std::vector< std::pair<int, int> > vertices_t;
typedef std::vector< std::string > types_t;
typedef std::vector< node_t > nodes_t;
typedef std::map<std::string, tree_node_t > tree_t;

std::string remove_file;
json remove_vertices;

void tokenize(std::string const &str, const char delim, std::vector<std::string> &out) {
  size_t start; size_t end = 0;

  while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
    end = str.find(delim, start);
    out.push_back(str.substr(start, end - start));
  }
}

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t\r") {
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}


vertices_t readVertices(std::string edges_file) {
    std::vector< std::pair<int, int> > vertices;

    std::ifstream file(edges_file);
    std::string str;
    std::vector< std::string > header;
    ptrdiff_t EndNodes_1_idx;
    ptrdiff_t EndNodes_2_idx;
    while (std::getline(file, str)) {
        // split by comma
        if (header.size() == 0) {
            // parse the first line as header
            tokenize(str, ',', header);
            EndNodes_1_idx = find(header.begin(), header.end(), "EndNodes_1") - header.begin();
            EndNodes_2_idx = find(header.begin(), header.end(), "EndNodes_2") - header.begin();
            if (EndNodes_1_idx >= header.size()) {
                fprintf(stderr, "Could not find the header entry EndNodes_1.\n");
                exit(-1);
            }
            if (EndNodes_2_idx >= header.size()) {
                fprintf(stderr, "Could not find the header entry EndNodes_2.\n");
                exit(-1);
            }
            continue;
        }
        std::vector< std::string > pieces;
        tokenize(str, ',', pieces);
        if (pieces.size() != header.size()) {
            fprintf(stderr, "strange line [%zu] found as: %s\n", vertices.size(), str.c_str());
        } else {
            vertices.push_back(std::make_pair<int, int>(std::stoi(pieces[EndNodes_1_idx]), std::stoi(pieces[EndNodes_2_idx])));
        }
    }

    return vertices;
}

std::pair<nodes_t, types_t> readNodes(std::string nodes_file) {
    std::vector< node_t > nodes;
    std::vector< std::string > types;

    std::ifstream file(nodes_file);
    std::string str;
    std::vector< std::string > header;
    ptrdiff_t avx_1_idx;
    ptrdiff_t avx_2_idx;
    ptrdiff_t avx_3_idx;
    ptrdiff_t belongroot_idx;

    while (std::getline(file, str)) {
        // remove \r from str
        str = trim(str);

        // split by comma
        if (header.size() == 0) {
            // parse the first line as header
            tokenize(str, ',', header);
            avx_1_idx = find(header.begin(), header.end(), "avx_1") - header.begin();
            avx_2_idx = find(header.begin(), header.end(), "avx_2") - header.begin();
            avx_3_idx = find(header.begin(), header.end(), "avx_3") - header.begin();
            belongroot_idx = find(header.begin(), header.end(), "belongroot") - header.begin();
            if (avx_1_idx >= header.size()) {
                fprintf(stderr, "Could not find the header entry avx_1.\n");
                exit(-1);
            }
            if (avx_2_idx >= header.size()) {
                fprintf(stderr, "Could not find the header entry avx_2.\n");
                exit(-1);
            }
            if (avx_3_idx >= header.size()) {
                fprintf(stderr, "Could not find the header entry avx_3.\n");
                exit(-1);
            }
            if (belongroot_idx >= header.size()) {
                fprintf(stderr, "Could not find the header entry belongroot.\n");
                exit(-1);
            }
            continue;
        }
        std::vector< std::string > pieces;
        tokenize(str, ',', pieces);
        if (pieces.size() != header.size()) {
            fprintf(stderr, "strange line [%zu] found as: %s\n", nodes.size(), str.c_str());
        } else {
            // we would like to map belongroot to float from string
            node_t entry = { std::stof(pieces[avx_1_idx]), std::stof(pieces[avx_2_idx]), std::stof(pieces[avx_3_idx]), pieces[belongroot_idx] };
            nodes.push_back(entry);
            types.push_back( pieces[belongroot_idx] );
        }
    }
    auto erg = std::make_pair(nodes, types);
    return erg;
}

bool inside(std::array<float, 6> bounds, node_t p) {
    if (p.avx_1 >= bounds[0] && p.avx_1 <= bounds[3] && p.avx_2 >= bounds[1] && p.avx_2 <= bounds[4] && p.avx_3 >= bounds[2] && p.avx_3 <= bounds[5])
      return true;
    return false;
}

std::array<float, 6> computeBounds(std::vector< node_t > nodes) {
    if (nodes.size() == 0) {
        fprintf(stderr, "Error: cannot compute bounds of nothing.\n");
        exit(-1);
    }
    float increase_by = 1.0/10.0; // in percent
    std::array<float, 6> bounds = {nodes[0].avx_1,nodes[0].avx_2,nodes[0].avx_3,nodes[0].avx_1,nodes[0].avx_2,nodes[0].avx_3};
    for (int i = 0; i < nodes.size(); i++) {
        if (nodes[i].avx_1 < bounds[0])
            bounds[0] = nodes[i].avx_1;
        if (nodes[i].avx_2 < bounds[1])
            bounds[1] = nodes[i].avx_2;
        if (nodes[i].avx_3 < bounds[2])
            bounds[2] = nodes[i].avx_3;
        if (nodes[i].avx_1 > bounds[3])
            bounds[3] = nodes[i].avx_1;
        if (nodes[i].avx_2 > bounds[4])
            bounds[4] = nodes[i].avx_2;
        if (nodes[i].avx_3 > bounds[5])
            bounds[5] = nodes[i].avx_3;
    }
    float size = std::min(std::min(bounds[3]-bounds[0],bounds[4]-bounds[1]),bounds[5]-bounds[2]);
    bounds[0] -= size*increase_by;
    bounds[1] -= size*increase_by;
    bounds[2] -= size*increase_by;
    bounds[3] += size*increase_by;
    bounds[4] += size*increase_by;
    bounds[5] += size*increase_by;

    return bounds;
}


tree_t getQuadTree(std::vector< node_t > nodes, tree_t tree, std::string startLevel ) {
    // if we have an empty tree
    if (tree.begin() == tree.end()) {
        // add a first node
        tree_node_t entry;
        for (int i = 0; i < nodes.size(); i++) {
            entry.children.push_back(i);
        }
        entry.bounds = computeBounds(nodes);
        tree.insert({"0", entry});
        startLevel = "0";
    }

    auto indices = tree[startLevel].children;
    auto bounds = tree[startLevel].bounds;
    std::vector<float> widths = {bounds[3]-bounds[0], bounds[4]-bounds[1], bounds[5]-bounds[2]};

    std::string newK = startLevel + "0";
    std::vector<int> newChildren; // array with indices
    std::array<float, 6> bb = { bounds[0], bounds[1], bounds[2], 0, 0, 0}; // start index
    bb[3] = bb[0] + widths[0]/2.0;
    bb[4] = bb[1] + widths[1]/2.0;
    bb[5] = bb[2] + widths[2]/2.0;
    for (int i = 0; i < indices.size(); i++) {
        auto p = nodes[indices[i]];
        if (inside(bb, p)) {
            newChildren.push_back(indices[i]);
        }
    }
    if (newChildren.size() > 0) {
        tree_node_t entry;
        entry.children = newChildren;
        entry.bounds = bb;
        tree.insert({newK, entry});
    }

    newK = startLevel + "1";
    newChildren.clear(); // array with indices
    bb = { bounds[0] + widths[0]/2, bounds[1], bounds[2], 0, 0, 0}; // start index
    bb[3] = bb[0] + widths[0]/2.0;
    bb[4] = bb[1] + widths[1]/2.0;
    bb[5] = bb[2] + widths[2]/2.0;
    for (int i = 0; i < indices.size(); i++) {
        auto p = nodes[indices[i]];
        if (inside(bb, p)) {
            newChildren.push_back(indices[i]);
        }
    }
    if (newChildren.size() > 0) {
        tree_node_t entry;
        entry.children = newChildren;
        entry.bounds = bb;
        tree.insert({newK, entry});
    }

    newK = startLevel + "2";
    newChildren.clear(); // array with indices
    bb = { bounds[0], bounds[1] + widths[1]/2, bounds[2], 0, 0, 0}; // start index
    bb[3] = bb[0] + widths[0]/2.0;
    bb[4] = bb[1] + widths[1]/2.0;
    bb[5] = bb[2] + widths[2]/2.0;
    for (int i = 0; i < indices.size(); i++) {
        auto p = nodes[indices[i]];
        if (inside(bb, p)) {
            newChildren.push_back(indices[i]);
        }
    }
    if (newChildren.size() > 0) {
        tree_node_t entry;
        entry.children = newChildren;
        entry.bounds = bb;
        tree.insert({newK, entry});
    }

    newK = startLevel + "3";
    newChildren.clear(); // array with indices
    bb = { bounds[0] + widths[0]/2, bounds[1] + widths[1]/2, bounds[2], 0, 0, 0}; // start index
    bb[3] = bb[0] + widths[0]/2.0;
    bb[4] = bb[1] + widths[1]/2.0;
    bb[5] = bb[2] + widths[2]/2.0;
    for (int i = 0; i < indices.size(); i++) {
        auto p = nodes[indices[i]];
        if (inside(bb, p)) {
            newChildren.push_back(indices[i]);
        }
    }
    if (newChildren.size() > 0) {
        tree_node_t entry;
        entry.children = newChildren;
        entry.bounds = bb;
        tree.insert({newK, entry});
    }

    newK = startLevel + "4";
    newChildren.clear(); // array with indices
    bb = { bounds[0], bounds[1], bounds[2] + widths[2]/2, 0, 0, 0}; // start index
    bb[3] = bb[0] + widths[0]/2.0;
    bb[4] = bb[1] + widths[1]/2.0;
    bb[5] = bb[2] + widths[2]/2.0;
    for (int i = 0; i < indices.size(); i++) {
        auto p = nodes[indices[i]];
        if (inside(bb, p)) {
            newChildren.push_back(indices[i]);
        }
    }
    if (newChildren.size() > 0) {
        tree_node_t entry;
        entry.children = newChildren;
        entry.bounds = bb;
        tree.insert({newK, entry});
    }

    newK = startLevel + "5";
    newChildren.clear(); // array with indices
    bb = { bounds[0] + widths[0]/2, bounds[1], bounds[2] + widths[2]/2, 0, 0, 0}; // start index
    bb[3] = bb[0] + widths[0]/2.0;
    bb[4] = bb[1] + widths[1]/2.0;
    bb[5] = bb[2] + widths[2]/2.0;
    for (int i = 0; i < indices.size(); i++) {
        auto p = nodes[indices[i]];
        if (inside(bb, p)) {
            newChildren.push_back(indices[i]);
        }
    }
    if (newChildren.size() > 0) {
        tree_node_t entry;
        entry.children = newChildren;
        entry.bounds = bb;
        tree.insert({newK, entry});
    }

    newK = startLevel + "6";
    newChildren.clear(); // array with indices
    bb = { bounds[0], bounds[1] + widths[1]/2, bounds[2] + widths[2]/2, 0, 0, 0}; // start index
    bb[3] = bb[0] + widths[0]/2.0;
    bb[4] = bb[1] + widths[1]/2.0;
    bb[5] = bb[2] + widths[2]/2.0;
    for (int i = 0; i < indices.size(); i++) {
        auto p = nodes[indices[i]];
        if (inside(bb, p)) {
            newChildren.push_back(indices[i]);
        }
    }
    if (newChildren.size() > 0) {
        tree_node_t entry;
        entry.children = newChildren;
        entry.bounds = bb;
        tree.insert({newK, entry});
    }

    newK = startLevel + "7";
    newChildren.clear(); // array with indices
    bb = { bounds[0] + widths[0]/2, bounds[1] + widths[1]/2, bounds[2] + widths[2]/2, 0, 0, 0}; // start index
    bb[3] = bb[0] + widths[0]/2.0;
    bb[4] = bb[1] + widths[1]/2.0;
    bb[5] = bb[2] + widths[2]/2.0;
    for (int i = 0; i < indices.size(); i++) {
        auto p = nodes[indices[i]];
        if (inside(bb, p)) {
            newChildren.push_back(indices[i]);
        }
    }
    if (newChildren.size() > 0) {
        tree_node_t entry;
        entry.children = newChildren;
        entry.bounds = bb;
        tree.insert({newK, entry});
    }

    return tree;
}

// no side-effects
std::map<int, int> getNumConnectionByPosition( std::vector< node_t > positions, std::vector< std::pair<int, int> > vertices) {
    std::map<int, int> numConnectionByPosition;
    for (int i = 0; i < positions.size(); i++) {
        numConnectionByPosition.insert({i, 0});
    }
    for (int i = 0; i < vertices.size(); i++) {
        numConnectionByPosition[vertices[i].first]++;
        numConnectionByPosition[vertices[i].second]++;
    }
    return numConnectionByPosition;
}

void shuffle(std::vector<int> *array, std::mt19937 rgen) {
    // seed the random number generator only once
    //std::random_device rdev;
    //std::mt19937 rgen(rdev());

    int currentIndex = array->size();
    while (currentIndex != 0) {
        std::uniform_int_distribution<int> idist(0,currentIndex-1); // including, including
        currentIndex--;
        //std::uniform_real_distribution<int> distribution(0,currentIndex-1);
        int randomIndex = (int)idist(rgen);
        int tmp = (*array)[currentIndex];
        (*array)[currentIndex] = (*array)[randomIndex];
        (*array)[randomIndex] = tmp;
    }
}

// estimate the expected time (number of iterations) to reach a specific target
class estimateExponential {
    public:
        std::vector< std::pair<float, float> > _data;
        float curveA = 0.0f;
        float curveB = 0.0f;
        void insert(float x, float y) { // store only the last 10 results
            _data.push_back(std::make_pair(x, y));
            if (_data.size() > 30) {
                _data.erase(_data.begin());
            }
        }
        void clear() {
            _data.resize(0);
        }
        void fit() {
            int n,i;
            n = _data.size();
            std::vector<float> Y(n);
            float sumx=0,sumy=0,sumxy=0,sumx2=0;
            float a,b,A;
            for (i=0; i <= n-1; i++) {
                Y[i] = log(_data[i].second);
            }
            for (i=0; i <= n-1; i++) {
                sumx += _data[i].first;
                sumx2 += _data[i].first * _data[i].first;
                sumy += Y[i];
                sumxy += _data[i].first * Y[i];
            }
            A = ((sumx2*sumy - sumx*sumxy)*1.0/(n*sumx2 - sumx*sumx)*1.0);
            b = ((n*sumxy - sumx*sumy)*1.0/(n*sumx2 - sumx*sumx)*1.0);
            a = exp(A);

            curveA = a;
            curveB = b;
            //printf("\n\n The curve is Y= %4.3fe^%4.3fX",a,b);
        }

        float estimateY(float x) {
            fit();
            return curveA* exp(curveB * x);
        }
        float estimateX(float y) {
            fit();
            return (std::log(y) - std::log(curveA)) / curveB;
        }
};


// diffusionSteps is a part of the number of vertices so 0.5 is half the vertice many iterations
std::vector< std::array<float, 3> > diffuse( types_t types, vertices_t vertices, std::map<int, int> numConnectionByPosition, float maxChange) {

    boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();

    std::vector < std::array<float, 3> > probs;
    probs.reserve(types.size());
    for (int i = 0; i < types.size(); i++) {
        if (types[i] == "" || types[i] == "none") {
            probs.push_back( std::array<float, 3>({0.0, 0.0, 0.0}));
        } else if (types[i] == "venous") {
            probs.push_back( std::array<float, 3>({0.0, 1.0, 1.0}));
        } else {
            probs.push_back( std::array<float, 3>({1.0, 0.0, 1.0}));
        }
    }
    float lr = 0.01;
    std::vector< std::array<float, 2> > probs_tmp(probs.size());
    for (int i = 0; i < probs.size(); i++) {
        probs_tmp[i] = std::array<float, 2>({probs[i][0], probs[i][1]});
    }

    std::random_device rdev;
    std::mt19937 rgen(rdev());

    std::vector<int> indices(vertices.size());
    //for (int i = 0; i < indices.size(); i++)
    //    indices[i] = i;
    //shuffle(&indices, rgen); // everything is slower if the indices are randomized in memory

    std::vector<float> sAs(vertices.size(), 0.0); // init with 0
    std::vector<float> sVs(vertices.size(), 0.0);

    double summed_change = 0.0; // change per iteration
    unsigned int num_valued = 0; // number of probs that are not 0/0
    estimateExponential est; // estimate the time for finishing

    int numIters = types.size(); // * diffusionSteps;
    for (int iter = 0; iter < numIters; iter++) {
        if (iter % 10 == 0) {
            int guessNumIters = 0;
            if (iter > 0) {
                guessNumIters = (int)est.estimateX(maxChange);
                //fprintf(stderr, "Estimate for number of iteration required is now: %d\n", guessNumIters);
            }
            boost::posix_time::ptime timeLocalEnd = boost::posix_time::microsec_clock::local_time();
            boost::posix_time::time_period tp(timeLocal, timeLocalEnd);
            std::string time_left = boost::posix_time::to_simple_string( (timeLocalEnd - timeLocal)/(float)iter * (float)(guessNumIters-iter) );
            fprintf(stdout, "\033[0G\033[K%d/%d [ETA: %s, Δ%.04f (%.02f%%)]", iter, numIters, time_left.substr(0,8).c_str(), summed_change, 100*((float)num_valued/probs.size()));
            fflush(stdout);
            if (iter > 10 && summed_change < maxChange)
                break; // stop the loop
        }
        for (int i = 0; i < indices.size(); i++)
            indices[i] = i;
        shuffle(&indices, rgen); // TODO: jumping around in memory is slow, better to sort before running

        float sA;
        float sV;
        int p1, p2;
        int i;

        // first compute the sA and sV for all pairs of vertices (can be done in parallel)
        tbb::parallel_for( tbb::blocked_range<int>(0,vertices.size()),
                       [&](tbb::blocked_range<int> r) {
            for (int i_idx=r.begin(); i_idx<r.end(); ++i_idx) {
                int i = indices[i_idx];
                int p1 = vertices[i].first-1;
                int p2 = vertices[i].second-1;

                // store half the value as needed in next loop
                sAs[i_idx] = 0.5 * lr * ((probs_tmp[p1][0] / numConnectionByPosition[p1]) + (probs_tmp[p2][0] / numConnectionByPosition[p2]));
                sVs[i_idx] = 0.5 * lr * ((probs_tmp[p1][1] / numConnectionByPosition[p1]) + (probs_tmp[p2][1] / numConnectionByPosition[p2]));
            }
        });

        /*for (int i_idx = 0; i_idx < vertices.size(); i_idx++) {
            i = indices[i_idx];
            p1 = vertices[i].first-1;
            p2 = vertices[i].second-1;

            sAs[i_idx] = lr * ((probs_tmp[p1][0] / numConnectionByPosition[p1]) + (probs_tmp[p2][0] / numConnectionByPosition[p2]));
            sVs[i_idx] = lr * ((probs_tmp[p1][1] / numConnectionByPosition[p1]) + (probs_tmp[p2][1] / numConnectionByPosition[p2]));
        }*/

        // TODO: parallelize the next loop
        // we could collect all sA per vertex and sum in parallel
        // We just need sufficient space to store the max number of connections by node and
        // initialize that space with 0.0.
        for (int i_idx = 0; i_idx < vertices.size(); i_idx++) {
            i = indices[i_idx];
            p1 = vertices[i].first-1;
            p2 = vertices[i].second-1;

            sA = sAs[i_idx]; // lr * ((probs_tmp[p1][0] / numConnectionByPosition[p1]) + (probs_tmp[p2][0] / numConnectionByPosition[p2]));
            sV = sVs[i_idx]; // lr * ((probs_tmp[p1][1] / numConnectionByPosition[p1]) + (probs_tmp[p2][1] / numConnectionByPosition[p2]));
            if (probs[p1][2] != 1 && sA > 0) { // not fixed
                probs[p1][0] += sA;
            }
            if (probs[p2][2] != 1 && sA > 0) { // not fixed
                probs[p2][0] += sA;
            }
            if (probs[p1][2] != 1 && sV > 0) { // not fixed
                probs[p1][1] += sV; // distribute to both sides
            }
            if (probs[p2][2] != 1 && sV > 0) { // not fixed
                probs[p2][1] += sV;
            }
        }

        // normalize probs sometimes
        //if (iter % 10 == 0) {
            /*tbb::parallel_for( tbb::blocked_range<int>(0,probs.size()),
                        [&](tbb::blocked_range<int> r) {
                for (int i=r.begin(); i<r.end(); ++i) {
                    float s = (probs[i][0] + probs[i][1])/10.0;
                    if (s > 0 && probs[i][2] != 1) { // not fixed
                        probs[i][0] /= s;
                        probs[i][1] /= s;
                    }
                }
            }); */

        tbb::parallel_for( tbb::blocked_range<int>(0,probs.size()),
                       [&](tbb::blocked_range<int> r) {
            for (int i=r.begin(); i<r.end(); ++i) {
                // make sure that probabilities sum up to 1 for each node
                float s = (probs[i][0] + probs[i][1])*10.0;
                if (s > 0 && probs[i][2] != 1) { // not fixed
                    probs[i][0] /= s;
                    probs[i][1] /= s;
                }
            }
        });
        //}

        /*tbb::parallel_for( tbb::blocked_range<int>(0,probs.size()),
                    [&](tbb::blocked_range<int> r) {
            for (int idx=r.begin(); idx<r.end(); ++idx) {
              probs_tmp[idx][0] = probs[idx][0];
              probs_tmp[idx][1] = probs[idx][1];
            }
        });*/

        summed_change = 0.0;
        num_valued = 0;
        for (int idx = 0; idx < probs.size(); idx++) {
            summed_change += std::abs(probs[idx][0]-probs_tmp[idx][0]) + std::abs(probs[idx][1]-probs_tmp[idx][1]);
            if (probs[idx][0] != 0 || probs[idx][1] != 0)
                num_valued++;
            probs_tmp[idx][0] = probs[idx][0];
            probs_tmp[idx][1] = probs[idx][1];
        }
        est.insert(iter, summed_change);
    }
    for (int i = 0; i < probs.size(); i++) {
        // make sure that probabilities sum up for each node
        float s = (probs[i][0] + probs[i][1]);
        if (s > 0 && probs[i][2] != 1) { // not fixed
            probs[i][0] /= s;
            probs[i][1] /= s;
        }
    }
    fprintf(stdout, "\n");

    return probs;
}

bool is_approximately_equal(float a, float b) {
  float epsilon = std::numeric_limits<float>::epsilon();
  float scale = std::max(abs(a), abs(b));
  return std::abs(a - b) <= scale * (2*epsilon);
}

// no side-effects
std::vector< std::string > sampleProbs( std::vector< std::array<float, 3> > probs, std::vector< std::string > types) {
    std::vector< std::string > ps(types.size());
    for (int idx = 0; idx < types.size(); idx++) {
        float diff = probs[idx][0] - probs[idx][1];
        if (is_approximately_equal(probs[idx][0], probs[idx][1])) {
            ps[idx] = ""; // undecided
        } else if (diff > 0) {
            ps[idx] = "arterial";
        } else {
            ps[idx] = "venous";
        }

        /*        if (diff > 2.0*eps) {
            ps[idx] = "arterial";
        } else if (std::abs(diff) <= eps) {
            ps[idx] = ""; // undecided
        } else {
            ps[idx] = "venous";
        }*/
    }
    return ps;
}

// side effect of changing tree[].probs
void computeOccupancy( std::map<std::string, tree_node_t > &tree, std::vector< std::string > types) {
    // add the occupancy score to the tree nodes
    std::vector< std::string > keys;
    std::map<std::string, tree_node_t>::const_iterator iter = tree.begin();
    while (iter != tree.end()) {
        keys.push_back(iter->first);
        ++iter;
    }
    std::array<float, 2> probs({0.0,0.0});
    std::map<std::string, std::array<float, 2> > old_probs;
    for (int i = 0; i < keys.size(); i++) {
        tree_node_t node = tree[keys[i]];
        probs[0] = 0.0;
        probs[1] = 0.0;
        for (int j = 0; j < node.children.size(); j++) {
            std::string typ = types[node.children[j]];
            if (typ == "arterial")
                probs[0]++;
            else if (typ == "venous")
                probs[1]++;
            //else
            //    fprintf(stderr, "computeOccupancy, unknown node type \"%s\"\n", typ.c_str());
        }
        float s = probs[0] + probs[1];
        if (s > 0) {
            probs[0] /= s;
            probs[1] /= s;
        }
        old_probs.insert( std::make_pair(keys[i], std::array<float, 2>({tree[keys[i]].probs[0],tree[keys[i]].probs[1]}) ));
        tree[keys[i]].probs[0] = probs[0];
        tree[keys[i]].probs[1] = probs[1];
    }
    // can we draw a bit what this looks like now?
    // Lets say the probs for the first couple of levels?
    std::sort(keys.begin(), keys.end(), [](const std::string& first, const std::string& second){
        return first.size() < second.size();
    });

    int l = 1;
    int ll = 0;
    for (int i = 0; i < keys.size(); i++) {
        if (keys[i].size() > 4)
            continue; // ignore
        if (l < keys[i].size()) {
            fprintf(stdout, " "); // add a space whenever we change the scale
            l = keys[i].size();
            if (keys[i].size() == 4) { // if we are at the last level start a new line
                fprintf(stdout, "\n");
                ll = i-1;
            }
        }
        float a = tree[keys[i]].probs[0]; // arterial
        float b = tree[keys[i]].probs[1]; // venous
        std::string col("124m"); // red
        if (b > a)
            col = std::string("33m"); // blue
        if (is_approximately_equal(a, b)) {
            col = std::string("244m");
        }
        fprintf(stdout, "\033[0m\033[38;5;%s■\033[0m", col.c_str());
        if (keys[i].size() == 4 && (i-ll) > 0 && ((i-ll) % 80) == 0) {
            fprintf(stdout, "\n");
        }
    }
    // if we compare the old probabilities with the new probabilities, how many switches do we observe?
    int num_switches = 0;
    for (int i = 0; i < keys.size(); i++) {
        float a = tree[keys[i]].probs[0] - tree[keys[i]].probs[1];
        float b = old_probs[keys[i]][0] - old_probs[keys[i]][1];
        // change in sign?
        if ( ((0.0f < a) - (a < 0.0f)) != ((0.0f < b) - (b < 0.0f)) )
            num_switches++;
    }
    fprintf(stdout, " %d switches\n", num_switches);
}

bool getEdgeCandidate(vertices_t vertices, int *edge2remove) {
    std::map<int, std::vector<int> > connectionByPoint;

   /* var connectionByPoint = {}; */
    for (int i = 0; i < vertices.size(); i++) {
        int p1 = vertices[i].first-1;
        int p2 = vertices[i].second-1;
        if (connectionByPoint.find(p1) == connectionByPoint.end())
            connectionByPoint.insert(std::make_pair(p1, std::vector<int>()));
        if (std::find(connectionByPoint[p1].begin(), connectionByPoint[p1].end(), p2) == std::end(connectionByPoint[p1]))
            connectionByPoint[p1].push_back(p2);
        if (connectionByPoint.find(p2) == connectionByPoint.end())
            connectionByPoint.insert(std::make_pair(p2, std::vector<int>()));
        if (std::find(connectionByPoint[p2].begin(), connectionByPoint[p2].end(), p1) == std::end(connectionByPoint[p2]))
            connectionByPoint[p2].push_back(p1);
    }

    // one short cut would be if one of the two points do not have other neighbors anymore,
    // in that case we know its now disconnected
    std::vector<int> edgeCandidates;
    for (int i = 0; i < vertices.size(); i++) {
        int p1 = vertices[i].first-1;
        int p2 = vertices[i].second-1;
        if (connectionByPoint[p1].size() != 1 && connectionByPoint[p2].size() != 1)
            edgeCandidates.push_back(i);
    }
    if (edgeCandidates.size() > 0) {
        // return a random edge candidate
        std::srand(std::time(nullptr));
        (*edge2remove) = edgeCandidates[std::rand() % edgeCandidates.size()];
        return true;
    }
    return false; // none found
}

// TODO: it would be better to compute the partitions and return those
// Would be fully connected if the edge2remove is removed.
bool isFullyConnected(vertices_t vertices, nodes_t nodes, int edge2remove) {
    std::map<int, std::vector<int> > connectionByPoint;

   /* var connectionByPoint = {}; */
    for (int i = 0; i < vertices.size(); i++) {
        if (edge2remove == i)
            continue; // ignore this edge
        int p1 = vertices[i].first-1;
        int p2 = vertices[i].second-1;
        if (connectionByPoint.find(p1) == connectionByPoint.end())
            connectionByPoint.insert(std::make_pair(p1, std::vector<int>()));
        if (std::find(connectionByPoint[p1].begin(), connectionByPoint[p1].end(), p2) == std::end(connectionByPoint[p1]))
            connectionByPoint[p1].push_back(p2);
        if (connectionByPoint.find(p2) == connectionByPoint.end())
            connectionByPoint.insert(std::make_pair(p2, std::vector<int>()));
        if (std::find(connectionByPoint[p2].begin(), connectionByPoint[p2].end(), p1) == std::end(connectionByPoint[p2]))
            connectionByPoint[p2].push_back(p1);
    }

    // one short cut would be if one of the two points do not have other neighbors anymore,
    // in that case we know its now disconnected
    int p1 = vertices[edge2remove].first-1;
    int p2 = vertices[edge2remove].second-1;
    if (connectionByPoint[p1].size() == 1 || connectionByPoint[p2].size() == 1)
        return false; // removing that edge would orphan one node == disconnect the graph

    // we have all points
    int numPoints = nodes.size();
    // do a tree traversal
    std::vector<int> queue;
    queue.reserve(nodes.size());
    queue.push_back(vertices[0].first-1); // start with one point
    std::set<int> visitedPoints;
    //visitedPoints.reserve(nodes.size());
    visitedPoints.insert(vertices[0].first-1);

    while(queue.size() > 0) {
        int node = queue.back();
        queue.pop_back();
        for (int i = 0; i < connectionByPoint[node].size(); i++) {
            int np = connectionByPoint[node][i];
            if (visitedPoints.find(np) == visitedPoints.end()) {
                queue.push_back(connectionByPoint[node][i]);
                // we visited this point now
                visitedPoints.insert(np);
            }
        }
    }
    if (visitedPoints.size() == numPoints)
        return true;
    return false;
}

std::pair<vertices_t, types_t> step( nodes_t nodes, vertices_t vertices, types_t types, tree_t tree, std::map< int, int> &numConnectionByPosition, float stop) {

    //auto numConnectionByPositionBefore = getNumConnectionByPosition(nodes, vertices); // updated vertices
    auto numConnectionByPositionBefore = numConnectionByPosition;

    //auto probs = diffuse(types, vertices, numConnectionByPositionBefore, 0.01); // need the same types here all the time (not updated but initials)

    //auto newTypes = sampleProbs(probs, types);

    //computeOccupancy(tree, newTypes); // fill in the tree box probabilities

    // compute the overall score for our octree
    float summed_occupancy_score_before = 0.0;
    {
        std::vector< std::string > keys;
        std::map<std::string, tree_node_t>::const_iterator iter = tree.begin();
        while (iter != tree.end()) {
            keys.push_back(iter->first);
            ++iter;
        }
        for (int i = 0; i < keys.size(); i++) {
            summed_occupancy_score_before += std::abs(tree[keys[i]].probs[0] - tree[keys[i]].probs[1]);
        }
    }

    int idx_2_remove = -1;
    int max_attempts = 100; // how often to we try to find a vertex that after removing keeps a fully connected graph
    int attempt = 0;
    while (idx_2_remove == -1) {
        if (!getEdgeCandidate(vertices, &idx_2_remove)) { // randomly draw an edge
            fprintf(stderr, "Removing any edge would split this graph... giving up.\n");
            // but write out the resulting types for all nodes
            exit(-1);
        }
        if (!isFullyConnected(vertices, nodes, idx_2_remove)) { // lets pick another one ++ But we want to have two graphs that are not connected with each other...
            fprintf(stdout, "removing vertex %d [%d<->%d] makes graph disconnected.. [%zu]\n", idx_2_remove, vertices[idx_2_remove].first-1, vertices[idx_2_remove].second-1, vertices.size());

            numConnectionByPosition = numConnectionByPositionBefore;
            idx_2_remove = -1;
        }
        if (attempt > max_attempts) {
            fprintf(stdout, "Info: Max number of attempts to find an edge candidate, giving up.\n");
            return std::make_pair(vertices, types);
        }
    }

    auto verticesNew = vertices;
    verticesNew.erase(verticesNew.begin() + idx_2_remove);

    auto numConnectionByPositionAfter = getNumConnectionByPosition(nodes, verticesNew);
    auto probs = diffuse(types, verticesNew, numConnectionByPositionAfter, stop);

    auto newTypes = sampleProbs(probs, types);
    // make a copy of the old probs in the tree first
    std::map<std::string, std::array<float, 2> > old_probs;
    std::map<std::string, tree_node_t>::const_iterator iter = tree.begin();
    while (iter != tree.end()) {
        old_probs.insert(std::make_pair(iter->first, std::array<float, 2>({ iter->second.probs[0], iter->second.probs[1] })));
        ++iter;
    }

    computeOccupancy(tree, newTypes); // overwrites probs in tree

    float summed_occupancy_score = 0.0;
    {
        std::vector< std::string > keys;
        std::map<std::string, tree_node_t>::const_iterator iter = tree.begin();
        while (iter != tree.end()) {
            keys.push_back(iter->first);
            ++iter;
        }
        for (int i = 0; i < keys.size(); i++) {
            summed_occupancy_score += std::abs(tree[keys[i]].probs[0] - tree[keys[i]].probs[1]);
        }
    }

    if (summed_occupancy_score <= summed_occupancy_score_before) {
        fprintf(stdout, "remove vertex %d as it makes our graph more or equally balanced.. [%zu, score: %.04f, %.04f]\n", idx_2_remove, verticesNew.size(), summed_occupancy_score, summed_occupancy_score_before-summed_occupancy_score);
        numConnectionByPosition = numConnectionByPositionAfter;
        remove_vertices.push_back(idx_2_remove);
        if (boost::filesystem::exists(remove_file)) {
          std::ofstream out(remove_file);
          out << remove_vertices;
        }
        return std::make_pair(verticesNew, newTypes);
    }

    fprintf(stdout, "removing vertex %d would make our graph less balanced, undo now [%zu]\n", idx_2_remove, verticesNew.size());
    // get our old probs back in the tree for the next iteration
    iter = tree.begin();
    while (iter != tree.end()) {
        tree[iter->first].probs[0] = old_probs[iter->first][0]; // .insert(std::make_pair(iter->first, std::array<float, 2>({ iter->second.probs[0], iter->second.probs[1] })));
        tree[iter->first].probs[1] = old_probs[iter->first][1]; // .insert(std::make_pair(iter->first, std::array<float, 2>({ iter->second.probs[0], iter->second.probs[1] })));
        ++iter;
    }

    return std::make_pair(vertices, types);
}

int main(int argc, char *argv[]) {
  setlocale(LC_NUMERIC, "en_US.utf-8");

  boost::posix_time::ptime timeLocal = boost::posix_time::microsec_clock::local_time();

  bool verbose = false;
  std::string output("");
  std::string nodes_file("");
  std::string edges_file("");
  remove_file = std::string("");
  int steps = 5;
  float stop = 0.3; // maximum change in a diffusion step

  po::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "prune: Run a diffusion quad-tree pruning on a given graph.")
      ("version,V", "Print the version number.")
      ("verbose,v", po::bool_switch(&verbose), "Print more verbose output during processing.")
      ("edges,e", po::value< std::string >(&edges_file), "The edges csv file.")
      ("nodes,n", po::value< std::string >(&nodes_file), "The nodes csv file.")
      ("output,o", po::value<std::string>(&output), "Path to output csv file.")
      ("stopping,c", po::value<float>(&stop), "When to assume diffusion solution has converged based on summed overall change [default 0.3].")
      ("remove_edges,r", po::value<std::string>(&remove_file), "Path to a json that contains edges that should be removed initially. This file will be updated continuously.")
      ("steps,s", po::value<int>(&steps), "Number of connections to remove.")
  ;
  // allow positional arguments to map to rawdata
  po::positional_options_description p;
  p.add("raw-data", -1);

  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(p).run(), vm);
    po::notify(vm);
  } catch(std::exception& e) {
        std::cout << e.what() << "\n";
        return 1;
  }

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  std::string versionString = std::string("0.0.1.") + boost::replace_all_copy(std::string(__DATE__), " ", ".");
  if (vm.count("version")) {
    fprintf(stdout, "version: %s\n", versionString.c_str());
    return 0;
  }
  if (edges_file.empty() || nodes_file.empty()) {
    std::cout << desc << std::endl;
    return 0;
  }

  if (!boost::filesystem::exists(edges_file)) {
    fprintf(stderr, "Error: file for edges (\"%s\") does not exist.\n", edges_file.c_str());
    exit(-1);
  }
  if (!boost::filesystem::exists(nodes_file)) {
    fprintf(stderr, "Error: file for nodes (\"%s\") does not exist.\n", nodes_file.c_str());
    exit(-1);
  }

  auto vertices = readVertices(edges_file);
  fprintf(stdout, " reading edges [%zu]\n", vertices.size());

  auto nodes_types = readNodes(nodes_file);
  auto nodes = nodes_types.first;
  auto types = nodes_types.second;
  fprintf(stdout, " reading nodes [%zu], types [%zu]\n", nodes.size(), types.size());

  // remove some vertices if we have a remove_file
  // We need to make this work concurrently with other programs running the same algorithm.
  // But this might be problematic. Two runs can remove conflicting vertices (disconnecting
  // a part of the graph). We would need a process like a blockchain. Programs with many
  // removals can add to the chain and extend it.
  if (boost::filesystem::exists(remove_file)) {
    std::ifstream ifs(remove_file);
    remove_vertices = json::parse(ifs);
    fprintf(stderr, "Remove known edges before processing [%zu-%zu=%zu]\n", vertices.size(), remove_vertices.size(), vertices.size()-remove_vertices.size());
    for (int i = 0; i < remove_vertices.size(); i++) {
        vertices.erase(vertices.begin() + remove_vertices[i]);
    }
  }

  // compute a QuadTree
  std::map<std::string, tree_node_t> tree; // start with an empty tree
  tree = getQuadTree( nodes, tree, "" );
  for (int level = 2; level <= 3; level++) {
    // what nodes do we have for level k?
    std::vector<std::string> keys;
    auto treeIter = tree.begin();
    while (treeIter != tree.end()) {
        if (treeIter->first.size() == level) { // extend all keys of length level
        keys.push_back(treeIter->first);
           //tree = getQuadTree(nodes, tree, treeIter->first);
        }
       ++treeIter;
    }
    for (int i = 0; i < keys.size(); i++) {
        tree = getQuadTree(nodes, tree, keys[i]);
    }
  }

  types_t newTypes;
  auto numConnectionByPosition = getNumConnectionByPosition( nodes, vertices);
  auto probs = diffuse(types, vertices, numConnectionByPosition, stop); // need the same types here all the time (not updated but initials)
  newTypes = sampleProbs(probs, types);
  computeOccupancy(tree, newTypes); // fill in the tree box probabilities

  // how many steps do we want to do?
  for (int s = 0; s < steps; s++) {
    auto erg = step(nodes, vertices, types, tree, numConnectionByPosition, stop);
    vertices = erg.first;
    newTypes = erg.second; // this does not work in a loop, its not the initial types anymore but the returned types
  }

  // we need to save the results (vertices and types)
  // use output
  std::ofstream out(output);
  out << "NodeID,type" << std::endl;
  for (int i = 1; i < newTypes.size(); i++) {
    out << i << "," << newTypes[i] << std::endl;
  }
  out.close();

  // if we have a remove vertices file, write out our updated list
  if (boost::filesystem::exists(remove_file)) {
    std::ofstream out(remove_file);
    out << remove_vertices;
  }

}
