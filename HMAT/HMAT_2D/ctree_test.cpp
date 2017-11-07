#include "include/point.hpp"
#include "include/ctree.hpp"
#include "include/node.hpp"
#include <iostream>

bool check(std::vector<int> dfs_points, std::vector<Point> tree_vector){
    int i= 0;
    for(std::vector<Point>::iterator iter = tree_vector.begin(); iter < tree_vector.end(); iter++, i++){
        if(iter->getId() != dfs_points[i]) {
            std::cout << iter->getId() << " " << dfs_points[i] << " " << i << std::endl;
            return false;
        }
    }
    return true;
}

std::vector<Point> DFS_traversing(Node* x, std::vector<Point>&  tree_vector){
    if(x->getPPoints().size()<=1){
        return x->getPPoints();
    }
    else{
        if (x->getTl_Child() != NULL) {
            std::vector<Point> t = DFS_traversing(x->getTl_Child(),tree_vector);
            for(std::vector<Point>::iterator it = t.begin(); it < t.end(); it++){
                tree_vector.push_back(*it);
            }
        }
        if (x->getTr_Child() != NULL) {
            std::vector<Point> t = DFS_traversing(x->getTr_Child(),tree_vector);
            for(std::vector<Point>::iterator it = t.begin(); it < t.end(); it++){
                tree_vector.push_back(*it);
            }
        }
        if (x->getBl_Child() != NULL) {
            std::vector<Point> t = DFS_traversing(x->getBl_Child(),tree_vector);
            for(std::vector<Point>::iterator it = t.begin(); it < t.end(); it++){
                tree_vector.push_back(*it);
            }
        }
        if (x->getBr_Child() != NULL) {
            std::vector<Point> t = DFS_traversing(x->getBr_Child(),tree_vector);
            for(std::vector<Point>::iterator it = t.begin(); it < t.end(); it++){
                tree_vector.push_back(*it);
            }
        }
        return x->getPPoints();
    }
}

int main() {
    unsigned n=16;  // number of points
    std::vector<int> x = {24,22,73,63,14,17,39,99,83,41,40,4,83,30,65,23}, y = {24,55,87,91,30,1,9,28,67,85,10,84,57,72,86,56}, v = {66,63,27,46,83,58,46,8,95,57,2,79,34,21,64,95};
    std::vector<Point> PPoints; // initalizing Polygon Points properties
    PPoints.reserve(n);
    for (int i=0; i<n; i++){
        Point p;
        p.setId(i);
        p.setX(x[i]);
        p.setY(y[i]);
        p.setV(v[i]);
        PPoints.push_back(p);
    }
    unsigned d = 2;
    cTree test_tree(PPoints, d);                    // creating the tree
    std::vector<Point> tree_vector, tree_vector_t;  // dummy vectors
    /*Point p;
    p.setId(-1);
    p.setX(0);
    p.setY(0);
    p.setV(0);
    tree_vector.push_back(p);*/
    tree_vector_t = DFS_traversing(test_tree.getRoot(),tree_vector);    // traversing the tree in DFS, so we know with what order should the node of the tree be printed given standard points
    std::vector<int> dfs_points = {11,9,13,9,13,15,1,1,15,1,15,1,15,1,15,1,9,11,13,15,3,14,3,14,3,14,2,2,3,14,8,12,8,12,2,3,8,12,14,4,0,0,4,5,10,6,6,10,6,10,6,10,0,4,5,6,10,7};
    bool correct = check(dfs_points,tree_vector);                       // function to check if the IDs of the points that are in the vector are the same with the IDs we are looking for
    if (correct) std::cout << "Correct" << std::endl;
    else std::cout << "Wrong" << std::endl;
}
