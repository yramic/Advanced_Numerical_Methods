#include "include/is_admissible.hpp"
#include "include/node.hpp"
#include <iostream>

int main()
{
    AdmissibilityH adm;   // defining an admissibility rule
    Node* x = new Node(); // defining two Nodes with their Bounding Boxes´ coordinates
    Node* y = new Node();
    x->setX1_b(7);
    x->setX2_b(38);
    x->setY1_b(62);
    x->setY2_b(90);
    y->setX1_b(11);
    y->setX2_b(17);
    y->setY1_b(30);
    y->setY2_b(46);
    double eta = 2;
    if(adm.is_admissible(x,y,eta)) { // checking the admissibility between the two bounding boxes
        std::cout << "Yes" << std::endl;
    } else {
        std::cout << "No"  << std::endl;
    }
}
