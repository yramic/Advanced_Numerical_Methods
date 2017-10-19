#include <iostream>
#include "include/is_admissible.hpp"
#include "include/node.hpp"
using namespace std;

int main(int argc, char *argv[])
{
        AdmissibilityH adm;
        Node* x = new Node();
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
        if (adm.is_admissible(x,y,eta)){
            cout << "Yes" << endl;
        }
        else {
            cout << "No" << endl;
        }
}
