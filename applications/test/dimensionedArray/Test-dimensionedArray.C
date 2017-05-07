#include "dimensionedArray.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    dimensionedArray D("D", dimLength, array::zeros(4));

    Info << D << endl;


    List<token> l;

    l.append(word("D"));
    l.append(token::BEGIN_SQR);
    l.append(0);
    l.append(1);
    l.append(0);
    l.append(0);
    l.append(0);
    l.append(0);
    l.append(0);
    l.append(token::END_SQR);
    l.append(token::BEGIN_LIST);
    l.append(2.5e-2);
    l.append(1.0e-2);
    l.append(0.9e-2);
    l.append(0.8e-2);
    l.append(0.24e-2);
    l.append(token::END_LIST);

    ITstream its("D", l);


    its >> D;

    Info << D << endl;


    return 0;
}
