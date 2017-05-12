#include "IFstream.H"
#include "dictionary.H"
#include "dimensionedArray.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    dimensionedArray D("D", dimLength, array::zeros(5));

    Info << D << endl;


    dictionary dict(IFstream("testDict")());
    ITstream its("D", dict.lookup("D"));
    its >> D;

    Info << D << endl;

    return 0;
}
