#include "IFstream.H"
#include "dictionary.H"
#include "dimensionedArray.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    dimensionedArray D("D", dimLength, array::zeros(5));
    dimensionedArray2 SigmaS("SigmaS", dimless/dimLength, array2::zeros(5, 5));

    Info << D << endl;
    Info << SigmaS << endl;


    dictionary dict(IFstream("testDict")());

    ITstream *its;

    its = new ITstream("D", dict.lookup("D"));
    *its >> D;
    Info << D << endl;

    its = new ITstream("SigmaS", dict.lookup("SigmaS"));
    *its >> SigmaS;
    Info << SigmaS << endl;

    return 0;
}
