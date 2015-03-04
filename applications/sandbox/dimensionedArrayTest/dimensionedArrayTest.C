#include "fvCFD.H"

#include "array.H"


int main(int argc, char *argv[])
{

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
  // l.append(5);
  l.append(token::BEGIN_LIST);
  l.append(2.5e-2);
  l.append(1.0e-2);
  l.append(0.9e-2);
  l.append(0.8e-2);
  l.append(0.24e-2);
  l.append(token::END_LIST);

  ITstream its("D", l);


  dimensioned<array> D
  (
      "D",
      dimLength,
      array::zeros(5)
  );

  // dimensioned<array> D
  // (
  //     "D",
  //     dimLength,
  //     its
  // );

  // array D(array::zero);

  its >> D;

  Info << D << endl;


  // Info << D << endl;
  // Info << D[0] << endl;
  // Info << D.component(1) << endl;
  // D[0] = 3;
  // Array<int> D1(D);
  // Info << D1[0] << endl;
  // double x;
  // D.component(x, 2);
  // Info << x << endl;
  // D.replace(2, 3.14);
  // Info << D << endl;
  // array A(5, 1);
  // Info << A << endl;
  // A /= 2;
  // Info << -A << endl;
  // Info << name(A) << endl;
  // Info << array::zero(5) << endl;
  // Info << (A + A + A) << endl;
  // Info << (3*A - A) << endl;
  // Info << (A && A) << endl;
  // Info << (A == A) << endl;
  // Info << (A == D) << endl;
  // // Info << pow(A, 0) << endl;
  // // Info << pow(A, 1) << endl;
  // // Info << pow(A, 2) << endl;
  // Info << magSqr(A) << endl;
  // Info << mag(A) << endl;
  // Info << cmptMag(A) << endl;
  // Info << minMod(A, D) << endl;
  // Info << dot(3, A) << endl;
  // Info << (A & A) << endl;
  // // Info << dot(A, D) << endl;


  return 0;
}
