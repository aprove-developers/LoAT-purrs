#include <iostream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

int main() {

#if 0
  symbol x("x");
  symbol y("y");
  symbol z("z");
  symbol t("t");

  for (unsigned i = 0; i < 14229; ++i)
    symbol s = symbol();
  
  symbol a;
  symbol b;

  ex e = (-a*pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a, -1)+b*pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)+pow(1+b*a-b-a,-1)-pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)+pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1))*z+(a*pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1)-b*pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)-pow(-2*b*pow(a,2)+b*a+b*pow(a,3)-pow(a,2)-pow(a,4)+2*pow(a,3),-1)*pow(a,3)+pow(b,3)*pow(-pow(b,3)*a-b*a+2*pow(b,2)*a-2*pow(b,3)+pow(b,2)+pow(b,4),-1)+pow(1+b*a-b-a,-1))*y+(b*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(b,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-4*b*a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+2*b*pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+2*pow(b,2)*a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(b,2)*pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)+pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1))*t+(-b*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(b,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(a,3)*pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1)+pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(b,3)*pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)+pow(-2*b*pow(a,2)+b*a+b*pow(a,3)-pow(a,2)-pow(a,4)+2*pow(a,3),-1)*pow(a,3)-pow(b,3)*pow(-pow(b,3)*a-b*a+2*pow(b,2)*a-2*pow(b,3)+pow(b,2)+pow(b,4),-1)+pow(1+b*a-b-a,-1)+4*b*a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-2*b*pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-2*pow(b,2)*a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(b,2)*pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)-pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1))*x;

  cout << "e = " << e << endl << endl;

  e = e.numer_denom();

  cout << "e = " << e << endl;

#else

  symbol x("x");
  symbol y("y");
  symbol z("z");
  symbol t("t");

  symbol a;
  symbol b;

  ex e = (-a*pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a, -1)+b*pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)+pow(1+b*a-b-a,-1)-pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)+pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1))*z+(a*pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1)-b*pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)-pow(-2*b*pow(a,2)+b*a+b*pow(a,3)-pow(a,2)-pow(a,4)+2*pow(a,3),-1)*pow(a,3)+pow(b,3)*pow(-pow(b,3)*a-b*a+2*pow(b,2)*a-2*pow(b,3)+pow(b,2)+pow(b,4),-1)+pow(1+b*a-b-a,-1))*y+(b*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(b,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-4*b*a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+2*b*pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+2*pow(b,2)*a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(b,2)*pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)+pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1))*t+(-b*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(b,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(a,3)*pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1)+pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(b,3)*pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)+pow(-2*b*pow(a,2)+b*a+b*pow(a,3)-pow(a,2)-pow(a,4)+2*pow(a,3),-1)*pow(a,3)-pow(b,3)*pow(-pow(b,3)*a-b*a+2*pow(b,2)*a-2*pow(b,3)+pow(b,2)+pow(b,4),-1)+pow(1+b*a-b-a,-1)+4*b*a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-2*b*pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-2*pow(b,2)*a*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(b,2)*pow(a,2)*pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)-pow(2*b*a-pow(b,2)*a+b-a+pow(b,3)-2*pow(b,2),-1)-pow(-1+2*b*pow(a,2)-4*b*a-pow(b,2)*pow(a,2)+2*pow(b,2)*a+2*b-pow(a,2)+2*a-pow(b,2),-1)+pow(b*pow(a,2)-2*b*a+b+2*pow(a,2)-pow(a,3)-a,-1))*x;

  cout << "e = " << e << endl << endl;

  e = e.numer_denom();

  cout << "e = " << e << endl;
#endif

  return 0;
}
