/*

 File: cinterval.cpp, 2002/03/21

 CoStLy (COmplex interval STandard functions LibrarY), Version 0.2

 Copyright (C) Markus Neher, markus.neher@math.uni-karlsruhe.de
               Ingo Eble,    ingoeble@web.de

 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "cinterval.h"
#include "error.h"

using std::ostream;

cinterval operator - (const cinterval& z)
{
  return cinterval( -z.re(), -z.im() );
}

void transform(double& x1, double& x2, double& y1, double& y2, int& q, const cinterval& z)
{
  x1 = inf( z.re() );
  x2 = sup( z.re() );
  y1 = inf( z.im() );
  y2 = sup( z.im() );

  if( (x1>=0 && y1>0) || (x1>0 && y2>0) ) //first quadrant with intersection of positive real axis
    {
      q = 1;
      return;
    }
  else if( (y1>=0 && x2<0) || (y1>0 && x1<0) ) //second quadrant with intersection of positive imaginary axis
    {
      q = 2;

      double tmp = x1;
      x1 = y1;
      y1 = -x2;
      x2 = y2;
      y2 = -tmp;

      return;
    }
  else if( (y2<0 && x2<=0) || (x2<0 && y1<0) ) //third quadrant with intersection of negative real axis
    {
      q = 3;

      double tmp = x1;
      x1 = -x2;
      x2 = -tmp;
      tmp = y1;
      y1 = -y2;
      y2 = -tmp;

      return;
    }
  else if( (y2<=0 && x1>0) || (y2<0 && x2>0) ) //fourth quadrant with intersection of negative imaginary axis
    {
      q = 4;

      double tmp = y1;
      y1 = x1;
      x1 = -y2;
      y2 = x2;
      x2 = -tmp;

      return;
    }
}

cinterval transform_back(const double& x1, const double& x2, const double& y1, const double& y2, int q)
{
  switch( q )
    {
    case 1: //first quadrant
      {
	return cinterval( Interval(x1,x2), Interval(y1,y2) );
      }
    case 2: //second quadrant
      {
	return cinterval( Interval(y1,y2), Interval(-x2,-x1) );
      }
    case 3: //third quadrant
      {
	return cinterval( Interval(-x2,-x1), Interval(-y2,-y1) );
      }
    case 4: //fourth quadrant
      {
	return cinterval( Interval(-y2,-y1), Interval(x1,x2) );
      }
    };
}

void real_part_function(const double& x1, const double& x2, const double& y1, const double& y2, double& inf, double& sup)
{
  /*

    Calculate the minimum and maximum of u(x,y)=x/(x^2+y^2), x in [x1,x2], y in [y1,y2].

  */

  double n,tmp;

  if( y1>0 )
    {
      if( y1>=x2 )
	{
	  MUL_UP( tmp, x1, x1 );
	  MUL_UP( n, y2, y2 );
	  ADD_UP_UPD( n, tmp ); // n = x1^2 + y2^2 upwards rounded
	  DIV_DOWN( inf, x1, n ); // inf = x1 / n downwards rounded

	  MUL_DOWN( tmp, x2, x2 );
	  MUL_DOWN( n, y1, y1 );
	  ADD_DOWN_UPD( n, tmp ); // n = x2^2 + y1^2 downwards rounded
	  DIV_UP( sup, x2, n ); // sup = z / n upwards rounded
	}
      else if( y2>=x1 )
	{
	  if( y1>=x1 )
	    {
	      MUL_DOWN( n, 2.0, y1 ); // n = 2.0 * y1 downwards rounded
	      DIV_UP( sup, 1.0, n ); // sup = z / n upwards rounded
	    }
	  else
	    {
	      MUL_DOWN( tmp, x1, x1 );
	      MUL_DOWN( n, y1, y1 );
	      ADD_DOWN_UPD( n, tmp ); // n = x1^2 + y1^2 downwards rounded
	      DIV_UP( sup, x1, n ); // sup = z / n upwards rounded
	    }

	  MUL_UP( tmp, x1, x1 );
	  MUL_UP( n, y2, y2 );
	  ADD_UP_UPD( n, tmp ); // n = x1^2 + y2^2 upwards rounded
	  DIV_DOWN( inf, x1, n ); // inf = z / n downwards rounded

	  if( y2<x2 )
	    {
	      MUL_UP( tmp, x2, x2 );
	      MUL_UP( n, y2, y2 );
	      ADD_UP_UPD( n, tmp ); // n = x2^2 + y2^2 upwards rounded
	      DIV_DOWN( tmp, x2, n ); // tmp = z / n downwards rounded

	      inf = (inf<tmp)?inf:tmp; // inf = min{ inf, tmp }
	    }
	}
      else
	{
	  MUL_UP( tmp, x2, x2 );
	  MUL_UP( n, y2, y2 );
	  ADD_UP_UPD( n, tmp ); // n = x2^2 + y2^2 upwards rounded
	  DIV_DOWN( inf, x2, n ); // inf = z / n downwards rounded

	  MUL_DOWN( tmp, x1, x1 );
	  MUL_DOWN( n, y1, y1 );
	  ADD_DOWN_UPD( n, tmp ); // n = x1^2 + y1^2 downwards rounded
	  DIV_UP( sup, x1, n ); // sup = z / n upwards rounded
	}
    }
  else
    {
      DIV_UP( sup, 1.0, x1 ); // sup = 1/x1 upwards rounded

      double ym = (-y1>y2)?y1:y2;
      double aym = (ym>0)?ym:-ym;
      double sqrym;
      MUL_UP( sqrym, ym, ym );

      if( aym<=x1 )
	{
	  MUL_UP( tmp, x2, x2 );
	  n = sqrym;
	  ADD_UP_UPD( n, tmp ); // n = x2^2 + ym^2 upwards rounded
	  DIV_DOWN( inf, x2, n ); // inf = z / n downwards rounded
	}
      else if( aym<x2 )
	{
	  MUL_UP( tmp, x2, x2 );
	  n = sqrym;
	  ADD_UP_UPD( n, tmp ); // n = x2^2 + ym^2 upwards rounded
	  DIV_DOWN( inf, x2, n ); // inf = z / n downwards rounded

	  MUL_UP( tmp, x1, x1 );
	  n = sqrym;
	  ADD_UP_UPD( n, tmp ); // n = x1^2 + ym^2 upwards rounded
	  DIV_DOWN( tmp, x1, n ); // tmp = z / n downwards rounded
	  
	  inf = (inf<tmp)?inf:tmp; // inf = min{ inf, tmp }
	}
      else
	{
	  MUL_UP( tmp, x1, x1 );
	  n = sqrym;
	  ADD_UP_UPD( n, tmp ); // n = x1^2 + ym^2 upwards rounded
	  DIV_DOWN( inf, x1, n ); // inf = z / n downwards rounded
	}
    }
}

void imag_part_function(const double& x1, const double& x2, const double& y1, const double& y2, double& inf, double& sup)
{
  /*

    Calculate the minimum and maximum of v(x,y)=-y/(x^2+y^2), x in [x1,x2], y in [y1,y2].

  */

  double n,tmp;

  if( y1>0 )
    {
      if( y2<=x1 )
	{
	  MUL_UP( tmp, x1, x1 );
	  MUL_UP( n, y2, y2 );
	  ADD_UP_UPD( n, tmp ); // n = x1^2 + y2^2 upwards rounded
	  DIV_DOWN( inf, -y2, n ); // inf = z / n downwards rounded

	  MUL_DOWN( tmp, x2, x2 );
	  MUL_DOWN( n, y1, y1 );
	  ADD_DOWN_UPD( n, tmp ); // n = x2^2 + y1^2 downwards rounded
	  DIV_UP( sup, -y1, n ); // sup = z / n upwards rounded
	}
      else if( y1<x2 )
	{
	  if( y1<=x1 )
	    {
	      MUL_UP( n, 2.0, x1 ); // n = 2.0 * x1 upwards rounded
	      DIV_DOWN( inf, -1.0, n ); // sup = z / n upwards rounded
	    }
	  else
	    {
	      MUL_UP( tmp, x1, x1 );
	      MUL_UP( n, y1, y1 );
	      ADD_UP_UPD( n, tmp ); // n = x1^2 + y1^2 upwards rounded
	      DIV_DOWN( inf, -y1, n ); // inf = z / n downwards rounded
	    }

	  MUL_DOWN( tmp, x2, x2 );
	  MUL_DOWN( n, y1, y1 );
	  ADD_DOWN_UPD( n, tmp ); // n = x2^2 + y1^2 upwards rounded
	  DIV_UP( sup, -y1, n ); // inf = z / n downwards rounded

	  if( y2>x2 )
	    {
	      MUL_DOWN( tmp, x2, x2 );
	      MUL_DOWN( n, y2, y2 );
	      ADD_DOWN_UPD( n, tmp ); // n = x2^2 + y2^2 upwards rounded
	      DIV_UP( tmp, -y2, n ); // tmp = z / n downwards rounded

	      sup = (sup>tmp)?sup:tmp; // sup = max{ sup, tmp }
	    }
	}
      else
	{
	  MUL_UP( tmp, x1, x1 );
	  MUL_UP( n, y1, y1 );
	  ADD_UP_UPD( n, tmp ); // n = x1^2 + y1^2 upwards rounded
	  DIV_DOWN( inf, -y1, n ); // inf = z / n downwards rounded

	  MUL_DOWN( tmp, x2, x2 );
	  MUL_DOWN( n, y2, y2 );
	  ADD_DOWN_UPD( n, tmp ); // n = x2^2 + y2^2 downwards rounded
	  DIV_UP( sup, -y2, n ); // sup = z / n upwards rounded
	}
    }
  else
    {
      if( -y1<=x1 )
	{
	  MUL_DOWN( tmp, x1, x1 );
	  MUL_DOWN( n, y1, y1 );
	  ADD_DOWN_UPD( n, tmp ); // n = x1^2 + y1^2 downwards rounded
	  DIV_UP( sup, -y1, n ); // sup = z / n upwards rounded
	}
      else
	{
	  MUL_DOWN( n, 2.0, x1 );
	  DIV_UP( sup, 1.0, n );
	}
      
      if( y2<=x1 )
	{
	  MUL_UP( tmp, x1, x1 );
	  MUL_UP( n, y2, y2 );
	  ADD_UP_UPD( n, tmp ); // n = x1^2 + y2^2 upwards rounded
	  DIV_DOWN( inf, -y2, n ); // inf = z / n downwards rounded
	}
      else
	{
	  MUL_UP( n, 2.0, x1 );
	  DIV_DOWN( inf, -1.0, n );
	}
    }
}

/*--------------------------------------------------------------+
 | complex interval division:                                   |
 |                                                              |
 |    Z1 /= Z2    <=>   Z1 = Z1 / Z2                            |
 |                                                              |
 | is implemented as                                            |
 |                                                              |
 |    Z1 *= 1/Z2                                                |
 |                                                              |
 | An enclosure for 1/Z2 is found by using the monotony         |
 | of the real- and imaginary part of 1/z.                      |
 |                                                              |
 |  1/z = x/(x^2+y^2) - i * y/(x^2+y^2)                         |
 +--------------------------------------------------------------*/

cinterval& cinterval::operator /= (const cinterval& z)
{
  if( in( 0.0, z.re() ) && in( 0.0, z.im() ) )
    throw division_by_zero();

  double x1,x2,y1,y2;
  int q;
  
  transform(x1,x2,y1,y2,q,z);

  double a1,a2,b1,b2;

  real_part_function(x1,x2,y1,y2,a1,a2);
  imag_part_function(x1,x2,y1,y2,b1,b2);

  return (*this) *= transform_back(a1,a2,b1,b2,q);
}

//

Interval abs(const cinterval& z)
{
  return sqrt(sqr(z.re())+sqr(z.im())); 
} 

#ifdef HAS_Complex
Complex  mid (const cinterval& z)
{
  return Complex( mid(z.re()), mid(z.im()) );
}

Complex  diam(const cinterval& z)
{
  return Complex( diam(z.re()), diam(z.im()) );
}
#endif

//Binary operators

cinterval operator + (const cinterval& z1, const cinterval& z2)
{
  cinterval w( z1 );
  return w += z2;
}

cinterval operator - (const cinterval& z1, const cinterval& z2)
{
  cinterval w( z1 );
  return w -= z2;
}

cinterval operator * (const cinterval& z1, const cinterval& z2)
{
  cinterval w( z1 );
  return w *= z2;
}

cinterval operator / (const cinterval& z1, const cinterval& z2)
{
  cinterval w( z1 );
  return w /= z2;
}

// Interval o cinterval, o \in { +,-,* }.

cinterval operator + (const Interval& i, const cinterval& z)
{
  return cinterval( i+z.re(), z.im() );
}

cinterval operator - (const Interval& i, const cinterval& z)
{
  return cinterval( i-z.re(), z.im() );
}

cinterval operator * (const Interval& i, const cinterval& z)
{
  return cinterval( i*z.re(), i*z.im() );
}

cinterval operator / (const Interval& i, const cinterval& z)
{
  cinterval w( i, Interval::ZERO() );
  return w /= z;
}

// cinterval o Interval, o \in { +,-,* }.

cinterval operator + (const cinterval& z, const Interval& i)
{
  return cinterval( z.re()+i, z.im() );
}

cinterval operator - (const cinterval& z, const Interval& i)
{
  return cinterval( z.re()-i, z.im() );
}

cinterval operator * (const cinterval& z, const Interval& i)
{
  return cinterval( z.re()*i, z.im()*i );
}

cinterval operator / (const cinterval& z, const Interval& i)
{
  return cinterval( z.re()/i, z.im()/i );
}

// double o cinterval, o \in { +,-,* }.

cinterval operator + (const double& d, const cinterval& z)
{
  return cinterval( d+z.re(), z.im() );
}

cinterval operator - (const double& d, const cinterval& z)
{
  return cinterval( d-z.re(), z.im() );
}

cinterval operator * (const double& d, const cinterval& z)
{
  return cinterval( d*z.re(), d*z.im() );
}

cinterval operator / (const double& d, const cinterval& z)
{
  cinterval w( Interval(d,d), Interval::ZERO() );
  return w /= z;
}

// cinterval o double, o \in { +,-,* }.

cinterval operator + (const cinterval& z, const double& d)
{
  return cinterval( z.re()+d, z.im() );
}

cinterval operator - (const cinterval& z, const double& d)
{
  return cinterval( z.re()-d, z.im() );
}

cinterval operator * (const cinterval& z, const double& d)
{
  return cinterval( z.re()*d, z.im()*d );
}

cinterval operator / (const cinterval& z, const double& d)
{
  return cinterval( z.re()/d, z.im()/d );
}

// output

ostream& operator << (ostream& os, const cinterval& z)
{
  os << "(" << z.re() << "," << z.im() << ")";
  return os;

}

//

bool operator <= (const double& d, const cinterval& z)
{
  return ( d <= z.re() && 0.0 <= z.im() );
}

//

cinterval operator & (const cinterval& z1, const cinterval& z2)
{
  return cinterval( intersect( z1.re(), z2.re() ), intersect( z1.im(), z2.im() ) );
}

/*

  End of File: cinterval.cpp

*/
