/*

 File: cinterval.cpp, 2002/12/06

 CoStLy (COmplex interval STandard functions LibrarY), Version 0.3

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
#ifdef FILIB_NAMESPACES
using filib::Double;
#endif

cinterval operator - (const cinterval& z)
{
  return cinterval( -z.re(), -z.im() );
}

int sign(const double& a)
{
  if( a==0.0 ) return 0;
  else         return ( a>0.0 )?1:-1;
}

double f(const double& a, const double& b, const double& c, const double& d, int round)
{
  double value;

  //calculates    (a*c+b*d) / (c*c+d*d)
  switch( round )
    {
    case -1: //downwards rounded
      {
	double z, n, tmp;
	MUL_DOWN( z, a, c ); // z = a*c, downwards rounded
	MUL_DOWN( tmp, b, d ); // tmp = b*d, downwards rounded
	ADD_DOWN_UPD( z, tmp ); // z = z + tmp, downwards rounded

	MUL_UP( n, c, c ); // n = c*c, upwards rounded
	MUL_UP( tmp, d, d ); // tmp = d*d, upwards rounded
	ADD_UP_UPD( n, tmp ); // n = n + tmp, upwards rounded

	DIV_DOWN( value, z, n );

	return value;
      }
    case  0: //nearest rounded
      {
      }
    case  1: //upwards rounded
      {
	double z, n, tmp;
	MUL_UP( z, a, c ); // z = a*c, upwards rounded
	MUL_UP( tmp, b, d ); // tmp = b*d, upwards rounded
	ADD_UP_UPD( z, tmp ); // z = z + tmp, upwards rounded

	MUL_DOWN( n, c, c ); // n = c*c, downwards rounded
	MUL_DOWN( tmp, d, d ); // tmp = d*d, downwards rounded
	ADD_DOWN_UPD( n, tmp ); // n = n + tmp, downwards rounded

	DIV_UP( value, z, n );
	
	return value;
      }
    }
}

bool p[5];

double MinMax(bool minimum, const double& a, const double& b, const double& y0, const Interval& X, int i, int j)
{
  /* 

     berechnet inneres Minimum bzw. Maximum von f = (ac+bd)/(cc+dd) 
     auf dem Interval X = [c.inf,c.sup] ( a,b,d=y0 fest ).         
     Falls minimum = TRUE wird das Minimum berechnet,              
     sonst das Maximum.                                             

  */

  double minmax;

  if( minimum ) minmax =  Double::MAX();
  else          minmax = -Double::MAX();

  if( inf(X)==sup(X) )
    {
      if( p[i] && p[j] ) minmax = f( a, b, inf(X), y0, 1-2*minimum );
      p[i] = false;
      p[j] = false;
    }
  else
    {
      if( a==0.0 )
	{
	  if( b==0.0 || y0==0.0 )
	    {
	      minmax = 0.0;
	      p[i]  = false;
	      p[j]  = false;
	    }
	  else
	    {
	      if( in( 0.0, X ) )
		{
		  if( minimum && ( sign(b)!=sign(y0) ) ) 
		    {
		      DIV_DOWN( minmax, b, y0 );
		      p[i]  = false;
		      p[j]  = false;
		    }
		  else
		    {
		      if( !minimum && ( sign(b)==sign(y0) ) )
			{
			  DIV_UP( minmax, b, y0 );
			  p[i]  = false;
			  p[j]  = false;
			}
		    }
		}
	    }
	}
      else  // a != 0.0 
	{
	  if( y0==0.0 )
	    {
	      if( minimum )
		{
		  if( a>0.0 ) { DIV_DOWN( minmax, a, sup(X) ); }
		  else        { DIV_DOWN( minmax, a, inf(X) ); }
		}
	      else 
		{
		  if( a>0.0 ) { DIV_UP( minmax, a, inf(X) ); }
		  else        { DIV_UP( minmax, a, sup(X) ); }
		}

	      p[i] = false;
	      p[j] = false;
	    }
	  else  // y0 != 0.0,  Ber. von Extremstelle und Minimum,Maximum
	    {
	      // IF NOTBEKANNT THEN
	      // Berechnung von t = sign(a) * ( |b/a| + sqrt( 1+|b/a|^2 ) ) 
	      
	      Interval T,Q( abs( Interval(b)/Interval(a) ) ); 
	      
	      if( a<0.0 ) T = -( Q + sqrt( 1+Q*Q ) );
	      else        T =  ( Q + sqrt( 1+Q*Q ) );
	      
	      // Fallunterscheidung zur Min-,Max- Best:
	      double x, ay0 = fabs(y0);
	      if( ( sign(b)==sign(y0) )== minimum ) 
		{
		  // Berechne Unterschranke für x = |y0| * t, t aus T.
		  MUL_DOWN( x, ay0, inf(T) );
		}
	      else
		{
		  // Berechne Unterschranke für x = |y0| / t, t aus T.
		  DIV_DOWN( x, ay0, sup(T) );
		}
	      if( minimum ) x = -x; // Unterschranke x wird im Falle minimum=TRUE zu einer Oberschranke
	      
	      if( in( x, X ) )
		{
		  // Berechne für minimum=TRUE eine Unterschranke für a / ( 2*x )
		  // andernfalls eine Oberschranke
		  if( minimum )
		    {
		      double tmp;
		      MUL_UP( tmp, 2.0, x );
		      DIV_DOWN( minmax, a, tmp );
		    }
		  else
		    {
		      double tmp;
		      MUL_DOWN( tmp, 2.0, x );
		      DIV_UP( minmax, a, tmp );
		    }
		  p[i] = false;
		  p[j] = false;
		}
	    }
	}
    }
  return minmax;
}

/*--------------------------------------------------------------+
 | complex interval division:                                   |
 |                                                              |
 |    Z1 /= Z2    <=>   Z1 = Z1 / Z2                            |
 |                                                              |
 | see: Lohner, R., Gudenberg, W.v., Complex Interval Division  |
 | with Maximum Accuracy, Proc.Comp.Arith., Ed. Hwang, 1985.    |
 |                                                              |
 +--------------------------------------------------------------*/

cinterval& cinterval::operator /= (const cinterval& z)
{
  if( in( 0.0, z.re() ) && in( 0.0, z.im() ) )
    throw division_by_zero();

  Interval
    rhsRe( z.re() ),
    rhsIm( z.im() );

  double 
    lhsReInf = inf( real_part ),
    lhsReSup = sup( real_part ),
    lhsImInf = inf( imag_part ),
    lhsImSup = sup( imag_part ),
    rhsReInf = inf( rhsRe     ),
    rhsReSup = sup( rhsRe     ),
    rhsImInf = inf( rhsIm     ),
    rhsImSup = sup( rhsIm     );

  // *** Berechnung des Realteils ***

  bool 
    a_repeat = ( rhsReInf < 0.0 ) && ( 0.0 < rhsReSup ),
    b_repeat = ( rhsImInf < 0.0 ) && ( 0.0 < rhsImSup );

  int rep;
  if( a_repeat || b_repeat ) rep = 2;
  else                       rep = 1;

  double a0,b0;
  if( rhsReInf >= 0.0 ) a0 = lhsReSup;
  else                  a0 = lhsReInf;
  if( rhsImInf >= 0.0 ) b0 = lhsImSup;
  else                  b0 = lhsImInf;

  double realteilSUP, realteilINF;
  realteilSUP = -Double::MAX();

  for(int j = 1; j <= rep; j++)
    {
      p[1] = p[2] = p[3] = p[4] = true;

      double 
	mm1 = MinMax( false, a0, b0, rhsImInf, rhsRe, 1, 2 ),
	mm2 = MinMax( false, a0, b0, rhsImSup, rhsRe, 3, 4 ),
	mm3 = MinMax( false, b0, a0, rhsReInf, rhsIm, 1, 3 ),
	mm4 = MinMax( false, b0, a0, rhsReSup, rhsIm, 2, 4 );

      double 
	max1 = ( mm3>mm4 )?mm3:mm4,
	max2 = ( mm1>mm2 )?mm1:mm2,
	max3 = ( max1>max2 )?max1:max2;

      realteilSUP = ( realteilSUP>max3 )?realteilSUP:max3;

      if( p[1] )
	{
	  double y = f( a0, b0, rhsReInf, rhsImInf, +1 );
	  realteilSUP = ( realteilSUP>y )?realteilSUP:y;
	}
      if( p[2] )
	{
	  double y = f( a0, b0, rhsReSup, rhsImInf, +1 );
	  realteilSUP = ( realteilSUP>y )?realteilSUP:y;
	}
      if( p[3] )
	{
	  double y = f( a0, b0, rhsReInf, rhsImSup, +1 );
	  realteilSUP = ( realteilSUP>y )?realteilSUP:y;
	}
      if( p[4] )
	{
	  double y = f( a0, b0, rhsReSup, rhsImSup, +1 );
	  realteilSUP = ( realteilSUP>y )?realteilSUP:y;
	}
      
      if( a_repeat ) a0 = lhsReSup;
      else if( b_repeat ) b0 = lhsImSup;
    }

  if( rhsReInf>=0.0 ) a0 = lhsReInf;
  else                a0 = lhsReSup;
  if( rhsImInf>=0.0 ) b0 = lhsImInf;
  else                b0 = lhsImSup;

  realteilINF = Double::MAX();

  for(int j = 1; j <= rep; j++ )
    {
      p[1] = p[2] = p[3] = p[4] = true;
      
      double 
	mm1 = MinMax( true, a0, b0, rhsImInf, rhsRe, 1, 2 ),
	mm2 = MinMax( true, a0, b0, rhsImSup, rhsRe, 3, 4 ),
	mm3 = MinMax( true, b0, a0, rhsReInf, rhsIm, 1, 3 ),
	mm4 = MinMax( true, b0, a0, rhsReSup, rhsIm, 2, 4 );

      double
	min1 = ( mm3<mm4 )?mm3:mm4,
	min2 = ( mm1<mm2 )?mm1:mm2,
	min3 = ( min1<min2 )?min1:min2;

      realteilINF = ( realteilINF<min3 )?realteilINF:min3;

      if( p[1] )
	{
	  double y = f( a0, b0, rhsReInf, rhsImInf, -1 );
	  realteilINF = ( realteilINF<y )?realteilINF:y;
	}
      if( p[2] )
	{
	  double y = f( a0, b0, rhsReSup, rhsImInf, -1 );
	  realteilINF = ( realteilINF<y )?realteilINF:y;
	}
      if( p[3] )
	{
	  double y = f( a0, b0, rhsReInf, rhsImSup, -1 );
	  realteilINF = ( realteilINF<y )?realteilINF:y;
	}
      if( p[4] )
	{
	  double y = f( a0, b0, rhsReSup, rhsImSup, -1 );
	  realteilINF = ( realteilINF<y )?realteilINF:y;
	}

      if( a_repeat ) a0 = lhsReInf;
      else if( b_repeat ) b0 = lhsImInf;
    }

  // *** Berechnung des Imaginaerteils: g(a, b, c, d) = f(b, -a, c, d) ***

  a_repeat = ( rhsImInf < 0.0 ) && ( 0.0 < rhsImSup );
  b_repeat = ( rhsReInf < 0.0 ) && ( 0.0 < rhsReSup );

  // IF a_repeat OR b_repeat THEN rep:= 2 ELSE rep:= 1;  STIMMT NOCH 

  if( rhsReInf >= 0.0 ) b0 = lhsImSup;
  else                  b0 = lhsImInf;
  if( rhsImInf >= 0.0 ) a0 = lhsReInf;
  else                  a0 = lhsReSup;

  double imagteilSUP, imagteilINF;
  imagteilSUP = -Double::MAX();

  for(int j = 1; j <= rep; j++)
    {
      p[1] = p[2] = p[3] = p[4] = true;

      double
	mm1 = MinMax( false,  b0, -a0, rhsImInf, rhsRe, 1, 2 ),
	mm2 = MinMax( false,  b0, -a0, rhsImSup, rhsRe, 3, 4 ),
	mm3 = MinMax( false, -a0,  b0, rhsReInf, rhsIm, 1, 3 ),
	mm4 = MinMax( false, -a0,  b0, rhsReSup, rhsIm, 2, 4 );

      double
	max1 = ( mm1>mm2 )?mm1:mm2,
	max2 = ( mm3>mm4 )?mm3:mm4,
	max3 = ( max1>max2 )?max1:max2;
      
      imagteilSUP = ( imagteilSUP>max3 )?imagteilSUP:max3;

      if( p[1] )
	{
	  double y = f( b0, -a0, rhsReInf, rhsImInf, +1 );
	  imagteilSUP = ( imagteilSUP>y )?imagteilSUP:y;
	}
      if( p[2] )
	{
	  double y = f( b0, -a0, rhsReSup, rhsImInf, +1 );
	  imagteilSUP = ( imagteilSUP>y )?imagteilSUP:y;
	}
      if( p[3] )
	{
	  double y = f( b0, -a0, rhsReInf, rhsImSup, +1 );
	  imagteilSUP = ( imagteilSUP>y )?imagteilSUP:y;
	}
      if( p[4] )
	{
	  double y = f( b0, -a0, rhsReSup, rhsImSup, +1 );
	  imagteilSUP = ( imagteilSUP>y )?imagteilSUP:y;
	}

      if( b_repeat ) b0 = lhsImSup;
      else if( a_repeat ) a0 = lhsReInf;
    }
 
  if( rhsReInf>=0.0 ) b0 = lhsImInf;
  else                b0 = lhsImSup;
  if( rhsImInf>=0.0 ) a0 = lhsReSup;
  else                a0 = lhsReInf; 
  
  imagteilINF = Double::MAX();
  
  for(int j = 1; j <= rep; j++)
    {
      p[1] = p[2] = p[3] = p[4] = true;
      
      double
	mm1 = MinMax( true,  b0, -a0, rhsImInf, rhsRe, 1, 2 ),
	mm2 = MinMax( true,  b0, -a0, rhsImSup, rhsRe, 3, 4 ),
	mm3 = MinMax( true, -a0,  b0, rhsReInf, rhsIm, 1, 3 ),
	mm4 = MinMax( true, -a0,  b0, rhsReSup, rhsIm, 2, 4 );
      
      double
	min1 = ( mm1<mm2 )?mm1:mm2,
	min2 = ( mm3<mm4 )?mm3:mm4,
	min3 = ( min1<min2 )?min1:min2;
      
      imagteilINF = ( imagteilINF<min3 )?imagteilINF:min3;
      
      if( p[1] )
	{
	  double y = f( b0, -a0, rhsReInf, rhsImInf, -1 );
	  imagteilINF = ( imagteilINF<y )?imagteilINF:y;
	}
      if( p[2] )
	{
	  double y = f( b0, -a0, rhsReSup, rhsImInf, -1 );
	  imagteilINF = ( imagteilINF<y )?imagteilINF:y;
	}
      if( p[3] )
	{
	  double y = f( b0, -a0, rhsReInf, rhsImSup, -1 );
	  imagteilINF = ( imagteilINF<y )?imagteilINF:y;
	}
      if( p[4] )
	{
	  double y = f( b0, -a0, rhsReSup, rhsImSup, -1 );
	  imagteilINF = ( imagteilINF<y )?imagteilINF:y;
	}
      
      if( b_repeat ) b0 = lhsImInf;
      else if( a_repeat ) a0 = lhsReSup;
    }

  real_part = Interval( realteilINF, realteilSUP );
  imag_part = Interval( imagteilINF, imagteilSUP );

  return *this;
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

std::ostream& operator << (std::ostream& os, const cinterval& z)
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
