/*

 File: cimath.cpp, 2002/12/06

 CoStLy (COmplex interval STandard functions LibrarY), Version 0.3

 Copyright (C) Markus Neher, markus.neher@math.uni-karlsruhe.de
               Ingo Eble,    ingoeble@web.de

 Implemetation in C++ with 

        filib++:     Copyright (C) 2002 Markus Neher 
                                        Ingo Eble

	C-XSC:	     Copyright (C) 2000 Markus Neher 
	                                Ingo Eble

 Documentation: See attached file CoStLy03.pdf


 Most of this implementation of complex interval standard functions 
 was taken from a Pascal-XSC modul written by Andreas Westphal
 ((C) 1999 Andreas Westphal, Walter Kraemer, Universitaet Karlsruhe).


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


/*
  Include header files
*/

#include "cimath.h" //Declaration of the complex functions

#ifdef FILIB_VERSION

#include "Interval.h" //fi_lib++ Header: Macro Version
typedef Interval interval;
typedef double flnumber;

#ifdef FILIB_NAMESPACES
using namespace filib;
#endif

#else //C-XSC Version

#include "rmath.hpp" //"real" standard functions 
#include "imath.hpp" //"interval" standard functions 
#include "dot.hpp"   //"dotprecision" standard functions 
typedef real flnumber;

#endif

#include "error.h"

/*
  Some useful values
*/

inline const interval& PI()
{ 
  static const interval pi( 
#ifdef FILIB_VERSION
			   interval::PI()
#else
			   acos(interval(-1.0))
#endif
			   );
  return pi;
}

inline const interval& HALFPI()
{ 
  static const interval hp( 
#ifdef FILIB_VERSION
			   interval::PI() / 2.0
#else //C-XSC Version
                           acos(interval(0.0))
#endif
			   );
  return hp;
}

inline const interval& QUARTERPI()
{
  static const interval qp( 
#ifdef FILIB_VERSION
			   interval::PI() / 4.0
#else //C-XSC Version
                           acos(interval(0.0)) / 2.0
#endif
			   );
  return qp;
}

inline const interval& SQRT_2()
{
  static const interval square_root_2( sqrt(interval(2.0,2.0)) );
  return square_root_2;
}

inline const interval& INV_SQRT_2()
{
  static const interval inv_square_root_2( SQRT_2()/2.0 );
  return inv_square_root_2;
}

inline const interval& ZERO_INTERVAL()
{
  static const interval zero_interval(
#ifdef FILIB_VERSION
				      interval::ZERO()
#else
				      0.0
#endif
				      );
  return zero_interval;
}

inline const interval& ONE_INTERVAL()
{
  static const interval one_interval(
#ifdef FILIB_VERSION
				      interval::ONE()
#else
				      1.0
#endif
				      );
  return one_interval;
}

/* ------------------------------------------------------------------------- */
/*                          Single-valued functions                          */
/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */
/*     Power operator  pow  is not listed here, since it relies on the       */
/*     (multi-valued) logarithm                                              */
/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*     The hyperbolic functions exp, sin, cos, sinh, cosh are separable:     */
/*     Their real and imaginary parts are products of real functions         */
/*                                                                           */
/*     With Re(z)=x, Im(z)=y :                                               */
/*                                                                           */
/*          exp   :   Re(exp(z)) = exp(x) * cos(y)                           */
/*                    Im(exp(z)) = exp(x) * sin(y)                           */
/*                                                                           */
/*          sin   :   Re(sin(z)) = sin(x) * cosh(y)                          */
/*                    Im(sin(x)) = cos(x) * sinh(y)                          */
/*                                                                           */
/*          cos   :   Re(cos(z)) = cos(x) * cosh(y)                          */
/*                    Im(sin(x)) = -sin(x) * sinh(y)                         */
/*                                                                           */
/*          sinh  :   Re(sinh(z)) = sinh(x) * cos(y)                         */
/*                    Im(sinh(z)) = cosh(x) * sin(y)                         */
/*                                                                           */
/*          cosh  :   Re(cosh(z)) = cosh(x) * cos(y)                         */
/*                    Im(cosh(z)) = sinh(x) * sin(y)                         */
/*                                                                           */
/* ------------------------------------------------------------------------- */
 
cinterval exp(const cinterval& z) 
{
#ifdef FILIB_VERSION
  interval A( exp( z.re() ) ), B( z.im() );   
#else //C-XSC Version
  interval A( exp( Re(z) ) ), B( Im(z) );  
#endif
  return cinterval( A*cos( B ) , A*sin( B ) );
}

cinterval cos(const cinterval& z) 
{
#ifdef FILIB_VERSION
  interval A( z.re() ), B( z.im() );   
#else //C-XSC Version
  interval A( Re(z) ), B( Im(z) );  
#endif 
  return cinterval( cos( A )*cosh( B ) , -sin( A )*sinh( B ) );
}

cinterval sin(const cinterval& z) 
{
#ifdef FILIB_VERSION
  interval A( z.re() ), B( z.im() );   
#else //C-XSC Version
  interval A( Re(z) ), B( Im(z) );  
#endif 
  return cinterval( sin( A )*cosh( B ) , cos( A )*sinh( B ) );
} 
    
cinterval cosh(const cinterval& z) 
{
#ifdef FILIB_VERSION
  interval A( z.re() ), B( z.im() );   
#else //C-XSC Version
  interval A( Re(z) ), B( Im(z) );  
#endif 
  return cinterval( cos( B )*cosh( A ) , sin( B )*sinh( A ) );
}

cinterval sinh(const cinterval& z) 
{
#ifdef FILIB_VERSION
  interval A( z.re() ), B( z.im() );   
#else //C-XSC Version
  interval A( Re(z) ), B( Im(z) );  
#endif 
  return cinterval( cos( B )*sinh( A ) , sin( B )*cosh( A ) );
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*     Tangent is NOT separable, naive evaluation may yield                  */
/*     range overestimation.                                                 */
/*                                                                           */
/*     Accurate range evaluation requires discussion of                      */
/*     several cases and uses the symmetry condition                         */
/*                                                                           */
/*     tan(z) = transp( tan( transp(z) ) )                                   */
/*                                                                           */
/* ------------------------------------------------------------------------- */

void ReAdd(const flnumber& hx,const flnumber& hy,interval& re_tan,bool& re_first,bool both) 
{
  if( both ) 
    {
      interval HX( hx ),HY( hy ),cosHX( cos( HX ) );
      if( re_first ) 
	{ //Naive evaluation of real part
	  re_tan = sin( HX ) * cosHX /( sqr( cosHX ) + sqr( sinh( HY )) );
	  re_first = false;
	}
      else 
	{ //Interval hull of former calculations
	  re_tan
#ifdef FILIB_VERSION
	    = hull( re_tan, sin( HX ) * cosHX /( sqr( cosHX ) + sqr( sinh( HY )) ) )
#else //C-XSC Version
	    |= sin( HX ) * cosHX /( sqr( cosHX ) + sqr( sinh( HY )) )
#endif
	    ;
	}
    }
}

void ReAdd(const interval& hx,const flnumber& hy,interval& re_tan,bool& re_first,bool both)
{
  if( both ) 
    {
      interval cosHX( cos( hx ) );
      if( re_first ) 
	{ //Naive evaluation of real part
	  re_tan = sin( hx ) * cosHX /( sqr( cosHX ) + sqr( sinh( interval(hy) )) );
	  re_first = false;
	}
      else 
	{ //Interval hull of former calculations
	  re_tan
#ifdef FILIB_VERSION
	    = hull( re_tan, sin( hx ) * cosHX /( sqr( cosHX ) + sqr( sinh( interval(hy) )) ) )
#else //C-XSC Version
	    |= sin( hx ) * cosHX /( sqr( cosHX ) + sqr( sinh( interval(hy) )) )
#endif
	    ;
	}
    }
}

void ImAdd(const flnumber& hx,const flnumber& hy,interval& im_tan,bool& im_first) 
{
  interval HX( hx ),HY( hy ),sinhHY( sinh( HY ) );
  if( im_first ) 
    { //Naive evaluation of imaginary part
      im_tan = sinhHY * cosh( HY ) /( sqr( cos( HX )) + sqr( sinhHY ) );
      im_first = false;
    }
  else 
    { //Interval hull of former calculations
      im_tan
#ifdef FILIB_VERSION
	= hull( im_tan, sinhHY * cosh( HY ) /( sqr( cos( HX )) + sqr( sinhHY ) ) )
#else //C-XSC Version
	|= sinhHY * cosh( HY ) /( sqr( cos( HX )) + sqr( sinhHY ) )
#endif
	;
    }
}

void ImAdd(const interval& hy,const flnumber& hx,interval& im_tan,bool& im_first) 
{
  interval sinhHY( sinh( hy ) );
  if( im_first ) 
    { //Naive evaluation of imaginary part
      im_tan = sinhHY * cosh( hy ) /( sqr( cos( interval( hx ) )) + sqr( sinhHY ) );
      im_first = false;
    }
  else 
    { //Interval hull of former calculations
      im_tan
#ifdef FILIB_VERSION
	= hull( im_tan, sinhHY * cosh( hy ) /( sqr( cos( interval( hx ) )) + sqr( sinhHY ) ) )
#else //C-XSC Version
	|= sinhHY * cosh( hy ) /( sqr( cos( interval( hx ) )) + sqr( sinhHY ) )
#endif
	;
    }
}

#ifndef FILIB_VERSION
bool disjoint(const interval& x,const interval& y) 
{
  flnumber ix( Inf(x) ),iy( Inf(y) ),sx( Sup(x) ), sy( Sup(y) );
  flnumber inf( ( ix > iy )? ix : iy );
  flnumber sup( ( sx < sy )? sx : sy );

  return ( inf > sup );
}
#endif

void htan(const interval& x,const interval& y,interval& re_tan,interval& im_tan,bool both) 
{
  /* If BOTH=TRUE, then enclosures of  Re(tan(x+i*y)) and                */
  /* Im(tan(x+i*y)) are calculated; otherwise only Im_tan is calculated. */
  /* y > 0 is assumed here (upper half plane).                           */

  bool Re_first( true ),
       Im_first( true );
  	
  /* Real part of arguments from left to right */
  
  interval hint(0.0),argx(0.0),argy(0.0);
#ifdef FILIB_VERSION
  if( inf(x) < (-inf(QUARTERPI())) ) 
    { //Real part of argument intersects [-pi/2,-pi/4] 

      if( (-inf(QUARTERPI())) < sup(x) ) hint = interval( inf(x), -inf(QUARTERPI()) );
      else                               hint = x; //hint = Intersection of x and [-pi/2,-pi/4] 
    
      ReAdd(inf(x),sup(y),re_tan,Re_first,both);     //possible maximum 
      ReAdd(sup(hint),sup(y),re_tan,Re_first,both);  //possible maximum 
      ImAdd(sup(hint),sup(y),im_tan,Im_first);       //possible minimum 
      ImAdd(sup(hint),inf(y),im_tan,Im_first);       //possible minimum 

      if( both ) 
	{ //minimum of real part in hint 
	  if( inf(y) > 0 ) 
	    {
	      argx = atan(coth(interval(inf(y)))); 
	      if( disjoint(argx , hint) ) 
		{ 
		  if( inf(x) < inf(argx) ) ReAdd(inf(x),inf(y),re_tan,Re_first,both);    //minimum 
		  else                     ReAdd(sup(hint),inf(y),re_tan,Re_first,both); //minimum 
		}
	      else ReAdd(argx,inf(y),re_tan,Re_first,both);
	    }
	  else ReAdd(inf(x),inf(y),re_tan,Re_first,both); //minimum 
	}
      //maximum of imaginary part in hint 

      if( inf(x) < -sup(QUARTERPI()) ) 
	{
	  argy = acoth(-tan(interval(inf(x)))); //maximum for fixed x and arbitrary y
	  if( disjoint(argy , y) ) 
	    { //maximum in one of the left corners 
	      if( inf(y)< sup(argy) ) ImAdd(inf(x),inf(y),im_tan,Im_first);
	      else                    ImAdd(inf(x),sup(y),im_tan,Im_first);
	    }
	  else ImAdd(argy,inf(x),im_tan,Im_first);
	}
      else ImAdd(y,inf(x),im_tan,Im_first);
    }

  hint = interval(-inf(QUARTERPI()) , 0);
  if( !disjoint( hint , x ) ) 
    { //Real part of arguments intersects [-pi/4,0] 
      hint = intersect(hint, x);
      ReAdd(sup(hint),sup(y),re_tan,Re_first,both); //maximum 
      ReAdd(inf(hint),inf(y),re_tan,Re_first,both); //minimum 
      ImAdd(inf(hint),sup(y),im_tan,Im_first); //maximum 
      ImAdd(sup(hint),inf(y),im_tan,Im_first); //minimum 
    }
  
  hint = interval( 0 ,inf(QUARTERPI()) );
  if( !disjoint( hint , x ) ) 
    { //Real part of argument intersects [0,pi/4] 
      hint = intersect(hint, x);
      ReAdd(sup(hint),inf(y),re_tan,Re_first,both); //maximum 
      ReAdd(inf(hint),sup(y),re_tan,Re_first,both); //minimum 
      ImAdd(sup(hint),sup(y),im_tan,Im_first); //maximum 
      ImAdd(inf(hint),inf(y),im_tan,Im_first); //minimum 
    }

  if( sup(x) > inf(QUARTERPI()) ) 
    { //Real part of argument intersects [pi/4,pi/2] 
      if( inf(x) < inf(QUARTERPI()) ) hint = interval(inf(QUARTERPI()),sup(x));
      else                           hint = x;
      if( both ) 
	{ //Maximum of real part in hint 
	  if( inf(y) == 0 ) ReAdd(sup(x),inf(y),re_tan,Re_first,both);
	  else 
	    {
	      argx = atan(coth(interval(inf(y)))); //Maximum for fixed argy, arbitrary argx 
	      if( disjoint(argx , x ) ) 
		{ //Maximum in one of the bottom corners
		  if(sup(argx)< inf(hint)) ReAdd(inf(hint),inf(y),re_tan,Re_first,both);
		  else                     ReAdd(sup(x),inf(y),re_tan,Re_first,both);
		}
	    }
	}
  
      //Maximum of imaginary part in hint 
      if( sup(x) < inf(QUARTERPI()) )  
	{
	  argy = acoth(tan(interval(sup(x)))); //Maximum for fixed x and arbitrary y 
	  if( disjoint(argy , y) ) 
	    { //Maximum in one of the right corners 
	      if(sup(y)< inf(argy)) ImAdd(sup(x),sup(y),im_tan,Im_first);
	      else                  ImAdd(sup(x),inf(y),im_tan,Im_first);
	    }
	}
      else ImAdd(y,sup(x),im_tan,Im_first);
      //Possible minima : left corners  
      ImAdd(inf(hint),sup(y),im_tan,Im_first);
      ImAdd(inf(hint),inf(y),im_tan,Im_first);
      //Possible minima: top corners 
      ReAdd(inf(hint),sup(y),re_tan,Re_first,both);
      ReAdd(sup(x),sup(y),re_tan,Re_first,both);  
    }  
#else //C-XSC Version
  if( Inf(x) < (-Inf(QUARTERPI())) ) 
    { //Real part of argument intersects [-pi/2,-pi/4] 
      if( (-Inf(QUARTERPI())) < Sup(x) ) hint = interval( Inf(x), -Inf(QUARTERPI()) );
      else                               hint = x; //hint = Intersection of x and [-pi/2,-pi/4] 
    
      ReAdd(Inf(x),Sup(y),re_tan,Re_first,both);     //possible maximum 
      ReAdd(Sup(hint),Sup(y),re_tan,Re_first,both);  //possible maximum 
      ImAdd(Sup(hint),Sup(y),im_tan,Im_first);       //possible minimum 
      ImAdd(Sup(hint),Inf(y),im_tan,Im_first);       //possible minimum 

      if( both ) 
	{ //minimum of real part in hint 
	  if( Inf(y) > 0 ) 
	    {
	      argx = atan(coth(interval(Inf(y)))); 
	      if( disjoint(argx , hint) ) 
		{ 
		  if( Inf(x) < Inf(argx) ) ReAdd(Inf(x),Inf(y),re_tan,Re_first,both);    //minimum 
		  else                     ReAdd(Sup(hint),Inf(y),re_tan,Re_first,both); //minimum 
		}
	      else ReAdd(argx,Inf(y),re_tan,Re_first,both);
	    }
	  else ReAdd(Inf(x),Inf(y),re_tan,Re_first,both); //minimum 
	}
      //maximum of imaginary part in hint 

      if( Inf(x) < -Sup(QUARTERPI()) ) 
	{
	  argy = acoth(-tan(interval(Inf(x)))); //maximum for fixed x and arbitrary y
	  if( disjoint(argy , y) ) 
	    { //maximum in one of the left corners 
	      if( Inf(y)< Sup(argy) ) ImAdd(Inf(x),Inf(y),im_tan,Im_first);
	      else                    ImAdd(Inf(x),Sup(y),im_tan,Im_first);
	    }
	  else ImAdd(argy,Inf(x),im_tan,Im_first);
	}
      else ImAdd(y,Inf(x),im_tan,Im_first);
    }

  hint = interval(-Inf(QUARTERPI()) , 0);
  if( !disjoint( hint , x ) ) 
    { //Real part of arguments intersects [-pi/4,0] 
      hint &= x;
      ReAdd(Sup(hint),Sup(y),re_tan,Re_first,both); //maximum 
      ReAdd(Inf(hint),Inf(y),re_tan,Re_first,both); //minimum 
      ImAdd(Inf(hint),Sup(y),im_tan,Im_first); //maximum 
      ImAdd(Sup(hint),Inf(y),im_tan,Im_first); //minimum 
    }
  
  hint = interval( 0 ,Inf(QUARTERPI()) );
  if( !disjoint( hint , x ) ) 
    { //Real part of argument intersects [0,pi/4] 
      hint &= x;
      ReAdd(Sup(hint),Inf(y),re_tan,Re_first,both); //maximum 
      ReAdd(Inf(hint),Sup(y),re_tan,Re_first,both); //minimum 
      ImAdd(Sup(hint),Sup(y),im_tan,Im_first); //maximum 
      ImAdd(Inf(hint),Inf(y),im_tan,Im_first); //minimum 
    }
  
  if( Sup(x) > (Inf(QUARTERPI())) ) 
    { //Real part of argument intersects [pi/4,pi/2] 
      if( Inf(x) < (Inf(HALFPI())) ) hint = interval(Inf(QUARTERPI()),Sup(x));
      else                           hint = x;
      if( both ) 
	{ //Maximum of real part in hint 
	  if( Inf(y) == 0 ) ReAdd(Sup(x),Inf(y),re_tan,Re_first,both);
	  else 
	    {
	      argx = atan(coth(interval(Inf(y)))); //Maximum for fixed argy, arbitrary argx 
	      if( disjoint(argx , x ) ) 
		{ //Maximum in one of the bottom corners
		  if(Sup(argx)< Inf(hint)) ReAdd(Inf(hint),Inf(y),re_tan,Re_first,both);
		  else                     ReAdd(Sup(x),Inf(y),re_tan,Re_first,both);
		}
	    }
	}
  
      //Maximum of imaginary part in hint 
      if( Sup(x) < (Inf(QUARTERPI())) )  
	{
	  argy = acoth(tan(interval(Sup(x)))); //Maximum for fixed x and arbitrary y 
	  if( disjoint(argy , y) ) 
	    { //Maximum in one of the right corners 
	      if(Sup(y)< Inf(argy)) ImAdd(Sup(x),Sup(y),im_tan,Im_first);
	      else                  ImAdd(Sup(x),Inf(y),im_tan,Im_first);
	    }
	}
      else ImAdd(y,Sup(x),im_tan,Im_first);
      //Possible minima : left corners  
      ImAdd(Inf(hint),Sup(y),im_tan,Im_first);
      ImAdd(Inf(hint),Inf(y),im_tan,Im_first);
      //Possible minima: top corners 
      ReAdd(Inf(hint),Sup(y),re_tan,Re_first,both);
      ReAdd(Sup(x),Sup(y),re_tan,Re_first,both);  
    }
#endif
}

void h_tan(const interval& x,const interval& y,interval& re_tan,interval& im_tan,int& Error) 
{
  interval hx(0.0),hx2(0.0);
  interval Imtan(0.0),Retan(0.0);
  bool BOTH( false );
  bool bisection( false ); // Pi/2 mod Pi  in x ? 
  
  Error = 0;
#ifdef FILIB_VERSION
  if ( x == ZERO_INTERVAL() ) 
    {
      re_tan = x;               //Re(tan) = 0
      im_tan = tanh(y);
    }
  else
    {
      if( (!disjoint(ZERO_INTERVAL(),y)) && (!disjoint(ZERO_INTERVAL(),cos(x))) ) Error = 1;
      else
	{
	  if ( y == ZERO_INTERVAL() )  
	    {
	      re_tan = tan(x);
	      im_tan = y;         //Im(tan) = 0
	    }
	  else
	    { 
	      //Since  z equivalent (z mod pi), use  
	      //equivalent argument in [ -pi/2 , pi/2 ] .  
	      flnumber r_int;

	      if( inf(x) > 0 ) modf( inf(x/PI()) + 0.5, &r_int );
	      else             modf( inf(x/PI()) - 0.5, &r_int );
	      
	      hx = x - r_int * PI();

	      if( inf(PI()) < (2*sup(hx)) )  
		{
		  bisection = true;
		  hx2 = interval( -sup(HALFPI()), sup(hx)-inf(HALFPI()) );
		  hx  = interval(  inf(hx)      , sup(HALFPI())         );
		} 

	      BOTH = true;
	      //Use  tan(transp(z))=transp(tan(z))  for cinterval z. 
	      //Only arguments in upper half plane are required! 
	      if( sup(y) <= 0 )  
		{
		  if( bisection )
		    {
		      htan(hx , -y, re_tan, im_tan, BOTH );
		      htan(hx2, -y, Retan , Imtan , BOTH );
		      im_tan = -hull( im_tan, Imtan );
		      re_tan =  hull( re_tan, Retan );
		    }
		  else
		    {
		      htan(hx,-y,re_tan,im_tan,BOTH);
		      im_tan = -im_tan;
		    }
		}
	      else
		{
		  if( inf(y) >= 0 )  
		    {
		      if( bisection )
			{
			  htan(hx,y,re_tan,im_tan,BOTH);
			  htan(hx2,y,Retan,Imtan,BOTH);
			  im_tan = hull( im_tan, Imtan );
			  re_tan = hull( re_tan, Retan );
			}
		      else htan(hx,y,re_tan,im_tan,BOTH);
		    }
		  else
		    {
		      if( -inf(y) < sup(y) )
			{
			  if( bisection )
			    {
			      htan(hx,interval(0,sup(y)),re_tan,im_tan,BOTH);
			      htan(hx2,interval(0,sup(y)),Retan,Imtan,BOTH);
			      im_tan = hull( im_tan, Imtan );
			      re_tan = hull( re_tan, Retan );
			      htan(hx,interval(0,-inf(y)),re_tan,Imtan,!BOTH);
			      im_tan = hull( im_tan, -Imtan ); 
			      htan(hx2,interval(0,-inf(y)),re_tan,Imtan,!BOTH);
			      im_tan = hull( im_tan, -Imtan );
			    }
			  else
			    {
			      htan(hx,interval(0,sup(y)),re_tan,im_tan,BOTH);
			      htan(hx,interval(0,-inf(y)),re_tan,Imtan,!BOTH);
			      im_tan = hull( im_tan, -Imtan );
			    }
			}
		      else
			{
			  if( bisection )
			    {
			      htan(hx,interval(0,-inf(y)),re_tan,Imtan,BOTH);
			      htan(hx,interval(0,sup(y)),re_tan,im_tan,!BOTH);
			      im_tan = hull( im_tan, -Imtan );
			      htan(hx2,interval(0,-inf(y)),Retan,Imtan,BOTH);
			      im_tan = hull( im_tan, Imtan );
			      re_tan = hull( re_tan, Retan );
			      htan(hx2,interval(0,sup(y)),re_tan,im_tan,!BOTH);
			      im_tan = hull( im_tan, -Imtan );
			    }
			  else
			    {
			      htan(hx,interval(0,-inf(y)),re_tan,Imtan,BOTH);
			      htan(hx,interval(0,sup(y)),re_tan,im_tan,!BOTH);
			      im_tan = hull( im_tan, -Imtan );
			    }
			}
		    }
		}
	    }
	}
    }
#else //C-XSC Version
  if ( x == ZERO_INTERVAL() )  
    {
      re_tan = x;               //Re(tan) = 0
      im_tan = tanh(y);
    }
  else
    {
      if( (!disjoint(ZERO_INTERVAL(),y)) && (!disjoint(ZERO_INTERVAL(),cos(x))) ) Error = 1;
      else
	{
	  if ( y == ZERO_INTERVAL() )  
	    {
	      re_tan = tan(x);
	      im_tan = y;         //Im(tan) = 0
	    }
	  else
	    { 
	      //Since  z equivalent (z mod pi), use  
	      //equivalent argument in [ -pi/2 , pi/2 ] .  

	      if( sign(Inf(x)) == 1 ) hx = x - int( _double( Inf(x/PI()) + 0.5 ) ) * PI();
	      else    		      hx = x - int( _double( Inf(x/PI()) - 0.5 ) ) * PI();

	      if( Inf(PI()) < (2*Sup(hx)) )  
		{
		  bisection = true;
		  hx2 = interval( -Sup(HALFPI()),Sup(hx)-Inf(HALFPI()) );
		  hx  = interval(  Inf(hx)      ,Sup(HALFPI()) );
		} 

	      BOTH = true;
	      //Use  tan(transp(z))=transp(tan(z))  for cinterval z. 
	      //Only arguments in upper half plane are required! 
	      if( Sup(y) <= 0 )  
		{
		  if( bisection )
		    {
		      htan(hx,-y,re_tan,im_tan,BOTH);
		      htan(hx2,-y,Retan,Imtan,BOTH);
		      im_tan = -( im_tan | Imtan);
		      re_tan |= Retan;
		    }
		  else
		    {
		      htan(hx,-y,re_tan,im_tan,BOTH);
		      im_tan = -im_tan;
		    }
		}
	      else
		{
		  if( Inf(y) >= 0 )  
		    {
		      if( bisection )
			{
			  htan(hx,y,re_tan,im_tan,BOTH);
			  htan(hx2,y,Retan,Imtan,BOTH);
			  im_tan |= Imtan;
			  re_tan |= Retan;
			}
		      else htan(hx,y,re_tan,im_tan,BOTH);
		    }
		  else
		    {
		      if( -Inf(y) < Sup(y) )
			{
			  if( bisection )
			    {
			      htan(hx,interval(0,Sup(y)),re_tan,im_tan,BOTH);
			      htan(hx2,interval(0,Sup(y)),Retan,im_tan,BOTH);
			      im_tan |= im_tan;
			      re_tan |= Retan;
			      htan(hx,interval(0,-Inf(y)),re_tan,Imtan,!BOTH);
			      im_tan |=  -Imtan; 
			      htan(hx2,interval(0,-Inf(y)),re_tan,Imtan,!BOTH);
			      im_tan |= -Imtan;
			    }
			  else
			    {
			      htan(hx,interval(0,Sup(y)),re_tan,im_tan,BOTH);
			      htan(hx,interval(0,-Inf(y)),re_tan,Imtan,!BOTH);
			      im_tan |= -Imtan;
			    }
			}
		      else
			{
			  if( bisection )
			    {
			      htan(hx,interval(0,-Inf(y)),re_tan,Imtan,BOTH);
			      htan(hx,interval(0,Sup(y)),re_tan,im_tan,!BOTH);
			      im_tan |= -Imtan;
			      htan(hx2,interval(0,-Inf(y)),Retan,Imtan,BOTH);
			      im_tan |=  Imtan;
			      re_tan |=  Retan;
			      htan(hx2,interval(0,Sup(y)),re_tan,im_tan,!BOTH);
			      im_tan |= -Imtan;
			    }
			  else
			    {
			      htan(hx,interval(0,-Inf(y)),re_tan,Imtan,BOTH);
			      htan(hx,interval(0,Sup(y)),re_tan,im_tan,!BOTH);
			      im_tan |= -Imtan;
			    }
			}
		    }
		}
	    }
	}
    }
#endif
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*      Computation of cot, tanh, coth based on tan.                         */
/*                                                                           */
/*      cot(z)  = tan( pi/2 - z )                                            */
/*      tanh(z) = transp( i * tan( transp( i * z ) ) )                       */
/*      coth(z) = i * cot( i * z )                                           */
/*              = i * tan( pi/2 - i * z )                                    */
/*                                                                           */
/* ------------------------------------------------------------------------- */

void h_cot(const interval& x,const interval& y,interval& re_cot,interval& im_cot,int& Error) 
{
  interval hx(0.0),hx2(0.0);
  interval Imcot(0.0),Recot(0.0);
  bool BOTH( false );
  bool bisection( false ); //0 mod Pi  in x ? 
  //If yes: bisect input interval
  
  Error =  0;
#ifdef FILIB_VERSION
  if( (!disjoint(ZERO_INTERVAL(),y))&&(!disjoint(ZERO_INTERVAL(),sin(x))) ) 
    {
      Error = 1;
    }
  else
    {
      if( y == ZERO_INTERVAL() ) 
	{
	  re_cot = cot(x);
	  im_cot = y;               //Im(cot) = 0 
	}
      else
	{
	  if( x == ZERO_INTERVAL() ) 
	    {
	      re_cot =  x;              //Re(cot) = 0 
	      im_cot =  coth(y);
	    }
	  else 
	    {
	      /* z equivalent (z mod pi), use equivalent          */
	      /* argument in [ 0 , Pi ] . Since                   */
	      /*          cot(z)=tan(Pi/2 - z)                    */
	      /* compute                                          */
	      /*    hz = Pi/2 - z  mod Pi  in  [ -Pi/2 , Pi/2 ]   */
	      /* and use procedure htan.                          */
	      if( (sup(x)-inf(x)) > inf(PI()) ) 
		{
		  hx = interval( -sup(HALFPI()), sup(HALFPI()) );
		}
	      else
		{
		  /* z equivalent (z mod pi), use               */
		  /* equivalent argument in [ -pi/2 , pi/2 ] .  */
		  
		  flnumber r_int;
		  
		  if( inf(x) < 0 ) 
		    {
		      modf( -inf(x/PI()), &r_int );
		      hx = -x + ( -1 - 2 * r_int ) * HALFPI(); 
		    }
		  else
		    {
		      modf(  inf(x/PI()), &r_int );
		      hx = -x + (  1 + 2 * r_int ) * HALFPI(); 
		    }
		  
		  if( sup(hx) > inf(HALFPI()) )  
		    {
		      bisection = true;
		      hx2 =  interval( -sup(HALFPI()), sup(hx)-inf(HALFPI()) );
		      hx  =  interval(  inf(hx)    , sup(HALFPI())         );
		    }
		  else bisection = false; 
		}
	   
	      BOTH =  true;
	      /* Use tan(transp(z))=transp(tan(z))  for cinterval z. */
	      /* Requires only arguments in upper half plane!        */
	      if( sup(y) <= 0) 
		{
		  htan(hx,-y,re_cot,im_cot,BOTH);
		  if( bisection ) 
		    {
		      htan(hx2,-y,Recot,Imcot,BOTH);
		      re_cot =  hull( re_cot, Recot );
		      im_cot = -hull( im_cot, Imcot );
		    }
		  else im_cot = -im_cot;
		}
	      else
		{
		  if ( inf(y) >= 0 ) 
		    {
		      if( bisection )
			{
			  htan(hx,y,re_cot,im_cot,BOTH);
			  htan(hx2,y,Recot,Imcot,BOTH);
			  re_cot = hull( re_cot, Recot );
			  im_cot = hull( im_cot, Imcot ); 
			}
		      else htan(hx,y,re_cot,im_cot,BOTH);
		    }
		  else
		    {
		      if( -inf(y) < sup(y) )
			{
			  htan(hx,interval(0,sup(y)),re_cot,im_cot,BOTH);
			  htan(hx,interval(0,-inf(y)),Recot,Imcot,!BOTH);
			  im_cot = hull( im_cot, -Imcot );
			  if( bisection )
			    {
			      htan(hx2,interval(0,sup(y)),Recot,Imcot,BOTH);
			      re_cot = hull( re_cot, Recot );
			      im_cot = hull( im_cot, Imcot );
			      htan(hx,interval(0,-inf(y)),Recot,Imcot,!BOTH);
			      im_cot = hull( im_cot, -Imcot );
			    }
			}
		      else
			{
			  htan(hx,interval(0,-inf(y)),re_cot,Imcot,BOTH);
			  htan(hx,interval(0,sup(y)),Recot,im_cot,!BOTH);
			  im_cot = hull( im_cot, -Imcot );
			  if( bisection )
			    {
			      htan(hx2,interval(0,-inf(y)),Recot,Imcot,BOTH);
			      re_cot = hull( re_cot, Recot );
			      im_cot = hull( im_cot, -Imcot ); 
			      htan(hx2,interval(0,sup(y)),Recot,Imcot,!BOTH);
			      im_cot = hull( im_cot, Imcot );
			    }
			}
		    }
		}
	    }
	}
    }
#else //C-XSC Version
  if( (!disjoint(ZERO_INTERVAL(),y))&&(!disjoint(ZERO_INTERVAL(),sin(x))) ) 
    {
      Error = 1;
    }
  else
    {
      if( y == ZERO_INTERVAL() ) 
	{
	  re_cot = cot(x);
	  im_cot = y;               //Im(cot) = 0 
	}
      else
	{
	  if( x == ZERO_INTERVAL() ) 
	    {
	      re_cot =  x;              //Re(cot) = 0 
	      im_cot =  coth(y);
	    }
	  else 
	    {
	      /* z equivalent (z mod pi), use equivalent          */
	      /* argument in [ 0 , Pi ] . Since                   */
	      /*          cot(z)=tan(Pi/2 - z)                    */
	      /* compute                                          */
	      /*    hz = Pi/2 - z  mod Pi  in  [ -Pi/2 , Pi/2 ]   */
	      /* and use procedure htan.                          */
	      if( (Sup(x)-Inf(x)) > Inf(PI()) ) 
		{
		  hx = interval(-Sup(HALFPI()),Sup(HALFPI()));
		}
	      else
		{
		  /* z equivalent (z mod pi), use               */
		  /* equivalent argument in [ -pi/2 , pi/2 ] .  */
		  if( Inf(x) < 0 ) hx = -x + ( -1 - 2 * int( _double( -Inf(x/PI()) ) ) ) * HALFPI(); 
		  else	           hx = -x + (  1 + 2 * int( _double(  Inf(x/PI()) ) ) ) * HALFPI();
		  
		  if( Sup(hx) > Inf(HALFPI()) )  
		    {
		      bisection = true;
		      hx2 = interval(-Sup(HALFPI()),Sup(hx)-Inf(HALFPI()) );
		      hx  = interval( Inf(hx)      ,Sup(HALFPI())         );
		    }
		  else bisection = false; 
		}
	   
	      BOTH =  true;
	      /* Use tan(transp(z))=transp(tan(z))  for cinterval z. */
	      /* Requires only arguments in upper half plane!        */
	      if( Sup(y) <= 0) 
		{
		  htan(hx,-y,re_cot,im_cot,BOTH);
		  if( bisection ) 
		    {
		      htan(hx2,-y,Recot,Imcot,BOTH);
		      re_cot |= Recot;
		      im_cot = -(im_cot | Imcot);
		    }
		  else im_cot =  -im_cot;
		}
	      else
		{
		  if ( Inf(y) >= 0 ) 
		    {
		      if( bisection )
			{
			  htan(hx,y,re_cot,im_cot,BOTH);
			  htan(hx2,y,Recot,Imcot,BOTH);
			  re_cot |= Recot;
			  im_cot |= Imcot; 
			}
		      else htan(hx,y,re_cot,im_cot,BOTH);
		    }
		  else
		    {
		      if( -Inf(y) < Sup(y) )
			{
			  htan(hx,interval(0,Sup(y)),re_cot,im_cot,BOTH);
			  htan(hx,interval(0,-Inf(y)),Recot,Imcot,!BOTH);
			  im_cot |= (-Imcot);
			  if( bisection )
			    {
			      htan(hx2,interval(0,Sup(y)),Recot,Imcot,BOTH);
			      re_cot |= Recot;
			      im_cot |= Imcot;
			      htan(hx,interval(0,-Inf(y)),Recot,Imcot,!BOTH);
			      im_cot |= (-Imcot);
			    }
			}
		      else
			{
			  htan(hx,interval(0,-Inf(y)),re_cot,Imcot,BOTH);
			  htan(hx,interval(0,Sup(y)),Recot,im_cot,!BOTH);
			  im_cot |= (-Imcot);
			  if( bisection )
			    {
			      htan(hx2,interval(0,-Inf(y)),Recot,Imcot,BOTH);
			      re_cot |= Recot;
			      im_cot |= (-Imcot); 
			      htan(hx2,interval(0,Sup(y)),Recot,Imcot,!BOTH);
			      im_cot |= Imcot;
			    }
			}
		    }
		}
	    }
	}
    }
#endif
}

cinterval tan(const cinterval& z) 
{
  interval re_tan(0.0),im_tan(0.0);
  int Error(0);
  
  h_tan(Re(z),Im(z),re_tan,im_tan,Error);

  if( Error == 0 )
    return cinterval(re_tan,im_tan);
  else
    {
      if( Error == 1 )//Pole of tangent contained in argument.
	throw function_not_defined();
      else//Real part of argument too far from [-pi/2,pi/2].
	throw function_not_defined();
    }
}

cinterval cot(const cinterval& z) 
{
  interval re_cot(0.0),im_cot(0.0);
  int Error(0);
  
  h_cot(Re(z),Im(z),re_cot,im_cot,Error);

  if( Error == 0 ) 
    return cinterval(re_cot,im_cot);
  else
    {
      if( Error == 1 )//Pole of cotangent contained in argument.
	throw function_not_defined();
      else//Real part of argument too far from [0,pi].
	throw function_not_defined();
    }
}

cinterval tanh(const cinterval& z) 
{
  interval re_tanh(0.0),im_tanh(0.0);
  int Error(0);
  
  h_tan(Im(z),Re(z),im_tanh,re_tanh,Error);

  if( Error == 0 )
      return cinterval(re_tanh,im_tanh);
  else
    {
      if( Error == 1 )//Pole of tangent hyperbolicus contained in argument.
	throw function_not_defined();
      else//Imaginary part of argument too far from [-pi/2,pi/2].
	throw function_not_defined();
    }
}

cinterval coth(const cinterval& z) 
{
  interval m_re_coth(0.0),im_coth(0.0);
  int Error(0);
    
  h_cot(-Im(z),Re(z),im_coth,m_re_coth,Error);

  if( Error == 0 )
    return cinterval(-m_re_coth,im_coth);
  else
    {
      if( Error == 1 )//Pole of cotangens hyperbolicus contained in argument.
	throw function_not_defined();
      else//Imaginary part of argument too far from [-pi/2,pi/2].
	throw function_not_defined();
    }
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*                          Multi-valued functions                           */
/*                                                                           */
/* ------------------------------------------------------------------------- */

interval arg(const cinterval& z) 
{
  interval hxl(0.0),hyl(0.0);
  interval hxu(0.0),hyu(0.0);
#ifdef FILIB_VERSION  
  if( ZERO_INTERVAL() == Re(z) ) 
    {
      if( ZERO_INTERVAL() == Im(z) ) throw function_not_defined();//arg not defined at (0,0).
      else                         //z purely imaginary  
	{
	  if( 0 < inf(Im(z)) )  return HALFPI();       //Im(z) positive
	  else 
	    if( 0 > sup(Im(z)) )  return -HALFPI(); //Im(z) negative
	    else return interval(-sup(HALFPI()),sup(HALFPI()));
	}
    }
  else 
    {    
      hxl = interval(inf(Re(z))); hxu = interval(sup(Re(z)));
      hyl = interval(inf(Im(z))); hyu = interval(sup(Im(z))); 
      if( ZERO_INTERVAL() <= Re(z) ) 
	{
	  if( ZERO_INTERVAL() <= Im(z) )            //all quadrants 
	    return 2.0 * interval(-sup(HALFPI()),sup(HALFPI()));
	  else 
	    {
	      if( inf(Im(z)) >= 0 )                //1. and 2. quadrant 
		return interval(inf(atan(hyl/hxu)),sup(atan(hyl/hxl)+PI()));//Minimum, maximum 
	      else                                   //3. and 4. quadrant 
		return interval(inf(atan(hyu/hxl)-PI()),sup(atan(hyu/hxu)));//Minimum, maximum 
	    }
	}
      else
	{ 
	  if( ZERO_INTERVAL() <= Im(z) ) 
	    {
	      if( inf(Re(z)) >= 0 )                //1. and 4. quadrant 
		{
		  if ( inf(Re(z)) == 0 ) 
		    return interval(-sup(HALFPI()),sup(HALFPI())); //Minimum, maximum 
		  else
		    return interval(inf(atan(hyl/hxl)),sup(atan(hyu/hxl))); //Minimum, maximum 
		}
	      else                              //2. and 3. quadrant (intersection) 
		{
		  if( sup(Re(z)) == 0 ) 
		    return interval(inf(HALFPI()),sup(3.0*HALFPI())); //Minimum, maximum 
		  else
		    return interval(inf(atan(hyu/hxu)),sup(atan(hyl/hxu))) + PI(); //Minimum, maximum 
		}
	    }
	  else                 //argument lies in single quadrant  
	    {
	      if( inf(Im(z)) >= 0 )                   //upper half plane
		{
		  if( inf(Re(z)) >= 0 )                     //1. quadrant 
		    {
		      if( inf(Re(z)) == 0 ) 
			return interval(inf(atan(hyl/hxu)),sup(HALFPI())); //Minimum, maximum 
		      else 
			return interval(inf(atan(hyl/hxu)),sup(atan(hyu/hxl))); //Minimum, maximum 
		    }
		  else                                            //2. quadrant 
		    {
		      if( sup(Re(z)) == 0 )  
			return interval(inf(HALFPI()),sup(atan(hyl/hxl)+PI())); //Maximum 
		      else
			return interval(inf(atan(hyu/hxu)),sup(atan(hyl/hxl))) + PI(); //Minimum, maximum 
		    }
		}
	      else                                          //lower half plane 
		{
		  if( inf(Re(z)) >= 0 )                     //4. quadrant 
		    {
		      if( inf(Re(z)) == 0 ) 
			return interval(-sup(HALFPI()),sup(atan(hyu/hxu))); //Minimum, maximum 
		      else
			return interval(inf(atan(hyl/hxl)),sup(atan(hyu/hxu))); //Minimum, maximum 
		    }
		  else                                            //3. quadrant 
		    {
		      if( sup(Re(z)) == 0 ) 
			return interval(-sup(HALFPI()),sup(atan(hxu/hyl)-PI())); //Minimum, maximum 
		      else
			return interval(inf(atan(hyu/hxl)),sup(atan(hyl/hxu))) - PI(); //Minimum, maximum 
		    }
		}
	    }
	}
    }
#else //C-XSC Version
  if( ZERO_INTERVAL() == Re(z) ) 
    {
      if( ZERO_INTERVAL() == Im(z) ) throw function_not_defined();//arg not defined at (0,0).
      else                         //z purely imaginary  
	{
	  if( 0 < Inf(Im(z)) )  return HALFPI();       //Im(z) positive
	  else 
	    if( 0 > Sup(Im(z)) )  return -HALFPI(); //Im(z) negative
	    else return interval(-Sup(HALFPI()),Sup(HALFPI()));
	}
    }
  else 
    {    
      hxl = interval(Inf(Re(z))); hxu = interval(Sup(Re(z)));
      hyl = interval(Inf(Im(z))); hyu = interval(Sup(Im(z))); 
      if( ZERO_INTERVAL() <= Re(z) ) 
	{
	  if( ZERO_INTERVAL() <= Im(z) )            //all quadrants 
	    return 2 * interval(-Sup(HALFPI()),Sup(HALFPI()));
	  else 
	    {
	      if( Inf(Im(z)) >= 0 )                //1. and 2. quadrant 
		return interval(Inf(atan(hyl/hxu)),Sup(atan(hyl/hxl)+PI()));//Minimum, maximum 
	      else                                   //3. and 4. quadrant 
		return interval(Inf(atan(hyu/hxl)-PI()),Sup(atan(hyu/hxu)));//Minimum, maximum 
	    }
	}
      else
	{ 
	  if( ZERO_INTERVAL() <= Im(z) ) 
	    {
	      if( Inf(Re(z)) >= 0 )                //1. and 4. quadrant 
		{
		  if ( Inf(Re(z)) == 0 ) 
		    return interval(-Sup(HALFPI()),Sup(HALFPI())); //Minimum, maximum 
		  else
		    return interval(Inf(atan(hyl/hxl)),Sup(atan(hyu/hxl))); //Minimum, maximum 
		}
	      else                              //2. and 3. quadrant (intersection) 
		{
		  if( Sup(Re(z)) == 0 ) 
		    return interval(Inf(HALFPI()),Sup(3*HALFPI())); //Minimum, maximum 
		  else
		    return interval(Inf(atan(hyu/hxu)),Sup(atan(hyl/hxu))) + PI(); //Minimum, maximum 
		}
	    }
	  else                 //argument lies in single quadrant  
	    {
	      if( Inf(Im(z)) >= 0 )                   //upper half plane
		{
		  if( Inf(Re(z)) >= 0 )                     //1. quadrant 
		    {
		      if( Inf(Re(z)) == 0 ) 
			return interval(Inf(atan(hyl/hxu)),Sup(HALFPI())); //Minimum, maximum 
		      else 
			return interval(Inf(atan(hyl/hxu)),Sup(atan(hyu/hxl))); //Minimum, maximum 
		    }
		  else                                            //2. quadrant 
		    {
		      if( Sup(Re(z)) == 0 )  
			return interval(Inf(HALFPI()),Sup(atan(hyl/hxl)+PI())); //Maximum 
		      else
			return interval(Inf(atan(hyu/hxu)),Sup(atan(hyl/hxl))) + PI(); //Minimum, maximum 
		    }
		}
	      else                                          //lower half plane 
		{
		  if( Inf(Re(z)) >= 0 )                     //4. quadrant 
		    {
		      if( Inf(Re(z)) == 0 ) 
			return interval(-Sup(HALFPI()),Sup(atan(hyu/hxu))); //Minimum, maximum 
		      else
			return interval(Inf(atan(hyl/hxl)),Sup(atan(hyu/hxu))); //Minimum, maximum 
		    }
		  else                                            //3. quadrant 
		    {
		      if( Sup(Re(z)) == 0 ) 
			return interval(-Sup(HALFPI()),Sup(atan(hxu/hyl)-PI())); //Minimum, maximum 
		      else
			return interval(Inf(atan(hyu/hxl)),Sup(atan(hyl/hxu))) - PI(); //Minimum, maximum 
		    }
		}
	    }
	}
    }
#endif
}  

inline interval abs(const interval& x,const interval& y) 
{
  return sqrt(sqr(x)+sqr(y));
} 

cinterval
#ifdef FILIB_VERSION
log
#else //C-XSC Version
ln 
#endif
(const cinterval& z) 
{
  flnumber re_min(0.0),im_min(0.0);
  flnumber re_max(0.0),im_max(0.0);
  interval re_ln(0.0);
#ifdef FILIB_VERSION
  flnumber srez( sup( Re(z) ) ),irez( inf( Re(z) ) ), simz( sup( Im(z) ) ),iimz( inf( Im(z) ) );

  if( srez < 0 ) re_min = -srez;
  else 
    if( irez < 0 ) re_min = 0;
    else re_min = irez;

  if( simz < 0 ) im_min = -simz;
  else 
    if( iimz < 0 ) im_min = 0;
    else im_min = iimz;

  if( re_min == 0 && im_min == 0 )//Singularity of logarithm in argument.
      throw function_not_defined();
  else 
    if( (irez <= 0)&&(iimz <= 0)&&(simz >= 0) )//Argument intersects negative real axis.
	throw function_not_defined();
    else
      {
	re_max = std::max( srez, -irez );
	im_max = std::max( simz, -iimz );

	re_ln = interval( inf( log( abs(interval(re_min),interval(im_min)) ) ),
			  sup( log( abs(interval(re_max),interval(im_max)) ) )  );

	return cinterval(re_ln,arg(z));
      }
#else //C-XSC Version
  flnumber srez( Sup( Re(z) ) ),irez( Inf( Re(z) ) ), simz( Sup( Im(z) ) ),iimz( Inf( Im(z) ) );

  if( srez < 0 ) re_min = -srez;
  else 
    if( irez < 0 ) re_min = 0;
    else re_min = irez;

  if( simz < 0 ) im_min = -simz;
  else 
    if( iimz < 0 ) im_min = 0;
    else im_min = iimz;
  
  if( re_min == 0 && im_min == 0 )//Singularity of logarithm in argument.
    throw function_not_defined();
  else 
    if( (irez <= 0)&&(iimz <= 0)&&(simz >= 0) )//Argument intersects negative real axis.
      throw function_not_defined();
    else
      { 
	re_max = max( srez, -irez );
	im_max = max( simz, -iimz );

	re_ln = interval( Inf( ln( abs(interval(re_min),interval(im_min)) ) ),
			  Sup( ln( abs(interval(re_max),interval(im_max)) ) )  );

	return cinterval(re_ln,arg(z));
      }
#endif
}

inline interval sign(const interval& x) 
{
#ifdef FILIB_VERSION
  flnumber ix = inf(x), sx = sup(x);
  flnumber inf, sup;

  if( ix > 0 ) inf = 1;
  else 
    if( ix ) inf = -1;
    else     inf =  0;

  if( sx > 0 ) sup = 1;
  else 
    if( sx ) sup = -1;
    else     sup =  0;

  return interval( inf, sup ); 
#else //C-XSC Version
  return interval(sign(Inf(x)),sign(Sup(x))); 
#endif
}

interval re_sqrt(const interval& x,const interval& y) 
{
  //Formulas for special quadrants 
#ifdef FILIB_VERSION
  if( sup(x) < 0 )
    {
      if( y == ZERO_INTERVAL() ) return ZERO_INTERVAL();
      else return INV_SQRT_2()*abs(y)/sqrt(abs(x,y)-x);
    }
  else 
    {
      if( y == ZERO_INTERVAL() ) return sqrt(x);
      else return INV_SQRT_2()*sqrt(abs(x,y)+x);
    }
#else //C-XSC Version
  if( Sup(x) < 0 )
    {
      if( y == ZERO_INTERVAL() ) return ZERO_INTERVAL();
      else return INV_SQRT_2()*abs(y)/sqrt(abs(x,y)-x);
    }
  else 
    {
      if( y == ZERO_INTERVAL() ) return sqrt(x);
      else return INV_SQRT_2()*sqrt(abs(x,y)+x);
    }
#endif
}

interval im_sqrt(const interval& x,const interval& y) 
{
  //Formulas for special quadrants 
#ifdef FILIB_VERSION 
  if( sup(x) < 0 )
    {
      if( y == ZERO_INTERVAL() ) return sqrt(-x); //arg(z) = pi
      else return INV_SQRT_2()*sign(y)*sqrt(abs(x,y)-x);
    }
  else
    {
      if( y == ZERO_INTERVAL() ) return ZERO_INTERVAL(); 
      else return INV_SQRT_2()*y/sqrt(abs(x,y)+x);
    }
#else //C-XSC Version
  if( Sup(x) < 0 )
    {
      if( y == ZERO_INTERVAL() ) return sqrt(-x); //arg(z) = pi
      else return INV_SQRT_2()*sign(y)*sqrt(abs(x,y)-x);
    }
  else
    {
      if( y == ZERO_INTERVAL() ) return ZERO_INTERVAL(); 
      else return INV_SQRT_2()*y/sqrt(abs(x,y)+x);
    }
#endif
}

cinterval sqrt(const cinterval& z) 
{
#ifdef FILIB_VERSION
  interval rwert(0.0),iwert(0.0);
  interval hxl( inf(Re(z)) ),hyl( inf(Im(z)) );
  interval hxu( sup(Re(z)) ),hyu( sup(Im(z)) );
  
  if( (sup(hxl) < 0.0) && (inf(hyl) < 0.0) && (inf(hyu) >= 0.0) ) 
    throw function_not_defined();//Argument intersects negative real axis

  if( sup(hyu) < 0.0 )  //lower half plane (without intersection)
    {
      rwert = re_sqrt(hxl,hyu); //Minimum 
      rwert = interval(inf(rwert),sup(re_sqrt(hxu,hyl))); //with maximum 
      iwert = im_sqrt(hxl,hyl); //Minimum 
      iwert = interval(inf(iwert),sup(im_sqrt(hxu,hyu))); //with maximum 
    } 
  else
    {
      if( inf(hyl) > 0.0 )  //upper half plane 
	{
	  rwert = re_sqrt(hxl,hyl); //Minimum 
	  rwert = interval(inf(rwert),sup(re_sqrt(hxu,hyu)));//with maximum 
	  iwert = im_sqrt(hxu,hyl); //Minimum 
	  iwert = interval(inf(iwert),sup(im_sqrt(hxl,hyu)));//with maximum 
	} 
      else  //Zero contained im imaginary part of argument! 
	{
	  //right half plane (no intersection) 
	  rwert = sqrt(hxl); //Minimum 
	  rwert = interval(inf(rwert),std::max(sup(re_sqrt(hxu,hyu)),sup(re_sqrt(hxu,hyl))));
	  iwert = im_sqrt(hxl,hyl); //Minimum 
	  iwert = interval(inf(iwert),sup(im_sqrt(hxl,hyu))); //with maximum 
	}
    }
  return cinterval(rwert,iwert);
#else //C-XSC Version
  interval rwert(0.0),iwert(0.0);
  interval hxl( InfRe(z) ),hyl( InfIm(z) );
  interval hxu( SupRe(z) ),hyu( SupIm(z) );

  if( (Sup(hxl) < 0.0) && (Sup(hyl) <= 0.0) && (Inf(hyu) >= 0.0) ) 
    throw function_not_defined();//Argument intersects negative real axis
  
  if( Sup(hyu) < 0 )  //lower half plane (without intersection)
    {
      rwert = re_sqrt(hxl,hyu); //Minimum 
      rwert = interval(Inf(rwert),Sup(re_sqrt(hxu,hyl))); //with maximum 
      iwert = im_sqrt(hxl,hyl); //Minimum 
      iwert = interval(Inf(iwert),Sup(im_sqrt(hxu,hyu))); //with maximum 
    } 
  else
    {
      if( Inf(hyl) >= 0 )  //upper half plane 
	{
	  rwert = re_sqrt(hxl,hyl); //Minimum 
	  rwert = interval(Inf(rwert),Sup(re_sqrt(hxu,hyu)));//with maximum 
	  iwert = im_sqrt(hxu,hyl); //Minimum 
	  iwert = interval(Inf(iwert),Sup(im_sqrt(hxl,hyu)));//with maximum 
	} 
       else  //Zero contained im imaginary part of argument! 
	{
	  //right half plane (no intersection) 
	  rwert = sqrt(hxl); //Minimum 
	  rwert = interval(Inf(rwert),max(Sup(re_sqrt(hxu,hyu)),Sup(re_sqrt(hxu,hyl))));
	  iwert = im_sqrt(hxl,hyl); //Minimum 
	  iwert = interval(Inf(iwert),Sup(im_sqrt(hxl,hyu))); //with maximum 
	} 
    }
  return cinterval(rwert,iwert);
#endif
}

std::list<cinterval> sqrt_all(const cinterval& z) 
{
#ifdef FILIB_VERSION
  interval rwert(0.0),iwert(0.0);
  interval hxl( inf(Re(z)) ),hyl( inf(Im(z)) );
  interval hxu( sup(Re(z)) ),hyu( sup(Im(z)) );
  
  if( sup(hyu) < 0.0 )  //lower half plane
    {
      rwert = re_sqrt(hxl,hyu); //Minimum 
      rwert = interval(inf(rwert),sup(re_sqrt(hxu,hyl))); //with maximum 
      iwert = im_sqrt(hxl,hyl); //Minimum 
      iwert = interval(inf(iwert),sup(im_sqrt(hxu,hyu))); //with maximum 
    } 
  else
    {
      if( inf(hyl) > 0.0 )  //upper half plane 
	{
	  rwert = re_sqrt(hxl,hyl); //Minimum 
	  rwert = interval(inf(rwert),sup(re_sqrt(hxu,hyu)));//with maximum 
	  iwert = im_sqrt(hxu,hyl); //Minimum 
	  iwert = interval(inf(iwert),sup(im_sqrt(hxl,hyu)));//with maximum 
	} 
      else  //Zero contained im imaginary part of argument! 
	{
	  //right half plane (no intersection) 
	  rwert = re_sqrt(hxl,ZERO_INTERVAL()); //Minimum 
	  rwert = interval(inf(rwert),std::max(sup(re_sqrt(hxu,hyu)),sup(re_sqrt(hxu,hyl))));
	  iwert = im_sqrt(hxl,hyl); //Minimum 
	  iwert = interval(inf(iwert),sup(im_sqrt(hxl,hyu))); //with maximum 
	}
    }
#else //C-XSC Version
  interval rwert(0.0),iwert(0.0);
  interval hxl( InfRe(z) ),hyl( InfIm(z) );
  interval hxu( SupRe(z) ),hyu( SupIm(z) );

  if( Sup(hyu) < 0 )  //lower half plane
    {
      rwert = re_sqrt(hxl,hyu); //Minimum 
      rwert = interval(Inf(rwert),Sup(re_sqrt(hxu,hyl))); //with maximum 
      iwert = im_sqrt(hxl,hyl); //Minimum 
      iwert = interval(Inf(iwert),Sup(im_sqrt(hxu,hyu))); //with maximum 
    } 
  else
    {
      if( Inf(hyl) >= 0 )  //upper half plane 
	{
	  rwert = re_sqrt(hxl,hyl); //Minimum 
	  rwert = interval(Inf(rwert),Sup(re_sqrt(hxu,hyu)));//with maximum 
	  iwert = im_sqrt(hxu,hyl); //Minimum 
	  iwert = interval(Inf(iwert),Sup(im_sqrt(hxl,hyu)));//with maximum 
	} 
       else  //Zero contained im imaginary part of argument! 
	{
	  //right half plane (no intersection) 
	  rwert = re_sqrt(hxl,ZERO_INTERVAL()); //Minimum 
	  rwert = interval(Inf(rwert),max(Sup(re_sqrt(hxu,hyu)),Sup(re_sqrt(hxu,hyl))));
	  iwert = im_sqrt(hxl,hyl); //Minimum 
	  iwert = interval(Inf(iwert),Sup(im_sqrt(hxl,hyu))); //with maximum 
	} 
    }
#endif
  
  cinterval w(rwert,iwert);
  std::list<cinterval> res;
  res.push_back( w); 
  res.push_back(-w);

  return res;
}

cinterval root(const cinterval& z,int n) 
{
  if( n==0 ) return
#ifdef FILIB_VERSION
	       cinterval( ONE_INTERVAL(), ZERO_INTERVAL() )
#else //C-XSC Version
	       cinterval( 1.0 )
#endif
	       ;
  if( n==1 ) return z;
  if( n==2 ) return sqrt( z );
  if( n>=3 ) return pow( z, 1.0/interval(n) );
}

/*
  For use in 'root_all'
*/

interval root(const interval& i,unsigned int n)
{
  flnumber inf_i(
#ifdef FILIB_VERSION
		 inf(i)
#else
		 Inf(i)
#endif
		 );

  if( inf_i < 0.0 ) throw function_not_defined();
  else if( inf_i > 0 ) return exp( 1.0/
#ifdef FILIB_VERSION
				   interval(n) * log(i)
#else
				   interval(int(n)) * ln (i)
#endif
				   );
  else
    {
      flnumber sup_res(
#ifdef FILIB_VERSION
		       sup(exp( 1.0/interval(n) * log(interval(sup(i))) ))
#else
		       Sup(exp( 1.0/interval(int(n)) * ln(interval(Sup(i))) ))
#endif
		       );

      return interval( 0, sup_res );
    }
}

std::list<cinterval> root_all(const cinterval& z,unsigned int n)
{
  std::list<cinterval> res;

  if( n == 0 )
#ifdef FILIB_VERSION
    {
      res.push_back( cinterval( ONE_INTERVAL(), ZERO_INTERVAL() ) );
      return res;
    }
#else //C-XSC Version
    {
      res.push_back( cinterval( 1.0 ) );
      return res;
    }
#endif
  else if( n == 1 )
    {
      res.push_back(z);
      return res;
    }
  else if( n == 2 ) return sqrt_all( z );
  else
    {
      interval root_abs_z(root(abs(z),n)), arg_z(arg(z));

      for(unsigned int i = 0; i < n; i++)
	{
	  interval arg_z_plus_2_i_pi_div_n( (arg_z+2*i*PI())/
#ifdef FILIB_VERSION
						     n
#else
						     int(n)
#endif
					   );

	  res.push_back( cinterval( root_abs_z * cos(arg_z_plus_2_i_pi_div_n) ,
				    root_abs_z * sin(arg_z_plus_2_i_pi_div_n)  ) );
	}

      return res;
    }
}

cinterval sqr(const cinterval& z) 
{
  interval 
    A( Re(z) ),
    B( Im(z) );

  return cinterval( sqr( A ) - sqr( B ) , 2.0 * A * B );
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*    arcsin from diploma thesis of Gabriele Buehler                         */
/*                                                                           */
/* ------------------------------------------------------------------------- */

inline interval g(const interval& s) 
{
  return sqrt( 1.0 + s ) - 1.0;
}

inline interval s_re(const interval& ix,const interval& iy) 
{
  return 2.0 * iy / ( ix*ix + iy*iy - 1.0 );
}

interval re_arcsin(const flnumber& x,const flnumber& y)
{
  //Real part 
  if( x == 0.0 ) 
    return interval( 0.0 );
  else
    {
      interval ix( x );

      if( y == 0.0 )
	{
	  if( x >= 1.0 )
	    return HALFPI();
	  else
	    return asin(ix);
	}
      else
	{
	  interval 
	    hilf1(0.0),
	    hilf2(0.0),
	    hilf3(0.0),
	    nenner(0.0),
	    zaehler(0.0),
	    sqrs_re(0.0),
	    iy( y );

	  if( (x*x - y*y) < 0.5 )
	    {
	      hilf1 = sqrt(sqr(ix+1.0) + iy*iy );
	      hilf2 = sqrt(sqr(ix-1.0) + iy*iy );
	      nenner = hilf1 + hilf2;
	      hilf3 = (2.0 * ix)/nenner;
	      return asin(hilf3);
	    }
	  else                    //x*x - y*y > 0.5 
	    {
	      hilf1 = ix*ix + iy*iy;
	      hilf2 = hilf1 - 1.0;
	      nenner = hilf1 + 1.0 + sqrt( sqr( hilf2 ) + 4.0 * iy*iy );
	      sqrs_re =sqr( s_re(ix,iy) );
#ifdef FILIB_VERSION
	      if( sup(hilf1) < 1.0 )
#else //C-XSC Version
	      if( Sup(hilf1) < 1.0 )
#endif
		{
		  zaehler = 2.0 - 2.0 * ix*ix - hilf2 * g(sqrs_re);
		  hilf3 = sqrt(zaehler / nenner);
		  return HALFPI() - asin(hilf3);
		}
	      else                //x*x - y*y > 1.0 
		{
		  zaehler = 2.0 * iy*iy + hilf2 * g(sqrs_re);
		  hilf3 = sqrt(zaehler / nenner);
		  return HALFPI() - asin(hilf3);
		}
	    }
	}
    }
}

interval im_arcsin(const flnumber& x,const flnumber& y)
{
  /* Interval computation of imaginary part of arcsin(z) */
  interval 
    hilf1(0.0),
    hilf2(0.0),
    hilf3(0.0),
    hilf4(0.0),
    hilf5(0.0),
    hilf6(0.0),
    ix( abs(interval(x)) ),
    iy( y );
  interval t(0.0),r(0.0);
  flnumber xc( fabs(x) );
  
  //Imaginary part 
  if( y == 0.0 )           //y = 0.0 
    {
      if( xc <= 1.0 )
	return interval( 0.0 );
      else                           //x> 1.0 
	{
	  if( xc > 1.1 )       //x> 1.1 
	    {
	      t = 0.5 * ( (ix + 1.0) + (ix - 1.0) );
	      return
#ifdef FILIB_VERSION
		log
#else //C-XSC Version
		ln
#endif
		(t + sqrt(sqr(t) - 1.0));
	    }
	  else                          //1< x < 1.1 
	    {
	      t = ix;
	      r = t - 1.0;
	      return
#ifdef FILIB_VERSION
		log
#else //C-XSC Version
		ln
#endif
		(1.0 + (r + sqrt(sqr(r) + 2.0 * r)));
	    }
	}
    }
  else                       //y <> 0.0 
    {
      hilf1 = ix*ix + iy*iy;
      hilf2 = hilf1 + 1.0;
      if( xc == 0.0 )         //x= 0.0 
	return
#ifdef FILIB_VERSION
	  log
#else //C-XSC Version
	  ln
#endif
	  (sqrt(1.0 + iy * iy) + iy);
      else                   //x <> 0.0 
        {
          t  = 0.5 * ( sqrt( hilf2 + 2.0 * ix) + sqrt( hilf2 - 2.0 * ix ));
#ifdef FILIB_VERSION
          if( sup(t) <= 1.1 )   //t <= 1.1 
#else //C-XSC Version
          if( Sup(t) <= 1.1 )   //t <= 1.1 
#endif
            {
              if( xc == 1.0 )   //x = 1.0 
                r =  g(iy * iy /4.0) + 0.5 * iy;
              else             //x <>1.0 
                {
                  hilf3 = iy/(ix + 1.0);
                  hilf4 = (ix + 1.0) * g(hilf3 * hilf3);
                  if( xc < 1.0 ) //x < 1.0 
                    {
                      hilf5  = iy/( 1.0 -ix);
                      hilf6  = (1.0 - ix) * g(hilf5*hilf5);
                      r = 0.5 * ( hilf4 + hilf6);
                    }
                  else              //x > 1.0 
                    r  = 0.5 * (hilf4 + (ix - 1.0) + sqrt(hilf2 -2.0 * ix));
                }
	      return
#ifdef FILIB_VERSION
		log
#else //C-XSC Version
		ln
#endif
		(1.0 + (r + sqrt(sqr(r) + 2.0 * r)));
            }
          else                    //t > 1.1 
	    return
#ifdef FILIB_VERSION
	      log
#else //C-XSC Version
	      ln
#endif
	      (t + sqrt(sqr(t) - 1.0));
        }
    }
}

inline interval fortsetz_asin(const flnumber& x,const flnumber& y)
{
  interval hilf( re_arcsin(x,y) );

  if( hilf ==
#ifdef FILIB_VERSION
       ZERO_INTERVAL()
#else //C-XSC Version
       0.0 
#endif
      )

    return PI();
  else
    return PI() - hilf;

}

inline void z(const cinterval& z,flnumber& re_l,flnumber& re_u,flnumber& im_l,flnumber& im_u)
{
#ifdef FILIB_VERSION
  re_l = inf(Re(z));
  re_u = sup(Re(z));
  im_l = inf(Im(z));
  im_u = sup(Im(z));
#else //C-XSC Version
  re_l = Inf(Re(z));
  re_u = Sup(Re(z));
  im_l = Inf(Im(z));
  im_u = Sup(Im(z));
#endif
}

inline void z(const interval& x,const interval& y,flnumber& x_l,flnumber& x_u,flnumber& y_l,flnumber& y_u)
{
#ifdef FILIB_VERSION
  x_l = inf(x);
  x_u = sup(x);
  y_l = inf(y);
  y_u = sup(y);
#else //C-XSC Version
  x_l = Inf(x);
  x_u = Sup(x);
  y_l = Inf(y);
  y_u = Sup(y);
#endif
}

interval real_asin(const cinterval& c) 
{
  flnumber xl(0.0),xu(0.0),yl(0.0),yu(0.0),maxy(0.0),max(0.0);
  flnumber null( 0.0 ),eins( 1.0 );
  bool re_spiegel( false );
  interval ergxl(0.0),ergxu(0.0),ergx(0.0);
  cinterval c1( c );
  
  z(c1,xl,xu,yl,yu);
  //Die Funktion z weisst den reellen Zahlen xl,xu,yl,yu die entsprechen 
  //Zahlen des komplexen Intervalls c1 zu 

  //Arguments in left half plane are reflected to right half plane
  if( xu <= null )
    {
      Re(c1) = -Re(c1);
      re_spiegel = true;
    }
  else re_spiegel = false;
  //Arguments in lower half plane are reflected to upper half plane
  if( yu <= null ) Im(c1) = -Im(c1);

  z(c1,xl,xu,yl,yu);
  if( -yl >= yu )   //Compute maximum of yl and yu (absolute values) 
    maxy = -yl;
  else
    maxy = yu;
  
  //Discuss special cases
  if( yl >= null )      //upper half plane 
    {
      if( xl >= null )    //1. Quadrant 
	{
	  ergxl = re_arcsin(xl,yu);
	  ergxu = re_arcsin(xu,yl);
	}
      else                 //0 in real part of c 
	{
	  ergxl = -re_arcsin(-xl,yl);
	  ergxu = re_arcsin(xu,yl);
	}
    }
  else                  //0 in imaginary part of c 
    {
      if( xl >= eins )    //Intersection without winding point inside argument
	{
	  ergxl = re_arcsin(xl,yl);
	  ergxu = fortsetz_asin(xl,yu);
	}
      else                  //xl < 1.0 
	{
	  if( xl < -eins ) //xl < lower winding point 
	    {
	      if( xu < eins )   //only winding point (-1,0) within argument
		{
		  ergxl = -HALFPI();
		  ergxu = re_arcsin(xu,null);
		}
	      else                //both winding points within argument
		{
		  ergxl = -HALFPI();
		  ergxu =  HALFPI();
		}
	    }
	  else   
	    {                //-1 <= xl < 1 
	      if( xl <= null )    //xl <= 0 
		{
		  if( xu <= eins ) //No intersection, but origin contained 
		    {
		      ergxl = -re_arcsin(-xl,null);
		      ergxu = re_arcsin(xu,null);
		    }
		  else               //winding point (1.0) and origin contained
		    {
		      ergxl = -re_arcsin(-xl,null);
		      ergxu =  HALFPI();
		    }
		}
	      else               //0 <= xl < 1 
		{
		  if( xu < eins ) //No WP, origin contained 
		    {
		      ergxl = re_arcsin(xl,maxy);
		      ergxu = re_arcsin(xu,null);
		    }
		  else            //WP (1,0) 
		    {
		      ergxl = re_arcsin(xl,max);
		      ergxu = HALFPI();
		    }
		}
	    }
	}
    }
#ifdef FILIB_VERSION
  ergx = interval(inf(ergxl),sup(ergxu));
#else //C-XSC Version
  ergx = interval(Inf(ergxl),Sup(ergxu));
#endif
  if( re_spiegel ) 
    return -ergx;
  else
    return ergx;
}

interval imag_asin(const cinterval& c) 
{
  flnumber xl(0.0),xu(0.0),yl(0.0),yu(0.0),maxx(0.0),maxy(0.0);
  flnumber null( 0.0 ),eins( 1.0 );
  bool im_spiegel( false );
  interval ergyl(0.0),ergyu(0.0),ergy(0.0);
  cinterval c1( c );

  z(c1,xl,xu,yl,yu);
  //Argument in left half plane reflected to right half plane
  if( xu <= null ) Re(c1) = -Re(c1);
  //Argument in lower half plane reflected to upper half plane
  if( yu <= null )
    {
      Im(c1) = -Im(c1);
      im_spiegel = true;
    }
  else  im_spiegel = false;

  z(c1,xl,xu,yl,yu);
  if( -xl >= yu )  //Maximum of xl and xu (absolute value)
    maxx = -xl;
  else
    maxx  = xu;

  if( -yl >= yu ) //Maximum of yl and yu (absolute value)
    maxy = -yl;
  else
    maxy = yu;
  
  //Discuss special cases
  if( yl >= null )            //upper half plane 
    {
      if( xl > null )           //1. Quadrant 
	{
	  ergyl = im_arcsin(xl,yl);
	  ergyu = im_arcsin(xu,yu);
	}
      else                   //0 in real part of c 
	{
	  ergyl = im_arcsin(null,yl);
	  ergyu = im_arcsin(maxx,yu);
	}
    }
  else                     //0 in imaginary part of c 
    {
      if( xl >= eins )       //intersection without WP in argument 
	{
	  ergyl = -im_arcsin(xu,maxy);
	  ergyu = im_arcsin(xl,null);
	}
      else                   //xl < 1.0 
	{
	  ergyl = -im_arcsin(maxx,-yl);
	  ergyu = im_arcsin(maxx,yu);
	}
    }
#ifdef FILIB_VERSION  
  ergy = interval(inf(ergyl),sup(ergyu));
#else //C-XSC Version
  ergy = interval(Inf(ergyl),Sup(ergyu));
#endif
  if( im_spiegel )
    return -ergy;
  else
    return ergy;
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*     Computation of  arccos, arsinh and arcosh  based on  arcsin           */
/*                                                                           */
/*     Arccos(z) = -/+ ( Arcsin(z) - pi/2 )                                  */
/*     Arsinh(z) = i * Arcsin( -i * z ) ( mod i*2*pi )                       */
/*     Arcosh(z) = -/+ i * Arccos(z)                                         */
/*                                                                           */
/*     Only principal values are computed.                                   */
/*                                                                           */
/* ------------------------------------------------------------------------- */

cinterval asin(const cinterval& c) 
{
#ifdef FILIB_VERSION
  if(  (inf(Im(c)) <= 0.0) && (0.0 <= sup(Im(c))) &&   //Zero in imaginary part and
       ((inf(Re(c)) < -1.0) || (sup(Re(c)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();
#else //C-XSC Version
  if(  (Inf(Im(c)) <= 0.0) && (0.0 <= Sup(Im(c))) &&   //Zero in imaginary part and
       ((Inf(Re(c)) < -1.0) || (Sup(Re(c)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();
#endif
  return cinterval(real_asin(c),imag_asin(c));
}

cinterval acos(const cinterval& c) 
{
#ifdef FILIB_VERSION
  if(  (inf(Im(c)) <= 0.0) && (0.0 <= sup(Im(c))) &&   //Zero in imaginary part and
       ((inf(Re(c)) < -1.0) || (sup(Re(c)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();
#else //C-XSC Version
  if(  (Inf(Im(c)) <= 0.0) && (0.0 <= Sup(Im(c))) &&   //Zero in imaginary part and
       ((Inf(Re(c)) < -1.0) || (Sup(Re(c)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();
#endif
  return cinterval(real_asin(c)-HALFPI(),imag_asin(c));
}

cinterval asinh(const cinterval& c) 
{
#ifdef FILIB_VERSION
  if(  (inf(Re(c)) <= 0.0) && (0.0 <= sup(Re(c))) &&   //Zero in real part and
       ((inf(Im(c)) < -1.0) || (sup(Im(c)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();
#else //C-XSC Version
  if(  (Inf(Re(c)) <= 0.0) && (0.0 <= Sup(Re(c))) &&   //Zero in real part and
       ((Inf(Im(c)) < -1.0) || (Sup(Im(c)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();
#endif
  cinterval hc( -Im(c), Re(c) ); //hc = i*c 
  
  return cinterval(imag_asin(hc),-real_asin(hc)); //arsinh(c) = -i*arcsin(i*c) 
}

cinterval acosh(const cinterval& c) 
{
#ifdef FILIB_VERSION
  if(  (inf(Im(c)) <= 0.0) && (0.0 <= sup(Im(c))) &&   //Zero in imaginary part and
       (inf(Re(c)) < 1.0)                           ) //real part intersects (-INFINITY,1)
    throw function_not_defined();
#else //C-XSC Version
  if(  (Inf(Im(c)) <= 0.0) && (0.0 <= Sup(Im(c))) &&   //Zero in imaginary part and
       (Inf(Re(c)) < 1.0)                           ) //real part intersects (-INFINITY,1)
    throw function_not_defined();
#endif
  cinterval hc( real_asin(c)-HALFPI(),imag_asin(c) );
  
  return cinterval(-Im(hc),Re(hc));               //arcosh = i * arccos(c) 
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*     Auxiliary functions for inverse tangent                               */
/*                                                                           */
/* ------------------------------------------------------------------------- */

interval qbetrag(const flnumber& lx,const flnumber& ux,const flnumber& ly,const flnumber& uy)
{
  flnumber infqbetr(0.0),supqbetr(0.0);
#ifdef FILIB_VERSION
  flnumber tmp1,tmp2;
  
  if( lx >= 0 ) 
    {
      if( ly > 0 ) 
	{
	  MUL_DOWN( tmp1, lx, lx );
	  MUL_DOWN( tmp2, ly, ly );
	  ADD_DOWN( infqbetr, tmp1, tmp2 );

	  MUL_UP( tmp1, ux, ux );
	  MUL_UP( tmp2, uy, uy );
	  ADD_UP( supqbetr, tmp1, tmp2 );
	}
      else
	if( uy < 0 )
	  {
	    MUL_DOWN( tmp1, lx, lx );
	    MUL_DOWN( tmp2, uy, uy );
	    ADD_DOWN( infqbetr, tmp1, tmp2 );
	    
	    MUL_UP( tmp1, ux, ux );
	    MUL_UP( tmp2, ly, ly );
	    ADD_UP( supqbetr, tmp1, tmp2 );
	  }
	else
	  {
	    MUL_DOWN( infqbetr, lx, lx );
	    if ( -ly < uy ) 
	      {
		MUL_UP( tmp1, ux, ux );
		MUL_UP( tmp2, uy, uy );
		ADD_UP( supqbetr, tmp1, tmp2 );
	      }
	    else 
	      {
		MUL_UP( tmp1, ux, ux );
		MUL_UP( tmp2, ly, ly );
		ADD_UP( supqbetr, tmp1, tmp2 );
	      }
	  }
      return interval(infqbetr,supqbetr);
    }
  else
    if( ux <= 0 ) return qbetrag( -ux,-lx,ly,uy );
    else
      {
	if( -lx < ux ) 
	  return qbetrag( 0.0,ux,ly,uy );
	else
	  return qbetrag( 0.0,-lx,ly,uy );
      }
#else //C-XSC Version
  dotprecision akku(0.0);
  
  if( lx >= 0 ) 
    {
      if( ly > 0 ) 
	{
	  akku = 0.0;
	  accumulate( akku , lx , lx ); accumulate( akku , ly , ly );
	  infqbetr = rnd( akku , RND_DOWN );
	  akku = 0.0;
	  accumulate( akku , ux , ux ); accumulate( akku , uy , uy );
	  supqbetr = rnd( akku , RND_UP );
	}
      else
	if( uy < 0 )
	  {
	    akku = 0.0;
	    accumulate( akku , lx , lx ); accumulate( akku , uy , uy );
	    infqbetr = rnd( akku , RND_DOWN );
	    akku = 0.0;
	    accumulate( akku , ux , ux ); accumulate( akku , ly , ly );
	    supqbetr = rnd( akku , RND_UP );
	  }
	else
	  {
	    
	    infqbetr = multdown( lx,lx );
	    if ( -ly < uy ) {
	      akku = 0.0;
	      accumulate( akku , ux , ux ); accumulate( akku , uy , uy );
	      supqbetr = rnd( akku , RND_UP ); 
	    }
	    else 
	      {
		akku = 0.0;
		accumulate( akku , ux , ux ); accumulate( akku , ly , ly );
		supqbetr = rnd( akku , RND_UP );
	      }
	  }
      return interval(infqbetr,supqbetr);
    }
  else
    if( ux <= 0 ) return qbetrag( -ux,-lx,ly,uy );
    else
      {
	if( -lx < ux ) 
	  return qbetrag( 0.0,ux,ly,uy );
	else
	  return qbetrag( 0.0,-lx,ly,uy );
      }
#endif
}

interval re_fun(const flnumber& hx,const flnumber& hy)
{
  interval qbetr( qbetrag(hx,hx,hy,hy) );
#ifdef FILIB_VERSION
  if( 1 > sup(qbetr) )
    return atan( interval(2*hx)/(1.0-qbetr) ) / 2.0;
  else
    {
      if( 1 < inf(qbetr) )
	return (atan( interval(2*hx)/(1.0-qbetr) ) - PI() ) / 2.0 ;
      else
	return -QUARTERPI();
    }
#else //C-XSC Version
  if( 1 > Sup(qbetr) )
    return atan( interval(2*hx)/(1-qbetr) ) / 2;
  else
    {
      if( 1 < Inf(qbetr) )
	return (atan( interval(2*hx)/(1-qbetr) ) - PI() )/2 ;
      else
	return -QUARTERPI();
    }
#endif
}

interval re_fun(const interval& hx,const flnumber& hy)
{
#ifdef FILIB_VERSION
  interval qbetr( qbetrag(inf(hx),sup(hx),hy,hy) );
  
  if( 1 > sup(qbetr) )
    return atan( 2.0*hx/(1.0-qbetr) ) / 2.0;
  else
    {
      if( 1 < inf(qbetr) ) 
	return (atan( 2.0*hx/(1.0-qbetr) ) - PI() ) / 2.0;
      else 
	return -QUARTERPI();
    }
#else //C-XSC Version
  interval qbetr( qbetrag(Inf(hx),Sup(hx),hy,hy) );
  
  if( 1 > Sup(qbetr) )
    return atan( 2*hx/(1-qbetr) ) / 2;
  else
    {
      if( 1 < Inf(qbetr) ) 
	return (atan( 2*hx/(1-qbetr) ) - PI() )/2;
      else 
	return -QUARTERPI();
    }
#endif
}

interval im_fun(const flnumber& hx,const flnumber& hy)
{
#ifdef FILIB_VERSION
  flnumber ly,uy;

  SUB_DOWN( ly, 1.0, hy );
  SUB_UP( uy, 1.0, hy );

  interval nenner( qbetrag(hx,hx,ly,uy) );

  if( 0 < inf(nenner) ) 
    return log( abs(1.0+4.0*interval(hy)/nenner) ) / 4.0;
  else//Argument too near to singularity.
    throw function_not_defined();
#else //C-XSC Version
  flnumber ly( subdown( 1 , hy ) ),uy( subup( 1 , hy ) );
  interval nenner( qbetrag(hx,hx,ly,uy) );

  if( 0 < Inf(nenner) ) 
    return ln( abs(1+4*interval(hy)/nenner) ) / 4;
  else//Argument too near to singularity.
    throw function_not_defined();
#endif
}

interval im_fun(const flnumber& hx,const interval& hy) 
{
#ifdef FILIB_VERSION
  flnumber ly,uy;

  SUB_DOWN( ly, 1.0, sup(hy) );
  SUB_UP( uy, 1.0, inf(hy) );

  interval nenner( qbetrag(hx,hx,ly,uy) );

  if( 0 < inf(nenner) )   
    return log(abs(1.0+4.0*hy/(sqr(interval(hx))+sqr(1.0-hy)))) / 4.0;
  else//Argument too near to singularity.
    throw function_not_defined();
#else //C-XSC Version
  flnumber ly( subdown( 1 , Sup(hy) ) ),uy( subup( 1 , Inf(hy) ) );
  interval nenner( qbetrag(hx,hx,ly,uy) );

  if( 0 < Inf(nenner) )   
    return ln(abs(1+4*hy/(sqr(interval(hx))+sqr(1-hy)))) / 4;
  else//Argument too near to singularity.
    throw function_not_defined();
#endif
}

void left_side(const interval& hx,const interval& hy,interval& re_a,interval& im_a) 
{
#ifdef FILIB_VERSION
  if( sup(sqr(interval(sup(hy))))<= inf( 1.0 + sqr(interval(sup(hx)))) )
    { //argument unter Extremalkurve y=sqrt(1+sqr(x)) 
      re_a = interval( inf(re_fun(inf(hx),sup(hy))),sup(re_fun(sup(hx),inf(hy))) ); //Minimum, maximum 
      im_a = interval( inf(im_fun(inf(hx),inf(hy))),sup(im_fun(sup(hx),sup(hy))) ); //Minimum, maximum 
    }
  else
    {
      if( inf(sqr(interval(inf(hy))))>= sup( 1.0 + sqr(interval(inf(hx)))) )
	{ //argument ueber Extremalkurve y=sqrt(1+sqr(x)) 
	  re_a = interval( inf(re_fun(sup(hx),sup(hy))),sup(re_fun(inf(hx),inf(hy))) ); //Minimum, maximum 
	  im_a = interval( inf(im_fun(inf(hx),sup(hy))),sup(im_fun(sup(hx),inf(hy))) ); //Minimum, maximum 
	}
      else //argument intersects extremal curve y=sqrt(1+sqr(x)) 
	{
	  re_a = re_fun(sup(hx),sup(hy));                //possible minimum 
	  re_a = hull( re_a, re_fun(inf(hx),sup(hy)) );  //possible minimum 
	  re_a = hull( re_a, re_fun(hx,inf(hy))      );  //maximum 
	  im_a = im_fun(inf(hx),inf(hy));                //possible minimum 
	  im_a = hull( im_a, im_fun(inf(hx),sup(hy)) );  //possible minimum 
	  im_a = hull( im_a, im_fun(sup(hx),hy)      );  //maximum 
	}
    }
#else //C-XSC Version
  if( Sup(sqr(interval(Sup(hy))))<= Inf( 1.0 + sqr(interval(Sup(hx)))) )
    { //argument unter Extremalkurve y=sqrt(1+sqr(x)) 
      re_a = interval( Inf(re_fun(Inf(hx),Sup(hy))),Sup(re_fun(Sup(hx),Inf(hy))) ); //Minimum, maximum 
      im_a = interval( Inf(im_fun(Inf(hx),Inf(hy))),Sup(im_fun(Sup(hx),Sup(hy))) ); //Minimum, maximum 
    }
  else
    {
      if( Inf(sqr(interval(Inf(hy))))>= Sup( 1.0 + sqr(interval(Inf(hx)))) )
	{ //argument ueber Extremalkurve y=sqrt(1+sqr(x)) 
	  re_a = interval( Inf(re_fun(Sup(hx),Sup(hy))),Sup(re_fun(Inf(hx),Inf(hy))) ); //Minimum, maximum 
	  im_a = interval( Inf(im_fun(Inf(hx),Sup(hy))),Sup(im_fun(Sup(hx),Inf(hy))) ); //Minimum, maximum 
	}
      else //argument intersects extremal curve y=sqrt(1+sqr(x)) 
	{
	  re_a = re_fun(Sup(hx),Sup(hy));   //possible minimum 
	  re_a |= re_fun(Inf(hx),Sup(hy));  //possible minimum 
	  re_a |= re_fun(hx,Inf(hy));       //maximum 
	  im_a = im_fun(Inf(hx),Inf(hy));   //possible minimum 
	  im_a |= im_fun(Inf(hx),Sup(hy));  //possible minimum 
	  im_a |= im_fun(Sup(hx),hy);       //maximum 
	}
    }
#endif
}

void right_side(const interval& hx,const interval& hy,interval& re_a,interval& im_a) 
{
  //Use transp(arctan(z))=arctan(transp(z)) 
  
  left_side(-hx,hy,re_a,im_a);                       
  re_a = -re_a;
}

void harctan(const interval& x,const interval& y,interval& re_arct,interval& im_arct) 
{
  interval h_re(0.0),h_im(0.0);
#ifdef FILIB_VERSION  
  if( inf(x) >= 0 ) //1. quadrant 
    right_side(x,y,re_arct,im_arct);
  else
    {
      if( sup(x) <= 0 ) //2. quadrant 
	left_side(x,y,re_arct,im_arct);
      else
	{
	  if( sup(y) < 1.0 )
	    {
	      re_arct = interval( inf(re_fun(inf(x),sup(y))),sup(re_fun(sup(x),sup(y))) );
	      im_arct = interval( inf(im_fun(std::max(sup(x),-inf(x)),inf(y))),sup(im_fun(0.0,sup(y))) );
	    }
	  else //Intersection in argument 
	    {
	      right_side(interval(0,sup(x)),y,re_arct,im_arct);
	      left_side(interval(inf(x),0),y,h_re,h_im);
	      re_arct = hull( re_arct, h_re + PI() );
	      im_arct = hull( im_arct, h_im                  ); 
	    }
	}
    }
#else //C-XSC Version
  if( Inf(x) >= 0 ) //1. quadrant 
    right_side(x,y,re_arct,im_arct);
  else
    {
      if( Sup(x) <= 0 ) //2. quadrant 
	left_side(x,y,re_arct,im_arct);
      else
	{
	  if( Sup(y) < 1.0 )
	    {
	      re_arct = interval( Inf(re_fun(Inf(x),Sup(y))),Sup(re_fun(Sup(x),Sup(y))) );
	      im_arct = interval( Inf(im_fun(max(Sup(x),-Inf(x)),Inf(y))),Sup(im_fun(0.0,Sup(y))) );
	    }
	  else //Intersection in argument 
	    {
	      right_side(interval(0,Sup(x)),y,re_arct,im_arct);
	      left_side(interval(Inf(x),0),y,h_re,h_im);
	      re_arct |= ( h_re + PI() );
	      im_arct |= h_im; 
	    }
	}
    }
#endif
}

void h_arctan(const interval& x,const interval& y,interval& re_arct,interval& im_arct,bool& singular) 
{
#ifdef FILIB_VERSION
  interval imarct(0.0),rearct(0.0);
  
  if( (inf(x) <= 0)&&(sup(x) >= 0)&&( ((inf(y) <= -1)&&(sup(y) >= -1))||((inf(y) <= 1)&&(sup(y) >= 1)) ) )
    singular = true;
  else
    {
      singular = false;
      if( y==ZERO_INTERVAL() ) 
	{
	  re_arct = atan(x);
	  im_arct = y;  //Im(arctan) = 0 
	}
      else
	{
	  if( (x == ZERO_INTERVAL())&&( inf(y) >- 1 )&&( sup(y) < 1 ) ) 
	    {
	      re_arct = x; //Re(arctan) = 0 
	      im_arct = atanh(y); //Only defined for  -1 < y < 1   
	    }
	  else
	    {
	      //Verwende Arctan(-z) = - Arctan(z) 
	      if( inf(y) >= 0 )
		harctan(x,y,re_arct,im_arct);
	      else  
		{
		  if( sup(y) <= 0 )
		    {
		      harctan(-x,-y,re_arct,im_arct); 
		      re_arct = - re_arct;
		      im_arct = - im_arct;
		    }
		  else
		    {
		      harctan(-x,interval(0,-inf(y)),rearct,imarct); 
		      harctan(x,interval(0,sup(y)),re_arct,im_arct);
		      re_arct = hull( re_arct, -rearct );
		      im_arct = hull( im_arct, -imarct );
		    } 
		}
	    }
	}
    }
#else //C-XSC Version
  cinterval z( ZERO_INTERVAL(), ZERO_INTERVAL() );
  interval imarct(0.0),rearct(0.0);
  
  if( (Inf(x) <= 0)&&(Sup(x) >= 0)&&( ((Inf(y) <= -1)&&(Sup(y) >= -1))||((Inf(y) <= 1)&&(Sup(y) >= 1)) ) )
    singular = true;
  else
    {
      singular = false;
      if( y==ZERO_INTERVAL() ) 
	{
	  re_arct = atan(x);
	  im_arct = y;  //Im(arctan) = 0 
	}
      else
	{
	  if( (x == ZERO_INTERVAL())&&( Inf(y) >- 1 )&&( Sup(y) < 1 ) ) 
	    {
	      re_arct = x; //Re(arctan) = 0 
	      im_arct = atanh(y); //Only defined for  -1 < y < 1   
	    }
	  else
	    {
	      //Verwende Arctan(-z) = - Arctan(z) 
	      if( Inf(y) >= 0 )
		harctan(x,y,re_arct,im_arct);
	      else  
		{
		  if( Sup(y) <= 0 )
		    {
		      harctan(-x,-y,re_arct,im_arct); 
		      re_arct = - re_arct;
		      im_arct = - im_arct;
		    }
		  else
		    {
		      harctan(-x,interval(0,-Inf(y)),rearct,imarct); 
		      harctan(x,interval(0,Sup(y)),re_arct,im_arct);
		      re_arct |= -rearct;
		      im_arct |= -imarct;
		    } 
		}
	    }
	}
    }
#endif
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*     Computation of  arccot, artanh and arcoth  based on arctan            */
/*                                                                           */
/*     arccot(z) = arctan(-z) +/-  pi/2 ,  z <> +/- i                        */
/*     artanh(z) = -i * arctan( i * z ) ,  z <> +/- 1                        */
/*     arcoth(z) = i * arccot( i * z )  ,  --- '' ---                        */
/*                                                                           */
/* ------------------------------------------------------------------------- */

cinterval atan(const cinterval& z) 
{
#ifdef FILIB_VERSION
  if(  (inf(Re(z)) <= 0.0) && (0.0 <= sup(Re(z))) &&   //Zero in real part and
       ((inf(Im(z)) < -1.0) || (sup(Im(z)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();
#else //C-XSC Version
  if(  (Inf(Re(z)) <= 0.0) && (0.0 <= Sup(Re(z))) &&   //Zero in real part and
       ((Inf(Im(z)) < -1.0) || (Sup(Im(z)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();
#endif
  interval im_arct(0.0),re_arct(0.0);
  bool singular( false );
  
  h_arctan(Re(z),Im(z),re_arct,im_arct,singular);

  return cinterval(re_arct,im_arct);
}

void h_arccot(const interval& x,const interval& y,interval& re_arcc,interval& im_arcc,bool& singular) 
{
  h_arctan(x,y,re_arcc,im_arcc,singular);

  if( !singular )
#ifdef FILIB_VERSION
    if( inf(re_arcc) >= 0.0 )
#else //C-XSC Version
    if( Inf(re_arcc) >= 0.0 )
#endif

      re_arcc = -re_arcc - HALFPI();
    else
      re_arcc =  HALFPI() - re_arcc;
}

cinterval acot(const cinterval& z) 
{
#ifdef FILIB_VERSION
  if(  (inf(Re(z)) <= 0.0) && (0.0 <= sup(Re(z))) &&   //Zero in real part and
       ((inf(Im(z)) < -1.0) || (sup(Im(z)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();
#else //C-XSC Version
  if(  (Inf(Re(z)) <= 0.0) && (0.0 <= Sup(Re(z))) &&   //Zero in real part and
       ((Inf(Im(z)) < -1.0) || (Sup(Im(z)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();
#endif
  interval im_arcc(0.0),re_arcc(0.0);
  bool singular( false );

  h_arccot(Re(z),Im(z),re_arcc,im_arcc,singular);
  
  return cinterval(re_arcc,im_arcc);
}

cinterval atanh(const cinterval& z) 
{
#ifdef FILIB_VERSION
  if(  (inf(Im(z)) <= 0.0) && (0.0 <= sup(Im(z))) &&   //Zero in imaginary part and
       ((inf(Re(z)) < -1.0) || (sup(Re(z)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();
#else //C-XSC Version
  if(  (Inf(Im(z)) <= 0.0) && (0.0 <= Sup(Im(z))) &&   //Zero in imaginary part and
       ((Inf(Re(z)) < -1.0) || (Sup(Re(z)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();
#endif
  interval im_art(0.0),m_re_art(0.0);
  bool singular( false );
  
  h_arctan(Im(z),-Re(z),im_art,m_re_art,singular);
  
  return cinterval(-m_re_art,im_art); //atanh(z) = -i * atan(i*z)
}

cinterval acoth(const cinterval& z) 
{
#ifdef FILIB_VERSION
  if(  (inf(Im(z)) <= 0.0) && (0.0 <= sup(Im(z))) &&   //Zero in imaginary part and
       (sup(Re(z)) > -1.0) && (inf(Re(z)) < 1.0 )   )  //real part intersects (-1,1)
    throw function_not_defined();
#else //C-XSC Version
  if(  (Inf(Im(z)) <= 0.0) && (0.0 <= Sup(Im(z))) &&   //Zero in imaginary part and
       (Sup(Re(z)) > -1.0) && (Inf(Re(z)) < 1.0 )   )  //real part intersects (-1,1)
    throw function_not_defined();
#endif
  interval im_arco(0.0),m_re_arco(0.0);
  bool singular( false );

  h_arccot(-Im(z),Re(z),im_arco,m_re_arco,singular); 

  if( singular ) throw function_not_defined();//Singularity of arcoth in argument.
  else return cinterval(-m_re_arco,im_arco); 
}

/* ------------------------------------------------------------------------- */
/*                                                                           */
/*     Computation of power function based on logarithm                      */
/*                                                                           */
/* ------------------------------------------------------------------------- */

cinterval power(const cinterval& z,int n) 
{
  if( n < 0 ) return 1.0/power( z , -n );
  else if( n == 0 ) return cinterval( ONE_INTERVAL(), ZERO_INTERVAL() );
  else if( n == 1 ) return z;
  else if( n == 2 ) return sqr( z );
  else
    {
#ifdef FILIB_VERSION
       if( (inf(Re(z)) > 0)||(inf(Im(z)) > 0)||(sup(Im(z)) < 0) ) 
	 return exp( n * log(z) );
#else  //C-XSC Version
       if( (Inf(Re(z)) > 0)||(Inf(Im(z)) > 0)||(Sup(Im(z)) < 0) ) 
	 return exp( n * ln(z) );
#endif
      else
	{//Zero in argument
	  cinterval w( z ),u(
#ifdef FILIB_VERSION
			     ONE_INTERVAL(), ZERO_INTERVAL()
#else //C-XSC Version
	                     1.0
#endif
			     );
	  while( n > 0 )
	    {
	      if( (n & 1) == 0 )//even exponent
		{
		  w = sqr(w);
		  n >>= 1;//Division by 2, realized with a right shift
		}
	      else
		{
		  u *= w;
		  n -= 1;
		}
	    }
	  return u;
	}
    }
}

cinterval pow(const cinterval& z,const interval& n) 
{
  interval re_z(Re(z)), im_z(Im(z));

#ifdef FILIB_VERSION
  if( ( inf(re_z) > 0 )||( inf(im_z) > 0 )||( sup(im_z) < 0 ) ) 
    return exp(n*log(z));
#else //C-XSC Version
  if( ( Inf(re_z) > 0 )||( Inf(im_z) > 0 )||( Sup(im_z) < 0 ) ) 
    return exp(n*ln(z));
#endif

  else throw function_not_defined();//Z contains zero.
}

std::list<cinterval> pow_all(const cinterval& z, const interval& n)
{
  interval abs_z(abs(z));

  if( 
#ifdef FILIB_VERSION
     inf(abs_z)
#else
     Inf(abs_z)
#endif
     > 0 )
    {
      interval sqrt_2( sqrt(interval(2.0)) ), dist_1(0.0), dist_2(0.0);
      interval exp_n_ln_abs_z(
#ifdef FILIB_VERSION
			      exp(n * log(abs_z))
#else
			      exp(n * ln (abs_z))
#endif
			      );
      dist_1 = 0.5 * sqrt_2 * 
#ifdef FILIB_VERSION
	inf(exp_n_ln_abs_z)
#else
	Inf(exp_n_ln_abs_z)
#endif
	;
      dist_2 = 
#ifdef FILIB_VERSION
	sup(exp_n_ln_abs_z)
#else
	Sup(exp_n_ln_abs_z)
#endif
	;

      std::list<cinterval> res;
      
#ifdef FILIB_VERSION
      res.push_back( cinterval( interval(  inf(dist_1),  sup(dist_2) ),
				interval( -sup(dist_1),  sup(dist_2) ) ) );
      res.push_back( cinterval( interval( -sup(dist_2),  sup(dist_1) ),
				interval(  inf(dist_1),  sup(dist_2) ) ) );
      res.push_back( cinterval( interval( -sup(dist_2), -inf(dist_1) ),
				interval( -sup(dist_2),  sup(dist_1) ) ) );
      res.push_back( cinterval( interval( -sup(dist_1),  sup(dist_2) ),
				interval( -sup(dist_2), -inf(dist_1) ) ) );
#else
      res.push_back( cinterval( interval(  Inf(dist_1),  Sup(dist_2) ), 
				interval( -Sup(dist_1),  Sup(dist_2) ) ) );
      res.push_back( cinterval( interval( -Sup(dist_2),  Sup(dist_1) ), 
				interval(  Inf(dist_1),  Sup(dist_2) ) ) );
      res.push_back( cinterval( interval( -Sup(dist_2), -Inf(dist_1) ), 
				interval( -Sup(dist_2),  Sup(dist_1) ) ) );
      res.push_back( cinterval( interval( -Sup(dist_1),  Sup(dist_2) ), 
				interval( -Sup(dist_2), -Inf(dist_1) ) ) );
#endif
  
      return res;
    }
  else //z contains zero
    {
      if(
#ifdef FILIB_VERSION
	 inf(n)
#else
	 Inf(n)
#endif
	 < 0 ) //return entire complex plane
	{
	  std::list<cinterval> res;
	  
#ifdef FILIB_VERSION
	  res.push_back( cinterval( interval( -Double::MAX(),Double::MAX() ),
				    interval( -Double::MAX(),Double::MAX() ) ) ); 
#else
	  res.push_back( cinterval( interval( -MaxReal,MaxReal ),
				    interval( -MaxReal,MaxReal ) ) );
#endif
	  return res;
	}
      else if(
#ifdef FILIB_VERSION
	      inf(n) == 0 && sup(n) == 0
#else
	      Inf(n) == 0 && Sup(n) == 0
#endif
	      )
	{
	  std::list<cinterval> res;
	  
	  res.push_back( cinterval( interval( 0.0, 1.0 ),
				    interval( ZERO_INTERVAL()   ) ) ); 

	  return res;
	}
      else
	{
	  interval exp_n_ln_sup_abs_z(
#ifdef FILIB_VERSION
				      exp(n * log(interval(sup(abs_z))))
#else
				      exp(n * ln (interval(Sup(abs_z))))
#endif
				      );
	  flnumber d_2 = 
#ifdef FILIB_VERSION
	    sup(exp_n_ln_sup_abs_z)
#else
	    Sup(exp_n_ln_sup_abs_z)
#endif
	    ;

	  std::list<cinterval> res;
	  
	  res.push_back( cinterval( interval( 0.0, d_2 ),
				    interval( 0.0, d_2 ) ) ); 
	  
	  return res;
	}
    }
}

cinterval pow(const cinterval& z,const cinterval& n) 
{
  interval re_z(Re(z)), im_z(Im(z));

#ifdef FILIB_VERSION
  if( ( inf(re_z) > 0 )||( inf(im_z) > 0 )||( sup(re_z) < 0 )||( sup(im_z) < 0 ) ) 
    return exp(n*log(z));
#else //C-XSC Version
  if( ( Inf(re_z) > 0 )||( Inf(im_z) > 0 )||( Sup(re_z) < 0 )||( Sup(im_z) < 0 ) ) 
    return exp(n*ln(z));
#endif

  else throw function_not_defined();//Base contains zero.
}

std::list<cinterval> pow_all(const cinterval& z, const cinterval& n)
{
  if( Im(n) == ZERO_INTERVAL() ) return pow_all( z, Re(n) );
  else 
    { // nonreal eponent: return entire complex plane
      std::list<cinterval> res;
      
#ifdef FILIB_VERSION
       res.push_back( cinterval( interval( -Double::MAX(),Double::MAX() ),
                                 interval( -Double::MAX(),Double::MAX() ) ) ); 
#else
       res.push_back( cinterval( interval( -MaxReal,MaxReal ),
                                 interval( -MaxReal,MaxReal ) ) );
#endif

      return res;
    }
}

/*

  End of File: cimath.cpp

*/
