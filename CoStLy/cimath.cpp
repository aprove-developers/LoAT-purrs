/*

 File: cimath.cpp, 2002/03/19

 CoStLy (COmplex interval STandard functions LibrarY), Version 0.2

 Copyright (C) Markus Neher, markus.neher@math.uni-karlsruhe.de
               Ingo Eble,    ingoeble@web.de

 This implementation of some complex standard functions for interval arithmetic
 was taken from a Pascal-XSC modul written by Andreas Westphal.
 Implemetation ...

 ... in C++ with 

            filib++:	     Copyright (C) 2002 Markus Neher 
                                                Ingo Eble

	    C-XSC:	     Copyright (C) 2000 Markus Neher 
                                                Ingo Eble

 ... in Pascal-XSC:          Copyright (C) 1999 Andreas Westphal
                                                Walter Kraemer

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

 
 Original header of the Pascal-xsc file written by Andreas Westphal:

 ************************************************************************** 
 ************************************************************************** 
 ***                                                                    *** 
 ***        Weitere Standardfunktionen fuer komplexe Intervalle         *** 
 ***                                                                    *** 
 ************************************************************************** 
 ************************************************************************** 
 ************************************************************************** 
 ***  Literatur :                                                       *** 
 ***                                                                    *** 
 ***   Braune, K.D. : 'Hochgenaue Standardfunktionen fuer reelle und    *** 
 ***     und komplexe Punkte und Intervalle in beliebigen Gleitpunkt-   *** 
 ***     rastern', Dissertation, Karlsruhe, 1987.                       *** 
 ***   Buehler, G. : 'Standardfunktionen fuer komplexe Intervalle im    *** 
 ***     64 Bit IEEE Datenformat', Diplomarbeit, Karlsruhe, 1993.       *** 
 ***   Kraemer, W. : 'Inverse Standardfunktionen fuer reelle und        *** 
 ***     komplexe Intervallargumente mit a priori Fehlerabschaetzungen  *** 
 ***     fuer beliebige Datenformate', Dissertation, Karlsruhe, 1987.   *** 
 ***                                                                    *** 
 ************************************************************************** 
 ***  Autor: Andreas Westphal                                           *** 
 ***  Betreuer: Walter Kraemer                                          ***  
 ***     Institut fuer Wissenschaftliches Rechnen  und                  ***   
 ***     Mathematische Modellbildung, Universitaet Karlsruhe (TH)       *** 
 ***                       Stand: Maerz 1999                            *** 
 ************************************************************************** 

*/


//Include header files
#include "cimath.h" //Declaration of the complex functions
#include "error.h"

#include "Interval.h" //fi_lib++ Header: Macro Version

typedef Interval interval;

inline const interval& HALFPI()
{ 
  static const interval hp( interval::PI() / 2.0 );
  return hp;
}

inline const interval& QUARTERPI()
{
  static const interval qp( interval::PI() / 4.0 );
  return qp;
}

/* ***************************************************************************/
/* ***************************************************************************/
/* ***                      Single-valued functions                      *** */
/* ***************************************************************************/
/* ***************************************************************************/


/* ***************************************************************************/
/* *** Power operator  pow  is not listed here, since it relies on the    ****/
/* *** (multi-valued) logarithm                                           ****/
/* ***************************************************************************/
 

/* ***************************************************************************/
/* *** The hyperbolic functions exp, sin, cos, sinh, cosh are separable:  ****/
/* *** Their real and imaginary parts are products of real functions      ****/
/* ***************************************************************************/
/* ***   With Re(z)=x, Im(z)=y :                                          ****/
/* ***                                                                    ****/
/* ***        exp   :   Re(exp(z)) = exp(x) * cos(y)                      ****/
/* ***                  Im(exp(z)) = exp(x) * sin(y)                      ****/
/* ***                                                                    ****/
/* ***        sin   :   Re(sin(z)) = sin(x) * cosh(y)                     ****/
/* ***                  Im(sin(x)) = cos(x) * sinh(y)                     ****/
/* ***                                                                    ****/
/* ***        cos   :   Re(cos(z)) = cos(x) * cosh(y)                     ****/
/* ***                  Im(sin(x)) = -sin(x) * sinh(y)                    ****/
/* ***                                                                    ****/
/* ***        sinh  :   Re(sinh(z)) = sinh(x) * cos(y)                    ****/
/* ***                  Im(sinh(z)) = cosh(x) * sin(y)                    ****/
/* ***                                                                    ****/
/* ***        cosh  :   Re(cosh(z)) = cosh(x) * cos(y)                    ****/
/* ***                  Im(cosh(z)) = sinh(x) * sin(y)                    ****/
/* ***                                                                    ****/
/* ***************************************************************************/
 
cinterval exp(const cinterval& z) 
{
  interval 
    A( exp( z.re() ) ), 
    B(      z.im()   );
  return cinterval( A*cos( B ) , A*sin( B ) );
}

cinterval cos(const cinterval& z) 
{
  interval 
    A( z.re() ), 
    B( z.im() );
  return cinterval( cos( A )*cosh( B ) , -sin( A )*sinh( B ) );
}

cinterval sin(const cinterval& z) 
{
  interval 
    A( z.re() ), 
    B( z.im() );
  return cinterval( sin( A )*cosh( B ) , cos( A )*sinh( B ) );
} 
    
cinterval cosh(const cinterval& z) 
{
  interval 
    A( z.re() ), 
    B( z.im() );
  return cinterval( cos( B )*cosh( A ) , sin( B )*sinh( A ) );
}

cinterval sinh(const cinterval& z) 
{
  interval 
    A( z.re() ), 
    B( z.im() );
  return cinterval( cos( B )*sinh( A ) , sin( B )*cosh( A ) );
}

/* *************************************************************/
/* *** Tangent is NOT separable, naive evaluation may yield ****/
/* *** range  overestimation.                               ****/
/* ***                                                      ****/
/* *** Accurate range evaluation requires discussion of     ****/
/* *** several cases and uses the symmetry condition        ****/
/* ***                                                      ****/
/* ***      tan(z) = transp( tan( transp(z) ) )             ****/
/* ***                                                      ****/
/* *************************************************************/

void ReAdd(const double& hx,const double& hy,interval& re_tan,bool& re_first,bool both) 
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
	  re_tan = hull( re_tan, sin( HX ) * cosHX /( sqr( cosHX ) + sqr( sinh( HY )) ) );
	}
    }
}

void ReAdd(const interval& hx,const double& hy,interval& re_tan,bool& re_first,bool both) 
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
	  re_tan = hull( re_tan, sin( hx ) * cosHX /( sqr( cosHX ) + sqr( sinh( interval(hy) )) ) );
	}
    }
}

void ImAdd(const double& hx,const double& hy,interval& im_tan,bool& im_first) 
{
  interval HX( hx ),HY( hy ),sinhHY( sinh( HY ) );
  if( im_first ) 
    { //Naive evaluation of imaginary part
      im_tan = sinhHY * cosh( HY ) /( sqr( cos( HX )) + sqr( sinhHY ) );
      im_first = false;
    }
  else 
    { //Interval hull of former calculations
      im_tan = hull( im_tan, sinhHY * cosh( HY ) /( sqr( cos( HX )) + sqr( sinhHY ) ) );
    }
}

void ImAdd(const interval& hy,const double& hx,interval& im_tan,bool& im_first) 
{
  interval sinhHY( sinh( hy ) );
  if( im_first ) 
    { //Naive evaluation of imaginary part
      im_tan = sinhHY * cosh( hy ) /( sqr( cos( interval( hx ) )) + sqr( sinhHY ) );
      im_first = false;
    }
  else 
    { //Interval hull of former calculations
      im_tan = hull( im_tan, sinhHY * cosh( hy ) /( sqr( cos( interval( hx ) )) + sqr( sinhHY ) ) );
    }
}

void htan(const interval& x,const interval& y,interval& re_tan,interval& im_tan,bool both) 
{
  /* If BOTH=TRUE, then enclosures of  Re(tan(x+i*y)) and                */
  /* Im(tan(x+i*y)) are calculated; otherwise only Im_tan is calculated. */
  /* y > 0 is assumed here (upper half plane).                           */

  bool Re_first( true ),
       Im_first( true );
  	
  /* Real part of arguments from left to right */
  
  interval hint(0.0),argx(0.0),argy(0.0);
  
  if( inf(x) < (-inf(QUARTERPI())) ) 
    { //Real part of argument intersects [-pi/2,-pi/4] 

      if( (-inf(QUARTERPI())) < sup(x) ) hint = interval( inf(x), -inf(QUARTERPI()) );
      else                            hint = x; //hint = Intersection of x and [-pi/2,-pi/4] 
    
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
}


void h_tan(const interval& x,const interval& y,interval& re_tan,interval& im_tan,int& Error) 
{
  interval hx(0.0),hx2(0.0);
  interval Imtan(0.0),Retan(0.0);
  bool BOTH( false );
  bool bisection( false ); // Pi/2 mod Pi  in x ? 
  
  Error = 0;
  if ( x == interval::ZERO() )  
    {
      re_tan = x;               //Re(tan) = 0
      im_tan = tanh(y);
    }
  else
    {
      if( (!disjoint(interval::ZERO(),y)) && (!disjoint(interval::ZERO(),cos(x))) ) Error = 1;
      else
	{
	  if ( y == interval::ZERO() )  
	    {
	      re_tan = tan(x);
	      im_tan = y;         //Im(tan) = 0
	    }
	  else
	    { 
	      //Since  z equivalent (z mod pi), use  
	      //equivalent argument in [ -pi/2 , pi/2 ] .

	      double r_int;

	      if( inf(x) > 0 ) modf( inf(x/interval::PI()) + 0.5, &r_int );
	      else             modf( inf(x/interval::PI()) - 0.5, &r_int );
	      
	      hx = x - r_int * interval::PI();
	      

	      if( inf(interval::PI()) < (2*sup(hx)) )  
		{
		  bisection = true;
		  hx2 = interval( -sup(HALFPI()), sup(hx)-inf(HALFPI()) );
		  hx  = interval(  inf(hx)    , sup(HALFPI())         );
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
}
	
/* ****************************************************************/
/* ***  Computation of cot, tanh, coth based on tan.           ****/
/* ****************************************************************/
/* ****************************************************************/
/* ***  cot(z)  = tan( pi/2 - z )                              ****/
/* ***  tanh(z) = transp( i * tan( transp( i * z ) ) )         ****/
/* ***  coth(z) = i * cot( i * z )                             ****/
/* ***          = i * tan( pi/2 - i * z )                      ****/
/* ****************************************************************/

void h_cot(const interval& x,const interval& y,interval& re_cot,interval& im_cot,int& Error) 
{
  interval hx(0.0),hx2(0.0);
  interval Imcot(0.0),Recot(0.0);
  bool BOTH( false );
  bool bisection( false ); //0 mod Pi  in x ? 
  //If yes: bisect input interval
  
  Error =  0;
  if( (!disjoint(interval::ZERO(),y))&&(!disjoint(interval::ZERO(),sin(x))) ) 
    {
      Error = 1;
    }
  else
    {
      if( y == interval::ZERO() ) 
	{
	  re_cot = cot(x);
	  im_cot = y;               //Im(cot) = 0 
	}
      else
	{
	  if( x == interval::ZERO() ) 
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
	      if( (sup(x)-inf(x)) > inf(interval::PI()) ) 
		{
		  hx = interval( -sup(HALFPI()), sup(HALFPI()) );
		}
	      else
		{
		  /* z equivalent (z mod pi), use               */
		  /* equivalent argument in [ -pi/2 , pi/2 ] .  */
		  
		  double r_int;
		  
		  if( inf(x) < 0 ) 
		    {
		      modf( -inf(x/interval::PI()), &r_int );
		      hx = -x + ( -1 - 2 * r_int ) * HALFPI(); 
		    }
		  else
		    {
		      modf(  inf(x/interval::PI()), &r_int );
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

/* ***************************************************************************/


/* ***************************************************************************/
/* ***************************************************************************/
/* ***                      Multi-valued functions                        ****/
/* ***************************************************************************/
/* ***************************************************************************/

interval arg(const cinterval& z) 
{
  interval wert(0.0); 
  interval hxl(0.0),hyl(0.0);
  interval hxu(0.0),hyu(0.0);
  
  if( interval::ZERO() == Re(z) ) 
    {
      if( interval::ZERO() == Im(z) ) throw function_not_defined();//arg not defined at (0,0).
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
      if( interval::ZERO() <= Re(z) ) 
	{
	  if( interval::ZERO() <= Im(z) )            //all quadrants 
	    return 2.0 * interval(-sup(HALFPI()),sup(HALFPI()));
	  else 
	    {
	      if( inf(Im(z)) >= 0 )                //1. and 2. quadrant 
		return interval(inf(atan(hyl/hxu)),sup(atan(hyl/hxl)+interval::PI()));//Minimum, maximum 
	      else                                   //3. and 4. quadrant 
		return interval(inf(atan(hyu/hxl)-interval::PI()),sup(atan(hyu/hxu)));//Minimum, maximum 
	    }
	}
      else
	{ 
	  if( interval::ZERO() <= Im(z) ) 
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
		    return interval(inf(atan(hyu/hxu)),sup(atan(hyl/hxu))) + interval::PI(); //Minimum, maximum 
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
			return interval(inf(HALFPI()),sup(atan(hyl/hxl)+interval::PI())); //Maximum 
		      else
			return interval(inf(atan(hyu/hxu)),sup(atan(hyl/hxl))) + interval::PI(); //Minimum, maximum 
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
			return interval(-sup(HALFPI()),sup(atan(hxu/hyl)-interval::PI())); //Minimum, maximum 
		      else
			return interval(inf(atan(hyu/hxl)),sup(atan(hyl/hxu))) - interval::PI(); //Minimum, maximum 
		    }
		}
	    }
	}
    }
}  

inline interval abs(const interval& x,const interval& y) 
{
  return sqrt(sqr(x)+sqr(y)); 
} 


cinterval ln(const cinterval& z) 
{
  double re_min(0.0),im_min(0.0);
  double re_max(0.0),im_max(0.0);
  interval re_ln(0.0);
  
  double srez( sup( Re(z) ) ),irez( inf( Re(z) ) ), simz( sup( Im(z) ) ),iimz( inf( Im(z) ) );

  if( srez < 0 )  re_min = -srez;
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
    if( (srez <= 0)&&(im_min==0) )//Argument intersects negative real axis.
      throw function_not_defined();
    else
      { 
	re_max = max( srez, -irez );
	im_max = max( simz, -iimz );
	re_ln = interval( inf( log( abs(interval(re_min),interval(im_min)) ) ),
			  sup( log( abs(interval(re_max),interval(im_max)) ) )  );
	return cinterval(re_ln,arg(z));
      }
}

inline interval sign(const interval& x) 
{
  double ix = inf(x), sx = sup(x);
  double inf, sup;

  if( ix > 0 ) inf = 1;
  else 
    if( ix ) inf = -1;
    else     inf =  0;

  if( sx > 0 ) sup = 1;
  else 
    if( sx ) sup = -1;
    else     sup =  0;

  return interval( inf, sup ); 
}

interval re_sqrt(const interval& x,const interval& y) 
{
  //Formulas for special quadrants 

  if( sup(x) < 0 )
    {
      if( y == interval::ZERO() ) return interval::ZERO();
      else return 1.0/sqrt(interval(2.0))*abs(y)/sqrt(abs(x,y)-x);
    }
  else 
    {
      if( y == interval::ZERO() ) return sqrt(x);
      else return 1.0/sqrt(interval(2.0))*sqrt(abs(x,y)+x);
    }
}

interval im_sqrt(const interval& x,const interval& y) 
{
  //Formulas for special quadrants 

  if( sup(x) < 0 )
    {
      if( y == interval::ZERO() ) return 1.0/sqrt(-x);
      else return 1.0/sqrt(interval(2.0))*sign(y)*sqrt(abs(x,y)-x);
    }
  else
    {
      if( y == interval::ZERO() ) return interval::ZERO(); 
      else return 1.0/sqrt(interval(2.0))*y/sqrt(abs(x,y)+x);
    }
}

cinterval sqrt(const cinterval& z) 
{
  interval rwert(0.0),iwert(0.0);
  interval hxl( inf(Re(z)) ),hyl( inf(Im(z)) );
  interval hxu( sup(Re(z)) ),hyu( sup(Im(z)) );
  
  if( (sup(hxl) < 0.0) && (sup(hyl) <= 0.0) && (inf(hyu) >= 0.0) ) 
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
	  rwert = interval(inf(rwert),max(sup(re_sqrt(hxu,hyu)),sup(re_sqrt(hxu,hyl))));
	  iwert = im_sqrt(hxl,hyl); //Minimum 
	  iwert = interval(inf(iwert),sup(im_sqrt(hxl,hyu))); //with maximum 
	}
    }
  return cinterval(rwert,iwert);
}

cinterval sqrt(const cinterval& z,int n) 
{
  if( n==0 ) return cinterval( interval::ONE(), interval::ZERO() );
  if( n==1 ) return z;
  if( n==2 ) return sqrt( z );
  if( n>=3 ) return pow( z, interval( 1.0/n ) );
}

cinterval sqr(const cinterval& z) 
{
  interval 
    A( Re(z) ),
    B( Im(z) );
  return cinterval( sqr( A ) - sqr( B ) , 2.0 * A * B );
}


/* ********************************************************/
/* *** arcsin from diploma thesis of Gabriele Buehler. ****/
/* ********************************************************/

inline interval g(const interval& s) 
{
  return sqrt( 1.0 + s ) - 1.0;
}

inline interval s_re(const interval& ix,const interval& iy) 
{
  return 2.0 * iy / ( ix*ix + iy*iy - 1.0 );
}

interval re_arcsin(const double& x,const double& y) 
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
	      if( sup(hilf1) < 1.0 )
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

interval im_arcsin(const double& x,const double& y) 
{
  /* Interval computation of imaginary part of arcsin(z) */
  interval 
    hilf1(0.0),
    hilf2(0.0),
    hilf3(0.0),
    hilf4(0.0),
    hilf5(0.0),
    hilf6(0.0),
    nenner(0.0),
    zaehler(0.0),
    ix( abs(interval(x)) ),
    iy( y );
  interval t(0.0),r(0.0);
  double xc( abs(x) );
  
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
	      return log(t + sqrt(sqr(t) - 1.0));
	    }
	  else                          //1< x < 1.1 
	    {
	      t = ix;
	      r = t - 1.0;
	      return log(1.0 + (r + sqrt(sqr(r) + 2.0 * r)));
	    }
	}
    }
  else                       //y <> 0.0 
    {
      hilf1 = ix*ix + iy*iy;
      hilf2 = hilf1 + 1.0;
      if( xc == 0.0 )         //x= 0.0 
        return log(sqrt(1.0 + iy * iy) + iy);
      else                   //x <> 0.0 
        {
          t  = 0.5 * ( sqrt( hilf2 + 2.0 * ix) + sqrt( hilf2 - 2.0 * ix ));
          if( sup(t) <= 1.1 )   //t <= 1.1 
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
              return log(1.0 + (r + sqrt(sqr(r) + 2.0 * r)));
            }
          else                    //t > 1.1 
            return log(t + sqrt(sqr(t) - 1.0));
        }
    }
}

inline interval fortsetz_asin(const double& x,const double& y) 
{
  interval hilf( re_arcsin(x,y) );

  if( hilf == interval::ZERO() ) 
    return interval::PI();
  else
    return interval::PI() - hilf;
}

inline void z(const cinterval& z,double& re_l,double& re_u,double& im_l,double& im_u) 
{
  re_l = inf(Re(z));
  re_u = sup(Re(z));
  im_l = inf(Im(z));
  im_u = sup(Im(z));
}

inline void z(const interval& x,const interval& y,double& x_l,double& x_u,double& y_l,double& y_u) 
{
  x_l = inf(x);
  x_u = sup(x);
  y_l = inf(y);
  y_u = sup(y);
}

interval real_asin(const cinterval& c) 
{
  double xl(0.0),xu(0.0),yl(0.0),yu(0.0),maxy(0.0),max(0.0);
  bool re_spiegel( false );
  double null( 0.0 ),eins( 1.0 );
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
  
  ergx = interval(inf(ergxl),sup(ergxu));
  if( re_spiegel ) 
    return -ergx;
  else
    return ergx;
}

interval imag_asin(const cinterval& c) 
{
  double xl(0.0),xu(0.0),yl(0.0),yu(0.0),maxx(0.0),maxy(0.0),max(0.0);
  bool im_spiegel( false );
  double null( 0.0 ),eins( 1.0 );
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
  
  ergy = interval(inf(ergyl),sup(ergyu));
  
  if( im_spiegel )
    return -ergy;
  else
    return ergy;
}


/* ********************************************************************/
/* *** Computation of  arccos, arsinh and arcosh  based on  arcsin ****/
/* ********************************************************************/
/* ***                                                             ****/
/* ***  Arccos(z) = -/+ ( Arcsin(z) - pi/2 )                       ****/
/* ***  Arsinh(z) = i * Arcsin( -i * z ) ( mod i*2*pi )            ****/
/* ***  Arcosh(z) = -/+ i * Arccos(z)                              ****/
/* ***                                                             ****/
/* ***  Only principal values are computed.                        ****/
/* ********************************************************************/

cinterval asin(const cinterval& c) 
{
  if(  (inf(Im(c)) <= 0.0) && (0.0 <= sup(Im(c))) &&   //Zero in imaginary part and
       ((inf(Re(c)) < -1.0) || (sup(Re(c)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();

  return cinterval(real_asin(c),imag_asin(c));
}

cinterval acos(const cinterval& c) 
{
  if(  (inf(Im(c)) <= 0.0) && (0.0 <= sup(Im(c))) &&   //Zero in imaginary part and
       ((inf(Re(c)) < -1.0) || (sup(Re(c)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();
  
  return cinterval(real_asin(c)-HALFPI(),imag_asin(c));
}

cinterval asinh(const cinterval& c) 
{
  if(  (inf(Re(c)) <= 0.0) && (0.0 <= sup(Re(c))) &&   //Zero in real part and
       ((inf(Im(c)) < -1.0) || (sup(Im(c)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();

  cinterval hc( -Im(c), Re(c) ); //hc = i*c 
  
  return cinterval(imag_asin(hc),-real_asin(hc)); //arsinh(c) = -i*arcsin(i*c) 
}

cinterval acosh(const cinterval& c) 
{
  if(  (inf(Im(c)) <= 0.0) && (0.0 <= sup(Im(c))) &&   //Zero in imaginary part and
       (inf(Re(c)) < 1.0)                           ) //real part intersects (-INFINITY,1)
    throw function_not_defined();

  cinterval hc(real_asin(c)-HALFPI(),imag_asin(c));
  
  return cinterval(-Im(hc),Re(hc));               //arcosh = i * arccos(c) 
}

/* ***************************************************************************/


/* *************************************************/
/* ***  Auxiliary functions for inverse tangent ****/
/* *************************************************/

interval qbetrag(const double& lx,const double& ux,const double& ly,const double& uy) 
{
  double infqbetr(0.0),supqbetr(0.0);
  double tmp1,tmp2;
  
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
}

interval re_fun(const double& hx,const double& hy) 
{
  interval qbetr( qbetrag(hx,hx,hy,hy) );
  
  if( 1 > sup(qbetr) )
    return atan( interval(2*hx)/(1.0-qbetr) ) / 2.0;
  else
    {
      if( 1 < inf(qbetr) )
	return (atan( interval(2*hx)/(1.0-qbetr) ) - interval::PI() ) / 2.0 ;
      else
	return -QUARTERPI();
    }
}

interval re_fun(const interval& hx,const double& hy) 
{
  interval qbetr( qbetrag(inf(hx),sup(hx),hy,hy) );
  
  if( 1 > sup(qbetr) )
    return atan( 2.0*hx/(1.0-qbetr) ) / 2.0;
  else
    {
      if( 1 < inf(qbetr) ) 
	return (atan( 2.0*hx/(1.0-qbetr) ) - interval::PI() ) / 2.0;
      else 
	return -QUARTERPI();
    }
}

interval im_fun(const double& hx,const double& hy) 
{
  double ly,uy;

  SUB_DOWN( ly, 1.0, hy );
  SUB_UP( uy, 1.0, hy );

  interval nenner( qbetrag(hx,hx,ly,uy) );

  if( 0 < inf(nenner) ) 
    return log( abs(1.0+4.0*interval(hy)/nenner) ) / 4.0;
  else//Argument too near to singularity.
    throw function_not_defined();
}

interval im_fun(const double& hx,const interval& hy) 
{
  double ly,uy;

  SUB_DOWN( ly, 1.0, sup(hy) );
  SUB_UP( uy, 1.0, inf(hy) );

  interval nenner( qbetrag(hx,hx,ly,uy) );

  if( 0 < inf(nenner) )   
    return log(abs(1.0+4.0*hy/(sqr(interval(hx))+sqr(1.0-hy)))) / 4.0;
  else//Argument too near to singularity.
    throw function_not_defined();
}

void left_side(const interval& hx,const interval& hy,interval& re_a,interval& im_a) 
{
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
	      im_arct = interval( inf(im_fun(max(sup(x),-inf(x)),inf(y))),sup(im_fun(0.0,sup(y))) );
	    }
	  else //Intersection in argument 
	    {
	      right_side(interval(0,sup(x)),y,re_arct,im_arct);
	      left_side(interval(inf(x),0),y,h_re,h_im);
	      re_arct = hull( re_arct, h_re + interval::PI() );
	      im_arct = hull( im_arct, h_im                  ); 
	    }
	}
    }
}


void h_arctan(const interval& x,const interval& y,interval& re_arct,interval& im_arct,bool& singular) 
{
  cinterval z( interval::ZERO(), interval::ZERO() );
  interval imarct(0.0),rearct(0.0);
  
  if( (inf(x) <= 0)&&(sup(x) >= 0)&&( ((inf(y) <= -1)&&(sup(y) >= -1))||((inf(y) <= 1)&&(sup(y) >= 1)) ) )
    singular = true;
  else
    {
      singular = false;
      if( y==interval::ZERO() ) 
	{
	  re_arct = atan(x);
	  im_arct = y;  //Im(arctan) = 0 
	}
      else
	{
	  if( (x == interval::ZERO())&&( inf(y) >- 1 )&&( sup(y) < 1 ) ) 
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
}

/* ****************************************************************** */
/* *** Computation of  arccot, artanh and arcoth  based on arctan *** */
/* ****************************************************************** */
/* ***                                                            *** */
/* ***  arccot(z) = arctan(-z) +/-  pi/2 ,  z <> +/- i            *** */
/* ***  artanh(z) = -i * arctan( i * z ) ,  z <> +/- 1            *** */
/* ***  arcoth(z) = i * arccot( i * z )  ,  --- '' ---            *** */
/* ***                                                            *** */
/* ****************************************************************** */

cinterval atan(const cinterval& z) 
{
  if(  (inf(Re(z)) <= 0.0) && (0.0 <= sup(Re(z))) &&   //Zero in real part and
       ((inf(Im(z)) < -1.0) || (sup(Im(z)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();

  interval im_arct(0.0),re_arct(0.0);
  bool singular( false );
  
  h_arctan(Re(z),Im(z),re_arct,im_arct,singular);

  return cinterval(re_arct,im_arct);
}

void h_arccot(const interval& x,const interval& y,interval& re_arcc,interval& im_arcc,bool& singular) 
{
  h_arctan(x,y,re_arcc,im_arcc,singular);

  if( !singular )
    if( inf(re_arcc) >= 0.0 )
      re_arcc = -re_arcc - HALFPI();
    else
      re_arcc = HALFPI() - re_arcc;
}

cinterval acot(const cinterval& z) 
{
  if(  (inf(Re(z)) <= 0.0) && (0.0 <= sup(Re(z))) &&   //Zero in real part and
       ((inf(Im(z)) < -1.0) || (sup(Im(z)) > 1.0 ))  ) //imaginary part intersects i*(-INFINITY,-1) and/or i*(1,+INFINITY)
    throw function_not_defined();

  interval im_arcc(0.0),re_arcc(0.0);
  bool singular( false );

  h_arccot(Re(z),Im(z),re_arcc,im_arcc,singular);
  
  return cinterval(re_arcc,im_arcc);
}

cinterval atanh(const cinterval& z) 
{
  if(  (inf(Im(z)) <= 0.0) && (0.0 <= sup(Im(z))) &&   //Zero in imaginary part and
       ((inf(Re(z)) < -1.0) || (sup(Re(z)) > 1.0 ))  ) //real part intersects (-INFINITY,-1) and/or (1,+INFINITY)
    throw function_not_defined();

  interval im_art(0.0),m_re_art(0.0);
  bool singular( false );
  
  h_arctan(Im(z),-Re(z),im_art,m_re_art,singular);

  return cinterval(-m_re_art,im_art); //atanh(z) = -i * atan(i*z)
}

cinterval acoth(const cinterval& z) 
{
  if(  (inf(Im(z)) <= 0.0) && (0.0 <= sup(Im(z))) &&   //Zero in imaginary part and
       (sup(Re(z)) > -1.0) && (inf(Re(z)) < 1.0 )   )  //real part intersects (-1,1)
    throw function_not_defined();

  interval im_arco(0.0),m_re_arco(0.0);
  bool singular( false );

  h_arccot(-Im(z),Re(z),im_arco,m_re_arco,singular); 

  if( singular ) throw function_not_defined();//Singularity of arcoth in argument.
  else return cinterval(-m_re_arco,im_arco); 
}

/* ********************************************************* */
/* *** Computation of power function based on logarithm  *** */
/* ********************************************************* */

cinterval power(const cinterval& bas,int n) 
{
  if( n < 0 ) return 1.0/power( bas , -n );
  else if( n == 0 ) return cinterval( interval::ONE(), interval::ZERO() );
  else if( n == 1 ) return bas;
  else if( n == 2 ) return sqr( bas );
  else
    {
      if( (inf(Re(bas)) > 0)||(inf(Im(bas)) > 0)||(sup(Im(bas)) < 0) ) 
	return exp( n * ln(bas) );
      else
	{//Zero in argument
	  cinterval z( bas ),u( interval::ONE(), interval::ZERO() );
	  while( n > 0 )
	    {
	      if( (n%2) == 0 )//even exponent
		{
		  z = sqr(z);
		  n >>= 1;//Division by 2, realized with a right shift
		}
	      else
		{
		  u *= z;
		  n -= 1;
		}
	    }
	  return u;
	}
    }
}

cinterval pow(const cinterval& bas,const interval& n) 
{
  if( ( inf(Re(bas)) > 0 )||( inf(Im(bas)) > 0 )||( sup(Re(bas)) < 0 )||( sup(Im(bas)) < 0 ) ) 
    return exp(n*ln(bas));
  else throw function_not_defined();//Base contains zero.
}

cinterval pow(const cinterval& bas,const cinterval& n) 
{
  if( ( inf(Re(bas)) > 0 )||( inf(Im(bas)) > 0 )||( sup(Re(bas)) < 0 )||( sup(Im(bas)) < 0 ) ) 
    return exp(n*ln(bas));
  else throw function_not_defined();//Base contains zero.
}

/*

  End of File: cimath.cpp

*/
