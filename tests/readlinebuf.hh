// Hacked (badly) by Roberto Bagnara <bagnara@cs.unipr.it>.

/******************************************************************************
 * $Revision$
 * $Date$
 * $Author$
 *
 * Contents: A streambuf which uses the GNU readline library for line I/O
 *           http://www.media.mit.edu/~vyzo/hacks/readlinebuf.html
 * (c) 2001 by Dimitris Vyzovitis [vyzo@media.mit.edu]
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA 02111-1307 USA
 *
 *****************************************************************************/

#ifndef _READLINEBUF_H_
#define _READLINEBUF_H_

#include <iostream>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <cstdio>

#ifdef HAVE_READLINE_READLINE_H
#include <readline/readline.h>
#else
extern "C" {
  char *readline(const char*);
  typedef int rl_command_func_t(int, int);
  int rl_bind_key(int, rl_command_func_t *);
  extern int rl_insert(int, int);
}
#endif

#ifdef HAVE_READLINE_HISTORY_H
#include <readline/history.h>
#else
extern "C" {
  void add_history(const char*);
}
#endif

#if (defined __GNUC__) && (__GNUC__ < 3)
#include <streambuf.h>
#else
#include <streambuf>
using std::streamsize;
using std::streambuf;
#endif

class readlinebuf : public streambuf {
public:
#if (defined __GNUC__) && (__GNUC__ < 3)
	typedef char char_type;
	typedef int int_type;
	typedef streampos pos_type;
	typedef streamoff off_type;
#endif
	static const int_type eof = EOF; // this is -1
	static const int_type not_eof = 0;

private:
	const char* prompt_;
	bool history_;
	char* line_;
	int low_;
	int high_;

	static char* my_readline (const char* prompt) {
	  char* line = readline(prompt);
	  if (line == NULL)
	    return line;
	  else {
	    size_t length = strlen(line);
	    char* new_line = (char*) malloc(length+2);
	    memcpy(new_line, line, length);
	    free((void*) line);
	    new_line[length++] = '\n';
	    new_line[length] = '\0';
	    return new_line;
	  }
	}

        static void my_add_history(char* line) {
	  size_t length = strlen(line);
	  line[length-1] = '\0';
	  add_history(line);
	  line[length-1] = '\n';
	}

protected:
		
	virtual int_type showmanyc() const { return high_ - low_; }
		
	virtual streamsize xsgetn( char_type* buf, streamsize n ) {
		int rd = n > (high_ - low_)? (high_ - low_) : n;
		memcpy( buf, line_ + low_, rd );
		low_ += rd;

		if ( rd < n ) {
			low_ = high_ = 0;
			free( line_ ); // free( NULL ) is a noop
			line_ = my_readline( prompt_ );
			if ( line_ ) {
				high_ = strlen( line_ );
				if ( history_ && high_ ) my_add_history( line_ );
				rd += xsgetn( buf + rd, n - rd );
			}
		}
			
		return rd; 
	}
		
	virtual int_type underflow() {
		if ( high_ == low_ ) {
			low_ = high_ = 0;
			free( line_ ); // free( NULL ) is a noop
			line_ = my_readline( prompt_ );
			if ( line_ ) {
				high_ = strlen( line_ );
				if ( history_ && high_ ) my_add_history( line_ );
			}
		}
			
		if ( low_ < high_ ) return line_[low_];
		else if ( line_ ) return underflow(); // empty line - retry
		else return eof;
	}
		
	virtual int_type uflow() {
		int_type c = underflow();
		if ( c != eof ) ++low_;
		return c;
	}
		
	virtual int_type pbackfail( int_type c = eof ) {
		if ( low_ > 0 )	--low_;
		else if ( c != eof ) {
			if ( high_ > 0 ) {
				char* nl = (char*)realloc( line_, high_ + 1 );
				if ( nl ) {
					line_ = (char*)memmove( nl + 1, line_, high_ );
					high_ += 1;
					line_[0] = char( c );
				} else return eof;
			} else {
				assert( !line_ );
				line_ = (char*)malloc( sizeof( char ) );
				*line_ = char( c );
				high_ = 1;
			}
		} else return eof;

		return not_eof;
	}
 		
public:
	readlinebuf( const char* prompt = NULL, bool history = true ) 
		: prompt_( prompt ), history_( history ),
		  line_( NULL ), low_( 0 ), high_( 0 ) {
		// Disable filename expansion.
		rl_bind_key( '\t', rl_insert );
		setbuf( 0, 0 );
	}
		
		
};

#endif
