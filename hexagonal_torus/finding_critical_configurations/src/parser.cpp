///////////////////////////////////////////////////////
//
//  parser.cpp
//
//  Purpose:  Implements generalized parsing tools, such as tokenizer parsers, 
//            and lexer
//       
///////////////////////////////////////////////////////
//
//  Copyright (c) 2013, Lawrence Livermore National Security, LLC.  Produced
//  at the Lawrence Livermore National Laboratory.  Written by Jeremy Mason,
//  reachable at jkylemason@gmail.com.
//  
//  CODE-636759. All rights reserved.
//  
//  This file is part of the Critical Configurations of Hard Disks on the 
//  Torus.  Please read LICENSE.txt for Our Notice and GNU General Public 
//  License information.
//  
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License (as published by
//  the Free Software Foundation) version 2, dated June 1991.
//  
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
//  conditions of the GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License along
//  with this program; if not, write to the Free Software Foundation, Inc.,
//  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//       
///////////////////////////////////////////////////////


#include "parser.h"


namespace MorseTheory
{

  //----------------------------------------
  // Tokenize - given the character string buffer sBuf and delimiters,
  // a set of tokens will be outputted as vectors of vectors of strings
  // (lines of strings)
  //----------------------------------------
  void Tokenize( vector<vector<string> > &vsTokens, const string &sBuf,  const string &sDelimiters)
  {
    typedef size_t usize;  
    usize iCurrentPos = 0;
    usize iDelimiterPos = 0;
    usize iEndLine = 0;
    
    while( iCurrentPos < sBuf.size() )
    {
      vector<string> vsCurrentLine;
      
      iEndLine = sBuf.find_first_of("\n", iEndLine + 1);
		
      if( iEndLine == string::npos )
        return;
      
      while( iCurrentPos < iEndLine )
      {
        string currentToken;
        
        // Find either space of tab or endline
        iDelimiterPos = sBuf.find_first_of(sDelimiters.c_str(), iDelimiterPos);
        iCurrentPos = sBuf.find_first_not_of(sDelimiters.c_str(), iCurrentPos);
		
        if( (iCurrentPos == string::npos) || (iDelimiterPos == string::npos) )
          break;
        
        for( usize i = iCurrentPos; i < iDelimiterPos; i++ )
          currentToken.append(sizeof(char), sBuf.at(i));
        
        if( currentToken.size() > 0 )
          vsCurrentLine.push_back(currentToken);
			
        iCurrentPos = iDelimiterPos + 1;
        iDelimiterPos ++;
      }
      vsTokens.push_back(vsCurrentLine);
    }
  }
  
}
