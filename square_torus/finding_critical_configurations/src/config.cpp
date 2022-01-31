////////////////////////////////////////////////////
//
//  config.h
//
//  Purpose:  Implements class declaration that interfaces with the
//            config files for the Morse Theory simulation
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


#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>
#include "config.h"


namespace MorseTheory
{

  //----------------------------------------
  // Read
  //----------------------------------------
  bool ConfigFile::Read( const string & ConfigFilename )
  {
    usize  nBufferSize = 0;
    char *pBuffer = ReadFileToBuf( nBufferSize, ConfigFilename);
    
    if( nBufferSize<= 0 || !pBuffer )
    {
      cerr << "ERROR: Unable to read config file" << endl; 
      return false;
    }
	
    if( !Parse( string( pBuffer, nBufferSize ) ) )
    {
      cerr << "ERROR: Unable to parse config file" << endl;
      return false;
    }
    
    if( pBuffer )
      delete[] pBuffer;
    
    return true;
  }
  
  //----------------------------------------
  //  Parser
  //----------------------------------------
  bool ConfigFile::Parse( const string & sBuf )
  {
    const char * KeywordStrings[] =
    {
      
      "#",        // comment
      
      "NDisks",
      "NTrials",
      
      "NThreads",
      "WIncrements",
      
      "TimeInterval",
      "BackupInterval"
    };
    
    enum EKeyword
    {
      eComment,     // comment
      
      eNDisks,
      eNTrials,
      
      eNThreads,
      eWIncrements,
      
      eTimeInterval,
      eBackupInterval,
      
      eNumKeywords  //  error check on number of keywords
    };

    vector< vector<string> > vsTokens;
    Tokenize( vsTokens, sBuf, " \t\n" );
  
    // check lists to see which variable is not initialized
    bool *vInitializationCheck;
    bool *vRequirementCheck;
    vInitializationCheck = new bool[ eNumKeywords ];
    vRequirementCheck    = new bool[ eNumKeywords ];

    for( usize i = 0; i < eNumKeywords; i ++ )
    {
      vInitializationCheck[i] = false;
      vRequirementCheck[i]    = true;
    }
    
    for( usize i = 0; i < vsTokens.size(); i ++ )
    {
      usize iFirstToken;
      
      if( vsTokens[i].size() == 0 )
      {
        iFirstToken = eComment;
      }
      else if( vsTokens[i][0].find_first_of( KeywordStrings[ eComment ] ) == 0 )
      {
        iFirstToken = eComment;
      }
      else
      {
        // Identify the keyword in the beginning of the line
        for( iFirstToken = 0; iFirstToken < eNumKeywords; iFirstToken ++ )
          if( strcmp( KeywordStrings[iFirstToken], vsTokens[i][0].c_str() ) == 0 )
            break;
      }
      vInitializationCheck[ iFirstToken ] = true;

      switch( iFirstToken ) // Look at first token of each line
      {
          
        case eComment:      // Comment
          break;
          
        case eNDisks:
          n_disks = std::atoi( vsTokens[i][1].c_str() );
          break;
        case eNTrials:
          n_trials = std::atoi( vsTokens[i][1].c_str() );
          break;
          
        case eNThreads:
          n_threads = std::atoi( vsTokens[i][1].c_str() );
          break;
        case eWIncrements:
          w_increments = std::atoi( vsTokens[i][1].c_str() );
          break;
          
        case eTimeInterval:
          time_interval = std::atoi( vsTokens[i][1].c_str() );
          break;
        case eBackupInterval:
          backup_interval = std::atoi( vsTokens[i][1].c_str() );
          break;
          
        default:
          cerr << "[ConfigFile] Error: syntax not recognized:  Line " << i
               << " Keyword: " << vsTokens[i][0] << endl;
          return false;
      }
    } // End loop over tokens
    
    bool Success = true;
    for( usize i = 0; i < eNumKeywords; i ++ )
    {
      if( ! vInitializationCheck[i] && vRequirementCheck[i] )
      {
        cerr << "[ConfigFile] Initialization Error: Missing parameter: "
             << KeywordStrings[i]  << " not optional." << endl;
        Success = false;
      }
    }

    if( !Success )
    {
      cerr << "[ConfigFile] Initialization Failed" << endl;
      return false;
    }
    
    delete [] vInitializationCheck;
    delete [] vRequirementCheck;
    
    return true;
  }


  //----------------------------------------
  //  Returns a Null terminated C-string
  //----------------------------------------
  char* ConfigFile::ReadFileToBuf( usize & nBufferSize, string filename )
  {
    FILE *pFile;
    char *pBuffer;

    if ( (pFile = fopen(filename.c_str(), "r" )) == NULL )
    {
      cerr << "[ReadFileToBuf] Cannot Open File: " << filename.c_str()  << endl;
      return NULL;
    }

    // Get the size of the file
    fseek( pFile, 0L, SEEK_END ); // Position to end of file
    nBufferSize = ftell( pFile ); // Get file length
    rewind( pFile );              // Back to start of file

    if( nBufferSize <= 0 )
      return NULL;
    
    // Read in the entire file and close the file handle
    pBuffer = new char[nBufferSize + 1];

    if ( fread( pBuffer, nBufferSize, 1, pFile ) <= 0 )
    {
      fclose( pFile );
      cerr << "[ReadFileToBuf] File read error" << endl;
      return NULL;
    }

    fclose( pFile );

    pBuffer[ nBufferSize ] = '\0'; // NULL termination of string
    return pBuffer;
  }

}
