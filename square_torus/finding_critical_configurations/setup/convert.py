#--------------------------------------------------------------------------
# 
#  Copyright (c) 2013, Lawrence Livermore National Security, LLC.  Produced
#  at the Lawrence Livermore National Laboratory.  Written by Jeremy Mason,
#  reachable at jkylemason@gmail.com.
#  
#  CODE-636759. All rights reserved.
#  
#  This file is part of the Critical Configurations of Hard Disks on the 
#  Torus.  Please read LICENSE.txt for Our Notice and GNU General Public 
#  License information.
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License (as published by
#  the Free Software Foundation) version 2, dated June 1991.
# 
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
#  conditions of the GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License along
#  with this program; if not, write to the Free Software Foundation, Inc.,
#  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#--------------------------------------------------------------------------

import os
import re

def radius(line):
    '''
    Correctly formats return value for radius.cpp
    '''
    output = re.search(r'^(\w+) = \[([^]]*)\]', line)
    entries = output.group(2).split(',')
    return 'arma::vec '+output.group(1)+';\n\t'+output.group(1)+' << '+' << '.join(entries)+';\n\t'+'return '+output.group(1)+';'
    
def energy(line):
    '''
    Correctly formats return value for energy.cpp
    '''
    output = re.search(r'^(\w+)', line)
    return 'double '+output.group(1)+';\n\t'+line+'\n\t'+'return '+output.group(1)+';'

def jacobian(line):
    '''
    Correctly formats return value for jacobian.cpp
    '''
    output = re.search(r'^(\w+) = \[([^]]*)\]', line)
    entries = output.group(2).split(',')
    return 'arma::vec '+output.group(1)+';\n\t'+output.group(1)+' << '+' << '.join(entries)+';\n\t'+'return '+output.group(1)+';'

def hessian(line):
    '''
    Correctly formats return value for hessian.cpp
    '''
    output = re.search(r'^(\w+) = [^[]*\[([^]]*)\][^[]*\[(\d+)\D+(\d+)\]', line)
    entries = output.group(2).split(',')
    rows = [' << '.join(row) for row in zip(*[iter(entries)]*int(output.group(3)))]
    return 'arma::mat '+output.group(1)+';\n\t'+output.group(1)+' << '+' << arma::endr << '.join(rows)+' << arma::endr;\n\t'+'return '+output.group(1)+';'

def convert(pathname):
    '''
    Reformats input MATLAB file as the output C++ file
    Requires an existing C++ file to modify in place
    '''
    (dirname, filename) = os.path.split(pathname)
    (shortname, extension) = os.path.splitext(filename)
    with open(pathname, encoding='utf-8') as in_file:
        in_string = in_file.read()
        
        # Replace powers with corresponding C functions
        in_string = re.sub(r'\b(\w+)\^2(?![.\d])', r'sqr(\1)', in_string)
        in_string = re.sub(r'\b(\w+)\^([.\d]+)', r'pow(\1,\2)', in_string)
        
        # Disposes of all lines that do not perform a computation
        in_string = re.sub(r'^(?:[^;]*\n)*', '', in_string)
        
        # Identifies all local variablese excluding return value
        vars = re.findall(r'([xyt]\d+) = ', in_string)
        
        # Extracts coordinates from arma::vec function argument
        a_pattern = re.compile(r'([xy]\d+) = in\d+\(:,(\d+)\);\n')
        from_arma = [a+' = coords('+str(int(b)-1)+');' for a, b in a_pattern.findall(in_string)]
        in_string = a_pattern.sub('', in_string)
        
        # Correctly formats return value
        lines = in_string.splitlines()
        lines[-1] = globals()[shortname](lines[-1])
        
        # Reassembles output string
        in_string = '\tdouble '+', '.join(vars)+';\n\t'+'\n\t'.join(from_arma)+'\n\t'+'\n\t'.join(lines)+'\n'
        
    with open(os.path.join(dirname, shortname+'.cpp'), mode='r+', encoding='utf-8') as out_file:
        out_string = out_file.read()
        out_file.seek(0)
        out_file.truncate(0)
        out_file.write(re.sub(r'(?<={\n)[^}]*(?=})', in_string, out_string))

if __name__ == '__main__':
    convert('radius.m')
    convert('energy.m')
    convert('jacobian.m')
    convert('hessian.m')
