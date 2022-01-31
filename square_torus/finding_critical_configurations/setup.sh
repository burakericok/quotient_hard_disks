#!/bin/bash
#
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

if [ $# -ne 1 ]
  then
    echo "script requires exactly one argument"
    exit 1
fi

if ! printf '%d' "$1" >/dev/null
  then
    echo "script requires an integer"
    exit 1
fi

NUMBER=$(($1))
STRING=$(printf "%02d" $NUMBER)

sed -E -i.bak "s:^(NDisks[^0-9]+)[0-9]+$:\1${NUMBER}:" options.txt
rm options.txt.bak 2>/dev/null
cp setup/$STRING/*.cpp src/
rm src/energy.*o src/hessian.*o src/jacobian.*o src/radius.*o 2>/dev/null
make parallel
