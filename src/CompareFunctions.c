/*
 * This file is part of the GWHL (Gravitational Wave Helper Lib).
 *
 * GWHL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GWHL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GWHL.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <inttypes.h>

int GWHL_compareInt(const void *a, const void *b){
  int aa = *((int*) a);
  int bb = *((int*) b);
  
  if(aa > bb) return 1;
  if(aa < bb) return -1;
  if(aa == bb) return 0;
  return 0;
}


int GWHL_compareLong(const void *a, const void *b){
  long aa = *((long*) a);
  long bb = *((long*) b);
  
  if(aa > bb) return 1;
  if(aa < bb) return -1;
  if(aa == bb) return 0;
  return 0;
}


int GWHL_compareDouble(const void *a, const void *b){
  double aa = *((double*) a);
  double bb = *((double*) b);
  
  if(aa > bb) return 1;
  if(aa < bb) return -1;
  if(aa == bb) return 0;
  return 0;
}


int GWHL_compareFloat(const void *a, const void *b){
  float aa = *((float*) a);
  float bb = *((float*) b);
  
  if(aa > bb) return 1;
  if(aa < bb) return -1;
  if(aa == bb) return 0;
  return 0;
}


int GWHL_compareUInt64(const void *a, const void *b){
  uint64_t aa = *((uint64_t*) a);
  uint64_t bb = *((uint64_t*) b);
  
  if(aa > bb) return 1;
  if(aa < bb) return -1;
  if(aa == bb) return 0;
  return 0;
}
