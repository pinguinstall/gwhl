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

#ifndef _____GWHL_STATSCONTAINER_
#define _____GWHL_STATSCONTAINER_

typedef struct GWHL_StatsContainer {
  double min;
  double max;
  double mean;
  double median;
  double q5;
  double q10;
  double q90;
  double q95;
  double skewness;
  double kurtosis;
  double stddev;
  double variance;
  unsigned long long int numPoints;
} GWHL_StatsContainer;

#endif
