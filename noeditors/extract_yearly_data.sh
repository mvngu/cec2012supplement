#!/usr/bin/env bash

###########################################################################
# Copyright (C) 2011 Minh Van Nguyen <nguyenminh2@gmail.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# http://www.gnu.org/licenses/
###########################################################################

# Extracting yearly data from the genetic programming dataset.
# NOTE: There are missing data for the years: 1951, 1952, 1953, 1954, 1955,
# 1956, 1957, 1959, 1960, 1961, 1962, 1964, 1965, 1966, 1967, 1968, 1971,
# 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1980, 1983
for year in 1950 \
    1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 \
    1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 \
    1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 \
    1981 1982 1983 1984 1985 1986 1987 1988 1989 1990 \
    1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 \
    2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 \
    2011; do
    grep "$year," gp-bib.csv > "$year".dat
done
