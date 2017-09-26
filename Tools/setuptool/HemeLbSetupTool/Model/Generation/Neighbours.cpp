// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "Neighbours.h"

const Neighbours& Neighbours::Get() {
  static Neighbours instance;
  return instance;
}

const std::vector<Index>& Neighbours::GetDisplacements() {
  return Get().neighbours;
}
const std::vector<int>& Neighbours::GetOpposites() {
  return Get().inverses;
}
 
Neighbours::Neighbours() : neighbours({
    {-1,-1,-1},
    {-1,-1, 0},
    {-1,-1,+1},
    {-1, 0,-1},
    {-1, 0, 0},
    {-1, 0,+1},
    {-1,+1,-1},
    {-1,+1, 0},
    {-1,+1,+1},
    { 0,-1,-1},
    { 0,-1, 0},
    { 0,-1,+1},
    { 0, 0,-1},
    // { 0, 0, 0},
    { 0, 0,+1},
    { 0,+1,-1},
    { 0,+1, 0},
    { 0,+1,+1},
    {+1,-1,-1},
    {+1,-1, 0},
    {+1,-1,+1},
    {+1, 0,-1},
    {+1, 0, 0},
    {+1, 0,+1},
    {+1,+1,-1},
    {+1,+1, 0},
    {+1,+1,+1}
  }),
	       inverses({25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0})
{
}
