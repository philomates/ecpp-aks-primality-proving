// Helper Functions


#include <stdio.h>

// check if number is in discriminant class number table
//
// To test correctness, drop this in main()
// assert(find(-427, hD2, 17));
// assert(find(-15, hD2, 17));
bool find(signed long int elem, int* array, int length)
{
  int i = 0;
  while(i < length)
  {
    if(array[i++] == elem)
      return true;
  }
  return false;
}
