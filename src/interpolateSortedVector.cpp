#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector interpolateSortedVectorRcpp_(NumericVector vector, int newSize)
{
  if(newSize < 0)
  {
    throw std::range_error("Inadmissible negative value of newSize!");
  }
  int size = vector.size();
  if(size == 0 || newSize == 0)
  {
    NumericVector newVector; //see http://stackoverflow.com/questions/30129968/how-to-initialize-numericvector-to-a-specific-size-after-it-has-been-declared
    return newVector;
  }
  NumericVector newVector = NumericVector(newSize);
  if(size == 1)
  {
    for(int i = 0; i != newSize; i++)
    {
      newVector[i] = vector[0];
    }
    return newVector;
  }
  if(newSize == 1)
  {
    double sum, mean;
    sum = 0.0;
    for(int j = 0; j != size; j++)
    {
      sum += vector[j];
    }
    mean = sum / size;
    newVector[0] = mean;
    return newVector;
  }
  newVector[0] = vector[0];
  double partialVectorIndex = 0.0;
  double indexIncrement = ((double)(size-1))/((double)(newSize-1));
  for(int i = 1; i != newSize-1; i++)
  {
    partialVectorIndex += indexIncrement;
    int ind = ((int) partialVectorIndex);
    double decimalPart = partialVectorIndex - ind;
    double vectorDiff = vector[ind+1] - vector[ind];
    newVector[i] = vector[ind] + decimalPart * vectorDiff; 
  }
  newVector[newSize - 1] = vector[size - 1];
  return newVector;
}