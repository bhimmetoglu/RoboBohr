/* 
  Functions for generating pair features
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void pairFeatures(int nSpecies, int natoms, int * elements, double * pos, double * featMat)
{
  int pairType[nSpecies*nSpecies];
  int nPairs = nSpecies*(nSpecies+1)/2; // Number of distinct pairs
  int dimFeat = natoms*(natoms-1)/2 * nPairs; // Dimension of the feature matrix
  int indexPair;
  double R1[3], R2[3];
  double dist;

  // Create pairType matrix
  int k = 0;
  for (int i = 0; i < nSpecies; i++){
    for (int j = i; j < nSpecies; j++){
      pairType[j + nSpecies*i] = k;
      k++;
    }
  }

  // Symmetrize
  for (int j = 0; j < nSpecies; j++){
    for (int i = j+1; i < nSpecies; i++){
      pairType[j + nSpecies*i] = pairType[i + nSpecies*j];
    }
  }

  /* TEST
  for (int i=0; i < nSpecies; i++){
    for (int j=0; j < nSpecies; j++){
      printf("%d ", pairType[j+nSpecies*i]);
      if (j == (nSpecies-1))
        printf("\n");
    }
  }

  */

  // Inititate featMat
  for (int i = 0; i < dimFeat; i++){
    featMat[i] = 0.0;
  }

  // Fill featMat
  k = 0;
  for (int iat = 0; iat < natoms; iat++){
    for (int jat = iat + 1; jat < natoms; jat++){ // I can parallelize this inner loop
      indexPair = pairType[elements[jat] + nSpecies*elements[iat]];

      // Distance between pairs
      for (int icart = 0; icart < 3; icart++){
        R1[icart] = pos[icart+3*iat];
        R2[icart] = pos[icart+3*jat];
      }
      //printf("%lf %lf %lf\n", R1[0], R1[1], R1[2]);
      dist = sqrt( (R1[0]-R2[0])*(R1[0]-R2[0]) + (R1[1]-R2[1])*(R1[1]-R2[1]) + (R1[2]-R2[2])*(R1[2]-R2[2]) );
      featMat[indexPair + nPairs*k] = dist;
      k++;
    } // -- jat -- //
  } // -- iat -- //

}
