/***********************************************************************
/
/  READ THE SPECTRUM TABLE
/
/  written by: Ji-hoon Kim
/  date:       February, 2010
/  modified1:
/
/  PURPOSE:
/    For RadiaitveTransferTraceSpectrum = TRUE, initialize the spectrum
/    table.  The table must be tabulated against column density (rows),
/    and should contain (1) the fraction of photons left in the ray, and
/    (2) the mean energy of the spectrum at a given column density.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#ifdef TRANSFER

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

int ReadRadiativeTransferSpectrumTable(float TemperatureUnits, float LengthUnits, 
				       float aUnits, float DensityUnits, float TimeUnits)
{

  FILE *fptr;
  char line[MAX_LINE_LENGTH];
  int i, nbins;
  float TemperatureRange[2];

  // open the file

  if ((fptr = fopen(RadiativeTransferTraceSpectrumTable, "r")) == NULL) {
    ENZO_VFAIL("Error opening spectrum table %s\n", 
	    RadiativeTransferTraceSpectrumTable)
  }

  // The second line gives the number of lines in this file

  fgets(line, MAX_LINE_LENGTH, fptr); // pass the first line
  fgets(line, MAX_LINE_LENGTH, fptr);
  if ((sscanf(line, "# %"ISYM, &nbins)) != 1) {
    ENZO_FAIL("Error reading number of bins (line 2)\n");
  }

  // Initialize

  RadiativeTransferSpectrumTable.NumberOfColumnDensityBins = nbins;
  RadiativeTransferSpectrumTable.columndensity_table = new float[nbins];
  RadiativeTransferSpectrumTable.fractionphotons_table[0] = new float[nbins];
  RadiativeTransferSpectrumTable.fractionphotons_table[1] = new float[nbins];
  RadiativeTransferSpectrumTable.fractionphotons_table[2] = new float[nbins];
  RadiativeTransferSpectrumTable.meanenergy_table = new float[nbins];

  // Read the spectrum table

  for (i = 0; i < nbins; i++) {
      if (fscanf(fptr, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %f",
		 &RadiativeTransferSpectrumTable.columndensity_table[i],
		 &RadiativeTransferSpectrumTable.fractionphotons_table[0][i],
		 &RadiativeTransferSpectrumTable.fractionphotons_table[1][i],
		 &RadiativeTransferSpectrumTable.fractionphotons_table[2][i],
		 &RadiativeTransferSpectrumTable.meanenergy_table[i])
	  != 5) {
	ENZO_VFAIL("Error reading RadiationData line %"ISYM"\n", i)

      }
  }

  // Calculate "pseudo cross-section" from table.  First search where
  // tau ~ 1.0, then estimate cross-section from that entry.

  int n = 0;
  float f0, f1, stopf;
  f0 = RadiativeTransferSpectrumTable.fractionphotons_table[0][0];
  stopf = f0 / M_E;
  do {
    f1 = RadiativeTransferSpectrumTable.fractionphotons_table[0][n++];
  } while (f1 > stopf);

  for (i = 0; i < 3; i++)
    RadiativeTransferSpectrumTable.pseudo_CrossSection[i] = 
      (RadiativeTransferSpectrumTable.fractionphotons_table[i][0] -
       RadiativeTransferSpectrumTable.fractionphotons_table[i][n]) /
      (RadiativeTransferSpectrumTable.columndensity_table[n] -
       RadiativeTransferSpectrumTable.columndensity_table[0]);

  fclose(fptr);


#ifdef OUTPUT_SPECTRUM_FOR_CHECK
  // Write out the spectrum table for check

  FILE *fptr2;

  fptr2 = fopen("spectrum_table.out", "w");
  for (i = 0; i < RadiativeTransferSpectrumTable.NumberOfColumnDensityBins; i++) 
    fprintf(fptr2, "%g  %"FSYM"  %"FSYM"  %"FSYM"  %f\n", 
	    RadiativeTransferSpectrumTable.columndensity_table[i],
	    RadiativeTransferSpectrumTable.fractionphotons_table[0][i],
	    RadiativeTransferSpectrumTable.fractionphotons_table[1][i],
	    RadiativeTransferSpectrumTable.fractionphotons_table[2][i],
	    RadiativeTransferSpectrumTable.meanenergy_table[i]);

  fclose(fptr2);
#endif

  return SUCCESS;

}

#endif
