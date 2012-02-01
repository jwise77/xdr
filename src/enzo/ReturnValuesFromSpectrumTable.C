/************************************************************************
/
/   RETURN VALUES AFTER READING (BLACK-BODY) SPECTRUM TABLE
/
/   written by: Ji-hoon Kim
/   date: February, 2010
/   modified1:
/
/   PURPOSE: For a given black body spectrum read in ReadSpectrumTable,
/            return the fraction of absorbed photons or the mean energy
/            of the spectrum at this column density
/
/   INPUTS:  type: type of absorber
/            ColumnDensity: column density so far
/            dColumnDensity: column density added in this cell
/
/   RETURNS: for mode = 0,1,2: the fraction of absorbed photons 
/                              by the absorber species HI, HeI, HeII
/            for mode = 3    : the mean energy of the spectrum
/
*************************************************************************/

#ifdef TRANSFER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "CosmologyParameters.h"

float ReturnEnergyFromSpectrumTable(float ColumnDensity, float dColumnDensity)
{

  int index_in, index_out;
  float mean_energy;
  float logC_start, logC_end, logC_step, logC_in, logC_out;

  int nbins  = RadiativeTransferSpectrumTable.NumberOfColumnDensityBins;
  logC_start = logf(RadiativeTransferSpectrumTable.columndensity_table[0]);
  logC_end   = logf(RadiativeTransferSpectrumTable.columndensity_table[nbins-1]);
  logC_step  = (logC_end - logC_start) / (nbins-1);

  /* calculate indices for ColumnDensity and ColumnDensity+dColumnDensity */

  if (ColumnDensity > 0) {
    logC_in   = logf(ColumnDensity);
    index_in  = min(nbins-1, max(0, int((logC_in  - logC_start)/logC_step)));
  } else {
    index_in = 0;
  }

  logC_out  = logf(ColumnDensity + dColumnDensity);
  index_out = min(nbins-1, max(0, int((logC_out - logC_start)/logC_step)));

  mean_energy = RadiativeTransferSpectrumTable.meanenergy_table[index_in];
  return mean_energy;

}

int ReturnValuesFromSpectrumTable(float ColumnDensity, float dColumnDensity,
				  float result[])
{

  int index_in, index_out, bin;
  float frac_in = 0.0, frac_out = 0.0, photon_fraction = 0.0;
  float mean_energy, change, tau;
  float logC_start, logC_end, logC_step, logC_in, logC_out;

  int nbins  = RadiativeTransferSpectrumTable.NumberOfColumnDensityBins;
  logC_start = logf(RadiativeTransferSpectrumTable.columndensity_table[0]);
  logC_end   = logf(RadiativeTransferSpectrumTable.columndensity_table[nbins-1]);
  logC_step  = (logC_end - logC_start) / (nbins-1);

  /* calculate indices for ColumnDensity and ColumnDensity+dColumnDensity */

  if (ColumnDensity > 0) {
    logC_in   = logf(ColumnDensity);
    index_in  = min(nbins-1, max(0, int((logC_in  - logC_start)/logC_step)));
  } else {
    index_in = 0;
  }

  logC_out  = logf(ColumnDensity + dColumnDensity);
  index_out = min(nbins-1, max(0, int((logC_out - logC_start)/logC_step)));

  for (bin = 0; bin < 3; bin++) {

    /* find photons left for ColumnDensity and ColumnDensity+dColumnDensity */

    frac_in  = RadiativeTransferSpectrumTable.fractionphotons_table[bin][index_in];
    frac_out = RadiativeTransferSpectrumTable.fractionphotons_table[bin][index_out];

    /* find fraction of photons absorbed by species w.r.t. the incoming number of photons;
       using the same logic for tau in Grid_WalkPhotonPackage, if "change" is smaller 
       than float precision, try linear interpolation to get absorbed photon fraction */

    /* tau = dColumnDensity * FindCrossSection(0, mean_energy), but we cannot use this!
       instead, we find the proportionality constant near ColumnDensity using
       frac_in = exp(-pseudo_CrossSection * ColumnDensity) */

    // JHW: Calculated in ReadRadiativeTransferSpectrumTable.C now
    
    //float pseudo_CrossSection = -log(frac_in) / 
    //  max(ColumnDensity, RadiativeTransferSpectrumTable.columndensity_table[0]);

    /* expf(-tau) = frac_out / frac_in, below is to avoid cases such as frac_in = 0.0 */

    change = frac_out / frac_in;
    tau = (isnan(change)) ? 
      dColumnDensity * RadiativeTransferSpectrumTable.pseudo_CrossSection[bin] : -logf(change);

    /* return photon_fraction */

    if (tau > 2.e1) 
      result[bin] = (1.0+BFLOAT_EPSILON);
    else if (tau > 1.e-4) 
      result[bin] = min(1.0f - frac_out / frac_in, 1.0f);
    else
      result[bin] = min(dColumnDensity * RadiativeTransferSpectrumTable.pseudo_CrossSection[bin], 1.0f);

  } // ENDFOR bin

//    if (index_in==0 && index_out==0)
//    fprintf(stderr, "RVFST: id_in = %d, id_out = %d, f_in =%f(%g), f_out = %f(%g), tau = %f, photon_f = %f\n", 
//	    index_in, index_out, frac_in, ColumnDensity, frac_out, ColumnDensity+dColumnDensity, tau, photon_fraction); 

  /* find mean energy */
  
  result[3] = RadiativeTransferSpectrumTable.meanenergy_table[index_in];

  return SUCCESS;

}

#endif

