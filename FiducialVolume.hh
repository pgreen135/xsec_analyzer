#pragma once
#include "TVector3.h"
// Definition of the fiducial volume for the analysis, together with a utility
// function template to check if a given point lies within it or not

// Boundaries of the neutrino vertex fiducial volume (cm)
// This is handled the same way for reco and in MC
double FV_X_MIN =   10.0;
double FV_X_MAX =  246.0;

double FV_Y_MIN = -101.0;
double FV_Y_MAX =  101.0;

double FV_Z_MIN =   10.0;
double FV_Z_MAX =  986.0;

// Use a template here so that this function can take float or double values as
// input
template <typename Number> bool point_inside_FV( Number x, Number y, Number z )
{
  bool x_inside_FV = ( FV_X_MIN < x ) && ( x < FV_X_MAX );
  bool y_inside_FV = ( FV_Y_MIN < y ) && ( y < FV_Y_MAX );
  bool z_inside_FV = ( FV_Z_MIN < z ) && ( z < FV_Z_MAX );
  return ( x_inside_FV && y_inside_FV && z_inside_FV );
}

inline bool point_inside_FV( const TVector3& pos ) {
  return point_inside_FV( pos.X(), pos.Y(), pos.Z() );
}

// Returns the number of Ar nuclei inside the fiducial volume
// PGreen: Changed to per nucleon, to allow direct comparison with other uboone Nue X-Secs
inline double num_Ar_targets_in_FV() {
  double volume = ( FV_X_MAX - FV_X_MIN ) * ( FV_Y_MAX - FV_Y_MIN )
    * ( FV_Z_MAX - FV_Z_MIN ); // cm^3
  constexpr double m_mol_Ar = 39.948; // g/mol
  constexpr double N_Avogadro = 6.02214076e23; // mol^(-1)
  constexpr double mass_density_LAr = 1.3836; // g/cm^3
  constexpr double n_nuc = 40; // Ar nucleons

  double num_Ar = volume * mass_density_LAr * N_Avogadro * n_nuc / m_mol_Ar;
  return num_Ar;
}

// Returns the total BNB muon neutrino flux (numu / cm^2) in the fiducial
// volume as a function of a given beam exposure (measured in
// protons-on-target)
// NOTE: This is currently approximated using the flux in the *active volume*.
// TODO: Revisit this approximation
inline double integrated_numu_flux_in_FV( double pot ) {
  // Obtained using the histogram hEnumu_cv (the central-value numu flux
  // in the MicroBooNE active volume as a function of neutrino energy)
  // stored in /pnfs/uboone/persistent/uboonebeam/bnb_gsimple
  // /bnb_gsimple_fluxes_01.09.2019_463_hist/. The ROOT commmands executed
  // were
  // root [3] hEnumu_cv->Scale( 1/(4997.*5e8)/(256.35*233.) )
  // root [4] hEnumu_cv->Integral()
  // (double) 7.3762291e-10
  // See the README file in that same folder for details.
  
  // normal
  /*
  double numu_per_cm2_per_POT_in_AV;
  if (useNuMI) {
    //numu_per_cm2_per_POT_in_AV = 1.844847e-11; // FHC NuMI Nue + Nuebar - Value from Krishan's thesis.
    numu_per_cm2_per_POT_in_AV = 1.84609e-11; // FHC NuMI Nue + Nuebar
    //numu_per_cm2_per_POT_in_AV = 1.71474e-11; // RHC NuMI Nue + Nuebar
  }
  else {
    numu_per_cm2_per_POT_in_AV = 7.3762291e-10;
  }
  double flux = pot * numu_per_cm2_per_POT_in_AV; // numu / cm^2
  */

  // combined
  // genie fake data
  //double flux = 1.844847e-11 * 12.22453e21 + 1.71474e-11 * 10.84891e21;
  // nuwro fake data
  double flux = 1.844847e-11 * 17.09056e21 + 1.71474e-11 * 8.508462e21;

  return flux;
}
