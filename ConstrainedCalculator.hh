#pragma once

// Standard library includes
#include <algorithm>
#include <stdexcept>

// STV analysis includes
#include "MatrixUtils.hh"
#include "SystematicsCalculator.hh"
#include "utils.hh"

// Calculates covariance matrices describing the uncertainty on the reco-space
// event counts. Distinguishes between "signal region" reco bins and "sideband"
// reco bins. For the former, covariances are calculated for both the total
// event count (signal + background) and background only. For the latter,
// covariances are calculated on just the total event count. Using the
// concept of "conditional covariance" (see MicroBooNE DocDB #32672), this
// class will calculate a constrained signal + background prediction
// together with its uncertainty. A constrained CV background prediction
// is also obtained that should be subtracted from the measured data points
// and the constrained signal + background prediction before unfolding.
class ConstrainedCalculator : public SystematicsCalculator {

  public:

    ConstrainedCalculator( const std::string& input_respmat_file_name,
      const std::string& syst_cfg_file_name = "",
      const std::string& respmat_tdirectoryfile_name = "" );

    virtual double evaluate_observable( const Universe& univ, int cm_bin,
      int flux_universe_index = -1 ) const override;

    virtual double evaluate_observable( const Universe& univ, int reco_bin,
      std::string event_category, int flux_universe_index = -1 ) const override;

    virtual double evaluate_mc_stat_covariance( const Universe& univ,
      int cm_bin_a, int cm_bin_b ) const override;

    virtual double evaluate_data_stat_covariance( int cm_bin_a, int cm_bin_b,
      bool use_ext ) const override;

    // This class uses a dimension for the covariance matrix that is different
    // than just the number of reco bins, so we need to override this virtual
    // function.
    virtual size_t get_covariance_matrix_size() const override;

    // Apply a sideband constraint before subtracting the EXT+MC background.
    // NOTE: this function assumes that the ordinary reco bins are all listed
    // before any sideband reco bins
    virtual MeasuredEvents get_measured_events(std::string type = "default") const override;

  protected:

    // Helper type used to decide how to handle the covariance calculation
    enum ConstrainedCalculatorBinType {
      // The bin is an "ordinary" one, and the observable of interest is the
      // total number of reconstructed events (signal + background)
      kOrdinaryRecoBinAll,

      // The bin is an "ordinary" one, and the observable of interest is the
      // number of reconstructed background events
      kOrdinaryRecoBinBkgd,

      // The bin is an "ordinary" one, and the observable of interest is the
      // number of reconstructed signal events
      kOrdinaryRecoBinSignal,

      // The bin is a "sideband" one, and the observable of interest is the
      // total number of reconstructed events (signal + background)
      kSidebandRecoBinAll
    };

    // Helper function for linking a covariance matrix bin number to
    // information about how to compute observables within it
    int get_reco_bin_and_type( int cm_bin,
      ConstrainedCalculatorBinType& bin_type ) const;

    // Number of "sideband" reco bins
    size_t num_sideband_reco_bins_ = 0u;
};

ConstrainedCalculator::ConstrainedCalculator(
  const std::string& input_respmat_file_name,
  const std::string& syst_cfg_file_name,
  const std::string& respmat_tdirectoryfile_name )
  : SystematicsCalculator( input_respmat_file_name,
  syst_cfg_file_name, respmat_tdirectoryfile_name )
{
  num_sideband_reco_bins_ = 0u;
  bool found_first_sideband_bin = false;
  for ( const auto& rbin : reco_bins_ ) {
    if ( rbin.type_ == kOrdinaryRecoBin ) {
      if ( found_first_sideband_bin ) throw std::runtime_error( "Ordinary"
        " reco bins must precede sideband ones in the UniverseMaker"
        " configuration file." );
    }
    else if ( rbin.type_ == kSidebandRecoBin ) {
      found_first_sideband_bin = true;
      ++num_sideband_reco_bins_;
    }
  }
}

double ConstrainedCalculator::evaluate_observable( const Universe& univ,
  int cm_bin, int flux_universe_index ) const
{
  // For the ConstrainedCalculator class, the observable of interest is the
  // total number of events (either signal + background or background only) in
  // the current bin in reco space
  double reco_bin_events = 0.;

  bool use_detVar_CV = this->is_detvar_universe( univ );
  bool use_altCV_CV = this->is_altcv_universe( univ );
  bool use_flugg_CV = this->is_flugg_universe( univ );

  // Get access to the CV universe. We need it regardless of the input universe
  // so that we can use it in the denominator of the smearceptance matrix
  // element. Note that we should use the detVarCV universe as the CV when the
  // input universe corresponds to a detector variation (or is the detVarCV
  // universe itself). Based on the check above, we assign a pointer to
  // either the regular or detVar CV here as appropriate.
  const Universe* cv_univ = nullptr;
  if ( use_detVar_CV ) {
    cv_univ = detvar_universes_.at( NFT::kDetVarMCCV ).get();
  }
  else if ( use_altCV_CV ) {
    // use altCVMC CV universe for altCV MC
    cv_univ = alt_cv_universes_.at( NFT::kAltCVMCGenieCV ).get();
  }
  else if ( use_flugg_CV ) {
    // use fluggMC CV universe for flugg MC
    cv_univ = flugg_universes_.at( NFT::kFluggMCCV ).get();
  }
  else {
    cv_univ = &this->cv_universe();
  }

  ConstrainedCalculatorBinType bin_type;
  int reco_bin = this->get_reco_bin_and_type( cm_bin, bin_type );

  // Look up the block index for the current reco bin. We will use this
  // to avoid double-counting when summing over true bins below.
  const auto& rbin = reco_bins_.at( reco_bin );
  int reco_block_index = rbin.block_index_;

  size_t num_true_bins = true_bins_.size();

  // We need to sum the contributions of the various true bins,
  // so loop over them while checking whether each one is associated
  // with either signal or background
  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    const auto& tbin = true_bins_.at( tb );

    if ( tbin.type_ == kSignalTrueBin ) {

      // Ignore signal true bins outside of the same block as the current
      // reco bin. This avoids issues with double-counting.
      int true_block_index = tbin.block_index_;
      if ( reco_block_index != true_block_index ) continue;

      // Get the CV event count for the current true bin
      double denom_CV = cv_univ->hist_true_->GetBinContent( tb + 1 );

      // For the systematic variation universes, we want to assess
      // uncertainties on the signal only through the smearceptance
      // matrix. We therefore compute the smearceptance matrix element
      // here and then apply it to the CV expected event count in
      // each true bin.
      // NOTE: ROOT histogram bin numbers are one-based (bin zero is always
      // the underflow bin). Our bin indices therefore need to be offset by
      // +1 in all cases here.
      double numer = univ.hist_2d_->GetBinContent( tb + 1, reco_bin + 1 );
      double denom = univ.hist_true_->GetBinContent( tb + 1 );

      // I plan to extract the flux-averaged cross sections in terms of the
      // *nominal* flux model (as opposed to the real flux). I therefore
      // vary the numerator of the smearceptance matrix for these while
      // keeping the denominator equal to the CV expectation under the
      // nominal flux model. This is the same strategy as is used in the
      // Wire-Cell CC inclusive analysis.
      if ( flux_universe_index >= 0 ) {
        denom = denom_CV;
      }

      // If the denominator is nonzero actually calculate the fraction.
      // Otherwise, just leave it zeroed out.
      // TODO: revisit this, think about MC statistical uncertainties
      // on the empty bins
      double smearcept = 0.;
      if ( denom > 0. ) smearcept = numer / denom;

      // Compute the expected signal events in this universe
      // by multiplying the varied smearceptance matrix element
      // by the unaltered CV prediction in the current true bin.
      double expected_signal = smearcept * denom_CV;

      // For the ordinary background reco bins, don't include any signal
      // contribution
      if ( bin_type == kOrdinaryRecoBinBkgd ) {
        expected_signal = 0.;
      }
      //// TODO: REVISIT THIS
      //// For the sideband bins, go ahead and vary the signal prediction
      //// just like the background one
      //else if ( bin_type == kSidebandRecoBinAll ) {
      //  expected_signal = numer;
      //}

      // Compute the expected signal events in the current reco bin
      // with the varied smearceptance matrix (and, for flux universes,
      // the varied integrated flux)
      reco_bin_events += expected_signal;
    }
    else if ( tbin.type_ == kBackgroundTrueBin ) {
      // For background events, we can use the same procedure regardless
      // of whether we're in the CV universe or not
      double background = univ.hist_2d_->GetBinContent( tb + 1, reco_bin + 1 );

      // For the ordinary signal reco bins, don't include any background
      // contribution
      if ( bin_type == kOrdinaryRecoBinSignal ) {
        background = 0.;
      }

      reco_bin_events += background;
    }
  } // true bins

  return reco_bin_events;
}

double ConstrainedCalculator::evaluate_observable( const Universe& univ,
  int cm_bin, std::string event_category, int flux_universe_index ) const
{
  // For the ConstrainedCalculator class, the observable of interest is the
  // total number of events (either signal + background or background only) in
  // the current bin in reco space
  double reco_bin_events = 0.;

  bool use_detVar_CV = this->is_detvar_universe( univ );
  bool use_altCV_CV = this->is_altcv_universe( univ );
  bool use_flugg_CV = this->is_flugg_universe( univ );

  // Get access to the CV universe. We need it regardless of the input universe
  // so that we can use it in the denominator of the smearceptance matrix
  // element. Note that we should use the detVarCV universe as the CV when the
  // input universe corresponds to a detector variation (or is the detVarCV
  // universe itself). Based on the check above, we assign a pointer to
  // either the regular or detVar CV here as appropriate.
  const Universe* cv_univ = nullptr;
  if ( use_detVar_CV ) {
    cv_univ = detvar_universes_.at( NFT::kDetVarMCCV ).get();
  }
  else if ( use_altCV_CV ) {
    // use altCVMC CV universe for altCV MC
    cv_univ = alt_cv_universes_.at( NFT::kAltCVMCGenieCV ).get();
  }
  else if ( use_flugg_CV ) {
    // use fluggMC CV universe for flugg MC
    cv_univ = flugg_universes_.at( NFT::kFluggMCCV ).get();
  }
  else {
    cv_univ = &this->cv_universe();
  }

  ConstrainedCalculatorBinType bin_type;
  int reco_bin = this->get_reco_bin_and_type( cm_bin, bin_type );

  // Look up the block index for the current reco bin. We will use this
  // to avoid double-counting when summing over true bins below.
  const auto& rbin = reco_bins_.at( reco_bin );
  int reco_block_index = rbin.block_index_;

  size_t num_true_bins = true_bins_.size();

  // We need to sum the contributions of the various true bins,
  // so loop over them while checking whether each one is associated
  // with either signal or background

  // Select only true bin of interest based on event category
  // Construct string to check against
  std::string category_string = "category == " + event_category;

  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    const auto& tbin = true_bins_.at( tb );

    if (tbin.signal_cuts_ != category_string) continue;

    if ( tbin.type_ == kSignalTrueBin ) {

      // Ignore signal true bins outside of the same block as the current
      // reco bin. This avoids issues with double-counting.
      int true_block_index = tbin.block_index_;
      if ( reco_block_index != true_block_index ) continue;

      // Get the CV event count for the current true bin
      double denom_CV = cv_univ->hist_true_->GetBinContent( tb + 1 );

      // For the systematic variation universes, we want to assess
      // uncertainties on the signal only through the smearceptance
      // matrix. We therefore compute the smearceptance matrix element
      // here and then apply it to the CV expected event count in
      // each true bin.
      // NOTE: ROOT histogram bin numbers are one-based (bin zero is always
      // the underflow bin). Our bin indices therefore need to be offset by
      // +1 in all cases here.
      double numer = univ.hist_2d_->GetBinContent( tb + 1, reco_bin + 1 );
      double denom = univ.hist_true_->GetBinContent( tb + 1 );

      // I plan to extract the flux-averaged cross sections in terms of the
      // *nominal* flux model (as opposed to the real flux). I therefore
      // vary the numerator of the smearceptance matrix for these while
      // keeping the denominator equal to the CV expectation under the
      // nominal flux model. This is the same strategy as is used in the
      // Wire-Cell CC inclusive analysis.
      if ( flux_universe_index >= 0 ) {
        denom = denom_CV;
      }

      // If the denominator is nonzero actually calculate the fraction.
      // Otherwise, just leave it zeroed out.
      // TODO: revisit this, think about MC statistical uncertainties
      // on the empty bins
      double smearcept = 0.;
      if ( denom > 0. ) smearcept = numer / denom;

      // Compute the expected signal events in this universe
      // by multiplying the varied smearceptance matrix element
      // by the unaltered CV prediction in the current true bin.
      double expected_signal = smearcept * denom_CV;

      // For the ordinary background reco bins, don't include any signal
      // contribution
      if ( bin_type == kOrdinaryRecoBinBkgd ) {
        expected_signal = 0.;
      }
      //// TODO: REVISIT THIS
      //// For the sideband bins, go ahead and vary the signal prediction
      //// just like the background one
      //else if ( bin_type == kSidebandRecoBinAll ) {
      //  expected_signal = numer;
      //}

      // Compute the expected signal events in the current reco bin
      // with the varied smearceptance matrix (and, for flux universes,
      // the varied integrated flux)
      reco_bin_events += expected_signal;
    }
    else if ( tbin.type_ == kBackgroundTrueBin ) {
      // For background events, we can use the same procedure regardless
      // of whether we're in the CV universe or not
      double background = univ.hist_2d_->GetBinContent( tb + 1, reco_bin + 1 );

      // For the ordinary signal reco bins, don't include any background
      // contribution
      if ( bin_type == kOrdinaryRecoBinSignal ) {
        background = 0.;
      }

      reco_bin_events += background;
    }
  } // true bins

  return reco_bin_events;
}

double ConstrainedCalculator::evaluate_mc_stat_covariance( const Universe& univ,
  int cm_bin_a, int cm_bin_b ) const
{
  ConstrainedCalculatorBinType bin_type_a, bin_type_b;
  int reco_bin_a = this->get_reco_bin_and_type( cm_bin_a, bin_type_a );
  int reco_bin_b = this->get_reco_bin_and_type( cm_bin_b, bin_type_b );

  TH2D* hist_stats_2d = nullptr;

  // set correct histogram
  if (bin_type_a == kOrdinaryRecoBinAll && bin_type_b == kOrdinaryRecoBinAll) {
    hist_stats_2d = univ.hist_reco2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinAll && bin_type_b == kOrdinaryRecoBinBkgd) {
    hist_stats_2d = univ.hist_reco_bkgd2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinAll && bin_type_b == kOrdinaryRecoBinSignal) {
    hist_stats_2d = univ.hist_reco_signal2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinAll && bin_type_b == kSidebandRecoBinAll) {
    return 0;
    //hist_stats_2d = univ.hist_reco2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinBkgd && bin_type_b == kOrdinaryRecoBinAll) {
    hist_stats_2d = univ.hist_reco_bkgd2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinBkgd && bin_type_b == kOrdinaryRecoBinBkgd) {
    hist_stats_2d = univ.hist_reco_bkgd2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinBkgd && bin_type_b == kOrdinaryRecoBinSignal) {
    return 0;
  }
  else if (bin_type_a == kOrdinaryRecoBinBkgd && bin_type_b == kSidebandRecoBinAll) {
    return 0;
    //hist_stats_2d = univ.hist_reco_bkgd2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinSignal && bin_type_b == kOrdinaryRecoBinAll) {
    hist_stats_2d = univ.hist_reco_signal2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinSignal && bin_type_b == kOrdinaryRecoBinBkgd) {
    return 0;
  }
  else if (bin_type_a == kOrdinaryRecoBinSignal && bin_type_b == kOrdinaryRecoBinSignal) {
    hist_stats_2d = univ.hist_reco_signal2d_.get();
  }
  else if (bin_type_a == kOrdinaryRecoBinSignal && bin_type_b == kSidebandRecoBinAll) {
    return 0;
    //hist_stats_2d = univ.hist_reco_signal2d_.get();
  }
  else if (bin_type_a == kSidebandRecoBinAll && bin_type_b == kOrdinaryRecoBinAll) {
    return 0;
    //hist_stats_2d = univ.hist_reco2d_.get();
  }
  else if (bin_type_a == kSidebandRecoBinAll && bin_type_b == kOrdinaryRecoBinBkgd) {
    return 0;
    //hist_stats_2d = univ.hist_reco_bkgd2d_.get();
  }
  else if (bin_type_a == kSidebandRecoBinAll && bin_type_b == kOrdinaryRecoBinSignal) {
    return 0;
    //hist_stats_2d = univ.hist_reco_signal2d_.get();
  }
  else if (bin_type_a == kSidebandRecoBinAll && bin_type_b == kSidebandRecoBinAll) {
    hist_stats_2d = univ.hist_reco2d_.get();
  }
  else {
    //std::cout << "Error [ConstrainedCalculator::evaluate_mc_stat_covariance]: undefined bin combination." << std::endl;
    std::cout << "bin_type_a = " << bin_type_a << ", bin_type_b = " << bin_type_b << std::endl;

    // Danger
    return 0.;
  }

  // get uncertainty for this combination of bins
  if (hist_stats_2d != nullptr) {
    // ROOT histograms use one-based bin indices, so I correct for that here
    double err = hist_stats_2d->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  else {
    //std:cout << "error, here" << std::endl;
    return 0.;
  }

  /*
  if ( bin_type_a != kOrdinaryRecoBinBkgd && bin_type_b != kOrdinaryRecoBinBkgd && bin_type_a != kOrdinaryRecoBinSignal && bin_type_b != kOrdinaryRecoBinSignal)
  {
    // ROOT histograms use one-based bin indices, so I correct for that here
    double err = univ.hist_reco2d_->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  else if (bin_type_a == kOrdinaryRecoBinBkgd && bin_type_b == kOrdinaryRecoBinBkgd) 
  {    
    double err = univ.hist_reco_bkgd2d_->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  else if (bin_type_a == kOrdinaryRecoBinSignal && bin_type_b == kOrdinaryRecoBinSignal) {
    double err = univ.hist_reco_signal2d_->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  /*
  else if (bin_type_a == kOrdinaryRecoBinBkgd && (bin_type_b == kOrdinaryRecoBinBkgd || bin_type_b == kOrdinaryRecoBinAll)) 
  {    
    double err = univ.hist_reco_bkgd2d_->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  else if ((bin_type_a == kOrdinaryRecoBinBkgd || bin_type_a == kOrdinaryRecoBinAll) && bin_type_b == kOrdinaryRecoBinBkgd) 
  {    
    double err = univ.hist_reco_bkgd2d_->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  else if (bin_type_a == kOrdinaryRecoBinSignal && (bin_type_b == kOrdinaryRecoBinSignal || bin_type_b == kOrdinaryRecoBinAll)) {
    double err = univ.hist_reco_signal2d_->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  else if ((bin_type_a == kOrdinaryRecoBinSignal || bin_type_a == kOrdinaryRecoBinAll) && bin_type_b == kOrdinaryRecoBinSignal) {
    double err = univ.hist_reco_signal2d_->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );
    double err2 = err * err;
    return err2;
  }
  
  else {
    return 0.;
  }
  */

  /*
  else {
    // TODO: add proper handling of MC stat correlations between "ordinary
    // background" bins
    if ( cm_bin_a != cm_bin_b ) return 0.;
  }
  */
  /*
  // Include only background bins when evaluating the MC stat uncertainty
  // for the "ordinary background" bin type
  size_t num_true_bins = true_bins_.size();
  double err2 = 0.;
  for ( size_t tb = 0u; tb < num_true_bins; ++tb ) {
    const auto& tbin = true_bins_.at( tb );
    if ( bin_type_a == kOrdinaryRecoBinBkgd && tbin.type_ == kBackgroundTrueBin ) {
      double bkgd_err = univ.hist_2d_->GetBinError( tb + 1, reco_bin_a + 1 );
      err2 += bkgd_err * bkgd_err;
    }
    else if ( bin_type_a == kOrdinaryRecoBinSignal && tbin.type_ == kSignalTrueBin ) {
      double sig_err = univ.hist_2d_->GetBinError( tb + 1, reco_bin_a + 1 );
      err2 += sig_err * sig_err;
    }
  }

  return err2;
  
  //return 0;
  */
}

double ConstrainedCalculator::evaluate_data_stat_covariance( int cm_bin_a,
  int cm_bin_b, bool use_ext ) const
{
  ConstrainedCalculatorBinType bin_type_a, bin_type_b;
  int reco_bin_a = this->get_reco_bin_and_type( cm_bin_a, bin_type_a );
  int reco_bin_b = this->get_reco_bin_and_type( cm_bin_b, bin_type_b );

  const TH2D* d_hist = nullptr;
  if ( use_ext ) d_hist = data_hists2d_.at( NFT::kExtBNB ).get(); // EXT data
  else d_hist = data_hists2d_.at( NFT::kOnBNB ).get(); // BNB data
  // ROOT histograms use one-based bin indices, so I correct for that here
  double err = d_hist->GetBinError( reco_bin_a + 1, reco_bin_b + 1 );

  bool orbb = (bin_type_a == kOrdinaryRecoBinBkgd || bin_type_b == kOrdinaryRecoBinBkgd || bin_type_a == kOrdinaryRecoBinSignal || bin_type_b == kOrdinaryRecoBinSignal);

  if ( orbb ) {
    // Don't evaluate the data statistical uncertainty for the "ordinary
    // background" bins unless we're considering EXT events
    if ( !use_ext ) return 0.;
  }

  double err2 = err * err;
  return err2;
}

size_t ConstrainedCalculator::get_covariance_matrix_size() const {
  size_t cm_size = 3u * num_ordinary_reco_bins_;
  //size_t cm_size = 2u * num_ordinary_reco_bins_;
  cm_size += num_sideband_reco_bins_;
  return cm_size;
}

// Note that zero-based bin indices are assumed by this function. The
// return value is the reco bin index to use for computing covariances.
// The bin_type variable will be loaded with information about the kind
// of bin that should be assumed by the ConstrainedCalculator class.
int ConstrainedCalculator::get_reco_bin_and_type( int cm_bin,
  ConstrainedCalculatorBinType& bin_type ) const
{
  if ( cm_bin < 0 ) throw std::runtime_error( "Negative bin index encountered"
    " in ConstrainedCalculator::get_reco_bin_and_type()" );

  int bin = 0;

  // Bin numbers within the first set of ordinary reco bins do not need
  // to be remapped
  if ( cm_bin < num_ordinary_reco_bins_ ) {
    bin = cm_bin;
    bin_type = kOrdinaryRecoBinAll;
  }
  // A second copy of the ordinary reco bins follows immediately after.
  else if ( cm_bin >= num_ordinary_reco_bins_ && cm_bin < 2*num_ordinary_reco_bins_) {
    bin = cm_bin - num_ordinary_reco_bins_;
    // total, signal 
    bin_type = kOrdinaryRecoBinSignal;
    // total, background
    //bin_type = kOrdinaryRecoBinBkgd;
  }
  // A second copy of the ordinary reco bins follows immediately after.
  // All bins after the copies are assumed to be sideband reco bins.
  else {
    bin = cm_bin - 2*num_ordinary_reco_bins_;
    if (cm_bin >= 2*num_ordinary_reco_bins_ && cm_bin < 3*num_ordinary_reco_bins_) bin_type = kOrdinaryRecoBinBkgd;
    else bin_type = kSidebandRecoBinAll;
    //bin = cm_bin - num_ordinary_reco_bins_;
    //bin_type = kSidebandRecoBinAll;
  }
  return bin;
}

MeasuredEvents ConstrainedCalculator::get_measured_events(std::string type) const
{
  const int two_times_ord_bins = 2*num_ordinary_reco_bins_;
  const int three_times_ord_bins = 3*num_ordinary_reco_bins_;
  const auto& cv_univ = this->cv_universe();

  // First create vectors of the measured event counts and central-value
  // prediction in each of the sideband bins. Note here and elsewhere in this
  // function that ROOT histogram bin indices are one-based to allow for
  // underflow. The TMatrixD element indices, on the other hand, are zero-based.
  TMatrixD sideband_data( num_sideband_reco_bins_, 1 );
  TMatrixD sideband_mc_plus_ext( num_sideband_reco_bins_, 1 );

  TH1D* d_hist = data_hists_.at( NFT::kOnBNB ).get(); // BNB data
  TH1D* ext_hist = data_hists_.at( NFT::kExtBNB ).get(); // EXT data
  for ( int s = 0; s < num_sideband_reco_bins_; ++s ) {
    // Zero-based reco bin index
    int r = s + num_ordinary_reco_bins_;

    // Switch to using the one-based TH1D index when retrieving these values
    double bnb_events = d_hist->GetBinContent( r + 1 );
    double ext_events = ext_hist->GetBinContent( r + 1 );
    double cv_mc_events = cv_univ.hist_reco_->GetBinContent( r + 1 );

    sideband_data( s, 0 ) = bnb_events;
    sideband_mc_plus_ext( s, 0 ) = cv_mc_events + ext_events;
  }

  // Create the vector of measured event counts in the ordinary reco bins
  TMatrixD ordinary_data( num_ordinary_reco_bins_, 1 );
  for ( int r = 0; r < num_ordinary_reco_bins_; ++r ) {
    // Switch to using the one-based TH1D index when retrieving these values
    double bnb_events = d_hist->GetBinContent( r + 1 );

    ordinary_data( r, 0 ) = bnb_events;
  }

  // Now create the vector which stores the unconstrained prediction for both
  // signal+background, signal-only and background-only bins
  TMatrixD cv_pred_vec( three_times_ord_bins, 1 );
  for ( int r = 0; r < three_times_ord_bins; ++r ) {
  //TMatrixD cv_pred_vec( two_times_ord_bins, 1 );
  //for ( int r = 0; r < two_times_ord_bins; ++r ) {
    // This will automatically handle the signal+background versus
    // background-only bin definitions correctly. Recall that this
    // function takes a zero-based index.
    double cv_mc_events = this->evaluate_observable( cv_univ, r );

    // Also get the EXT event count for the reco bin of interest
    ConstrainedCalculatorBinType dummy_bin_type;
    int reco_bin_index = this->get_reco_bin_and_type( r, dummy_bin_type );

    // We need to use a one-based bin index to retrieve this value from the TH1D
    double ext_events = ext_hist->GetBinContent( reco_bin_index + 1 );

    cv_pred_vec( r, 0 ) = cv_mc_events + ext_events;
  }

  // All we need now to apply the sideband constraint are submatrices of
  // the total covariance matrix. Build it and pull out the blocks that
  // we need.
  auto cov_map_ptr = this->get_covariances();
  auto tot_cov_mat = cov_map_ptr->at( "total" ).get_matrix();

  // Zero-based indices for the covariance matrix elements describing the
  // ordinary reco bins (ob) and sideband reco bins (sb)
  
  int first_ob_cm_idx = 0;
  int last_ob_cm_idx = three_times_ord_bins - 1;
  int first_sb_cm_idx = three_times_ord_bins;
  int last_sb_cm_idx = three_times_ord_bins + num_sideband_reco_bins_ - 1;
  
  /*
  int first_ob_cm_idx = 0;
  int last_ob_cm_idx = two_times_ord_bins - 1;
  int first_sb_cm_idx = two_times_ord_bins;
  int last_sb_cm_idx = two_times_ord_bins + num_sideband_reco_bins_ - 1;
  */

  // Covariance matrix block that describes the ordinary reco bins
  TMatrixD ordinary_cov_mat = tot_cov_mat->GetSub( first_ob_cm_idx,
    last_ob_cm_idx, first_ob_cm_idx, last_ob_cm_idx );

  // Block that describes the sideband reco bins
  TMatrixD sideband_cov_mat = tot_cov_mat->GetSub( first_sb_cm_idx,
    last_sb_cm_idx, first_sb_cm_idx, last_sb_cm_idx );

  // Block that describes correlations between the sideband and ordinary bins.
  // This version uses rows for sideband bins and columns for ordinary bins.
  TMatrixD s_o_cov_mat = tot_cov_mat->GetSub( first_sb_cm_idx, last_sb_cm_idx,
    first_ob_cm_idx, last_ob_cm_idx );

  // Invert the sideband covariance matrix in preparation for applying
  // the sideband constraint
  auto inverse_sideband_cov_mat = invert_matrix( sideband_cov_mat );
 
  // We're ready. Apply the sideband constraint to the prediction vector first.
  TMatrixD sideband_data_mc_diff( sideband_data,
    TMatrixD::EMatrixCreatorsOp2::kMinus, sideband_mc_plus_ext );
  TMatrixD temp1( *inverse_sideband_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kMult, sideband_data_mc_diff );
  TMatrixD add_to( s_o_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kTransposeMult, temp1 );

  TMatrixD constr_cv_pred_vec( cv_pred_vec,
    TMatrixD::EMatrixCreatorsOp2::kPlus, add_to );

  // Now get the corresponding updated covariance matrix
  TMatrixD temp2( *inverse_sideband_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kMult, s_o_cov_mat );
  TMatrixD subtract_from( s_o_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kTransposeMult, temp2 );

  TMatrixD constr_ordinary_cov_mat( ordinary_cov_mat,
    TMatrixD::EMatrixCreatorsOp2::kMinus, subtract_from );

  /*
  // Pull out the block of the constrained covariance matrix that describes
  // the uncertainty on signal and background
  // total -- background case
  // constrained
  auto* constr_sig_plus_bkgd_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, num_ordinary_reco_bins_ - 1 ));
  auto* constr_bkgd_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ) );
  auto* constr_first_offdiag_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ));
  auto* constr_second_offdiag_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, num_ordinary_reco_bins_ - 1 ));
  auto* constr_sig_cov_mat = new TMatrixD( *constr_sig_plus_bkgd_cov_mat - *constr_first_offdiag_cov_mat - *constr_second_offdiag_cov_mat + *constr_bkgd_cov_mat);
  // unconstrained
  auto* unconstr_sig_plus_bkgd_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, num_ordinary_reco_bins_ - 1 ) );
  auto* unconstr_bkgd_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ) );
  auto* unconstr_first_offdiag_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ));
  auto* unconstr_second_offdiag_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, num_ordinary_reco_bins_ - 1 ));
  auto* unconstr_sig_cov_mat = new TMatrixD( *unconstr_sig_plus_bkgd_cov_mat - *unconstr_first_offdiag_cov_mat - *unconstr_second_offdiag_cov_mat + *unconstr_bkgd_cov_mat);
  
  // build the constrainted background + unconstrained signal convariance matrix
  auto* sig_plus_bkgd_cov_mat = new TMatrixD( *unconstr_sig_cov_mat + *constr_first_offdiag_cov_mat + *constr_second_offdiag_cov_mat +  *constr_bkgd_cov_mat);

  // unconstrained signal + background
  auto* unconstr_mc_plus_ext_vec = new TMatrixD(cv_pred_vec.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, 0 ));
  // unconstrained background
  auto* unconstr_bkgd_vec = new TMatrixD( cv_pred_vec.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, 0 ));
  // unconstrained signal
  auto* unconstr_sig_vec = new TMatrixD( *unconstr_mc_plus_ext_vec, TMatrixD::EMatrixCreatorsOp2::kMinus, *unconstr_bkgd_vec );
  
  // constrained signal + background
  auto* constr_mc_plus_ext_vec = new TMatrixD(constr_cv_pred_vec.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, 0 ));
  // constrained background
  auto* reco_bkgd_vec = new TMatrixD(constr_cv_pred_vec.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, 0 ));
  
  // unconstrained signal + constrained background
  auto* mc_plus_ext_vec = new TMatrixD( *unconstr_sig_vec, TMatrixD::EMatrixCreatorsOp2::kPlus, *reco_bkgd_vec );
 
  // --- testing ---
  // fully constrained 
  //sig_plus_bkgd_cov_mat = constr_sig_plus_bkgd_cov_mat;
  //mc_plus_ext_vec = constr_mc_plus_ext_vec;
  */
  
  // total -- signal case
  /*
  // constrained
  auto* constr_sig_plus_bkgd_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, num_ordinary_reco_bins_ - 1 ));
  auto* constr_sig_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ) );
  auto* constr_first_offdiag_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ));
  auto* constr_second_offdiag_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, num_ordinary_reco_bins_ - 1 ));
  auto* constr_bkgd_cov_mat = new TMatrixD( *constr_sig_plus_bkgd_cov_mat - *constr_first_offdiag_cov_mat - *constr_second_offdiag_cov_mat + *constr_sig_cov_mat);
  // unconstrained
  auto* unconstr_sig_plus_bkgd_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, num_ordinary_reco_bins_ - 1 ) );
  auto* unconstr_sig_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ) );
  auto* unconstr_first_offdiag_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ));
  auto* unconstr_second_offdiag_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, num_ordinary_reco_bins_ - 1 ));
  auto* unconstr_bkgd_cov_mat = new TMatrixD( *unconstr_sig_plus_bkgd_cov_mat - *unconstr_first_offdiag_cov_mat - *unconstr_second_offdiag_cov_mat + *unconstr_sig_cov_mat);

  // build the constrainted background + unconstrained signal convariance matrix
  auto* sig_plus_bkgd_cov_mat = new TMatrixD( *constr_bkgd_cov_mat + *constr_first_offdiag_cov_mat + *constr_second_offdiag_cov_mat + *unconstr_sig_cov_mat);

  // unconstrained signal + background
  auto* unconstr_mc_plus_ext_vec = new TMatrixD(cv_pred_vec.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, 0 ));
  // unconstrained signal
  auto* unconstr_sig_vec = new TMatrixD( cv_pred_vec.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, 0 ));
  // unconstrained background
  auto* unconstr_bkgd_vec = new TMatrixD( *unconstr_mc_plus_ext_vec, TMatrixD::EMatrixCreatorsOp2::kMinus, *unconstr_sig_vec );
  
  // constrained signal + background
  auto* constr_mc_plus_ext_vec = new TMatrixD(constr_cv_pred_vec.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, 0 ));
  // constrained signal
  auto* constr_sig_vec = new TMatrixD(constr_cv_pred_vec.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, 0 ));
  // constrained background
  auto* reco_bkgd_vec = new TMatrixD( *constr_mc_plus_ext_vec, TMatrixD::EMatrixCreatorsOp2::kMinus, *constr_sig_vec );
  
  // unconstrained signal + constrained background 
  auto* mc_plus_ext_vec = new TMatrixD( *unconstr_sig_vec, TMatrixD::EMatrixCreatorsOp2::kPlus, *reco_bkgd_vec );
  */

  
  // total -- signal -- background case
  // constrained
  auto* constr_sig_plus_bkgd_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, num_ordinary_reco_bins_ - 1 ));
  auto* constr_sig_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ) );
  auto* constr_bkgd_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( two_times_ord_bins, three_times_ord_bins - 1, two_times_ord_bins, three_times_ord_bins - 1 ) );
  // and constrained correlations between total and signal 
  auto* constr_first_offdiag_totalsignal_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ));
  auto* constr_second_offdiag_totalsignal_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, num_ordinary_reco_bins_ - 1 ));
  // constrainted correlations between signal and background
  auto* constr_first_offdiag_signalbackground_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, two_times_ord_bins, three_times_ord_bins - 1 ));
  auto* constr_second_offdiag_signalbackground_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( two_times_ord_bins, three_times_ord_bins - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ));
  // between total and background
  auto* constr_first_offdiag_totalbackground_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, two_times_ord_bins, three_times_ord_bins - 1 ));
  auto* constr_second_offdiag_totalbackground_cov_mat = new TMatrixD(constr_ordinary_cov_mat.GetSub( two_times_ord_bins, three_times_ord_bins - 1, 0, num_ordinary_reco_bins_ - 1 ));

  // pull out blocks that describe unconstrained uncertainty on signal and background
  // unconstrained
  auto* unconstr_sig_plus_bkgd_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, num_ordinary_reco_bins_ - 1 ) );
  auto* unconstr_sig_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ) );
  auto* unconstr_bkgd_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( two_times_ord_bins, three_times_ord_bins - 1, two_times_ord_bins, three_times_ord_bins - 1 ) );
  // and unconstrained correlations between total and signal 
  auto* unconstr_first_offdiag_totalsignal_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ));
  auto* unconstr_second_offdiag_totalsignal_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, num_ordinary_reco_bins_ - 1 ));
  // unconstrainted correlations between signal and background
  auto* unconstr_first_offdiag_signalbackground_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, two_times_ord_bins, three_times_ord_bins - 1 ));
  auto* unconstr_second_offdiag_signalbackground_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( two_times_ord_bins, three_times_ord_bins - 1, num_ordinary_reco_bins_, two_times_ord_bins - 1 ));
  // between total and background
  auto* unconstr_first_offdiag_totalbackground_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( 0, num_ordinary_reco_bins_ - 1, two_times_ord_bins, three_times_ord_bins - 1 ));
  auto* unconstr_second_offdiag_totalbackground_cov_mat = new TMatrixD(ordinary_cov_mat.GetSub( two_times_ord_bins, three_times_ord_bins - 1, 0, num_ordinary_reco_bins_ - 1 ));

  // build the constrainted background + unconstrained signal convariance matrix
  // using total-signal
  auto* sig_plus_bkgd_cov_mat = new TMatrixD( *unconstr_sig_cov_mat + *constr_first_offdiag_totalsignal_cov_mat + *constr_second_offdiag_totalsignal_cov_mat + *constr_bkgd_cov_mat );
  // using total-background
  //auto* sig_plus_bkgd_cov_mat = new TMatrixD( *unconstr_sig_cov_mat + *constr_first_offdiag_totalbackground_cov_mat + *constr_second_offdiag_totalbackground_cov_mat + *constr_bkgd_cov_mat );
  // using signal-background (not right!)
 // auto* sig_plus_bkgd_cov_mat = new TMatrixD( *unconstr_sig_cov_mat + *constr_first_offdiag_signalbackground_cov_mat + *constr_second_offdiag_signalbackground_cov_mat + *constr_bkgd_cov_mat );  
  
  // fully constrained
  //auto* sig_plus_bkgd_cov_mat = constr_sig_plus_bkgd_cov_mat;


  // Get the constrained background prediction column vector
  auto* reco_bkgd_vec = new TMatrixD(constr_cv_pred_vec.GetSub( two_times_ord_bins, three_times_ord_bins - 1, 0, 0 ));
  
  // Get the unconstrained signal prediction column vector
  auto* unconstr_sig_vec = new TMatrixD(cv_pred_vec.GetSub( num_ordinary_reco_bins_, two_times_ord_bins - 1, 0, 0 ));

  // Get the unconstrained background prediction column vector
  auto* unconstr_bkgd_vec = new TMatrixD(cv_pred_vec.GetSub( two_times_ord_bins, three_times_ord_bins - 1, 0, 0 ));

  // unconstrained signal + background
  auto* unconstr_mc_plus_ext_vec = new TMatrixD(cv_pred_vec.GetSub( 0, num_ordinary_reco_bins_ - 1, 0, 0 ));

  // constrained background + unconstrained signal
  auto* mc_plus_ext_vec = new TMatrixD( *unconstr_sig_vec, TMatrixD::EMatrixCreatorsOp2::kPlus, *reco_bkgd_vec );
  

  // build total
  auto* test_total = new TMatrixD( *unconstr_sig_cov_mat + *unconstr_first_offdiag_signalbackground_cov_mat + *unconstr_second_offdiag_signalbackground_cov_mat + *unconstr_bkgd_cov_mat );
  std::cout << "test total: " << std::sqrt( test_total->E2Norm() ) << std::endl;

  // Plot the matrix as a colz plot
  TCanvas* cCorrTotal = new TCanvas;
  const auto total_covmat_tmatrixd = util::TH2DToTMatrixD(*test_total);
  //auto total_corrmat = util::CovarianceMatrixToCorrelationMatrix( total_covmat_tmatrixd );

  int nX = total_covmat_tmatrixd.GetNrows();
  int nY = total_covmat_tmatrixd.GetNcols();
  auto total_corrmat_hist = util::TMatrixDToTH2D(total_covmat_tmatrixd, "", "", 0, nX, 0, nY);
  
  //total_corrmat_hist.GetZaxis()->SetRangeUser(-0.5, 1.0);
  total_corrmat_hist.SetTitle("Test Total");
  total_corrmat_hist.Draw("COLZ");
  gStyle->SetOptStat(0); // Add this line to remove the stats box

  TCanvas* cCorrTotal2 = new TCanvas;
  const auto total_covmat_tmatrixd2 = util::TH2DToTMatrixD(*unconstr_sig_plus_bkgd_cov_mat);
  //auto total_corrmat = util::CovarianceMatrixToCorrelationMatrix( total_covmat_tmatrixd );

  int nX2 = total_covmat_tmatrixd2.GetNrows();
  int nY2 = total_covmat_tmatrixd2.GetNcols();
  auto total_corrmat_hist2 = util::TMatrixDToTH2D(total_covmat_tmatrixd2, "", "", 0, nX2, 0, nY2);
  
  //total_corrmat_hist2.GetZaxis()->SetRangeUser(-0.5, 1.0);
  total_corrmat_hist2.SetTitle("unconstr_sig_plus_bkgd_cov_mat");
  total_corrmat_hist2.Draw("COLZ");
  gStyle->SetOptStat(0); // Add this line to remove the stats box



  std::cout << "Total unconstrained: " << std::sqrt( unconstr_sig_plus_bkgd_cov_mat->E2Norm() ) << ", Total constrained: " << std::sqrt( constr_sig_plus_bkgd_cov_mat->E2Norm() ) << ", Signal Unconstrained + Background Constrained: " << std::sqrt( sig_plus_bkgd_cov_mat->E2Norm() ) << std::endl;
  std::cout << "Signal unconstrained: " << std::sqrt( unconstr_sig_cov_mat->E2Norm() ) << ", Signal constrained: " << std::sqrt( constr_sig_cov_mat->E2Norm() ) << std::endl;
  std::cout << "Background unconstrained: " << std::sqrt( unconstr_bkgd_cov_mat->E2Norm() ) << ", Background constrained: " << std::sqrt( constr_bkgd_cov_mat->E2Norm() ) << std::endl;

  //std::cout << "Unconstrained off-diagonal 1: " << std::sqrt( unconstr_first_offdiag_cov_mat->E2Norm() ) << ", Unconstrained off-diagonal 2: " << std::sqrt( unconstr_second_offdiag_cov_mat->E2Norm() ) << std::endl;
  //std::cout << "Constrained off-diagonal 1: " << std::sqrt( constr_first_offdiag_cov_mat->E2Norm() ) << ", Constrained off-diagonal 2: " << std::sqrt( constr_second_offdiag_cov_mat->E2Norm() ) << std::endl;
  
  std::cout << "Constrained off-diagonal total-signal 1: " << std::sqrt( constr_first_offdiag_totalsignal_cov_mat->E2Norm() ) << ", Constrained off-diagonal total-signal 2: " << std::sqrt( constr_second_offdiag_totalsignal_cov_mat->E2Norm() ) << std::endl;
  std::cout << "Constrained off-diagonal total-background 1: " << std::sqrt( constr_first_offdiag_totalbackground_cov_mat->E2Norm() ) << ", Constrained off-diagonal total-background 2: " << std::sqrt( constr_second_offdiag_totalbackground_cov_mat->E2Norm() ) << std::endl;
  std::cout << "Constrained off-diagonal signal-background 1: " << std::sqrt( constr_first_offdiag_signalbackground_cov_mat->E2Norm() ) << ", Constrained off-diagonal signal-background 2: " << std::sqrt( constr_second_offdiag_signalbackground_cov_mat->E2Norm() ) << std::endl;
  std::cout << "Unconstrained off-diagonal total-signal 1: " << std::sqrt( unconstr_first_offdiag_totalsignal_cov_mat->E2Norm() ) << ", Unconstrained off-diagonal total-signal 2: " << std::sqrt( unconstr_second_offdiag_totalsignal_cov_mat->E2Norm() ) << std::endl;
  std::cout << "Unconstrained off-diagonal total-background 1: " << std::sqrt( unconstr_first_offdiag_totalbackground_cov_mat->E2Norm() ) << ", Unconstrained off-diagonal total-background 2: " << std::sqrt( unconstr_second_offdiag_totalbackground_cov_mat->E2Norm() ) << std::endl;
  std::cout << "Unconstrained off-diagonal signal-background 1: " << std::sqrt( unconstr_first_offdiag_signalbackground_cov_mat->E2Norm() ) << ", Unconstrained off-diagonal signal-background 2: " << std::sqrt( unconstr_second_offdiag_signalbackground_cov_mat->E2Norm() ) << std::endl;
  

  // Get the ordinary reco bin data column vector after subtracting the
  // constrained background prediction
  auto* reco_data_minus_bkgd = new TMatrixD( ordinary_data, TMatrixD::EMatrixCreatorsOp2::kMinus, *reco_bkgd_vec );
  // uncontrained background
  auto* unconstr_reco_data_minus_bkgd = new TMatrixD( ordinary_data, TMatrixD::EMatrixCreatorsOp2::kMinus, *unconstr_bkgd_vec );

  if (type == "default") {
    MeasuredEvents result( reco_data_minus_bkgd, reco_bkgd_vec,
      mc_plus_ext_vec, sig_plus_bkgd_cov_mat );

    return result;
  }
  if (type == "unconstr sig+bkgd") {
    MeasuredEvents result( unconstr_reco_data_minus_bkgd, unconstr_bkgd_vec,
      unconstr_mc_plus_ext_vec, unconstr_sig_plus_bkgd_cov_mat );

    return result;
  }
  else if (type == "constr bkgd") {
    MeasuredEvents result( reco_data_minus_bkgd, reco_bkgd_vec,
      reco_bkgd_vec, constr_bkgd_cov_mat );

    return result;
  }
  else if (type == "unconstr bkgd") {
    MeasuredEvents result( unconstr_reco_data_minus_bkgd, unconstr_bkgd_vec,
      unconstr_bkgd_vec, unconstr_bkgd_cov_mat );

    return result;
  }
  else {
    throw std::runtime_error( "ConstrainedCalculator::get_measured_events: Invalid type." );
  }
}
