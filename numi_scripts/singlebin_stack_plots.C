// Standard library includes
#include <algorithm>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"

// STV analysis includes
#include "../FilePropertiesManager.hh"
#include "../MCC9SystematicsCalculator.hh"
#include "../PlotUtils.hh"
#include "../SliceBinning.hh"
#include "../SliceHistogram.hh"

using NFT = NtupleFileType;

void singlebin_stack_plots() {

  auto* syst_ptr = new MCC9SystematicsCalculator(
  
    // Retrained BDTs
    "/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_fhc_singlebin_withdata.root",
    // Fake Data
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_fhc_singlebin_NuWroFakeData_nueOnly_withPPFX.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_fhc_singleBin_NuWroFakeData_NuWroGenieUncertainty.root", 
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_fhc_singlebin_GenieFakeData_fakeDataWeights.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_rhc_singlebin_GenieFakeData_fakeDataWeights.root",
    // MC only
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_fhc_singlebin_withoutData.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_rhc_singlebin_withoutData.root",
    // no data
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_fhc_noData_withNuWroUncertainty.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_rhc_noData_withNuWroUncertainty.root",
    // Data
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_fhc_singlebin_withData.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_rhc_singlebin_withData.root", 
    "../systcalc.conf" );
  auto& syst = *syst_ptr;

  // Get access to the relevant histograms owned by the SystematicsCalculator
  // object. These contain the reco bin counts that we need to populate the
  // slices below.
  TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
  TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();  
  TH2D* category_hist = syst.cv_universe().hist_categ_.get();

  // Total MC+EXT prediction in reco bin space. Start by getting EXT.

  TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
  reco_mc_plus_ext_hist->SetDirectory( nullptr );

  // Add in the CV MC prediction
  reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );
  
  // Keys are covariance matrix types, values are CovMatrix objects that
  // represent the corresponding matrices
  auto* matrix_map_ptr = syst.get_covariances().release();
  auto& matrix_map = *matrix_map_ptr;

  auto* sb_ptr = new SliceBinning( "../singlebin_slice_config.txt" );
  auto& sb = *sb_ptr;

  const auto& slice = sb.slices_.at( 0 ); // only considering single slice

  // We now have all of the reco bin space histograms that we need as input.
  // Use them to make new histograms in slice space.
  SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
    *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

  SliceHistogram* slice_ext = SliceHistogram::make_slice_histogram(
    *reco_ext_hist, slice, &matrix_map.at("EXTstats") );

  SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
    *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

  auto chi2_result = slice_bnb->get_chi2( *slice_mc_plus_ext );
  std::cout << "\u03C7\u00b2 = "
    << chi2_result.chi2_ << '/' << chi2_result.num_bins_ << " bins,"
    << " p-value = " << chi2_result.p_value_ << '\n';

  // Build a stack of categorized central-value MC predictions plus the
  // extBNB contribution in slice space
  const auto& eci = EventCategoryInterpreter::Instance();
  eci.set_ext_histogram_style( slice_ext->hist_.get() );

  THStack* slice_pred_stack = new THStack( "mc+ext", "" );
  
  slice_pred_stack->Add( slice_ext->hist_.get() ); // extBNB

  const auto& cat_map = eci.label_map();

  // Legend
  TLegend *leg = new TLegend(0.1,0.88,0.9,0.99);
  // leg->SetTextFont(132);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(6);

  // Go in reverse so that signal ends up on top. Note that this index is
  // one-based to match the ROOT histograms
  
  int cat_bin_index = cat_map.size();
  
  for ( auto iter = cat_map.crbegin(); iter != cat_map.crend(); ++iter )
  {
    EventCategory cat = iter->first;

    TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
      cat_bin_index, cat_bin_index );
    temp_mc_hist->SetDirectory( nullptr );

    SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
      *temp_mc_hist, slice  );

    eci.set_mc_histogram_style( cat, temp_slice_mc->hist_.get() );

    slice_pred_stack->Add( temp_slice_mc->hist_.get() );

    std::string cat_col_prefix = "MC" + std::to_string( cat );

    std::cout << "Event catagory = " << cat << ": " << eci.label(cat) << ", entries: " << temp_slice_mc->hist_.get()->Integral() << std::endl;

    --cat_bin_index;
  }
  
  // Second loop to construct legend in desired order
  cat_bin_index = 1;
  for ( auto iter = cat_map.begin(); iter != cat_map.end(); ++iter )
  {
    EventCategory cat = iter->first;

    TH1D* temp_mc_hist = category_hist->ProjectionY( "temp_mc_hist",
      cat_bin_index, cat_bin_index );
    temp_mc_hist->SetDirectory( nullptr );

    SliceHistogram* temp_slice_mc = SliceHistogram::make_slice_histogram(
      *temp_mc_hist, slice  );

    eci.set_mc_histogram_style( cat, temp_slice_mc->hist_.get() );

    leg->AddEntry(temp_slice_mc->hist_.get(), eci.label(cat).c_str(), "f");

    ++cat_bin_index;
  }

  leg->AddEntry(slice_ext->hist_.get(), "EXT", "f");
  

  TCanvas* c1 = new TCanvas;
  c1->SetTopMargin(0.14);

  slice_bnb->hist_->SetLineColor( kBlack );
  slice_bnb->hist_->SetLineWidth( 3 );
  slice_bnb->hist_->SetMarkerStyle( kFullCircle );
  slice_bnb->hist_->SetMarkerSize( 0.8 );
  slice_bnb->hist_->SetStats( false );
  double ymax = std::max( slice_bnb->hist_->GetMaximum(),
    slice_mc_plus_ext->hist_->GetMaximum() ) * 1.07;
  slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., 100); //// ymax );
  slice_bnb->hist_->SetTitle("");

  slice_bnb->hist_->Draw( "e" );

  slice_pred_stack->Draw( "hist same" );

  slice_mc_plus_ext->hist_->SetLineWidth( 3 );

  slice_mc_plus_ext->hist_->SetMarkerColor(kBlack);
  slice_mc_plus_ext->hist_->SetLineColor(kBlack);
  slice_mc_plus_ext->hist_->DrawCopy( "hist same" );

  slice_mc_plus_ext->hist_->SetFillColor(kBlack);
  slice_mc_plus_ext->hist_->SetFillStyle(3004);
  slice_mc_plus_ext->hist_->Draw("e2 same");

  slice_bnb->hist_->Draw( "same e" );

  // Create the label text with the value
  TLatex label;
  label.SetTextAlign(12); // Set text alignment (left-aligned)
  label.SetNDC(); // Set position in normalized coordinates
  char labelText1[100];
  char labelText2[100];
  char labelText3[100];
  char labelText4[100];
  sprintf(labelText1, "RHC");
  //sprintf(labelText2, "7.766e+20 POT");
  sprintf(labelText2, "11.082e+20 POT");
  sprintf(labelText3, "#chi^{2} = %.2f / %d Bins", chi2_result.chi2_, chi2_result.num_bins_);
  sprintf(labelText4, "p-value = %.2f", chi2_result.p_value_);
  label.SetTextSize(0.04);
  label.DrawLatex(0.7, 0.80, labelText1);
  label.DrawLatex(0.7, 0.75, labelText2);
  //label.DrawLatex(0.7, 0.70, labelText3);
  //label.DrawLatex(0.25, 0.65, labelText4);


  // draw legend
  leg->Draw("Same");

  gPad->RedrawAxis();

  std::cout << "Data events: " << slice_bnb->hist_.get()->Integral() << std::endl;
  std::cout << "EXT events: " << slice_ext->hist_.get()->Integral() << std::endl;
  std::cout << "MC+EXT events: " << slice_mc_plus_ext->hist_->Integral() << std::endl;

  // Get the binning and axis labels for the current slice by cloning the
  // (empty) histogram owned by the Slice object
  TH1* slice_hist = dynamic_cast< TH1* >(
    slice.hist_->Clone("slice_hist") );

  slice_hist->SetDirectory( nullptr );

  // Keys are labels, values are fractional uncertainty histograms
  auto* fr_unc_hists = new std::map< std::string, TH1* >();
  auto& frac_uncertainty_hists = *fr_unc_hists;

  // Show fractional uncertainties computed using these covariance matrices
  // in the ROOT plot. All configured fractional uncertainties will be
  // included in the output pgfplots file regardless of whether they appear
  // in this vector.
  const std::vector< std::string > cov_mat_keys = { "total",
    "detVar_total", "flux", "flux_beamline", "reint", "xsec_total", "POT", "numTargets", "dirtNorm",
    "MCstats", "EXTstats", "NuWroGenie"
  };
  // show detvars
  //const std::vector< std::string > cov_mat_keys = {"total", "detVar_total",
  //  "detVarLYdown", "detVarLYrayl", "detVarLYatten", "detVarRecomb2", "detVarSCE", "detVarWMAngleXZ", "detVarWMAngleYZ",
  //  "detVarWMX", "detVarWMYZ", "detVarNumu" 
  //};  
  // show beamline uncertainties
  //const std::vector< std::string > cov_mat_keys = {"total", "flux_beamline",
  //"flux_Horn_2kA", "flux_Horn1_x_3mm", "flux_Horn1_y_3mm",
  //"flux_Beam_spot_1_1mm", "flux_Beam_spot_1_5mm", "flux_Horn2_x_3mm", "flux_Horn2_y_3mm",
  //"flux_Horns_0mm_water", "flux_Horns_2mm_water", "flux_Beam_shift_x_1mm", 
  //"flux_Beam_shift_y_1mm", "flux_Target_z_7mm" }; 

  // Loop over the various systematic uncertainties
  int color = 0;
  for ( const auto& pair : matrix_map ) {

    const auto& key = pair.first;
    const auto& cov_matrix = pair.second;

    SliceHistogram* slice_for_syst = SliceHistogram::make_slice_histogram(
      *reco_mc_plus_ext_hist, slice, &cov_matrix );

    // The SliceHistogram object already set the bin errors appropriately
    // based on the slice covariance matrix. Just change the bin contents
    // for the current histogram to be fractional uncertainties. Also set
    // the "uncertainties on the uncertainties" to zero.
    // TODO: revisit this last bit, possibly assign bin errors here
    for ( const auto& bin_pair : slice.bin_map_ ) {
      int global_bin_idx = bin_pair.first;
      double y = slice_for_syst->hist_->GetBinContent( global_bin_idx );
      double err = slice_for_syst->hist_->GetBinError( global_bin_idx );
      double frac = 0.;
      if ( y > 0. ) frac = err / y;

      slice_for_syst->hist_->SetBinContent( global_bin_idx, frac );
      slice_for_syst->hist_->SetBinError( global_bin_idx, 0. );
    }

    // Check whether the current covariance matrix name is present in
    // the vector defined above this loop. If it isn't, don't bother to
    // plot it, and just move on to the next one.
    auto cbegin = cov_mat_keys.cbegin();
    auto cend = cov_mat_keys.cend();
    auto iter = std::find( cbegin, cend, key );
    if ( iter == cend ) continue;

    frac_uncertainty_hists[ key ] = slice_for_syst->hist_.get();

    if ( color <= 9 ) ++color;
    if ( color == 5 ) ++color;
    if ( color >= 10 ) color += 10;

    slice_for_syst->hist_->SetLineColor( color );
    slice_for_syst->hist_->SetLineWidth( 3 );
  }

  TCanvas* c2 = new TCanvas;
  TLegend* lg2 = new TLegend( 0.2, 0.7, 0.9, 0.9 );
  lg2->SetTextAlign(12);
  lg2->SetNColumns(3);

  auto* total_frac_err_hist = frac_uncertainty_hists.at( "total" );
  total_frac_err_hist->SetStats( false );
  total_frac_err_hist->GetYaxis()->SetRangeUser( 0., 0.5);
  //total_frac_err_hist->GetMaximum() * 1.05 );
  total_frac_err_hist->SetLineColor( kBlack );
  total_frac_err_hist->SetLineWidth( 3 );
  total_frac_err_hist->SetTitle("FHC");
  total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");
  total_frac_err_hist->Draw( "hist" );

  stringstream ss; ss.precision(1);
  ss << std::fixed << "Total";
  ss << ": " << total_frac_err_hist->GetBinContent( 1 )*100. << "%";
  string label_str = ss.str();

  lg2->AddEntry( total_frac_err_hist, label_str.c_str(), "l" );

  for ( auto& pair : frac_uncertainty_hists ) {
    const auto& name = pair.first;
    TH1* hist = pair.second;
    // We already plotted the "total" one above
    if ( name == "total" ) continue;

    // refine labels
    std::string name_alt = ""; 
    if ( name == "EXTstats") name_alt = "EXT Stat";
    if ( name == "MCstats") name_alt = "MC Stat";
    if ( name == "detVar_total") name_alt = "DetVar Total";
    if ( name == "flux") name_alt = "Flux Hadronic";
    if ( name == "flux_beamline") name_alt = "Flux Beamline";
    if ( name == "numTargets") name_alt = "Num Targets";
    if ( name == "xsec_total") name_alt = "X-Sec Total";
    if ( name == "reint") name_alt = "Reinteractions";
    if ( name == "dirtNorm") name_alt = "Dirt Norm";

    stringstream ss; ss.precision(1);
    ss << std::fixed; 
    if (name_alt != "") ss << name_alt; 
    else ss << name;
    ss << ": " << hist->GetBinContent( 1 )*100. << "%";
    string label_str = ss.str();
    lg2->AddEntry( hist, label_str.c_str(), "l" );
    hist->Draw( "same hist" );

    std::cout << name << " frac err in bin #1 = "
      << hist->GetBinContent( 1 )*100. << "%\n";
  }

  lg2->Draw( "same" );

  std::cout << "Total frac error in bin #1 = "
    << total_frac_err_hist->GetBinContent( 1 )*100. << "%\n";
}

int main() {
  singlebin_stack_plots();
  return 0;
}
