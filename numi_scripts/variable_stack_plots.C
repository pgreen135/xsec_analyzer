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
#include "../utils.hh"

using NFT = NtupleFileType;

void scale_by_bin_width(SliceHistogram* pSlice)
{
    int num_slice_bins = pSlice->hist_->GetNbinsX();
    TMatrixD trans_mat( num_slice_bins, num_slice_bins );
    for ( int b = 0; b < num_slice_bins; ++b ) {
      const auto width = pSlice->hist_->GetBinWidth( b + 1 );
      // width *= other_var_width;
      trans_mat( b, b ) = 1 / width;
    }
    pSlice->transform(trans_mat);
}

void variable_stack_plots() {

  bool normaliseByBinWidth = true;

  auto* syst_ptr = new MCC9SystematicsCalculator(
    
    // Alternative Flux
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run3_original_flux.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run3_new_flux.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run3_reweighted_flux.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run3_new_flux_withData.root",

    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_test.root",


    // Full
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_sideband_full.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_sideband_full_noProton.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_extraLooseProtonSideband_noPionAngle.root",
    
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData.root",
    


    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_runDependentDetVars.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_runDependentAltCV.root",
    
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_fluggUncertainty.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_newDetVars.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withoutData.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_GenieFakeData.root",  // all slices
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_NuWroFakeData.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_NuWroFakeData_full_withNuWroGenie.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_FluggFakeData.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_FluggFakeData_FluggUncertainty.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_FluggFakeData_FluggUncertainty_alternate.root",

    // Sideband Test
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_sideband_total_looseProton.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_sideband_electron_angle_looseProtonTest.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_looseProtonSideband.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_strictProtonSideband.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_pi0Sideband.root",

    // AltFlux Test
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_standardFlux_test.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_altFlux_test.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_altFlux_asFakeData.root",

    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_withData_sideband_full.root",

    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_full_withData_sideband_Katrina.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run1_withData_sideband_Katrina.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run2_withData_sideband_Katrina.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run3_withData_sideband_Katrina.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run4_withData_sideband_Katrina.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run5_withData_sideband_Katrina.root",

    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_NuWroFakeData_sideband_full.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_NuWroFakeData.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_NuWroFakeData_sideband_full_noPi0.root",
    
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_GenieFakeData.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_combined_nuecc1pi_GenieFakeData_sideband_full.root",


    // Reweighted Flux [Some need re-running]
    // Total
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_withData.root", [OLD]
    // FHC
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_fhc_withData_withNuWroGenie.root",
    // RHC
    "/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_rhc_withData_withNuWroGenie.root",

    // Sideband
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_withData_sideband_withNuWroGenie.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_withData_sideband_withDelta2NpiFix.root",

    // NuWro Fake Data
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_NuWroFakeData.root", [OLD]
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_NuWroFakeData_sideband_alternate_withNuWroGenie.root",
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_NuWroFakeData_sideband_alternate_75pc_scaling.root",

    // Genie Fake Data
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_GenieFakeData.root", [OLD]

    // Reweight fake data (run 3 only)
    //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_ReweightFakeData.root", [OLD]


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

  auto* sb_ptr = new SliceBinning( "../nuecc1pi_slice_config.txt" );
  //auto* sb_ptr = new SliceBinning( "../pionangle_slice_config.txt" );
  //auto* sb_ptr = new SliceBinning( "../nuecc1pi_sideband_1_slice_config.txt" );
  //auto* sb_ptr = new SliceBinning( "../nuecc1pi_sideband_2_slice_config.txt" );
  

  //auto* sb_ptr = new SliceBinning( "../nuecc1pi_slice_config_oldflux.txt" );
 
  //auto* sb_ptr = new SliceBinning( "../nuecc1pi_sideband_2_Katrina_slice_config.txt" );

  auto& sb = *sb_ptr;

  unsigned int sl_idx = 4;
  const auto& slice = sb.slices_.at( sl_idx ); // only considering single slice at a time

  // -- confusion matrix -- 
  const auto cv_hist_2d = syst.get_cv_hist_2d();

  size_t start = std::numeric_limits<size_t>::max();
  size_t stop = 0;
  for (const auto& entry : slice.bin_map_) {
    const auto& set = entry.second;
    if (set.size() != 1) {
        throw std::runtime_error("Error: set in bin_map_ has more or less than 1 entry");
    }
    start = std::min(start, *set.begin());
    stop = std::max(stop, *set.rbegin());
  }

  TMatrixD cv_hist_2d_slice(stop - start + 1, stop - start + 1);
  for (int i = start; i <= stop; i++) {
      for (int j = start; j <= stop; j++) {
          cv_hist_2d_slice(i - start, j - start) = cv_hist_2d->GetBinContent(i+1, j+1);
      }
  }

  const auto cv_confusion_mat = util::CountsToConfusionMatrix(cv_hist_2d_slice, "row");  

  // Create a TH2D histogram
  const int num_bins = slice.hist_->GetNbinsX();
  TH2D *cv_confusion_hist = new TH2D(("cv_confusion slice "+std::to_string(sl_idx)).c_str(),"",
                          num_bins, 0, num_bins,
                          num_bins, 0, num_bins);

  // Fill histogram with values from TMatrixD
  for (int i = 0; i < cv_confusion_mat.GetNrows(); i++) {
      for (int j = 0; j < cv_confusion_mat.GetNcols(); j++) {
          cv_confusion_hist->SetBinContent(i+1, j+1, cv_confusion_mat(i, j));
      }
  }

  TCanvas* c_conf_1 = new TCanvas(("c_conf_1 slice "+std::to_string(sl_idx)).c_str(), "c confusion matrix", 800, 600);
  cv_confusion_hist->SetStats(false);
  cv_confusion_hist->SetMinimum(0);
  cv_confusion_hist->SetMaximum(1.0);
  // gStyle->SetPalette(util::CreateWhiteToBlueColorPalette(20));
  gStyle->SetPalette();

  // Set x and y axis labels to bin numbers
  for (int i = 1; i <= cv_confusion_hist->GetNbinsX(); i++) {
      cv_confusion_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
  }
  for (int i = 1; i <= cv_confusion_hist->GetNbinsY(); i++) {
      cv_confusion_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
  }
  // Increase the font size of the axis labels
  cv_confusion_hist->GetXaxis()->SetLabelSize(0.05);
  cv_confusion_hist->GetYaxis()->SetLabelSize(0.05);

  cv_confusion_hist->GetXaxis()->SetTitle("Bin Index");
  cv_confusion_hist->GetYaxis()->SetTitle("Bin Index");
  cv_confusion_hist->SetTitle("Confusion Matrix: Opening Angle");

  cv_confusion_hist->Draw("COLZ");

  // Fill histogram with slice data
  for (int i = 0; i < num_bins; i++) {
      for (int j = 0; j < num_bins; j++) {
          double bin_content = cv_hist_2d_slice(i, j);
          TLatex* latex = new TLatex(cv_confusion_hist->GetXaxis()->GetBinCenter(i+1), cv_confusion_hist->GetYaxis()->GetBinCenter(j+1), Form("#splitline{%.1f%%}{(%.1f)}", 100*cv_confusion_hist->GetBinContent(i+1, j+1), bin_content));
          latex->SetTextFont(42);
          latex->SetTextSize(0.02);
          latex->SetTextAlign(22);
          latex->Draw();
      }
  }

  c_conf_1->SaveAs(("cv_confusion_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + ".pdf").c_str());

  // -- confusion matrix -- 


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

  if (normaliseByBinWidth) scale_by_bin_width(slice_ext);
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

    if (normaliseByBinWidth) scale_by_bin_width(temp_slice_mc);

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
  
  TCanvas* c1 = new TCanvas("", "", 1080, 1080);
  TPad *upperPad = new TPad("upperPad", "Upper Pad", 0.01, 0.25, 0.99, 0.99);
  TPad *lowerPad = new TPad("lowerPad", "Lower Pad", 0.01, 0.01, 0.99, 0.24);
  upperPad->Draw();
  lowerPad->Draw();

  upperPad->cd();  // Switch to the upper pad
  gPad->SetBottomMargin(0.0125);
  gPad->SetTopMargin(0.14);

  slice_bnb->hist_->SetLineColor( kBlack );
  slice_bnb->hist_->SetLineWidth( 3 );
  slice_bnb->hist_->SetMarkerStyle( kFullCircle );
  slice_bnb->hist_->SetMarkerSize( 0.8 );
  slice_bnb->hist_->SetStats( false );

  if (normaliseByBinWidth) scale_by_bin_width(slice_bnb);
  if (normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext);

  double ymax = std::max( slice_bnb->hist_->GetMaximum(),
    slice_mc_plus_ext->hist_->GetMaximum() ) * 1.8;
  slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., ymax );
  //slice_bnb->hist_->GetYaxis()->SetRangeUser( 0., 20 );
  slice_bnb->hist_->SetTitle("");
  
  slice_bnb->hist_->GetXaxis()->SetLabelOffset(999); // Hide X-axis labels
  slice_bnb->hist_->GetXaxis()->SetTitleOffset(999); // Hide X-axis labels
  slice_bnb->hist_->GetXaxis()->SetTickLength(0.01);

  slice_bnb->hist_->Draw( "e" );

  slice_pred_stack->Draw( "hist same" );

  slice_mc_plus_ext->hist_->SetLineWidth( 3 ); // 3

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
  char labelText5[100];
  sprintf(labelText1, "FHC + RHC");
  //sprintf(labelText2, "7.77e+20 POT"); // FHC
  //sprintf(labelText2, "11.08e+20 POT"); // RHC
  sprintf(labelText2, "1.885e+21 POT"); // FHC + RHC
  
  //sprintf(labelText1, "Run 5 FHC");
  //sprintf(labelText2, "2.7973e+20 POT");
  //sprintf(labelText2, "3.859e+20 POT");
  //sprintf(labelText2, "5.003e+20 POT");
  //sprintf(labelText2, "4.958e+20 POT");
  //sprintf(labelText2, "2.231e+20 POT");
  //sprintf(labelText2, "1.3845e+21 POT");
  //sprintf(labelText2, "NuWro Fake Data");
  //sprintf(labelText2, "Genie Fake Data");
  //sprintf(labelText2, "New Flux");
  sprintf(labelText3, "#chi^{2} = %.2f / %d Bins", chi2_result.chi2_, chi2_result.num_bins_);
  sprintf(labelText4, "p-value = %.2f", chi2_result.p_value_);
  sprintf(labelText5, "%.1f #sigma", TMath::Sqrt( TMath::ChisquareQuantile( 1 - chi2_result.p_value_, 1 ) ));
  label.SetTextSize(0.04);
  //label.DrawLatex(0.65, 0.80, labelText1);
  //label.DrawLatex(0.65, 0.75, labelText2);
  //label.DrawLatex(0.65, 0.70, labelText3);
  //label.DrawLatex(0.65, 0.65, labelText4);
  //label.DrawLatex(0.65, 0.60, labelText5);

  label.DrawLatex(0.175, 0.80, labelText1);
  label.DrawLatex(0.175, 0.75, labelText2);
  label.DrawLatex(0.175, 0.70, labelText3);
  label.DrawLatex(0.175, 0.65, labelText4);
  //label.DrawLatex(0.175, 0.60, labelText5);

  // draw legend
  leg->Draw("Same");

  gPad->RedrawAxis();

  std::cout << "Data events: " << slice_bnb->hist_.get()->Integral() << std::endl;
  std::cout << "EXT events: " << slice_ext->hist_.get()->Integral() << std::endl;
  std::cout << "MC+EXT events: " << slice_mc_plus_ext->hist_->Integral() << std::endl;

  // Create and draw the ratio plot
  TH1D *h_ratio = (TH1D*)slice_bnb->hist_.get()->Clone("h_ratio");
  TH1D *h_ratio_error = (TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_ratio_error");
  TH1D *h_ratio_values = (TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_ratio_values");
  
  for(int i=0 ; i <= h_ratio_values->GetNbinsX() ; i++){
    h_ratio_values->SetBinError(i, 0);
  }
  
  //h_ratio->Divide(slice_mc_plus_ext->hist_.get());
  h_ratio->Divide(h_ratio_values);

  h_ratio_error->Divide(h_ratio_values);

  lowerPad->cd();  // Switch to the lower pad
  gPad->SetBottomMargin(0.35);

  h_ratio->GetYaxis()->SetRangeUser( 0, 2 );
  h_ratio->GetYaxis()->SetTitle("Ratio");
  h_ratio->GetYaxis()->SetTitleSize(0.12);
  h_ratio->GetYaxis()->SetTitleOffset(0.325);
  h_ratio->GetYaxis()->SetLabelSize(0.12);
  h_ratio->GetYaxis()->SetNdivisions(505);

  h_ratio->SetBit(TH1::kNoTitle);

  h_ratio->GetXaxis()->SetLabelSize(0.12);  // Adjust X-axis label size
  h_ratio->GetXaxis()->SetLabelOffset(0.01);
  h_ratio->GetXaxis()->SetTitleSize(0.12);
  h_ratio->GetXaxis()->SetTitleOffset(1.0);
  h_ratio->GetXaxis()->SetTickLength(0.03);
  h_ratio->Draw();
  h_ratio_error->Draw("e2 same");


  std::cout << "Bin 1 error: " << h_ratio->GetBinError(1) << std::endl;

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
    "MCstats", "EXTstats", "BNBstats", "NuWroGenie"
  };
  // show detvars
  /*
  const std::vector< std::string > cov_mat_keys = {"total", "detVar_total",
    "detVarLYdown", "detVarLYrayl", "detVarLYatten", "detVarRecomb2", "detVarSCE", "detVarWMAngleXZ", "detVarWMAngleYZ",
    "detVarWMX", "detVarWMYZ", "detVarNumu" 
  };
  */
  /*    
  // show beamline uncertainties
  const std::vector< std::string > cov_mat_keys = {"total", "flux_beamline",
  "flux_Horn_2kA", "flux_Horn1_x_3mm", "flux_Horn1_y_3mm",
  "flux_Beam_spot_1_1mm", "flux_Beam_spot_1_5mm", "flux_Horn2_x_3mm", "flux_Horn2_y_3mm",
  "flux_Horns_0mm_water", "flux_Horns_2mm_water", "flux_Beam_shift_x_1mm", 
  "flux_Beam_shift_y_1mm", "flux_Target_z_7mm" }; 
  */
  // show xsec uncertainties
  
  //const std::vector< std::string > cov_mat_keys = {"total", "xsec_multi", "xsec_unisim", "xsec_Theta_Delta2Npi" };
  

  // Loop over the various systematic uncertainties
  int color = 0;
  for ( const auto& pair : matrix_map ) {

    const auto& key = pair.first;
    const auto& cov_matrix = pair.second;

    std::cout << key << std::endl;

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

    std::cout << key << std::endl;

    if ( color <= 9 ) ++color;
    if ( color == 5 ) ++color;
    if ( color >= 10 ) color += 10;

    slice_for_syst->hist_->SetLineColor( color );
    slice_for_syst->hist_->SetLineWidth( 4 );

    if (key == "BNBstats") {
      slice_for_syst->hist_->SetLineStyle( 9 );
      slice_for_syst->hist_->SetLineWidth( 5 );
      slice_for_syst->hist_->SetLineColor( kBlack );
    }
    if (key == "NuWroGenie") slice_for_syst->hist_->SetLineColor( kOrange );
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
  total_frac_err_hist->SetLineWidth( 5 );
  total_frac_err_hist->SetTitle("Opening Angle");
  total_frac_err_hist->GetYaxis()->SetTitle("Fractional Uncertainty");
  total_frac_err_hist->Draw( "hist" );

  stringstream ss; ss.precision(1);
  ss << std::fixed << "Total";
  //ss << ": " << total_frac_err_hist->GetBinContent( 1 )*100. << "%";
  string label_str = ss.str();

  lg2->AddEntry( total_frac_err_hist, label_str.c_str(), "l" );

  std::cout << "HERE" << std::endl;

  for ( auto& pair : frac_uncertainty_hists ) {
    const auto& name = pair.first;
    TH1* hist = pair.second;
    // We already plotted the "total" one above
    if ( name == "total" ) continue;

    // refine labels
    std::string name_alt = "";
    if ( name == "BNBstats") name_alt = "NuMI Stat"; 
    if ( name == "EXTstats") name_alt = "EXT Stat";
    if ( name == "MCstats") name_alt = "MC Stat";
    if ( name == "detVar_total") name_alt = "DetVar Total";
    if ( name == "flux") name_alt = "Flux Hadronic";
    if ( name == "flux_beamline") name_alt = "Flux Beamline";
    if ( name == "flux_Flugg") name_alt = "Flux Flugg";
    if ( name == "numTargets") name_alt = "Num Targets";
    if ( name == "xsec_total") name_alt = "X-Sec Total";
    if ( name == "reint") name_alt = "Reinteractions";
    if ( name == "dirtNorm") name_alt = "Dirt Norm";

    stringstream ss; ss.precision(1);
    ss << std::fixed; 
    if (name_alt != "") ss << name_alt; 
    else ss << name;
    //ss << ": " << hist->GetBinContent( 1 )*100. << "%";
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
  variable_stack_plots();
  return 0;
}
