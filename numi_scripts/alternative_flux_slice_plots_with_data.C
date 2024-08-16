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

void scale_by_bin_width(SliceHistogram *pSlice) {
    int num_slice_bins = pSlice->hist_->GetNbinsX();
    TMatrixD trans_mat(num_slice_bins, num_slice_bins);
    for (int b = 0; b < num_slice_bins; ++b)
    {
        const auto width = pSlice->hist_->GetBinWidth(b + 1);
        // width *= other_var_width;
        trans_mat(b, b) = 1 / width;
    }
    pSlice->transform(trans_mat);
}

// ---------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------

void slice_plots() {

    bool normaliseByBinWidth = true;

    // ******************************************************************************************************************
    // Configure the input files ****************************************************************************************
    // ******************************************************************************************************************

    const bool using_fake_data = false;

    //auto* sb_ptr = new SliceBinning( "../nuecc1pi_slice_config.txt" );
    auto* sb_ptr = new SliceBinning( "../nueccinclusive_slice_config.txt" );
    auto& sb = *sb_ptr;

    // original flux
    auto *syst_ptr = new MCC9SystematicsCalculator(
        //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run3_original_flux.root",
        "/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nueccinclusive_run3_original_flux.root",
        "../systcalc.conf" );
    
    auto& syst = *syst_ptr;

    // get the original flux reco ext & mc + ext histograms
    TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
    TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
    TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
    reco_mc_plus_ext_hist->SetDirectory( nullptr );
    reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );

    // get the map of covariance matrices
    const auto* matrix_map_ptr = syst.get_covariances().release();
    const auto& matrix_map = *matrix_map_ptr;
    
    // alternative flux
    auto *syst_ptr_alt = new MCC9SystematicsCalculator(
        //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nuecc1pi_run3_new_flux.root",
        "/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output/univmake_output_nueccinclusive_run3_new_flux.root",
        "../systcalc.conf" );
    
    auto& syst_alt = *syst_ptr_alt;

    // get the alternative flux reco ext & mc + ext histograms
    TH1D* reco_ext_hist_alt =  syst_alt.data_hists_.at( NFT::kExtBNB ).get();
    TH1D* reco_mc_plus_ext_hist_alt = dynamic_cast< TH1D* >(reco_ext_hist_alt->Clone("reco_mc_plus_ext_hist_alt") );
    reco_mc_plus_ext_hist_alt->SetDirectory( nullptr );
    reco_mc_plus_ext_hist_alt->Add( syst_alt.cv_universe().hist_reco_.get() );

    // get the map of covariance matrices
    const auto* matrix_map_ptr_alt = syst_alt.get_covariances().release();
    const auto& matrix_map_alt = *matrix_map_ptr_alt;

    // ******************************************************************************************************************
    // Create slice plots ***********************************************************************************************
    // ******************************************************************************************************************

    size_t sl_idx = 0;  // can replace with loop here
    
    const auto &slice = sb.slices_.at(sl_idx);

    // get histograms
    SliceHistogram* slice_bnb = SliceHistogram::make_slice_histogram(
            *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );

    SliceHistogram* slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
            *reco_mc_plus_ext_hist, slice, &matrix_map.at("total") );

    SliceHistogram* slice_mc_plus_ext_alt = SliceHistogram::make_slice_histogram(
            *reco_mc_plus_ext_hist_alt, slice, &matrix_map_alt.at("total") );

    // Get chi2
    const auto chi2 = slice_bnb->get_chi2(*slice_mc_plus_ext);
    const auto chi2_alt = slice_bnb->get_chi2(*slice_mc_plus_ext_alt);

    if(normaliseByBinWidth) scale_by_bin_width(slice_bnb);
    if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext);
    if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext_alt);

    TCanvas* c = new TCanvas("", "", 1080, 1080);
    TPad *upperPad = new TPad("upperPad", "Upper Pad", 0.01, 0.25, 0.99, 0.99);
    TPad *lowerPad = new TPad("lowerPad", "Lower Pad", 0.01, 0.01, 0.99, 0.24);
    upperPad->Draw();
    lowerPad->Draw();

    upperPad->cd();  // Switch to the upper pad
    gPad->SetBottomMargin(0.0125);
    gPad->SetTopMargin(0.14);

    // set the color and line thickness of the histograms
    slice_bnb->hist_->SetLineColor(kBlack);
    slice_bnb->hist_->SetLineWidth(3);
    slice_bnb->hist_->SetMarkerStyle( kFullCircle );
    slice_bnb->hist_->SetMarkerSize( 0.8 );
    slice_bnb->hist_->SetStats( false );

    slice_mc_plus_ext->hist_->SetLineColor(kBlue);
    slice_mc_plus_ext->hist_->SetLineWidth(3);

    slice_mc_plus_ext_alt->hist_->SetLineColor(kRed);
    slice_mc_plus_ext_alt->hist_->SetLineWidth(3);

    slice_mc_plus_ext->hist_->SetMinimum(0.0);
    slice_mc_plus_ext_alt->hist_->SetMinimum(0.0);

    slice_mc_plus_ext->hist_->SetTitle("");

    slice_mc_plus_ext->hist_->GetYaxis()->SetRangeUser(0, slice_mc_plus_ext->hist_->GetMaximum() * 1.8);
    slice_mc_plus_ext->hist_->GetYaxis()->SetTitleSize(0.045);
    //slice_mc_plus_ext->hist_->GetYaxis()->SetTitleOffset(0.325);

    slice_mc_plus_ext->hist_->GetXaxis()->SetLabelOffset(999); // Hide X-axis labels
    slice_mc_plus_ext->hist_->GetXaxis()->SetTitleOffset(999); // Hide X-axis labels
    slice_mc_plus_ext->hist_->GetXaxis()->SetTickLength(0.01);

    // draw histograms
    slice_mc_plus_ext->hist_->DrawCopy( "hist" );
    slice_mc_plus_ext->hist_->SetFillColor(kBlue);
    slice_mc_plus_ext->hist_->SetFillStyle(3004);
    slice_mc_plus_ext->hist_->Draw("e2 same");

    slice_mc_plus_ext_alt->hist_->DrawCopy("hist same");
    slice_mc_plus_ext_alt->hist_->SetFillColor(kRed);
    slice_mc_plus_ext_alt->hist_->SetFillStyle(3005);
    slice_mc_plus_ext_alt->hist_->Draw("e2 same");

    slice_bnb->hist_->Draw("e same");

    // add legend
    TLegend* legendSlice = new TLegend(0.1,0.7,0.9,0.86);

    // Convert chi2 values to strings
    std::stringstream chi2SS; 
    chi2SS << std::setprecision(3) << "#chi^{2} = " << chi2.chi2_ << " / " << std::setprecision(1) << chi2.num_bins_ << " Bins, p = " << std::setprecision(2) <<  chi2.p_value_;
    std::string chi2Str = chi2SS.str();

    std::stringstream chi2AltSS; 
    chi2AltSS << std::setprecision(3) << "#chi^{2} = " << chi2_alt.chi2_ << " / " << std::setprecision(1) << chi2_alt.num_bins_ << " Bins, p = " << std::setprecision(2) << chi2_alt.p_value_;
    std::string chi2AltStr = chi2AltSS.str();

    legendSlice->AddEntry(slice_mc_plus_ext->hist_.get(), ("Original Flux: " + chi2Str).c_str(), "l");
    legendSlice->AddEntry(slice_mc_plus_ext_alt->hist_.get(), ("New Flux: " + chi2AltStr).c_str(), "l");
    legendSlice->Draw();

    // create and draw labels
    TLatex label;
    label.SetTextAlign(12); // Set text alignment (left-aligned)
    label.SetNDC(); // Set position in normalized coordinates
    char labelText1[100];
    char labelText2[100];
    sprintf(labelText1, "Run 3 RHC");
    sprintf(labelText2, "5.003e+20 POT");
    label.SetTextSize(0.04);
    label.DrawLatex(0.65, 0.61, labelText1);
    label.DrawLatex(0.65, 0.55, labelText2);
    //label.DrawLatex(0.175, 0.61, labelText1);
    //label.DrawLatex(0.175, 0.55, labelText2);

    // draw ratio
    lowerPad->cd();  // Switch to the lower pad
    gPad->SetBottomMargin(0.35);

    // Original flux
    TH1D *h_ratio = (TH1D*)slice_bnb->hist_.get()->Clone("h_ratio");
    TH1D *h_ratio_error = (TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_ratio_error");
    TH1D *h_ratio_values = (TH1D*)slice_mc_plus_ext->hist_.get()->Clone("h_ratio_values");

    for(int i=0 ; i <= h_ratio_values->GetNbinsX() ; i++){
        h_ratio_values->SetBinError(i, 0);
    }

    h_ratio->Divide(h_ratio_values);
    h_ratio_error->Divide(h_ratio_values);

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

    h_ratio->SetLineColor(kBlue);
    h_ratio->SetLineWidth(1);
    h_ratio->SetMarkerColor(kBlue);
    h_ratio->SetMarkerSize( 0.4 );

    h_ratio->Draw();
    h_ratio_error->Draw("e2 same");

    // New Flux
    TH1D *h_ratio_new = (TH1D*)slice_bnb->hist_.get()->Clone("h_ratio");
    TH1D *h_ratio_new_error = (TH1D*)slice_mc_plus_ext_alt->hist_.get()->Clone("h_ratio_error");
    TH1D *h_ratio_new_values = (TH1D*)slice_mc_plus_ext_alt->hist_.get()->Clone("h_ratio_values");

    for(int i=0 ; i <= h_ratio_new_values->GetNbinsX() ; i++){
        h_ratio_new_values->SetBinError(i, 0);
    }

    h_ratio_new->Divide(h_ratio_new_values);
    h_ratio_new_error->Divide(h_ratio_new_values);

    h_ratio_new->SetLineColor(kRed);
    h_ratio_new->SetLineWidth(1);
    h_ratio_new->SetMarkerColor(kRed);
    h_ratio_new->SetMarkerSize( 0.4 );

    h_ratio_new->Draw("same");
    h_ratio_new_error->Draw("e2 same");

    // draw line
    /*
    TLine *l = new TLine(h_ratio->GetXaxis()->GetXmin(), 1, h_ratio->GetXaxis()->GetXmax(), 1);
    l->SetLineColor(kBlack);
    l->SetLineStyle(9);
    l->SetLineWidth(2);
    l->Draw("same");
    */

}

int alternative_flux_slice_plots_with_data() {
    slice_plots();
    std::cout << "------------All Done------------" << std::endl;
    return 0;
}