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
#include "../ConstrainedCalculator.hh"
#include "../PlotUtils.hh"
#include "../SliceBinning.hh"
#include "../SliceHistogram.hh"

#include "../utils.hh"

using NFT = NtupleFileType;

void scale_by_bin_width(SliceHistogram *pSlice)
{
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

void slice_plots()
{
    bool normaliseByBinWidth = true;
    bool normaliseByUnconstrained = true;

    // ******************************************************************************************************************
    // Configure the input files ****************************************************************************************
    // ******************************************************************************************************************

    const bool using_fake_data = true;

    auto* sb_ptr = new SliceBinning( "../nuecc1pi_slice_config.txt" );
    //auto* sb_ptr = new SliceBinning( "../nuecc1pi_slice_config_noOpeningAngle.txt" );
    auto& sb = *sb_ptr;

    // ******************************************************************************************************************
    // Get the unconstrained objects ************************************************************************************
    // ******************************************************************************************************************
   // Get the object constaining all the selection and uncertainty information
    auto *syst_ptr = new MCC9SystematicsCalculator(

        // Data
        //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_withData_sideband.root",
        
        // NuWro
        //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_NuWroFakeData_sideband_alternate.root",
        "/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_NuWroFakeData_sideband_alternate_75pc_scaling.root",

        // Genie
        //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_GenieFakeData_sideband.root",

        "../systcalc.conf");
    auto &syst = *syst_ptr;

    
    // Get the measured events from the systematics calculator
    const auto meas = syst.get_measured_events();

    /*
    const auto &reco_signal = meas.reco_signal_; // Background-subtracted data event counts in the ordinary reco bins
    const auto &reco_bkgd = meas.reco_bkgd_; // Background that was subtracted from each reco bin to form the signal measurement
    const auto &reco_mc_plus_ext = meas.reco_mc_plus_ext_; // Total MC+EXT prediction in each reco bin
    const auto &reco_covmat = meas.cov_matrix_; // Covariance matrix for MC+EXT prediction
    const auto *reco_signal_plus_bkgd = new TMatrixD(*reco_signal, TMatrixD::EMatrixCreatorsOp2::kPlus, *reco_bkgd);

    // Get the unconstrained reco ext & mc + ext histograms
    TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();
    TH1D* reco_mc_plus_ext_hist = dynamic_cast< TH1D* >(
    reco_ext_hist->Clone("reco_mc_plus_ext_hist") );
    reco_mc_plus_ext_hist->SetDirectory( nullptr );
    reco_mc_plus_ext_hist->Add( syst.cv_universe().hist_reco_.get() );
    */
    

    // Get the map of covariance matrices
    const auto* matrix_map_ptr = syst.get_covariances().release();
    const auto& matrix_map = *matrix_map_ptr;

    // ******************************************************************************************************************
    // Get the constrained objects **************************************************************************************
    // ******************************************************************************************************************
    // Get the object constaining all the selection and uncertainty information for the constrained calculator
    auto* syst_ptr_constr = new ConstrainedCalculator(
        
        // Data
        //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_withData_sideband.root",
        
        // NuWro
        //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_NuWroFakeData_sideband_alternate.root",
        "/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_NuWroFakeData_sideband_alternate_75pc_scaling.root",

        // Genie
        //"/Users/patrick/Documents/MicroBooNE/CrossSections/NuePiXSec_Analysis/XSecAnalyzer/univmake_output_reweightedPPFX/univmake_output_nuecc1pi_combined_GenieFakeData_sideband.root",

        "../systcalc.conf" );
    auto &syst_constr = *syst_ptr_constr;

    // Get the measured events from the constrained calculator
    // Total unconstrained 
    const auto meas_unconstr = syst_constr.get_measured_events("unconstr sig+bkgd");
    const auto &reco_mc_plus_ext = meas_unconstr.reco_mc_plus_ext_;
    const auto &reco_covmat = meas_unconstr.cov_matrix_;

    // Total constrained (default)
    const auto meas_constr = syst_constr.get_measured_events();
    const auto &reco_signal_constr = meas_constr.reco_signal_;
    const auto &reco_bkgd_constr = meas_constr.reco_bkgd_;
    const auto &reco_mc_plus_ext_constr = meas_constr.reco_mc_plus_ext_;
    const auto &reco_covmat_constr = meas_constr.cov_matrix_;
    //const auto *reco_signal_plus_bkgd_constr = new TMatrixD(*reco_signal_constr, TMatrixD::EMatrixCreatorsOp2::kPlus, *reco_bkgd_constr);

    // Background Unconstrained
    const auto meas_bkgd_unconstr = syst_constr.get_measured_events("unconstr bkgd");
    const auto &reco_bkgd_unconstr = meas_bkgd_unconstr.reco_mc_plus_ext_;
    const auto &reco_bkdg_covmat_unconstr = meas_bkgd_unconstr.cov_matrix_;

    // Background Constrained
    const auto meas_bkgd_constr = syst_constr.get_measured_events("constr bkgd");
    //const auto &reco_bkgd_constr = meas_bkgd_constr.reco_mc_plus_ext_;
    const auto &reco_bkdg_covmat_constr = meas_bkgd_constr.cov_matrix_;

    // Get the map of covariance matrices (identical to unconstrained map)
    const auto* matrix_map_constr_ptr = syst_constr.get_covariances().release();
    const auto& matrix_map_constr = *matrix_map_constr_ptr;


    // ******************************************************************************************************************
    // Get the total covariance matrices for the unconstrained **********************************************************
    // ******************************************************************************************************************
    // Get the map of covariance matrices (identical to unconstrained map)
    //const auto* matrix_map_ptr = syst_constr.get_covariances().release();
    //const auto& matrix_map = *matrix_map_ptr;
    const auto total_covmat = matrix_map.at("total").cov_matrix_.get();
     
    // ******************************************************************************************************************
    // Get the total covariance matrices for the unconstrained case as a correlation matrix *****************************
    // ******************************************************************************************************************
    // Plot the matrix as a colz plot
    TCanvas* cCorrTotal = new TCanvas;
    const auto total_covmat_tmatrixd = util::TH2DToTMatrixD(*total_covmat);
    auto total_corrmat = util::CovarianceMatrixToCorrelationMatrix( total_covmat_tmatrixd );

    int nX = total_corrmat.GetNrows();
    int nY = total_corrmat.GetNcols();
    auto total_corrmat_hist = util::TMatrixDToTH2D(total_corrmat, "", "", 0, nX, 0, nY);
    
    total_corrmat_hist.GetZaxis()->SetRangeUser(-0.5, 1.0);
    total_corrmat_hist.Draw("COLZ");
    gStyle->SetOptStat(0); // Add this line to remove the stats box

    // Draw dotted lines to separate the sideband and the signal region
    // use the size of the TMatrixD total_covmat_tmatrixd and draw a vertical and a horizontal line between half the bins
    const auto nSignalBins = total_covmat->GetXaxis()->GetNbins() / 3;  

    TLine *l0a = new TLine(0, nSignalBins, nSignalBins, nSignalBins);
    l0a->Draw();
    l0a->SetLineStyle(2);
    TLine *l0b = new TLine(2*nSignalBins, nSignalBins, 3*nSignalBins, nSignalBins);
    l0b->Draw();
    TLine *l0c = new TLine(0, 2*nSignalBins, nSignalBins, 2*nSignalBins);
    l0c->Draw();
    
    TLine *l0d = new TLine(nSignalBins, 2*nSignalBins, nSignalBins, 3*nSignalBins);
    l0d->Draw();
    TLine *l0e = new TLine(2*nSignalBins, 0, 2*nSignalBins, nSignalBins);
    l0e->Draw();

    TLine *l1 = new TLine(0, nSignalBins, 2*nSignalBins, nSignalBins);
    l1->SetLineStyle(2);
    l1->Draw();
    TLine *l2 = new TLine(nSignalBins, 0, nSignalBins, 2*nSignalBins);
    l2->SetLineStyle(2);
    l2->Draw();
    TLine *l3 = new TLine(nSignalBins, 2*nSignalBins, 3*nSignalBins, 2*nSignalBins);
    l3->SetLineStyle(2);
    l3->Draw();
    TLine *l4 = new TLine(2*nSignalBins, nSignalBins, 2*nSignalBins, 3*nSignalBins);
    l4->SetLineStyle(2);
    l4->Draw();

    // Add labels to the x-axis
    TLatex *ltxX1 = new TLatex();
    ltxX1->SetTextSize(0.03);
    ltxX1->SetTextAlign(22);
    ltxX1->DrawLatex(0.5*nSignalBins, -3.0, "Signal bins");
    TLatex *ltxX2 = new TLatex();
    ltxX2->SetTextSize(0.03);
    ltxX2->SetTextAlign(22);
    ltxX2->DrawLatex(1.5*nSignalBins, -3.0, "#pi^{0} Sideband bins");
    TLatex *ltxX3 = new TLatex();
    ltxX3->SetTextSize(0.03);
    ltxX3->SetTextAlign(22);
    ltxX3->DrawLatex(2.5*nSignalBins, -3.0, "#nu_{e} Np Sideband bins");


    // Add labels to the y-axis
    TLatex *ltxY1 = new TLatex();
    ltxY1->SetTextSize(0.03);
    ltxY1->SetTextAlign(22);
    ltxY1->SetTextAngle(90); // Rotate the text by 90 degrees
    ltxY1->DrawLatex(-3.0, 0.5*nSignalBins, "Signal bins");
    TLatex *ltxY2 = new TLatex();
    ltxY2->SetTextSize(0.03);
    ltxY2->SetTextAlign(22);
    ltxY2->SetTextAngle(90); // Rotate the text by 90 degrees
    ltxY2->DrawLatex(-3.0, 1.5*nSignalBins, "#pi^{0} Sideband bins");
    TLatex *ltxY3 = new TLatex();
    ltxY3->SetTextSize(0.03);
    ltxY3->SetTextAlign(22);
    ltxY3->SetTextAngle(90); // Rotate the text by 90 degrees
    ltxY3->DrawLatex(-3.0, 2.5*nSignalBins, "#nu_{e} Np Sideband bins");
    
    /*
    // Draw dotted lines to separate the sideband and the signal region
    // use the size of the TMatrixD total_covmat_tmatrixd and draw a vertical and a horizontal line between half the bins
    const auto nSignalBins = total_covmat->GetXaxis()->GetNbins() / 2;
    TLine *l_0_1 = new TLine(0, nSignalBins, 2*nSignalBins, nSignalBins);
    l_0_1->SetLineStyle(2);
    l_0_1->Draw();
    TLine *l_0_2 = new TLine(nSignalBins, 0, nSignalBins, 2*nSignalBins);
    l_0_2->SetLineStyle(2);
    l_0_2->Draw();

    // Add labels to the x-axis
    TLatex *ltx_0_X1 = new TLatex();
    ltx_0_X1->SetTextSize(0.03);
    ltx_0_X1->SetTextAlign(22);
    ltx_0_X1->DrawLatex(0.5*nSignalBins, -0.12*nSignalBins, "Signal reco bins");
    TLatex *ltx_0_X2 = new TLatex();
    ltx_0_X2->SetTextSize(0.03);
    ltx_0_X2->SetTextAlign(22);
    ltx_0_X2->DrawLatex(1.5*nSignalBins, -0.12*nSignalBins, "Sideband reco bins");

    // Add labels to the y-axis
    TLatex *ltx_0_Y1 = new TLatex();
    ltx_0_Y1->SetTextSize(0.03);
    ltx_0_Y1->SetTextAlign(22);
    ltx_0_Y1->SetTextAngle(90); // Rotate the text by 90 degrees
    ltx_0_Y1->DrawLatex(-0.12*nSignalBins, 0.5*nSignalBins, "Signal reco bins");
    TLatex *ltx_0_Y2 = new TLatex();
    ltx_0_Y2->SetTextSize(0.03);
    ltx_0_Y2->SetTextAlign(22);
    ltx_0_Y2->SetTextAngle(90); // Rotate the text by 90 degrees
    ltx_0_Y2->DrawLatex(-0.12*nSignalBins, 1.5*nSignalBins, "Sideband reco bins");
    */


    // Add a title to the plot
    TPaveText *ptTotalCorr = new TPaveText(0.1, 0.94, 0.9, 0.98, "brNDC");
    ptTotalCorr->SetBorderSize(0);
    ptTotalCorr->AddText("Unconstrained Total Correlation Matrix");
    ptTotalCorr->Draw();

    std::string out_pdf_name_corr_total = "plots/corr_matrix_total_unconstrained.pdf";
    cCorrTotal->SaveAs(out_pdf_name_corr_total.c_str());
    

    // ******************************************************************************************************************
    // Loop over slices *************************************************************************************************
    // ******************************************************************************************************************
    for (size_t sl_idx = 0u; sl_idx < sb.slices_.size(); ++sl_idx)
    {
        const auto &slice = sb.slices_.at(sl_idx);

        // Make histograms in slice space.
        // Signal + Background
        SliceHistogram *slice_mc_plus_ext = SliceHistogram::make_slice_histogram(
            *reco_mc_plus_ext, slice, reco_covmat.get());

        //SliceHistogram *slice_reco_signal_plus_bkgd = SliceHistogram::make_slice_histogram(
        //   *reco_signal_plus_bkgd_constr, slice, matrix_map.at("BNBstats").get_matrix().get());

        TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();
        SliceHistogram* slice_reco_signal_plus_bkgd = SliceHistogram::make_slice_histogram(
            *reco_bnb_hist, slice, &matrix_map.at("BNBstats") );
      
        // Make constrained versions of the histograms
        SliceHistogram *slice_mc_plus_ext_constr= SliceHistogram::make_slice_histogram(
            *reco_mc_plus_ext_constr, slice, reco_covmat_constr.get());

        // Background only
        SliceHistogram *slice_bkgd_unconstr= SliceHistogram::make_slice_histogram(
            *reco_bkgd_unconstr, slice, reco_bkdg_covmat_unconstr.get());
       
        SliceHistogram *slice_bkgd_constr= SliceHistogram::make_slice_histogram(
            *reco_bkgd_constr, slice, reco_bkdg_covmat_constr.get());

        // sum covariance matrix
        TMatrixD* unconstrained_covmat = reco_covmat.get();  
        TMatrixD* constrained_covmat = reco_covmat_constr.get();

        // extract this slice
        int num_slice_bins = slice.bin_map_.size();

        TH2D* covmat_hist_unconstr = new TH2D( "covmat_hist_unconstr", "covariance; slice bin;"
          " slice bin; covariance", num_slice_bins, 0., num_slice_bins,
          num_slice_bins, 0., num_slice_bins );
        covmat_hist_unconstr->SetDirectory( nullptr );
        covmat_hist_unconstr->SetStats( false );

        TH2D* covmat_hist_constr = new TH2D( "covmat_hist_constr", "covariance; slice bin;"
          " slice bin; covariance", num_slice_bins, 0., num_slice_bins,
          num_slice_bins, 0., num_slice_bins );
        covmat_hist_constr->SetDirectory( nullptr );
        covmat_hist_constr->SetStats( false );

        // We're ready. Populate the new covariance matrix using the elements
        // of the one for the reco bin space
        double unconstrained_sum_sq = 0.;
        double constrained_sum_sq = 0.;

        for ( const auto& pair_a : slice.bin_map_ ) {
          // Global slice bin index
          int sb_a = pair_a.first;
          // Set of reco bins that correspond to slice bin sb_a
          const auto& rb_set_a = pair_a.second;
          for ( const auto& pair_b : slice.bin_map_ ) {
            int sb_b = pair_b.first;
            const auto& rb_set_b = pair_b.second;

            double cov = 0.;
            double cov_constr = 0.;
            for ( const auto& rb_m : rb_set_a ) {
              for ( const auto& rb_n : rb_set_b ) {
                // The TMatrixD object uses zero-based indices
                cov += unconstrained_covmat->operator()( rb_m, rb_n );
                cov_constr += constrained_covmat->operator()( rb_m, rb_n );
              } // reco bin index m
            } // reco bin index n
            covmat_hist_unconstr->SetBinContent( sb_a, sb_b, cov );
            covmat_hist_constr->SetBinContent( sb_a, sb_b, cov_constr );

            unconstrained_sum_sq += std::pow(cov,2);
            constrained_sum_sq += std::pow(cov_constr,2);
          } // slice bin index b
        } // slice bin index a

        TCanvas* c_covmat_unconstr = new TCanvas("", "", 1080, 1080);
        covmat_hist_unconstr->Draw("COLZ");

        // write to file
        std::string out_pdf_name3 = "plots/" + std::to_string( sl_idx ) + "_covmat_unconstr_";
        if ( sl_idx < 10 ) out_pdf_name3 += "0";
        out_pdf_name3 += std::to_string( sl_idx )+".pdf";
        c_covmat_unconstr->SaveAs(out_pdf_name3.c_str());
        delete c_covmat_unconstr;

        TCanvas* c_covmat_constr = new TCanvas("", "", 1080, 1080);
        covmat_hist_constr->Draw("COLZ");

        std::string out_pdf_name4 = "plots/" + std::to_string( sl_idx ) + "_covmat_constr_";
        if ( sl_idx < 10 ) out_pdf_name4 += "0";
        out_pdf_name4 += std::to_string( sl_idx )+".pdf";
        c_covmat_constr->SaveAs(out_pdf_name4.c_str());
        delete c_covmat_constr;

        std::cout << "Slice " << sl_idx << ": unconstrained covmat sum = " << std::sqrt(unconstrained_sum_sq) 
            << ", constrained covmat sum = " << std::sqrt(constrained_sum_sq) << std::endl;

        /*
        // Create plot for the unconstrained covariance matrix
        TCanvas* cCov = new TCanvas;
        auto data_corrmat = util::CovarianceMatrixToCorrelationMatrix( *reco_covmat );
        
        // Plot the TMatrixD data_covmat as a colz plot
        int nX_uncons = data_corrmat.GetNrows();
        int nY_uncons = data_corrmat.GetNcols();
        auto data_corrmat_hist = util::TMatrixDToTH2D(data_corrmat, "", "", 0, nX_uncons, 0, nY_uncons);
        
        data_corrmat_hist.GetZaxis()->SetRangeUser(-0.5, 1.0);
        data_corrmat_hist.Draw("COLZ");

        gStyle->SetOptStat(0); // Add this line to remove the stats box

        // Add a title to the plot
        TPaveText *pt = new TPaveText(0.1, 0.94, 0.9, 0.98, "brNDC");
        pt->AddText("Unconstrained Correlation Matrix");
        pt->Draw();

        std::string out_pdf_name_cov = "plots/corr_matrix_";
        if ( sl_idx < 10 ) out_pdf_name_cov += "0";
        out_pdf_name_cov += std::to_string( sl_idx )+"_unconstrained.pdf";
        cCov->SaveAs(out_pdf_name_cov.c_str());

        //  Create plot for the constrained covariance matrix
        TCanvas* cCovConstr = new TCanvas;
        auto data_corrmat_constr = util::CovarianceMatrixToCorrelationMatrix( *reco_covmat_constr );
        
        // Plot the TMatrixD data_covmat as a colz plot   
        int nX_cons = data_corrmat_constr.GetNrows();
        int nY_cons = data_corrmat_constr.GetNcols();
        auto data_corrmat_constr_hist = util::TMatrixDToTH2D(data_corrmat_constr, "", "", 0, nX_cons, 0, nY_cons);
        
        data_corrmat_constr_hist.GetZaxis()->SetRangeUser(-0.5, 1.0);
        data_corrmat_constr_hist.Draw("COLZ");
        
        gStyle->SetOptStat(0); // Add this line to remove the stats box

        // Add a title to the plot
        TPaveText *ptConstr = new TPaveText(0.1, 0.94, 0.9, 0.98, "brNDC");
        ptConstr->AddText("Constrained Correlation Matrix");
        ptConstr->Draw();

        std::string out_pdf_name_cov_constr = "plots/corr_matrix_";
        if ( sl_idx < 10 ) out_pdf_name_cov_constr += "0";
        out_pdf_name_cov_constr += std::to_string( sl_idx )+"_constrained.pdf";
        cCovConstr->SaveAs(out_pdf_name_cov_constr.c_str());
        */

        // Get chi2
        const auto chi2 = slice_reco_signal_plus_bkgd->get_chi2(*slice_mc_plus_ext);
        const auto chi2Const = slice_reco_signal_plus_bkgd->get_chi2(*slice_mc_plus_ext_constr);

        if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext);
        if(normaliseByBinWidth) scale_by_bin_width(slice_reco_signal_plus_bkgd);
        if(normaliseByBinWidth) scale_by_bin_width(slice_mc_plus_ext_constr);
        //if(normaliseByBinWidth) scale_by_bin_width(slice_reco_signal_plus_bkgd_constr);

        if(normaliseByBinWidth) scale_by_bin_width(slice_bkgd_unconstr);
        if(normaliseByBinWidth) scale_by_bin_width(slice_bkgd_constr);

        TCanvas* c = new TCanvas("", "", 1080, 1080);
        TPad *upperPad = new TPad("upperPad", "Upper Pad", 0.01, 0.25, 0.99, 0.99);
        TPad *lowerPad = new TPad("lowerPad", "Lower Pad", 0.01, 0.01, 0.99, 0.24);
        upperPad->Draw();
        lowerPad->Draw();

        upperPad->cd();  // Switch to the upper pad
        gPad->SetBottomMargin(0.0125);
        gPad->SetTopMargin(0.14);
        
        // Set the color and line thickness of the histograms
        slice_mc_plus_ext->hist_->SetLineColor(kBlue);
        slice_mc_plus_ext->hist_->SetLineWidth(1);
        slice_reco_signal_plus_bkgd->hist_->SetLineColor(kBlack);
        slice_reco_signal_plus_bkgd->hist_->SetLineWidth(3);

        slice_mc_plus_ext_constr->hist_->SetLineColor(kRed);
        slice_mc_plus_ext_constr->hist_->SetLineWidth(1);
        
        // Set the minimum value of the y-axis to 0
        slice_mc_plus_ext->hist_->SetMinimum(0.0);
        slice_reco_signal_plus_bkgd->hist_->SetMinimum(0.0);
        slice_mc_plus_ext_constr->hist_->SetMinimum(0.0);

        slice_reco_signal_plus_bkgd->hist_->GetYaxis()->SetRangeUser(0, slice_reco_signal_plus_bkgd->hist_->GetMaximum() * 1.5);

        slice_reco_signal_plus_bkgd->hist_->GetXaxis()->SetLabelOffset(999); // Hide X-axis labels
        slice_reco_signal_plus_bkgd->hist_->GetXaxis()->SetTitleOffset(999); // Hide X-axis labels
        slice_reco_signal_plus_bkgd->hist_->GetXaxis()->SetTickLength(0.01);

        // Draw the histograms
        slice_reco_signal_plus_bkgd->hist_->Draw("E");

        //slice_mc_plus_ext->hist_->Draw("E hist same");
        slice_mc_plus_ext->hist_->DrawCopy( "hist same" );
        slice_mc_plus_ext->hist_->SetFillColor(kBlue);
        slice_mc_plus_ext->hist_->SetFillStyle(3004);
        slice_mc_plus_ext->hist_->Draw("e2 same");

        slice_mc_plus_ext_constr->hist_->DrawCopy("hist same");
        slice_mc_plus_ext_constr->hist_->SetFillColor(kRed);
        slice_mc_plus_ext_constr->hist_->SetFillStyle(3005);
        slice_mc_plus_ext_constr->hist_->Draw("e2 same");
        
        // Convert chi2 values to strings
        std::stringstream chi2SS; 
        chi2SS << std::setprecision(3) << " - Chi2: " << chi2.chi2_;
        if (slice_mc_plus_ext->hist_->GetNbinsX() > 1) chi2SS << ", p-value: " << std::setprecision(2) << chi2.p_value_;
        std::string chi2Str = chi2SS.str();

        std::stringstream chi2ConstSS; 
        chi2ConstSS << std::setprecision(3) << " - Chi2: " << chi2Const.chi2_;
        if (slice_mc_plus_ext->hist_->GetNbinsX() > 1) chi2ConstSS << ", p-value: " << std::setprecision(2) << chi2Const.p_value_;
        std::string chi2ConstStr = chi2ConstSS.str();

        // Print fractional uncertainties in each bin
        for (int bin = 1; bin <= slice_mc_plus_ext->hist_->GetNbinsX(); bin++) {
            std::cout << "Bin " << bin << " -- ";
            // Unconstrained
            std::cout << "Unconstrained: " << slice_mc_plus_ext->hist_->GetBinError(bin) / slice_mc_plus_ext->hist_->GetBinContent(bin);
            // Constrained
            std::cout << ", Constrained: " << slice_mc_plus_ext_constr->hist_->GetBinError(bin) / slice_mc_plus_ext_constr->hist_->GetBinContent(bin) << std::endl;
        }
        
        // Add a legend
        TLegend* legendSlice = new TLegend(0.1,0.70,0.9,0.9);
        legendSlice->AddEntry(slice_mc_plus_ext->hist_.get(), ("MC+EXT (unconstr)" + chi2Str).c_str(), "l");
        legendSlice->AddEntry(slice_mc_plus_ext_constr->hist_.get(), ("MC+EXT (constr)" + chi2ConstStr).c_str(), "l");
        legendSlice->AddEntry(slice_reco_signal_plus_bkgd->hist_.get(), "NuWro Fake Data", "l");
        legendSlice->Draw();

        
        // create and draw ratio plot
        gPad->RedrawAxis();
        
        lowerPad->cd();  // Switch to the lower pad
        gPad->SetBottomMargin(0.35);

        // Unconstrained
        TH1D *h_ratio = (TH1D*)slice_reco_signal_plus_bkgd->hist_.get()->Clone("h_ratio");
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

        h_ratio->Draw();
        h_ratio_error->Draw("e2 same");

        // Constrained
        TH1D *h_ratio_constr = (TH1D*)slice_reco_signal_plus_bkgd->hist_.get()->Clone("h_ratio");
        TH1D *h_ratio_constr_error = (TH1D*)slice_mc_plus_ext_constr->hist_.get()->Clone("h_ratio_error");
        TH1D *h_ratio_constr_values = (TH1D*)slice_mc_plus_ext_constr->hist_.get()->Clone("h_ratio_values");

        for(int i=0 ; i <= h_ratio_constr_values->GetNbinsX() ; i++){
            h_ratio_constr_values->SetBinError(i, 0);
        }

        h_ratio_constr->Divide(h_ratio_constr_values);
        h_ratio_constr_error->Divide(h_ratio_constr_values);

        h_ratio_constr->SetLineColor(kRed);
        h_ratio_constr->SetLineWidth(1);

        h_ratio_constr->Draw("same");
        h_ratio_constr_error->Draw("e2 same");

        // write to file
        std::string out_pdf_name = "plots/constrained_plots_slice_";
        if ( sl_idx < 10 ) out_pdf_name += "0";
        out_pdf_name += std::to_string( sl_idx )+".pdf";
        c->SaveAs(out_pdf_name.c_str());
        //std::cout << "########################## Saved " << out_pdf_name << " ##########################" << std::endl;
        delete c;

        // Background only plots
        TCanvas* c2 = new TCanvas("", "", 1080, 1080);

        // normalise by unconstrained
        if (normaliseByUnconstrained) {
            TH1D *h_background_values = (TH1D*)slice_bkgd_unconstr->hist_.get()->Clone("h_background_values");
            for(int i=0 ; i <= h_background_values->GetNbinsX() ; i++){
                h_background_values->SetBinError(i, 0);
            }

            slice_bkgd_constr->hist_->Divide(h_background_values);
            slice_bkgd_unconstr->hist_->Divide(h_background_values);  
        }

        // Set the color and line thickness of the histograms
        slice_bkgd_unconstr->hist_->SetLineColor(kBlue);
        slice_bkgd_unconstr->hist_->SetLineWidth(1);
    
        slice_bkgd_constr->hist_->SetLineColor(kRed);
        slice_bkgd_constr->hist_->SetLineWidth(1);
        
        // Set the minimum value of the y-axis to 0
        slice_bkgd_unconstr->hist_->SetMinimum(0.0);
        slice_bkgd_constr->hist_->SetMinimum(0.0);

        slice_bkgd_unconstr->hist_->GetYaxis()->SetRangeUser(0, slice_bkgd_unconstr->hist_->GetMaximum() * 2);

        // Draw the histograms
        slice_bkgd_unconstr->hist_->DrawCopy( "hist" );
        slice_bkgd_unconstr->hist_->SetFillColor(kBlue);
        slice_bkgd_unconstr->hist_->SetFillStyle(3004);
        slice_bkgd_unconstr->hist_->Draw("e2 same");

        slice_bkgd_constr->hist_->DrawCopy("hist same");
        slice_bkgd_constr->hist_->SetFillColor(kRed);
        slice_bkgd_constr->hist_->SetFillStyle(3005);
        slice_bkgd_constr->hist_->Draw("e2 same");        

        // Print fractional uncertainties in each bin
        for (int bin = 1; bin <= slice_bkgd_unconstr->hist_->GetNbinsX(); bin++) {
            std::cout << "Bin " << bin << " -- ";
            // Unconstrained
            std::cout << "Unconstrained: " << slice_bkgd_unconstr->hist_->GetBinError(bin) / slice_bkgd_unconstr->hist_->GetBinContent(bin);
            // Constrained
            std::cout << ", Constrained: " << slice_bkgd_constr->hist_->GetBinError(bin) / slice_bkgd_constr->hist_->GetBinContent(bin) << std::endl;
        }
        
        // Add a legend
        TLegend* legendSlice2 = new TLegend(0.5,0.70,0.9,0.9);
        legendSlice2->AddEntry(slice_bkgd_unconstr->hist_.get(), "MC+EXT (unconstr)", "l");
        legendSlice2->AddEntry(slice_bkgd_constr->hist_.get(), "MC+EXT (constr)", "l");
        legendSlice2->Draw();

        // write to file
        std::string out_pdf_name2 = "plots/constrained_plots_background_slice_";
        if ( sl_idx < 10 ) out_pdf_name2 += "0";
        out_pdf_name2 += std::to_string( sl_idx )+".pdf";
        c2->SaveAs(out_pdf_name2.c_str());
        //std::cout << "########################## Saved " << out_pdf_name2 << " ##########################" << std::endl;
        delete c2;

    } // End of slice loop

    std::cout << "Finished Slice Loop" << std::endl;

    // ******************************************************************************************************************
    // Plot the histograms for the entire reco space ********************************************************************
    // ******************************************************************************************************************
    // Get the total covariance matrices for the unconstrained and constrained cases
    auto* cov_mat = matrix_map.at( "total" ).cov_matrix_.get();
    auto* cov_mat_constr = matrix_map_constr.at( "total" ).cov_matrix_.get();

    // Get the number of ordinary and sideband reco bins
    int num_ordinary_reco_bins = 0;
    int num_sideband_reco_bins = 0;
    for ( int b = 0; b < syst_constr.reco_bins_.size(); ++b ) 
    {
        const auto& rbin = syst_constr.reco_bins_.at( b );
        if ( rbin.type_ == kSidebandRecoBin ) ++num_sideband_reco_bins;
        else ++num_ordinary_reco_bins;
    }

    // Get the number of true signal bins
    int num_true_signal_bins = 0;
    for ( int t = 0; t < syst_constr.true_bins_.size(); ++t )
    {
        const auto& tbin = syst_constr.true_bins_.at( t );
        if ( tbin.type_ == kSignalTrueBin ) ++num_true_signal_bins;
    }

    TH1D* reco_data_hist = dynamic_cast< TH1D* >(
        syst_constr.data_hists_.at( NFT::kOnBNB )->Clone( "reco_data_hist" )
    );

    TH1D* reco_ext_hist = syst_constr.data_hists_.at( NFT::kExtBNB ).get();
    const auto& cv_univ = syst_constr.cv_universe();
    int num_reco_bins = reco_data_hist->GetNbinsX();

    // Clone the reco data hist twice. We will fill the clones with the CV
    // MC+EXT prediction and the constrained one
    TH1D* reco_mc_and_ext_hist = dynamic_cast< TH1D* >(
    reco_data_hist->Clone( "reco_mc_and_ext_hist" )
    );

    reco_mc_and_ext_hist->Reset();
    reco_mc_and_ext_hist->Add( reco_ext_hist );
    reco_mc_and_ext_hist->Add( cv_univ.hist_reco_.get() );

    TH1D* reco_constrained_hist = dynamic_cast< TH1D* >(
    reco_data_hist->Clone( "reco_constrained_hist" )
    );
    reco_constrained_hist->Reset();

    // reco_mc_and_ext_hist_no_constr serves as a safety check that
    // reco_mc_and_ext_hist uncertainties do not change when set the same way as reco_constrained_hist
    TH1D* reco_mc_and_ext_hist_no_constr = dynamic_cast< TH1D* >(
        reco_mc_and_ext_hist->Clone( "reco_mc_and_ext_hist_no_constr" )
    );

    // Get the post-constraint event counts and covariance matrix in the signal region
    for ( int rb = 0; rb < num_reco_bins; ++rb ) {
        if(cov_mat_constr->GetBinContent(rb + 1, rb + 1) < 0 || cov_mat->GetBinContent(rb + 1, rb + 1) < 0)
            throw std::runtime_error("Negative diagonal covariance matrix element");

        // const double err_constr = std::sqrt(
        //     std::max( 0., cov_mat_constr->GetBinContent(rb + 1, rb + 1) )
        // );
        // reco_mc_and_ext_hist->SetBinError( rb + 1, err_constr );

        const double err = std::sqrt(
            std::max( 0., cov_mat->GetBinContent(rb + 1, rb + 1) )
        );
        reco_mc_and_ext_hist_no_constr->SetBinError( rb + 1, err );

        if ( rb >= num_ordinary_reco_bins )
        {
            double data_evts = reco_data_hist->GetBinContent( rb + 1 );
            reco_constrained_hist->SetBinContent( rb + 1, data_evts );
            reco_constrained_hist->SetBinError( rb + 1, 0. );
            // reco_mc_and_ext_hist_no_constr->SetBinContent( rb + 1, data_evts );
            // reco_mc_and_ext_hist_no_constr->SetBinError( rb + 1, 0 );
        }
        else 
        {
            double constr_pred = meas_constr.reco_mc_plus_ext_->operator()( rb, 0 );
            double constr_err = std::sqrt(
            std::max( 0., meas_constr.cov_matrix_->operator()(rb, rb) )
            );

            reco_constrained_hist->SetBinContent( rb + 1, constr_pred );
            reco_constrained_hist->SetBinError( rb + 1, constr_err );

            double pred = meas_unconstr.reco_mc_plus_ext_->operator()( rb, 0 );
            double err = std::sqrt(
            std::max( 0., meas_unconstr.cov_matrix_->operator()(rb, rb) )
            );
            reco_mc_and_ext_hist_no_constr->SetBinContent( rb + 1, pred );
            reco_mc_and_ext_hist_no_constr->SetBinError( rb + 1, err );
        }
    }

    // ******************************************************************************************************************
    // Plot the histograms for sideband region **************************************************************************
    // ******************************************************************************************************************
    TCanvas* c2 = new TCanvas;

    reco_data_hist->SetLineColor( kBlack );
    reco_data_hist->SetLineWidth( 1 );

    //reco_mc_and_ext_hist->SetLineColor( kRed );
    //->SetLineStyle( 2 );
    //reco_mc_and_ext_hist->SetLineWidth( 3 );

    reco_constrained_hist->SetLineColor( kBlue );
    reco_constrained_hist->SetLineStyle( 9 );
    reco_constrained_hist->SetLineWidth( 3 );

    reco_mc_and_ext_hist_no_constr->SetLineColor( kRed );
    reco_mc_and_ext_hist_no_constr->SetLineStyle( 2 );
    reco_mc_and_ext_hist_no_constr->SetLineWidth( 3 );

    // Set the x-axis range to the first num_ordinary_reco_bins bins
    reco_data_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 3*num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 3*num_ordinary_reco_bins);
    reco_constrained_hist->GetXaxis()->SetRange(num_ordinary_reco_bins+1, 3*num_ordinary_reco_bins);

    reco_data_hist->GetYaxis()->SetRangeUser(0, reco_mc_and_ext_hist_no_constr->GetMaximum() * 1.2);

    // Draw the histograms
    reco_data_hist->Draw( "e" );
    reco_mc_and_ext_hist_no_constr->Draw( "same hist e" );
    //reco_mc_and_ext_hist_no_constr->Draw( "same hist e" );
    reco_constrained_hist->Draw( "same hist e" );
    reco_data_hist->Draw( "same e" );

    TLegend* lg2 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    lg2->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data", "l" );
    lg2->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT", "l" );
    lg2->AddEntry( reco_constrained_hist, "post-constraint", "l" );
    lg2->Draw( "same" );
    c2->SaveAs("plots/constrained_plots_reco_sideband.pdf");
    std::cout << "########################## Saved " << "plots/constrained_plots_reco_sideband.pdf" << " ##########################" << std::endl;


    // ******************************************************************************************************************
    // Plot the histograms for signal region ****************************************************************************
    // ******************************************************************************************************************
    TCanvas* c3 = new TCanvas;

    reco_data_hist->SetLineColor( kBlack );
    reco_data_hist->SetLineWidth( 1 );

    reco_mc_and_ext_hist_no_constr->SetLineColor( kRed );
    reco_mc_and_ext_hist_no_constr->SetLineStyle( 2 );
    reco_mc_and_ext_hist_no_constr->SetLineWidth( 3 );

    reco_constrained_hist->SetLineColor( kBlue );
    reco_constrained_hist->SetLineStyle( 9 );
    reco_constrained_hist->SetLineWidth( 2 );

    // Set the x-axis range to the first num_ordinary_reco_bins bins
    reco_data_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_mc_and_ext_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_mc_and_ext_hist_no_constr->GetXaxis()->SetRange(1, num_ordinary_reco_bins);
    reco_constrained_hist->GetXaxis()->SetRange(1, num_ordinary_reco_bins);

    // Draw the histograms
    reco_data_hist->Draw( "e" );
    reco_mc_and_ext_hist_no_constr->Draw( "same hist e" );
    // reco_mc_and_ext_hist_no_constr->Draw( "same hist e" );
    reco_constrained_hist->Draw( "same hist e" );
    reco_data_hist->Draw( "same e" );

    TLegend* lg3 = new TLegend( 0.15, 0.7, 0.3, 0.85 );
    lg3->AddEntry( reco_data_hist, using_fake_data ? "fake data" : "data", "l" );
    lg3->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT", "l" );
    // lg3->AddEntry( reco_mc_and_ext_hist_no_constr, "uB tune + EXT (safety check)", "l" );
    lg3->AddEntry( reco_constrained_hist, "uB tune + EXT (constrained)", "l" );
    lg3->Draw( "same" );

    c3->SaveAs("plots/constrained_plots_reco.pdf");
    std::cout << "########################## Saved " << "plots/constrained_plots_reco.pdf" << " ##########################" << std::endl;

    //delete reco_signal_plus_bkgd_constr, syst_ptr_constr, sb_ptr;
}

int constrained_background_slice_plots()
{
    slice_plots();
    std::cout << "------------All Done------------" << std::endl;
    return 0;
}
