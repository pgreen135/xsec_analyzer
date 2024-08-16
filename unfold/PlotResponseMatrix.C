void PlotResponseMatrix() {


	
	// Open the ROOT file
    TFile* file = TFile::Open("OutputFileMatrix.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << "OutputFileMatrix.root" << std::endl;
        return;
    }

    // Retrieve the TMatrixD object from the file
    TMatrixD* matrix = (TMatrixD*)file->Get("ResponseMatrix");


    // Create a TH2D histogram to hold the matrix data
    TH2D* hist = new TH2D("", "", 5, 0, 5, 5, 0, 5);

    int slice = 0;

    // Fill the histogram with the matrix data
    for (Int_t i = 0 + slice * 5; i < 5 + slice * 5; ++i) {
        for (Int_t j = 0 + slice * 5; j < 5 + slice * 5; ++j) {
            hist->SetBinContent(j + 1, i + 1, 5*(*matrix)(i, j));
        }
    }

    hist->SetStats(false);
    hist->SetMinimum(0);
    hist->SetMaximum(1.0);
  
    // Draw the histogram
    TCanvas* canvas = new TCanvas("c1", "", 1080, 1080);
    //hist->Draw("COLZ");


    // Set x and y axis labels to bin numbers
    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
    }
    for (int i = 1; i <= hist->GetNbinsY(); i++) {   
        hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
    }
    // Increase the font size of the axis labels
    hist->GetXaxis()->SetLabelSize(0.05);
    hist->GetYaxis()->SetLabelSize(0.05);

    hist->Draw("COLZ");
    // cv_hist_2d_slice.Draw("same TEXT");

    // Fill histogram with slice data
    for (int i = 0; i < hist->GetNbinsX(); i++) {
        for (int j = 0; j < hist->GetNbinsY(); j++) {
            double bin_content = hist->GetBinContent(i+1, j+1);
            TLatex* latex = new TLatex(hist->GetXaxis()->GetBinCenter(i+1), hist->GetYaxis()->GetBinCenter(j+1), Form("%.1f%%", 100*bin_content));
            latex->SetTextFont(42);
            latex->SetTextSize(0.02);
            latex->SetTextAlign(22);
            latex->Draw();
        }
    }

    // Save the canvas to a file
    canvas->SaveAs("ResponseMatrix.pdf");

    // Clean up
    delete canvas;
    file->Close();


/*
	TCanvas* c_conf_1 = new TCanvas(("c_conf_1 slice "+std::to_string(sl_idx)).c_str(), "c confusion matrix", 800, 600);
        cv_confusion_hist->SetStats(false);
        cv_confusion_hist->SetMinimum(0);
        cv_confusion_hist->SetMaximum(1.0);
        // gStyle->SetPalette(util::CreateWhiteToBlueColorPalette(20));
        gStyle->SetPalette();

        // Set x and y axis labels to bin numbers
        for (int i = 1; i <= cv_confusion_hist->GetNbinsX(); i++) {
            // const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 1 || i == 6));
            // const auto binLabel = overFlow ? Form("%d*", i) : Form("%d", i);
            cv_confusion_hist->GetXaxis()->SetBinLabel(i, Form("%d", i));
        }
        for (int i = 1; i <= cv_confusion_hist->GetNbinsY(); i++) {
            // const auto overFlow = (sl_idx == 2 && i ==6 ) || (sl_idx == 5 && (i == 1 || i == 6));
            // const auto binLabel = overFlow ? Form("%d*", i) : Form("%d", i);
            cv_confusion_hist->GetYaxis()->SetBinLabel(i, Form("%d", i));
        }
        // Increase the font size of the axis labels
        cv_confusion_hist->GetXaxis()->SetLabelSize(0.05);
        cv_confusion_hist->GetYaxis()->SetLabelSize(0.05);

        cv_confusion_hist->Draw("colz");
        // cv_hist_2d_slice.Draw("same TEXT");

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

        c_conf_1->SaveAs(("plots/cv_confusion_matrix_slice_" + std::string(sl_idx < 10 ? "0" : "") + std::to_string(sl_idx) + "_nuwro" + postfix + ".pdf").c_str());

        */
}