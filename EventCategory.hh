#pragma once

#include <string>
#include <map>

#include "TH1.h"

// Enum used to label event categories of interest for analysis plots
enum EventCategory {

  kCCNue1piXp = 0,
  kCCNueNpi = 1,
  kCCNuepizero = 2,
  kCCNueNp = 3,
  kCCNueOther = 4,
  kCCNumupizero = 5,
  kCCNumuOther = 6,
  kNCpizero = 7,
  kNCOther = 8,
  kOutFV= 9,
  kOutOfCryo = 10,

};

// Singleton class that helps manipulate EventCategory enum values
class EventCategoryInterpreter {

  public:

    // This is a singleton class, so we'll delete the copy constructor
    // the move constructor, and the assignment operators
    EventCategoryInterpreter( const EventCategoryInterpreter& ) = delete;
    EventCategoryInterpreter( EventCategoryInterpreter&& ) = delete;
    EventCategoryInterpreter& operator=( const EventCategoryInterpreter& )
      = delete;
    EventCategoryInterpreter& operator=( EventCategoryInterpreter&& )
      = delete;

    // Get a const reference to the singleton instance of the
    // EventCategoryInterpreter
    inline static const EventCategoryInterpreter& Instance() {

      // Create the EventCategoryInterpreter object using a static variable.
      // This ensures that the singleton instance is only created once.
      static std::unique_ptr<EventCategoryInterpreter>
        the_instance( new EventCategoryInterpreter() );

      // Return a reference to the singleton instance
      return *the_instance;
    }

    inline const std::map< EventCategory, std::string >& label_map() const
      { return event_category_to_label_map_; }

    inline std::string label( EventCategory ec ) const
      { return event_category_to_label_map_.at( ec ); }

    inline int color_code( EventCategory ec ) const
      { return event_category_to_color_map_.at( ec ); }

    inline void set_mc_histogram_style( EventCategory ec, TH1* mc_hist ) const
    {
      int color = color_code( ec );
      int style = event_category_to_style_map_.at( ec );
      mc_hist->SetFillColor( color );
      mc_hist->SetLineColor( color );
      mc_hist->SetFillStyle (style );
      mc_hist->SetStats( false );
    }

    inline void set_ext_histogram_style( TH1* ext_hist ) const {
      ext_hist->SetFillColor( 28 );
      ext_hist->SetLineColor( 28 );
      ext_hist->SetLineWidth( 2 );
      ext_hist->SetFillStyle( 3005 );
      ext_hist->SetStats( false );
    }

    inline void set_bnb_data_histogram_style( TH1* bnb_hist ) const {

      bnb_hist->SetLineColor( kBlack );
      bnb_hist->SetLineWidth( 3 );
      bnb_hist->SetMarkerStyle( kFullCircle );
      bnb_hist->SetMarkerSize( 0.8 );
      bnb_hist->SetStats( false );

      bnb_hist->GetXaxis()->SetTitleOffset( 0.0 );
      bnb_hist->GetXaxis()->SetTitleSize( 0.0 );
      bnb_hist->GetYaxis()->SetTitleSize( 0.05 );
      bnb_hist->GetYaxis()->CenterTitle( true );
      bnb_hist->GetXaxis()->SetLabelSize( 0.0 );

      // This prevents the first y-axis label label (0) to be clipped by the
      // ratio plot
      bnb_hist->SetMinimum( 1e-3 );
    }

    inline void set_stat_err_histogram_style( TH1* stat_err_hist ) const {
      stat_err_hist->SetFillColor( kBlack );
      stat_err_hist->SetLineColor( kBlack );
      stat_err_hist->SetLineWidth( 2 );
      stat_err_hist->SetFillStyle( 3004 );
    }

  private:

    EventCategoryInterpreter() {}

    std::map< EventCategory, std::string > event_category_to_label_map_ = {
      { kCCNue1piXp, "#nu_{e} CC 1#pi Xp"},
      { kCCNueNpi, "#nu_{e} CC N#pi"},
      { kCCNuepizero, "#nu_{e} CC #pi^{0}"},
      { kCCNueNp, "#nu_{e} CC Np"},
      { kCCNueOther, "#nu_{e} CC Other" },
      { kCCNumupizero, "#nu_{#mu} CC #pi^{0}"},
      { kCCNumuOther, "#nu_{#mu} CC Other"},
      { kNCpizero, "NC #pi^{0}"},
      { kNCOther, "NC Other"},
      { kOutFV, "Out of FV"},
      { kOutOfCryo, "Out of Cryo"}
    };

    std::map< EventCategory, int > event_category_to_color_map_ = {

      { kCCNue1piXp, kCyan+3 },
      { kCCNueNpi, kCyan-3 },
      { kCCNuepizero, kRed+1 }, 
      { kCCNueNp, kRed+2 },
      { kCCNueOther, kRed },
      { kCCNumupizero, kBlue-3 },
      { kCCNumuOther, kBlue-6 },
      { kNCpizero, kMagenta+3 },
      { kNCOther, kMagenta+1 },
      { kOutFV, kGray+1 },
      { kOutOfCryo, kRed-3 }
    };

    std::map< EventCategory, int > event_category_to_style_map_ = {

      { kCCNue1piXp, 1001 },
      { kCCNueNpi, 1001 },
      { kCCNuepizero, 1001 }, 
      { kCCNueNp, 1001 },
      { kCCNueOther, 1001 },
      { kCCNumupizero, 1001 },
      { kCCNumuOther, 1001 },
      { kNCpizero, 1001 },
      { kNCOther, 1001 },
      { kOutFV, 1001 },
      { kOutOfCryo, 3004 }
    };
};
