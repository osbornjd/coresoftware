#include "InttRawHitQA.h"

#include <qautils/QAHistManagerDef.h>
#include <qautils/QAUtil.h>

using namespace std;

//____________________________________________________________________________..
InttRawHitQA::InttRawHitQA(const std::string &name):
  SubsysReco(name)
{

}

//____________________________________________________________________________..
InttRawHitQA::~InttRawHitQA()
{

}

vector < InttRawHit* > InttRawHitQA::GetHits()
{
  
  vector < InttRawHit* > hits;
  auto raw_hit_num = node_inttrawhit_map_->get_nhits();
  for (unsigned int i = 0; i < raw_hit_num; i++)
    {
      auto hit = node_inttrawhit_map_->get_hit(i);
      hits.push_back( hit );
    }
  
  return hits;
}

void InttRawHitQA::ProcessHists()
{

  //////////////////////////////////////////////////////////////////
  // Normalization using the last event counter                   //
  //////////////////////////////////////////////////////////////////
  for( int felix=0; felix<InttQa::kFelix_num; felix++ )
    {
      hist_fee_chip_chan_[ felix ]			->Scale( 1.0 / event_counter_by_myself_ );
      hist_fee_bco_full_event_counter_[ felix ]		->Scale( 1.0 / event_counter_by_myself_ );
      hist_fee_bco_full_event_counter_diff_[ felix ]	->Scale( 1.0 / event_counter_by_myself_ );
      hist_event_counter_[ felix ]			->Scale( 1.0 / event_counter_by_myself_ );
      hist_event_counter_diff_[ felix ]			->Scale( 1.0 / event_counter_by_myself_ );
    }
  
  hist_nhit_		->Scale( 1.0 / event_counter_by_myself_ ); 
  hist_nhit_south_	->Scale( 1.0 / event_counter_by_myself_ ); 
  hist_nhit_north_	->Scale( 1.0 / event_counter_by_myself_ ); 
  hist_adc_		->Scale( 1.0 / event_counter_by_myself_ );
  hist_bco_		->Scale( 1.0 / event_counter_by_myself_ ); 
  hist_bco_full_	->Scale( 1.0 / event_counter_by_myself_ );

  //////////////////////////////////////////////////////////////////
  // pid dist                                                     //
  //////////////////////////////////////////////////////////////////
  for( int i=0; i<InttQa::kFelix_num; i++ )
    {
      hist_pid_->SetBinContent( i+1, hist_fee_chip_chan_[i]->GetEntries() );
    }
  hist_pid_->Scale( 1.0 / event_counter_by_myself_ );

  for (int felix = 0; felix < InttQa::kFelix_num; felix++)
    {
      for (int ladder = 0; ladder < InttQa::kFee_num; ladder++)
        {
          int xbin = ladder + 1;
          for (int ybin = 1; ybin <= hist_fee_chip_chan_[felix]->GetNbinsY(); ybin++)
            {
              for (int zbin = 1; zbin <= hist_fee_chip_chan_[felix]->GetNbinsZ(); zbin++)
                {
                  double content = hist_fee_chip_chan_[felix]->GetBinContent(xbin, ybin, zbin);
                  double x = hist_fee_chip_chan_[felix]->GetYaxis()->GetBinCenter(ybin);
                  double y = hist_fee_chip_chan_[felix]->GetZaxis()->GetBinCenter(zbin);

                  hist_hitmap_[felix][ladder]->Fill(x, y, content);
                }
            }
        }
    }

}

int InttRawHitQA::Init(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawHitQA::InitRun(PHCompositeNode *topNode)
{
  createHistos();

  /////////////////////////////////////////////////////////////////////////
  // INTT raw hit
  string node_name_inttrawhit = "INTTRAWHIT";
  node_inttrawhit_map_ =
    findNode::getClass<InttRawHitContainer>(topNode, node_name_inttrawhit);
  
  if (!node_inttrawhit_map_)
    {
      cerr << PHWHERE << node_name_inttrawhit << " node is missing." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for( int felix=0; felix<InttQa::kFelix_num;felix++ )
    {
      string name = getHistoPrefix() + "intt" + to_string( felix );
      hist_fee_chip_chan_[ felix ] = dynamic_cast<TH3D *>(hm->getHisto(name.c_str()));

      string name_bco_event = name + "_ladder_bco_full_event_counter";
      hist_fee_bco_full_event_counter_[ felix ] = dynamic_cast<TH3D *>(hm->getHisto(name_bco_event.c_str()));

      string name_bco_event_diff = name + "_ladder_bco_full_event_counter_diff";
      hist_fee_bco_full_event_counter_diff_[ felix ] = dynamic_cast<TH3D *>(hm->getHisto(name_bco_event_diff.c_str()));

      string name_event_counter = name + "_event_counter";
      hist_event_counter_[ felix ] = dynamic_cast<TH1D *>(hm->getHisto(name_event_counter.c_str()));
 
      string name_event_counter_diff = name + "_event_counter_diff";
      hist_event_counter_diff_[ felix ] = dynamic_cast<TH1D *>(hm->getHisto(name_event_counter_diff.c_str()));
    }

  for( int felix=0; felix<InttQa::kFelix_num;felix++ )
    {
      for( int ladder=0; ladder<InttQa::kFee_num; ladder++ )
	{
	  string name = getHistoPrefix() + "intt" + to_string( felix ) + "_" + to_string( ladder );
          hist_hitmap_[ felix ][ ladder ] = dynamic_cast<TProfile2D *>(hm->getHisto(name.c_str()));
	}
    }

  hist_nhit_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "nhit").c_str()));
 
  hist_nhit_south_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "nhit_south").c_str()));
  hist_nhit_north_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "nhit_north").c_str()));
  hist_pid_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "pid").c_str()));
  hist_adc_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "adc").c_str()));
  hist_bco_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "bco").c_str()));
  hist_bco_full_ = dynamic_cast<TH1D *>(hm->getHisto(std::string(getHistoPrefix() + "bco_full").c_str()));
  
  return Fun4AllReturnCodes::EVENT_OK;
}

void InttRawHitQA::createHistos()
{
  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  for( int felix=0; felix<InttQa::kFelix_num;felix++ )
    {
      string name = getHistoPrefix() + "intt" + to_string( felix );
      string title = name + ";FELIX CH;Chip;Channel;Entries/event";
      {
        auto h = new TH3D( name.c_str(), title.c_str(),
		    InttQa::kFee_num, 0, InttQa::kFee_num,
		    InttQa::kChip_num, 1, InttQa::kChip_num+1,
		    InttQa::kChan_num, 0, InttQa::kChan_num );
        hm->registerHisto(h);
      }

      string name_bco_event = name + "_ladder_bco_full_event_counter";
      string title_bco_event = name + ";FELIX_CH;BCO full;Event Counter;Entries/event";
      {
        auto h = new TH3D( name_bco_event.c_str(), title_bco_event.c_str(),
		    InttQa::kFee_num, 0, InttQa::kFee_num,
		    100, 0, TMath::Power(2, 40),
		    1e4, 0, 1e7 );
        hm->registerHisto(h);
      }
     
      string name_bco_event_diff = name + "_ladder_bco_full_event_counter_diff";
      string title_bco_event_diff = name + ";FELIX_CH;#Delta BCO full;#Delta Event Counter;Entries/event";
      int max = 1000;
      {
        auto h = new TH3D( name_bco_event_diff.c_str(), title_bco_event_diff.c_str(),
		    InttQa::kFee_num, 0, InttQa::kFee_num,
		    2 * max / 100, -max, max,
		    2 * max / 100, -max, max );
        hm->registerHisto(h);
      }

      string name_event_counter = name + "_event_counter";
      string title_event_counter = name_event_counter + ";Event Counter;Entries/event";
      {
        auto h = new TH1D(name_event_counter.c_str(), title_event_counter.c_str(), 1e4, 0, 1e7);
        hm->registerHisto(h);
      }
 
      string name_event_counter_diff = name + "_event_counter_diff";
      string title_event_counter_diff = name_event_counter_diff + ";Event Counter;Entries/event";
      {
        auto h = new TH1D(name_event_counter_diff.c_str(), title_event_counter_diff.c_str(), 2 * max / 100, -max, max);
        hm->registerHisto(h);
      }
    }

  for( int felix=0; felix<InttQa::kFelix_num;felix++ )
    {

      for( int ladder=0; ladder<InttQa::kFee_num; ladder++ )
	{
	  string name = getHistoPrefix() + "intt" + to_string( felix ) + "_" + to_string( ladder );
	  string title = name + ";Chip;Channel;Entries/event";
          auto h = new TProfile2D( name.c_str(), title.c_str(),
				InttQa::kChip_num, 1, InttQa::kChip_num,
				InttQa::kChan_num, 0, InttQa::kChan_num );
          hm->registerHisto(h);
	}
    }

  {
    auto h = new TH1D( std::string(getHistoPrefix() + "nhit").c_str(), "#INTTRAWHIT per event;#hit;Entries/event", 1e4, 0, 1e4 );
    InttQa::HistConfig( h );
    hm->registerHisto(h);
  }

  {
    auto h = new TH1D( std::string(getHistoPrefix() + "nhit_south").c_str(), "#INTTRAWHIT South;event;#hit/event", 1e4, 0, 1e7 );
    InttQa::HistConfig( h );
    hm->registerHisto(h);
  }
  
  {
    auto h = new TH1D( std::string(getHistoPrefix() + "nhit_north").c_str(), "#INTTRAWHIT North;event;#hit/event", 1e4, 0, 1e7 );
    InttQa::HistConfig( h );
    hm->registerHisto(h);
  }

  {
    auto h = new TH1D ( std::string(getHistoPrefix() + "pid").c_str(), "Packet ID distribution;pid;Entries/event", InttQa::kFelix_num, InttQa::kFirst_pid, InttQa::kFirst_pid + InttQa::kFelix_num );
    InttQa::HistConfig( h );
    hm->registerHisto(h);
  }
  
  {
    auto h = new TH1D( std::string(getHistoPrefix() + "adc").c_str(), "ADC distribution;ADC;Entries/event", 8, 0, 8 );
    InttQa::HistConfig( h );
    hm->registerHisto(h);
  }
  
  {
    auto h = new TH1D( std::string(getHistoPrefix() + "bco").c_str(), "BCO distribution;BCO;Entries/event", InttQa::kBco_max+10, -5, InttQa::kBco_max+5 );
    InttQa::HistConfig( h );
    hm->registerHisto(h);
  }
  
  {
    auto h = new TH1D( std::string(getHistoPrefix() + "bco_full").c_str(), "BCO full distribution;BCO full;Entries/event", 100, 0, TMath::Power( 2, 40 ) );
    InttQa::HistConfig( h );
    hm->registerHisto(h);
  }

}

int InttRawHitQA::process_event(PHCompositeNode *)
{
  auto hits = this->GetHits();
  
  auto raw_hit_num = hits.size();
  hist_nhit_->Fill( raw_hit_num );
  
  // if no raw hit is found, skip this event
  if( raw_hit_num == 0 )
    return Fun4AllReturnCodes::EVENT_OK;

  event_counter_by_myself_++;
  
  //////////////////////////////////////////////////////////////////
  // processes for each event                                     //
  //////////////////////////////////////////////////////////////////
  uint64_t bco_full = (node_inttrawhit_map_->get_hit( 0 )->get_bco());
  hist_bco_full_->Fill( bco_full );
  
  //////////////////////////////////////////////////////////////////
  // primary raw hit sweep to get some reference values           //
  //////////////////////////////////////////////////////////////////
  uint32_t event_counter_ref = node_inttrawhit_map_->get_hit(0)->get_event_counter();
  
  //////////////////////////////////////////////////////////////////
  // processes for each raw hit                                   //
  //////////////////////////////////////////////////////////////////
  // loop over all raw hits
  bool found = false;
  for (unsigned int i = 0; i < hits.size(); i++)
  {
    auto hit = hits[i];

    int felix		= hit->get_packetid() - InttQa::kFirst_pid;
    int felix_ch	= hit->get_fee();

    // uint16_t InttRawHit::get_chip_id
    int chip		= hit->get_chip_id();    
    if(chip  > InttQa::kChip_num)
      chip = chip - InttQa::kChip_num;
    
    int chan		= hit->get_channel_id();
    auto adc		= hit->get_adc();
    auto bco		= hit->get_FPHX_BCO();
    int event_counter	= hit->get_event_counter(); // uint32_t IttRawHit::get_event_counter()
    if( is_first_event_ == true )
      {
	event_counter = 0;
      }
    else if( event_counter - previous_event_counter_ > 1000 )
      {
	event_counter = -1; // it means bad
      }
    else
      {
	last_event_counter_ = event_counter;
      }

    int bco_diff = 0;
    if( (bco_full & 0x7f )  > bco )
      bco_diff = int(bco_full & 0x7f ) - bco;
    else
      bco_diff = int(bco_full & 0x7f ) + (128 - bco);
        
    //////////////////////////////////////////////////////////////////
    // Filling hists                                                //
    //////////////////////////////////////////////////////////////////

    hist_fee_chip_chan_[felix]->Fill( felix_ch, chip, chan );

    hist_adc_->Fill(hit->get_adc());
    hist_bco_->Fill(hit->get_FPHX_BCO());

    hist_fee_bco_full_event_counter_[felix]
        ->Fill(felix_ch, // chip, chan );
               hit->get_bco(),
               event_counter);

    hist_fee_bco_full_event_counter_diff_[felix]
        ->Fill(felix_ch,
               hit->get_bco() - bco_full,
               event_counter - event_counter_ref);

    hist_event_counter_[felix]
        ->Fill(event_counter);

    hist_event_counter_diff_[felix]
        ->Fill(event_counter - event_counter_ref);

    if (felix < 4)
      hist_nhit_south_->Fill(event_counter);
    else
      hist_nhit_north_->Fill(event_counter);

  }

  is_first_event_ = false;

  if( last_event_counter_ - previous_event_counter_ < 1000 )
    previous_event_counter_ = last_event_counter_; // in the case of reasonable event counter
  else
    previous_event_counter_ = -1; // in the case of a crazy event counter
  
  //cout << "-------------------------------------------------" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawHitQA::ResetEvent(PHCompositeNode *)
{
  // Intitialize for Clone hit counter
  return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawHitQA::EndRun(const int)
{
  this->ProcessHists();

  auto hm = QAHistManagerDef::getHistoManager();
  assert(hm);

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttRawHitQA::End(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

std::string InttRawHitQA::getHistoPrefix() const { return std::string("h_") + Name() + std::string("_"); }

int InttRawHitQA::Reset(PHCompositeNode *)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
