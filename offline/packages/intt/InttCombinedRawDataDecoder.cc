#include "InttCombinedRawDataDecoder.h"
#include "InttMapping.h"

#include <trackbase/InttDefs.h>
#include <trackbase/TrkrDefs.h>  // for hitkey, hitsetkey
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitv2.h>

#include <ffarawobjects/Gl1RawHit.h>
#include <ffarawobjects/InttRawHit.h>
#include <ffarawobjects/InttRawHitContainer.h>

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>  // for PHIODataNode
#include <phool/PHNodeIterator.h>
#include <phool/getClass.h>
#include <phool/phool.h>  // for PHWHERE

#include <TSystem.h>

#include <cstdlib>   // for exit
#include <iostream>  // for operator<<, endl, bas...
#include <map>       // for _Rb_tree_iterator

InttCombinedRawDataDecoder::InttCombinedRawDataDecoder(std::string const& name)
  : SubsysReco(name)
{
  // Do nothing
  // Consider calling LoadHotChannelMapRemote()
}

int InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)
{
  if (!topNode)
  {
    std::cout << "InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)" << std::endl;
    std::cout << "\tCould not retrieve topNode; doing nothing" << std::endl;
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator dst_itr(topNode);
  PHCompositeNode* dst_node = dynamic_cast<PHCompositeNode*>(dst_itr.findFirst("PHCompositeNode", "DST"));
  if (!dst_node)
  {
    if (Verbosity())
    {
      std::cout << "InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)" << std::endl;
    }
    if (Verbosity())
    {
      std::cout << "\tCould not retrieve dst_node; doing nothing" << std::endl;
    }
    exit(1);
    gSystem->Exit(1);

    return 1;
  }

  PHNodeIterator trkr_itr(dst_node);
  PHCompositeNode* trkr_node = dynamic_cast<PHCompositeNode*>(trkr_itr.findFirst("PHCompositeNode", "TRKR"));
  if (!trkr_node)
  {
    trkr_node = new PHCompositeNode("TRKR");
    dst_node->addNode(trkr_node);
  }

  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    if (Verbosity())
    {
      std::cout << "InttCombinedRawDataDecoder::InitRun(PHCompositeNode* topNode)" << std::endl;
    }
    if (Verbosity())
    {
      std::cout << "\tMaking TrkrHitSetContainer" << std::endl;
    }

    trkr_hit_set_container = new TrkrHitSetContainerv1;
    PHIODataNode<PHObject>* new_node = new PHIODataNode<PHObject>(trkr_hit_set_container, "TRKR_HITSET", "PHObject");
    trkr_node->addNode(new_node);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)
{
  TrkrHitSetContainer* trkr_hit_set_container = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!trkr_hit_set_container)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"TRKR_HITSET\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;
    gSystem->Exit(1);
    exit(1);

    return Fun4AllReturnCodes::DISCARDEVENT;
  }

  InttRawHitContainer* inttcont = findNode::getClass<InttRawHitContainer>(topNode, m_InttRawNodeName);
  if (!inttcont)
  {
    std::cout << PHWHERE << std::endl;
    std::cout << "InttCombinedRawDataDecoder::process_event(PHCompositeNode* topNode)" << std::endl;
    std::cout << "Could not get \"" << m_InttRawNodeName << "\" from Node Tree" << std::endl;
    std::cout << "Exiting" << std::endl;

    gSystem->Exit(1);
    exit(1);
  }
  auto gl1 = findNode::getClass<Gl1RawHit>(topNode, "GL1RAWHIT");
  if (!gl1)
  {
    std::cout << PHWHERE << " no gl1 container, exiting" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  uint64_t gl1rawhitbco = gl1->get_bco();
  // get the last 40 bits by bit shifting left then right to match
  // to the mvtx bco
  auto lbshift = gl1rawhitbco << 24;
  auto gl1bco = lbshift >> 24;

  TrkrDefs::hitsetkey hit_set_key = 0;
  TrkrDefs::hitkey hit_key = 0;
  TrkrHitSetContainer::Iterator hit_set_container_itr;
  TrkrHit* hit = nullptr;

  InttNameSpace::RawData_s raw;
  InttNameSpace::Offline_s ofl;
  for (unsigned int i = 0; i < inttcont->get_nhits(); i++)
  {
    InttRawHit* intthit = inttcont->get_hit(i);
    // uint64_t gtm_bco = intthit->get_bco();

    InttNameSpace::RawFromHit(raw, intthit);
    // raw.felix_server = InttNameSpace::FelixFromPacket(intthit->get_packetid());
    // raw.felix_channel = intthit->get_fee();
    // raw.chip = (intthit->get_chip_id() + 25) % 26;
    // raw.channel = intthit->get_channel_id();

    int adc = intthit->get_adc();
    // amp = intthit->get_amplitude();
    // int bco = intthit->get_FPHX_BCO();

    if(m_HotChannelSet.find(raw) != m_HotChannelSet.end())
    {
      continue;
    }

    ofl = InttNameSpace::ToOffline(raw);
    hit_key = InttDefs::genHitKey(ofl.strip_y, ofl.strip_x);  // col, row <trackbase/InttDefs.h>
    hit_set_key = InttDefs::genHitSetKey(ofl.layer, ofl.ladder_z, ofl.ladder_phi, intthit->get_bco() - gl1bco);
    hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
    hit = hit_set_container_itr->second->getHit(hit_key);

    if (hit)
    {
      continue;
    }

    hit = new TrkrHitv2;
    hit->setAdc(adc);
    hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int InttCombinedRawDataDecoder::LoadHotChannelMapLocal(std::string const& filename)
{
  if(filename.empty())
  {
    std::cout << "int InttCombinedRawDataDecoder::LoadHotChannelMapLocal(std::string const& filename)" << std::endl;
    std::cout << "\tArgument 'filename' is empty string" << std::endl;
    return 1;
  }
  CDBTTree cdbttree(filename);
  // need to checkt for error exception
  cdbttree.LoadCalibrations();

  m_HotChannelSet.clear();
  Long64_t N = cdbttree.GetSingleIntValue("size");
  for(Long64_t n = 0; n < N; ++n)
  {
    m_HotChannelSet.insert((struct InttNameSpace::RawData_s){
        .felix_server  = cdbttree.GetIntValue(n, "felix_server"),
        .felix_channel = cdbttree.GetIntValue(n, "felix_channel"),
        .chip          = cdbttree.GetIntValue(n, "chip"),
        .channel       = cdbttree.GetIntValue(n, "channel")
    });
  }

  return 0;
}

int InttCombinedRawDataDecoder::LoadHotChannelMapRemote(std::string const& name)
{
  if(name.empty())
  {
    std::cout << "int InttCombinedRawDataDecoder::LoadHotChannelMapRemote(std::string const& name)" << std::endl;
    std::cout << "\tArgument 'name' is empty string" << std::endl;
    return 1;
  }
  std::string database = CDBInterface::instance()->getUrl(name);
  CDBTTree cdbttree(database);
  cdbttree.LoadCalibrations();

  m_HotChannelSet.clear();
  Long64_t N = cdbttree.GetSingleIntValue("size");
  for(Long64_t n = 0; n < N; ++n)
  {
    m_HotChannelSet.insert((struct InttNameSpace::RawData_s){
        .felix_server  = cdbttree.GetIntValue(n, "felix_server"),
        .felix_channel = cdbttree.GetIntValue(n, "felix_channel"),
        .chip          = cdbttree.GetIntValue(n, "chip"),
        .channel       = cdbttree.GetIntValue(n, "channel")
    });
  }

  return 0;
}

/*
        Packet* p = evt->getPacket(itr->first);
        if(!p)continue;

        int N = p->iValue(0, "NR_HITS");
        full_bco = p->lValue(0, "BCO");

        if(Verbosity() > 20)std::cout << N << std::endl;

        for(int n = 0; n < N; ++n)
        {
        rawdata = InttNameSpace::RawFromPacket(itr->second, n, p);

        adc = p->iValue(n, "ADC");
        //amp = p->iValue(n, "AMPLITUE");
        bco = p->iValue(n, "FPHX_BCO");

        offline = InttNameSpace::ToOffline(rawdata);

        hit_key = InttDefs::genHitKey(offline.strip_y, offline.strip_x); //col, row <trackbase/InttDefs.h>
        hit_set_key = InttNameSpace::genHitSetKey(offline.layer, offline.ladder_z, offline.ladder_phi, bco);

        hit_set_container_itr = trkr_hit_set_container->findOrAddHitSet(hit_set_key);
        hit = hit_set_container_itr->second->getHit(hit_key);
        if(hit)continue;

        hit = new TrkrHitv2;
        hit->setAdc(adc);
        hit_set_container_itr->second->addHitSpecificKey(hit_key, hit);
        }

        delete p;
        }
        }
        if(Verbosity() > 20)
        {
        std::cout << std::endl;
        std::cout << "Identify():" << std::endl;
        trkr_hit_set_container->identify();
        std::cout << std::endl;
        }
*/
