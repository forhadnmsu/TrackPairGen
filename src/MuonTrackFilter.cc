#include <iomanip>
#include <algorithm>
#include <fun4all/Fun4AllReturnCodes.h>
#include <interface_main/SQTrackVector.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHIODataNode.h>
#include <phool/getClass.h>
#include <geom_svc/GeomSvc.h>
#include "MuonTrackFilter.h"
#include <TLorentzVector.h>

using namespace std;

MuonTrackFilter::MuonTrackFilter(const string& name)
  : SubsysReco   (name)
  , m_vec_trk    (0)
{
  ;
}

MuonTrackFilter::~MuonTrackFilter()
{
  ;
}

int MuonTrackFilter::Init(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

int MuonTrackFilter::InitRun(PHCompositeNode* topNode)
{
  m_vec_trk  = findNode::getClass<SQTrackVector >(topNode, "SQTruthTrackVector");
  if(!m_vec_trk) {
    cerr << Name() << ":  Failed at getting SQTruthTrackVector.  Abort." << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int MuonTrackFilter::process_event(PHCompositeNode* topNode)
{
	TLorentzVector mup, mum;
	bool has_positive = false, has_negative = false;

	// First pass: Check if we have at least one mu+ and one mu-
	for (unsigned int ii = 0; ii < m_vec_trk->size(); ii++) 
	{
		SQTrack* trk = m_vec_trk->at(ii);
		int charge = trk->get_charge();

		if (charge > 0) has_positive = true;
		if (charge < 0) has_negative = true;

		if (has_positive && has_negative) 
			break;  // No need to continue checking
	}
	if (!has_positive || !has_negative) 
	{
		return Fun4AllReturnCodes::EVENT_OK;
	}


	for (unsigned int ii = 0; ii < m_vec_trk->size(); ii++) 
	{
		SQTrack* trk = m_vec_trk->at(ii);
		int charge = trk->get_charge();
		TLorentzVector mom_vtx = trk->get_mom_vtx();

		if (charge > 0) 
		{
			mup = mom_vtx;
		} 
		else if (charge < 0) 
		{
			mum = mom_vtx;
		}
	}
	double angle = mup.Vect().Angle(mum.Vect()) * (180.0 / M_PI);
	//cout << "opening angle of the dimuon: "<< angle << endl;


	return (angle > angle_min && angle < angle_max)
		? Fun4AllReturnCodes::EVENT_OK : Fun4AllReturnCodes::ABORTEVENT;
}

int MuonTrackFilter::End(PHCompositeNode* topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}
