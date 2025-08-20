 #include "Fun4AllRUSEventInputManager.h"
 #include <interface_main/SQHit_v1.h>
 #include <interface_main/SQTrack_v1.h>
 #include <interface_main/SQHitVector_v1.h>
 #include <interface_main/SQTrackVector_v1.h>
 #include <interface_main/SQEvent_v1.h>
 #include <interface_main/SQRun_v1.h>
 #include <interface_main/SQSpill_v2.h>
 #include <interface_main/SQSpillMap_v1.h>
 #include <interface_main/SQStringMap.h>
 #include <interface_main/SQScaler_v1.h>
 #include <interface_main/SQSlowCont_v1.h>
 #include <fun4all/Fun4AllServer.h>
 #include <fun4all/Fun4AllSyncManager.h>
 #include <fun4all/Fun4AllReturnCodes.h>
 #include <fun4all/Fun4AllUtils.h>
 #include <fun4all/PHTFileServer.h>
 #include <ffaobjects/RunHeader.h>
 #include <ffaobjects/SyncObjectv2.h>
 #include <phool/getClass.h>
 #include <phool/PHCompositeNode.h>
 #include <phool/PHDataNode.h>
 #include <phool/recoConsts.h>
 #include <cstdlib>
 #include <memory>
 #include <TFile.h>
 #include <TTree.h>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

//Constructor
Fun4AllRUSEventInputManager::Fun4AllRUSEventInputManager(const std::string& name, const std::string& topnodename):
   Fun4AllInputManager(name, ""),
   segment(-999),
   isopen(0),
   events_total(0),
   events_thisfile(0),
   topNodeName(topnodename),
   _tree_name("save"),
   run_header(nullptr),
   spill_map(nullptr),
   event_header(nullptr),
   hit_vec(nullptr),
   trig_hit_vec(nullptr),
   true_track_vec(nullptr),
   _fin(nullptr),
   _tin(nullptr),
   mc_mode(true)
{
   Fun4AllServer* se = Fun4AllServer::instance();
   topNode = se->topNode(topNodeName.c_str());
   PHNodeIterator iter(topNode);

    PHCompositeNode* runNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode) {
    runNode = new PHCompositeNode("RUN");
    topNode->addNode(runNode);
  }

  PHCompositeNode* eventNode = static_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));
  if (!eventNode) {
    eventNode = new PHCompositeNode("DST");
    topNode->addNode(eventNode);
  }

     run_header = new SQRun_v1();
     PHIODataNode<PHObject>* runHeaderNode = new PHIODataNode<PHObject>(run_header, "SQRun", "PHObject");
     runNode->addNode(runHeaderNode);

     spill_map = new SQSpillMap_v1();
     PHIODataNode<PHObject>* spillNode = new PHIODataNode<PHObject>(spill_map, "SQSpillMap", "PHObject");
     runNode->addNode(spillNode);

     event_header = new SQEvent_v1();
     PHIODataNode<PHObject>* eventHeaderNode = new PHIODataNode<PHObject>(event_header, "SQEvent", "PHObject");
     eventNode->addNode(eventHeaderNode);

     hit_vec = new SQHitVector_v1();
     PHIODataNode<PHObject>* hitNode = new PHIODataNode<PHObject>(hit_vec, "SQHitVector", "PHObject");
     eventNode->addNode(hitNode);

     trig_hit_vec = new SQHitVector_v1();
     PHIODataNode<PHObject>* triggerhitNode = new PHIODataNode<PHObject>(trig_hit_vec, "SQTriggerHitVector", "PHObject");
     eventNode->addNode(triggerhitNode);

    
if (mc_mode == true) {
    std::cout << "mc_mode is true. Initializing true_track_vec." << std::endl;
    true_track_vec = new SQTrackVector_v1();
    PHIODataNode<PHObject>* true_track_node = new PHIODataNode<PHObject>(true_track_vec, "SQTruthTrackVector", "PHObject");
    eventNode->addNode(true_track_node);
} else {
    std::cerr << "mc_mode is false. true_track_vec will not be initialized." << std::endl;
}


   syncobject = new SyncObjectv2();
}

Fun4AllRUSEventInputManager::~Fun4AllRUSEventInputManager()
 {
   if(isopen)
   {
     fileclose();
   }
   delete syncobject;
 }

void Fun4AllRUSEventInputManager::VectToE1039() {
	// Example: Initialize default triggers
	cout << "run number"<< runID << endl;
	event_header->set_run_id(runID);
	SQSpill* spill = spill_map->get(spillID);
	if(!spill) {
		spill = new SQSpill_v2();
		spill->set_spill_id(spillID);
		spill->set_run_id(runID);

		spill_map->insert(spill);
		run_header->set_n_spill(spill_map->size());
	}
	event_header->set_spill_id(spillID);
	event_header->set_event_id(eventID);
	event_header->set_data_quality(0);
	for(int i = -16; i < 16; ++i) {event_header->set_qie_rf_intensity(i, rfIntensities[i+16]);} // need to be added
	//for(int i = -16; i < 16; ++i) {cout << "intensity: "<< rfIntensities[i+16] <<endl;} // need to be added
	
	if (event_header) {
		// Apply the FPGA triggers to the event header
		event_header->set_trigger(SQEvent::MATRIX1, fpgaTriggers[0]);
		event_header->set_trigger(SQEvent::MATRIX2, fpgaTriggers[1]);
		event_header->set_trigger(SQEvent::MATRIX3, fpgaTriggers[2]);
		event_header->set_trigger(SQEvent::MATRIX4, fpgaTriggers[3]);
		event_header->set_trigger(SQEvent::MATRIX5, fpgaTriggers[4]);
		// Apply the NIM triggers to the event header
		event_header->set_trigger(SQEvent::NIM1, nimTriggers[0]);
		event_header->set_trigger(SQEvent::NIM2, nimTriggers[1]);
		event_header->set_trigger(SQEvent::NIM3, nimTriggers[2]);
		event_header->set_trigger(SQEvent::NIM4, nimTriggers[3]);
		event_header->set_trigger(SQEvent::NIM5, nimTriggers[4]);
	}
	//we need to fill from all hits vector
	for (size_t i = 0; i < elementIDs->size(); ++i) {
		SQHit* hit = new SQHit_v1();
		hit->set_hit_id(i);
		hit->set_detector_id(detectorIDs->at(i));
		hit->set_element_id(elementIDs->at(i));
		hit->set_tdc_time(tdcTimes->at(i));
		hit->set_drift_distance(driftDistances->at(i));
		hit->set_in_time(hitsInTime->at(i));
		hit_vec->push_back(hit);
	}
	//we need to fill from the trig-hits vector (need to work on this later)
	for (size_t i = 0; i < triggerElementIDs->size(); ++i) {
		SQHit* hit = new SQHit_v1();
		hit->set_hit_id(i);
		hit->set_detector_id(triggerDetectorIDs->at(i));
		hit->set_element_id(triggerElementIDs->at(i));
		hit->set_tdc_time(triggerTdcTimes->at(i));
		hit->set_drift_distance(triggerDriftDistances->at(i));
		hit->set_in_time(triggerHitsInTime->at(i));
		trig_hit_vec->push_back(hit);
	}

	if(mc_mode==true){
		for (size_t i = 0; i < mc_track_id->size(); ++i) {
			SQTrack* trk = new SQTrack_v1();

			if (verbosity > 90) {
				cout << "charge: "<< mc_track_charges->at(i) << endl;
				cout << "id: "<< mc_track_id->at(i) << endl;
				cout << "pos x: "<< mc_pos_vtx_x->at(i) << endl;
				cout << "pos y: "<< mc_pos_vtx_y->at(i) << endl;
				cout << "pos z: "<< mc_pos_vtx_z->at(i) << endl;
				cout << "mom Px: "<<  mc_mom_vtx_px->at(i) << endl;
				cout << "mom Py: "<<  mc_mom_vtx_py->at(i) << endl;
				cout << "mom Pz: "<<  mc_mom_vtx_pz->at(i) << endl;
			}
			trk->set_track_id(mc_track_id->at(i));
			trk->set_charge(mc_track_charges->at(i));
			TVector3 vertex_position(mc_pos_vtx_x->at(i), mc_pos_vtx_y->at(i), mc_pos_vtx_z->at(i)); 
			TVector3 vertex_position_st1(mc_pos_st1_x->at(i), mc_pos_st1_y->at(i), mc_pos_st1_z->at(i)); 
			TVector3 vertex_position_st3(mc_pos_st3_x->at(i), mc_pos_st3_y->at(i), mc_pos_st3_z->at(i)); 

			trk->set_pos_vtx(vertex_position);
			trk->set_pos_vtx(vertex_position_st1);
			trk->set_pos_vtx(vertex_position_st3);

			TLorentzVector momentum_at_vertex;
			TLorentzVector momentum_at_st1;
			TLorentzVector momentum_at_st3;
			momentum_at_vertex.SetXYZM(
					mc_mom_vtx_px->at(i),
					mc_mom_vtx_py->at(i),
					mc_mom_vtx_pz->at(i),
					0.105658 // Muon mass in GeV/c^2
					);

			 momentum_at_st1.SetXYZM(
                                        mc_mom_st1_px->at(i),
                                        mc_mom_st1_py->at(i),
                                        mc_mom_st1_pz->at(i),
                                        0.105658 // Muon mass in GeV/c^2
                                        );  

			 momentum_at_st3.SetXYZM(
                                        mc_mom_st3_px->at(i),
                                        mc_mom_st3_py->at(i),
                                        mc_mom_st3_pz->at(i),
                                        0.105658 // Muon mass in GeV/c^2
                                        );  

			trk->set_mom_vtx(momentum_at_vertex);
			trk->set_mom_vtx(momentum_at_st1);
			trk->set_mom_vtx(momentum_at_st3);
			true_track_vec->push_back(trk);
		}
	}
}

int Fun4AllRUSEventInputManager::fileopen(const std::string &filenam) {
    if (isopen) {
        std::cout << "Closing currently open file "
                  << filename
                  << " and opening " << filenam << std::endl;
        fileclose();
    }   
    filename = filenam;

    // Print the absolute path for debugging
    std::cout << "Attempting to open file: " << std::filesystem::absolute(filenam) << std::endl;

    if (!std::filesystem::exists(filenam)) {
        std::cerr << "!!ERROR!! File does not exist: " << filenam << std::endl;
        return -1; 
    }   

    if (verbosity > 0) {
        std::cout << ThisName << ": opening file " << filename.c_str() << std::endl;
    }   

    events_thisfile = 0;

    _fin = TFile::Open(filenam.c_str(), "READ"); // Open the file dynamically
    if (!_fin || _fin->IsZombie()) {
        std::cerr << "!!ERROR!! Failed to open file " << filenam << std::endl;
        return -1;
    }

    _tin = (TTree*) _fin->Get(_tree_name.c_str());
    if (!_tin) {
        std::cerr << "!!ERROR!! Tree " << _tree_name << " not found in file " << filenam << std::endl;
        return -1; 
    }
_tin->SetBranchAddress("eventID", &eventID);    
_tin->SetBranchAddress("runID", &runID);    
_tin->SetBranchAddress("spillID", &spillID);    
_tin->SetBranchAddress("fpgaTriggers", fpgaTriggers);
_tin->SetBranchAddress("nimTriggers", nimTriggers);
_tin->SetBranchAddress("rfIntensities", rfIntensities);

_tin->SetBranchAddress("detectorIDs", &detectorIDs);    
_tin->SetBranchAddress("elementIDs", &elementIDs);    
_tin->SetBranchAddress("driftDistances", &driftDistances);    
_tin->SetBranchAddress("tdcTimes", &tdcTimes);    
_tin->SetBranchAddress("hitsInTime", &hitsInTime);    

_tin->SetBranchAddress("triggerDetectorIDs", &triggerDetectorIDs);    
_tin->SetBranchAddress("triggerDriftDistances", &triggerDriftDistances);    
_tin->SetBranchAddress("triggerElementIDs", &triggerElementIDs);    
_tin->SetBranchAddress("triggerTdcTimes", &triggerTdcTimes);    
_tin->SetBranchAddress("triggerHitsInTime", &triggerHitsInTime);

if(mc_mode){
	_tin->SetBranchAddress("mc_track_charges", &mc_track_charges);    
	_tin->SetBranchAddress("mc_track_id", &mc_track_id);    
	_tin->SetBranchAddress("mc_pos_vtx_x", &mc_pos_vtx_x);    
	_tin->SetBranchAddress("mc_pos_vtx_y", &mc_pos_vtx_y);    
	_tin->SetBranchAddress("mc_pos_vtx_z", &mc_pos_vtx_z);
	_tin->SetBranchAddress("mc_mom_vtx_px", &mc_mom_vtx_px);          
        _tin->SetBranchAddress("mc_mom_vtx_py", &mc_mom_vtx_py);      
        _tin->SetBranchAddress("mc_mom_vtx_pz", &mc_mom_vtx_pz);

	_tin->SetBranchAddress("mc_pos_st1_x", &mc_pos_st1_x);    
        _tin->SetBranchAddress("mc_pos_st1_y", &mc_pos_st1_y);    
        _tin->SetBranchAddress("mc_pos_st1_z", &mc_pos_st1_z);
        _tin->SetBranchAddress("mc_mom_st1_px", &mc_mom_st1_px);    
        _tin->SetBranchAddress("mc_mom_st1_py", &mc_mom_st1_py);    
        _tin->SetBranchAddress("mc_mom_st1_pz", &mc_mom_st1_pz);

	_tin->SetBranchAddress("mc_pos_st3_x", &mc_pos_st3_x);    
        _tin->SetBranchAddress("mc_pos_st3_y", &mc_pos_st3_y);    
        _tin->SetBranchAddress("mc_pos_st3_z", &mc_pos_st3_z);
        _tin->SetBranchAddress("mc_mom_st3_px", &mc_mom_st3_px);    
        _tin->SetBranchAddress("mc_mom_st3_py", &mc_mom_st3_py);    
        _tin->SetBranchAddress("mc_mom_st3_pz", &mc_mom_st3_pz);
}

    segment = 0;
    isopen = 1;
    AddToFileOpened(filenam); // Add file to the list of opened files
    return 0;
}

int Fun4AllRUSEventInputManager::run(const int nevents) {
    readagain:
    if (!isopen) {
        if (filelist.empty()) {
            if (verbosity > 0) {
                std::cout << Name() << ": No Input file open" << std::endl;
            }
            return -1;
        } else {
            if (OpenNextFile()) {
                std::cout << Name() << ": No Input file from filelist opened" << std::endl;
                return -1;
            }
        }
    }

    if (verbosity > 3) {
        std::cout << "Getting Event from " << Name() << std::endl;
    }

    if (events_thisfile >= _tin->GetEntries()) {
        fileclose();
        goto readagain;
    }

    _tin->GetEntry(events_thisfile);
    events_thisfile++;
    events_total++;
    std::cout << "EventID: "<< eventID<<std::endl;
    std::cout << "Fpga1 trigger: "<<fpgaTriggers[0] <<std::endl;

   SetRunNumber                (runID);
   mySyncManager->PrdfEvents   (events_thisfile);
   mySyncManager->SegmentNumber(spillID);
   mySyncManager->CurrentEvent (eventID);
  
   syncobject->RunNumber       (runID);
   syncobject->EventCounter    (events_thisfile);
   syncobject->SegmentNumber   (spillID);
   syncobject->EventNumber     (eventID);
   std::cout << "deugg 1: " <<std::endl;
      VectToE1039();
   std::cout << "deugg 2: " <<std::endl;
    if (RejectEvent() != Fun4AllReturnCodes::EVENT_OK) {
        ResetEvent();
        goto readagain;
    }
    return 0;
}

int Fun4AllRUSEventInputManager::ResetEvent()
{
    syncobject->Reset();
    return 0;
}

int Fun4AllRUSEventInputManager::fileclose()
 {
         if (!isopen) {
                 cout << Name() << ": fileclose: No Input file open" << endl;
                 return -1;
         }
  
         _fin->Close();
         isopen = 0;
  
         // if we have a file list, move next entry to top of the list
         // or repeat the same entry again
         if (!filelist.empty()) {
                 if (repeat) {
                         filelist.push_back(*(filelist.begin()));
                         if (repeat > 0) {
                                 repeat--;
                         }
                 }
                 filelist.pop_front();
         }
  
         return 0;
 }
  
 void
 Fun4AllRUSEventInputManager::Print(const string &what) const
 {
   Fun4AllInputManager::Print(what);
   return ;
 }
  
 int
 Fun4AllRUSEventInputManager::OpenNextFile()
 {
   while (!filelist.empty())
     {
       list<string>::const_iterator iter = filelist.begin();
       if (verbosity)
         {
           cout << PHWHERE << " opening next file: " << *iter << endl;
         }
       if (fileopen(*iter))
         {
           cout << PHWHERE << " could not open file: " << *iter << endl;
           filelist.pop_front();
         }
       else
         {
           return 0;
         }
  
     }
   return -1;
 }
 
 int Fun4AllRUSEventInputManager::PushBackEvents(const int i)
 {
   cerr << "!!ERROR!!  PushBackEvents():  Not implemented yet." << endl;
   // PushBackEvents is supposedly pushing events back on the stack which works
   // easily with root trees (just grab a different entry) but hard in these HepMC ASCII files.
   // A special case is when the synchronization fails and we need to only push back a single
   // event. In this case we save the evt pointer as save_evt which is used in the run method
   // instead of getting the next event.
 //  if (i > 0)
 //    {
 //      if (i == 1 && evt) // check on evt pointer makes sure it is not done from the cmd line
 //      {
 //        save_evt = evt;
 //        return 0;
 //      }
 //      cout << PHWHERE << ThisName
 //           << " Fun4AllRUSEventInputManager cannot push back " << i << " events into file"
 //           << endl;
 //      return -1;
 //    }
 //  if (!parser->coda)
 //    {
 //      cout << PHWHERE << ThisName
 //         << " no file open" << endl;
 //      return -1;
 //    }
   // Skipping events is implemented as
   // pushing a negative number of events on the stack, so in order to implement
   // the skipping of events we read -i events.
 //  int nevents = -i; // negative number of events to push back -> skip num events
 //  int errorflag = 0;
 //  while (nevents > 0 && ! errorflag)
 //    {
 //                      int * data_ptr = nullptr;
 //                      unsigned int coda_id = 0;
 //                      if(parser->coda->NextCodaEvent(coda_id, data_ptr))
 //                              evt = new EVIO_Event(data_ptr);
 //      if (! evt)
 //      {
 //        cout << "Error after skipping " << i - nevents 
 //             << " file exhausted?" << endl;
 //        errorflag = -1;
 //          fileclose();
 //      }
 //      else
 //      {
 //        if (verbosity > 3)
 //          {
 //                      //TODO implement this
 //            //cout << "Skipping evt no: " << evt->getEvtSequence() << endl;
 //          }
 //      }
 //      delete evt;
 //      nevents--;
 //    }
 //  return errorflag;
   return -1;
 }
  
 int
 Fun4AllRUSEventInputManager::GetSyncObject(SyncObject **mastersync)
 {
   // here we copy the sync object from the current file to the
   // location pointed to by mastersync. If mastersync is a 0 pointer
   // the syncobject is cloned. If mastersync allready exists the content
   // of syncobject is copied
   if (!(*mastersync))
     {
       if (syncobject) *mastersync = syncobject->clone();
     }
   else
     {
       *(*mastersync) = *syncobject; // copy syncobject content
     }
   return Fun4AllReturnCodes::SYNC_OK;
 }
  
 int
 Fun4AllRUSEventInputManager::SyncIt(const SyncObject *mastersync)
 {
   if (!mastersync)
     {
       cout << PHWHERE << Name() << " No MasterSync object, cannot perform synchronization" << endl;
       cout << "Most likely your first file does not contain a SyncObject and the file" << endl;
       cout << "opened by the Fun4AllDstInputManager with Name " << Name() << " has one" << endl;
       cout << "Change your macro and use the file opened by this input manager as first input" << endl;
       cout << "and you will be okay. Fun4All will not process the current configuration" << endl << endl;
       return Fun4AllReturnCodes::SYNC_FAIL;
     }
   int iret = syncobject->Different(mastersync);
   if (iret)
     {
       cout << "big problem" << endl;
       exit(1);
     }
   return Fun4AllReturnCodes::SYNC_OK;
 }
