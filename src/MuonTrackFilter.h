#ifndef MUON_TRACK_FILTER_H
#define MUON_TRACK_FILTER_H

#include <fun4all/SubsysReco.h>

class SQTrackVector;  

class MuonTrackFilter : public SubsysReco {
public:
    MuonTrackFilter(const std::string& name = "MuonTrackFilter");
    virtual ~MuonTrackFilter();

    int Init(PHCompositeNode *topNode);
    int InitRun(PHCompositeNode *topNode);
    int process_event(PHCompositeNode *topNode);
    int End(PHCompositeNode *topNode);

    void SetAngleThreshold(double min_angle, double max_angle) {
        angle_min = min_angle;
        angle_max = max_angle;
    }

private:
    double angle_min;
    double angle_max;
    SQTrackVector* m_vec_trk;  

};

#endif // MUON_TRACK_FILTER_H

