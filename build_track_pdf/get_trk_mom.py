import ROOT
import numpy as np

file_exp = ROOT.TFile("data_in/combined_exp.root", "READ")
tree_exp = file_exp.Get("tree")

data_list = []
MUON_MASS = 0.105658
total_events = 0

ranges = {
    'dimuon_px': (-5, 5),
    'dimuon_py': (-5, 5),
    'dimuon_pz': (10, 120),
    'mass': (0, 10)
}

momentum_ranges = {
    'mu_pos_px': (-6, 6),
    'mu_pos_py': (-6, 6),
    'mu_pos_pz': (1, 120),
    'mu_neg_px': (-6, 6),
    'mu_neg_py': (-6, 6),
    'mu_neg_pz': (1, 120)
}

for event in tree_exp:
    total_events += 1
    n_dimuons = len(event.rec_dimuon_px_pos_tgt)
    if n_dimuons == 0:
        continue

    for i in range(n_dimuons):
        vz1 = event.rec_dimuon_z_pos_vtx[i]
        vz2 = event.rec_dimuon_z_neg_vtx[i]

        # Apply vertex cuts
        if abs(vz1 - vz2) > 200 or not (-690 < vz1 < 500 and -690 < vz2 < 500):
            continue

        # Extract momentum components from target dimuons 
        mu_pos_px = event.rec_dimuon_px_pos_tgt[i] # for mu+
        mu_pos_py = event.rec_dimuon_py_pos_tgt[i]
        mu_pos_pz = event.rec_dimuon_pz_pos_tgt[i]
        mu_neg_px = event.rec_dimuon_px_neg_tgt[i] #for mu-
        mu_neg_py = event.rec_dimuon_py_neg_tgt[i]
        mu_neg_pz = event.rec_dimuon_pz_neg_tgt[i]

        # Calculate dimuon 4-vector and mass
        track_pos = ROOT.TLorentzVector()
        track_neg = ROOT.TLorentzVector()
        track_pos.SetXYZM(mu_pos_px, mu_pos_py, mu_pos_pz, MUON_MASS)
        track_neg.SetXYZM(mu_neg_px, mu_neg_py, mu_neg_pz, MUON_MASS)
        dimu = track_pos + track_neg
        mass = dimu.M()
        dimuon_px = dimu.Px()
        dimuon_py = dimu.Py()
        dimuon_pz = dimu.Pz()

        #j/psi removal
        if (2.9 < mass < 3.4 and 65 < dimuon_pz < 75):
            continue

        # Apply dimuon momentum and mass cuts
        # removing unphysical dimuons
        if not (ranges['dimuon_px'][0] <= dimuon_px <= ranges['dimuon_px'][1] and
                ranges['dimuon_py'][0] <= dimuon_py <= ranges['dimuon_py'][1] and
                ranges['dimuon_pz'][0] <= dimuon_pz <= ranges['dimuon_pz'][1] and
                ranges['mass'][0] <= mass <= ranges['mass'][1]):
            continue

        #keep the muons within the experimental accpetance ranges
        if not (momentum_ranges['mu_pos_px'][0] <= mu_pos_px <= momentum_ranges['mu_pos_px'][1] and
                momentum_ranges['mu_pos_py'][0] <= mu_pos_py <= momentum_ranges['mu_pos_py'][1] and
                momentum_ranges['mu_pos_pz'][0] <= mu_pos_pz <= momentum_ranges['mu_pos_pz'][1] and
                momentum_ranges['mu_neg_px'][0] <= mu_neg_px <= momentum_ranges['mu_neg_px'][1] and
                momentum_ranges['mu_neg_py'][0] <= mu_neg_py <= momentum_ranges['mu_neg_py'][1] and
                momentum_ranges['mu_neg_pz'][0] <= mu_neg_pz <= momentum_ranges['mu_neg_pz'][1]):
            continue

        momentum_data = [mu_pos_px, mu_pos_py, mu_pos_pz, mu_neg_px, mu_neg_py, mu_neg_pz]
        if any(np.isnan(x) or np.isinf(x) for x in momentum_data):
            continue
        data_list.append(momentum_data)
data_array = np.array(data_list, dtype=np.float64)
npy_file_path = "data_in/muon_kinematics.npy"
np.save(npy_file_path, data_array)
print(f"Saved {data_array.shape[0]} events to {npy_file_path} with shape {data_array.shape}")
file_exp.Close()
