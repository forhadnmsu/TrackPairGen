import ROOT
import array
import math

input_file = "data_in/muon_momenta.root"
output_file = "../mc_gen_comb/NormMomentumMap.root"

f_in = ROOT.TFile.Open(input_file)
if not f_in or f_in.IsZombie():
    raise RuntimeError(f"Cannot open input file: {input_file}")

tree = f_in.Get("tree")
if not tree:
    raise RuntimeError("Cannot find TTree 'muon_tree' in input file")

nbins = 50  # Number of bins per dimension (adjust as needed)
ranges = [
    (-6, 6),   # mu_pos_px
    (-4, 4),   # mu_pos_py
    (10, 100), # mu_pos_pz
    (-6, 6),   # mu_neg_px
    (-4, 4),   # mu_neg_py
    (10, 100)  # mu_neg_pz
]
dim = 6
nbins_array = array.array('i', [nbins] * dim)
xmin = array.array('d', [ranges[i][0] for i in range(dim)])  
xmax = array.array('d', [ranges[i][1] for i in range(dim)])  

h6d = ROOT.THnSparseD("MomentumMap", "6D Momentum Probability Map", dim, nbins_array, xmin, xmax)

h6d.GetAxis(0).SetTitle("mu_pos_px [GeV/c]")
h6d.GetAxis(1).SetTitle("mu_pos_py [GeV/c]")
h6d.GetAxis(2).SetTitle("mu_pos_pz [GeV/c]")
h6d.GetAxis(3).SetTitle("mu_neg_px [GeV/c]")
h6d.GetAxis(4).SetTitle("mu_neg_py [GeV/c]")
h6d.GetAxis(5).SetTitle("mu_neg_pz [GeV/c]")

pt_min = 0.0
pt_max = 5.0
theta_mu_max = 30.0  

n_entries = tree.GetEntries()
point = array.array('d', [0.0] * 6)  #array for [mu_pos_px, mu_pos_py, mu_pos_pz, mu_neg_px, mu_neg_py, mu_neg_pz]

for i in range(n_entries):
    tree.GetEntry(i)
    point[0] = tree.mu_pos_px
    point[1] = tree.mu_pos_py
    point[2] = tree.mu_pos_pz
    point[3] = tree.mu_neg_px
    point[4] = tree.mu_neg_py
    point[5] = tree.mu_neg_pz

    pt_pos = math.sqrt(point[0]**2 + point[1]**2)
    pt_neg = math.sqrt(point[3]**2 + point[4]**2)
    theta_pos = math.atan2(pt_pos, point[2]) * 180.0 / math.pi
    theta_neg = math.atan2(pt_neg, point[5]) * 180.0 / math.pi

    if (pt_pos >= pt_min and pt_pos <= pt_max and theta_pos <= theta_mu_max and
        pt_neg >= pt_min and pt_neg <= pt_max and theta_neg <= theta_mu_max):
        h6d.Fill(point)

integral = h6d.Integral(True)  # Integrate over all bins
if integral > 0:
    h6d.Scale(1.0 / integral)
else:
    print("Warning: Histogram integral is zero, cannot normalize")
f_out = ROOT.TFile.Open(output_file, "RECREATE", "", 9)  # Maximum compression
h6d.Write()
f_out.Close()
f_in.Close()

