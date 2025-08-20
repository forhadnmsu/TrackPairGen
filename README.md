# Combinatoric Generation and Reconstruction with RUS
This module allows you to perform combinatoric event generation and reconstruction using the **`Fun4Sim.C`** macro.

## Step 1: Run the Python Script and Apply Mass-Momentum Selection
Befrore you run the Fun4All macro, you need the appropiate track PDF. If you want to generate track momenta from unform distributions then, you can ignore this step.
This script uses the following event selections in the experimental commisoning data from SpinQuest, where the selections are already applied so that only target likely events are survibed.

Remove events in the range:
$2.9 < M_{\mu^+\mu^-} < 3.4$ and $65 < p_z^{\mu^+\mu^-} < 75$

**Notation:**
- $M_{\mu^+\mu^-}$ — invariant mass of the dimuon
- $p_z^{\mu^+\mu^-}$ — $P_z$ of the dimuon


In the first step, run the Python scripts to muons track momenta PDF, and you need the these python library: 

```bash
cd build_track_pdf
python3 get_trk_mom.py
python3 train_nf.py
python makePDF.py 
```
### Required Python Libraries

You need the following Python libraries in addition to the standard Python libraries:

- numpy
- torch
- matplotlib
- tqdm
- scikit-learn
- ROOT or uproot

## Step 2: Construct the PDF using normalizing flow and pass it to Fun4All: Steps to gerenrate and Reconstruct Combinatoric Simulations

1. **Clone the repository**:
    ```bash
    git clone https://github.com/forhadnmsu/TrackPairGen
    ```

2. **Compile the script** (if you haven’t done so already):
    ```bash
    cd TrackPairGen
    source setup.sh      
    cmake-this
    make-this
    ```

3. **Run a test job locally** and check the output files:
    ```bash
    cd mc_gen_comb 
    root -b 'Fun4Sim.C(5)'
    ```

4. **If everything looks correct**, run the job on the Rivanna HPC with a few events:
    ```bash
    ./gridrun.sh test 1 10 200
    ```

