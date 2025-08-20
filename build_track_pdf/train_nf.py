import numpy as np
import torch
import torch.nn as nn
import torch.distributions as D
import matplotlib.pyplot as plt
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
import uproot

torch.manual_seed(42)
np.random.seed(42)

feature_ranges = {
    'mu_pos_px': (-6, 6),
    'mu_pos_py': (-6, 6),
    'mu_pos_pz': (1, 120),
    'mu_neg_px': (-6, 6),
    'mu_neg_py': (-6, 6),
    'mu_neg_pz': (1, 120)
}

feature_map = {
    'mu_pos_px': 'mu_pos_px',
    'mu_pos_py': 'mu_pos_py',
    'mu_pos_pz': 'mu_pos_pz',
    'mu_neg_px': 'mu_neg_px',
    'mu_neg_py': 'mu_neg_py',
    'mu_neg_pz': 'mu_neg_pz'
}

data_path = "muon_kinematics.npy"
exp_data = np.load(data_path)
exp_6d = exp_data[:, [0, 1, 2, 3, 4, 5]]  # Adjust indices based on actual data structure
if np.any(np.isnan(exp_6d)) or np.any(np.isinf(exp_6d)):
    raise ValueError("Input data contains NaN or inf values")
scaler = StandardScaler()
exp_scaled = scaler.fit_transform(exp_6d)
exp_scaled = np.clip(exp_scaled, -10, 10)  # Clip to [-10, 10]

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
exp_tensor = torch.tensor(exp_scaled, dtype=torch.float32, device=device)
print(f"Loaded experimental data shape: {exp_tensor.shape}") 

class AffineCoupling(nn.Module):
    def __init__(self, dim, hidden_dim):
        super().__init__()
        self.dim = dim
        self.net = nn.Sequential(
            nn.Linear(dim // 2, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 2 * (dim // 2))
        )

    def forward(self, x):
        x1, x2 = x.chunk(2, dim=-1)
        theta = self.net(x1)
        s, t = theta.chunk(2, dim=-1)
        s = torch.clamp(s, min=-7, max=7)
        y2 = x2 * torch.exp(s + 1e-6) + t
        y = torch.cat([x1, y2], dim=-1)
        log_det = s.sum(dim=-1)
        return y, log_det

    def inverse(self, y):
        y1, y2 = y.chunk(2, dim=-1)
        theta = self.net(y1)
        s, t = theta.chunk(2, dim=-1)
        s = torch.clamp(s, min=-7, max=7)
        x2 = (y2 - t) * torch.exp(-(s + 1e-6))
        x = torch.cat([y1, x2], dim=-1)
        log_det = -s.sum(dim=-1)
        return x, log_det

class Permute(nn.Module):
    def __init__(self, permutation):
        super().__init__()
        self.register_buffer("permutation", permutation)
        self.register_buffer("inverse_permutation", torch.argsort(permutation))

    def forward(self, x):
        return x[:, self.permutation]

    def inv(self, x):
        return x[:, self.inverse_permutation]

def create_flow(dim, num_layers=12, hidden_dim=128):
    transforms = []
    for _ in range(num_layers):
        transforms.append(AffineCoupling(dim, hidden_dim).to(device))
        transforms.append(Permute(torch.randperm(dim)).to(device))
    return transforms

def forward(x, transforms):
    log_det = torch.zeros(x.shape[0], device=x.device)
    z = x
    for transform in transforms:
        if isinstance(transform, AffineCoupling):
            z, ld = transform(z)
            log_det += ld
        else:
            z = transform(z)
    return z, log_det

def inverse(z, transforms):
    log_det_total = torch.zeros(z.shape[0], device=z.device)
    x = z
    for transform in reversed(transforms):
        if isinstance(transform, AffineCoupling):
            x, log_det = transform.inverse(x)
        else:
            x = transform.inv(x)
            log_det = 0
        log_det_total += log_det
    return x, log_det_total

def log_prob(x, base_dist, transforms):
    z, log_det = forward(x, transforms)
    return base_dist.log_prob(z) + log_det

def train_flow(transforms, base_dist, data_tensor, epochs=20000, batch_size=512, lr=0.0001, grad_clip=1.0):
    params = [p for transform in transforms if isinstance(transform, AffineCoupling) for p in transform.parameters()]
    optimizer = torch.optim.Adam(params, lr=lr)
    losses = []
    num_samples = data_tensor.shape[0]
    for epoch in tqdm(range(epochs), desc="Training"):
        indices = torch.randint(0, num_samples, (batch_size,))
        x_batch = data_tensor[indices]
        loss = -log_prob(x_batch, base_dist, transforms).mean()
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(params, grad_clip)
        optimizer.step()
        losses.append(loss.item())
    return losses

def generate_samples(n_samples, base_dist, transforms, device):
    z = base_dist.sample((n_samples,))
    samples, _ = inverse(z, transforms)
    return samples.detach().cpu().numpy()

dim = 6  # Six muon  momentum features
n_generated = 500000
transforms = create_flow(dim)
base_dist = D.Independent(D.Normal(torch.zeros(dim, device=device), torch.ones(dim, device=device)), 1)
train_losses = train_flow(transforms, base_dist, exp_tensor)
generated_samples_scaled = generate_samples(n_generated, base_dist, transforms, device)
generated_samples = scaler.inverse_transform(generated_samples_scaled)

feature_names = [
    'mu_pos_px', 'mu_pos_py', 'mu_pos_pz',
    'mu_neg_px', 'mu_neg_py', 'mu_neg_pz'
]

range_table = []
for i, name in enumerate(feature_names):
    exp_min = exp_6d[:, i].min()
    exp_max = exp_6d[:, i].max()
    gen_min = generated_samples[:, i].min()
    gen_max = generated_samples[:, i].max()
    range_key = feature_map[name]
    expected_range = feature_ranges.get(range_key, (None, None))
    range_table.append({
        'Feature': name,
        'Exp Min': exp_min,
        'Exp Max': exp_max,
        'Gen Min': gen_min,
        'Gen Max': gen_max,
        'Expected Min': expected_range[0],
        'Expected Max': expected_range[1]
    })

mask = np.ones(n_generated, dtype=bool)

for i, feature in enumerate(feature_names):
    min_val, max_val = feature_ranges[feature]
    mask &= (generated_samples[:, i] >= min_val) & (generated_samples[:, i] <= max_val)

filtered_data = generated_samples[mask]
data_dict = {feature: filtered_data[:, i].astype(np.float32) for i, feature in enumerate(feature_names)}

root_file_path = "data_in/data_in/muon_momenta.root"
with uproot.recreate(root_file_path) as root_file:
    root_file["tree"] = data_dict
