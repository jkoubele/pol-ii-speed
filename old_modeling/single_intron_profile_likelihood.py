from pathlib import Path

import numpy as np
import pandas as pd
import torch
from matplotlib import pyplot as plt
from scipy import stats
import sys


from modeling.load_dataset import load_dataset_from_results_folder

sys.path.insert(0, str(Path("/cellfile/datapublic/jkoubele/dev-pol-ii-analysis")))

#%%

results_folder = Path('/cellfile/datapublic/jkoubele/drosophila_mutants/results/')
gene_data_list, dataset_metadata = load_dataset_from_results_folder(
    results_folder=results_folder,
    gene_names_file_name='protein_coding_genes.csv',
    log_output_folder=Path('/cellfile/datapublic/jkoubele/data_pol_ii/mouse_myocardium/results/model_results/'))

# %%
introns_df = pd.read_csv(results_folder / 'extracted_introns/introns.bed',
                         sep='\t',
                         names=['chromosome', 'start', 'end', 'name', 'score', 'strand']).set_index('name')


# gene_of_interest = 'FBgn0019985'
# intron_of_interest = 'FBgn0019985_9'
# gene_data = [x for x in gene_data_list if x.gene_name == gene_of_interest][0]

# intron_index = gene_data.intron_names.index(intron_of_interest)
gene_data = gene_data_list[0]
intron_index = 0

intron_data = gene_data.coverage[:, intron_index, :]
intron_name = gene_data.intron_names[intron_index]
intron_row = introns_df.loc[intron_name]

intron_number = intron_name.rsplit('_', maxsplit=1)[1]
intron_length = intron_row.end - intron_row.start

print(f"Intron length: {intron_length}")
coverage = intron_data.sum(axis=0).numpy()

coverage_locations = np.linspace(0.0, 1.0, len(coverage))
plt.bar(coverage_locations, coverage, width=1.05 / len(coverage_locations))

plt.xlabel("Position in intron (5' to 3')")
plt.ylabel("Coverage (# reads)")
plt.title(
    f"Intron {intron_number} of gene {gene_data.gene_name} \n coverage: {int(coverage.sum())} reads, length: {intron_length} bp")

plt.show()


# %%
def compute_loss_per_phi(coverage: torch.tensor, num_phi_to_evaluate: int = 501) -> tuple[np.ndarray, np.ndarray]:
    num_position_coverage = coverage.shape[0]
    locations = torch.linspace(start=1 / (2 * num_position_coverage),
                               end=1 - 1 / (2 * num_position_coverage),
                               steps=num_position_coverage)

    phi_linspace = torch.linspace(0.0, 1.0, num_phi_to_evaluate).unsqueeze(1)  # Shape: (num_phi_to_evaluate, 1)

    location_term = (1 - 2 * locations).unsqueeze(0)  # Shape (num_positions, 1)

    loss_per_location = -torch.log(1 + phi_linspace * location_term)  # Shape (num_phi_to_evaluate, num_positions)

    loss_per_phi = (loss_per_location * coverage).sum(dim=1)  # shape (num_phi_to_evaluate,)

    return phi_linspace.numpy(), loss_per_phi.numpy()


#
phi_linspace, loss_per_phi = compute_loss_per_phi(coverage)

best_phi_index = np.argmin(loss_per_phi)
best_phi = phi_linspace[best_phi_index].item()
best_loss = loss_per_phi[best_phi_index]

chi2_critical_value = stats.chi2.ppf(0.95, df=1)
test_statistics = 2 * (loss_per_phi - best_loss)
phi_in_confidence_interval = test_statistics <= chi2_critical_value
indices = np.where(test_statistics <= chi2_critical_value)[0]

confidence_interval = (phi_linspace[indices[0]].item(),
                       phi_linspace[indices[-1]].item())

# %%

# Components
rect = 1 - best_phi
triangle = 2 * best_phi * (1 - coverage_locations)
y = rect + triangle

# Plot
plt.plot(coverage_locations, y, color='black', label='Total density')
plt.fill_between(coverage_locations, rect, y, alpha=0.25, color='purple',
                 label='Reads from currently transcribed introns' if best_phi > 0 else None)
plt.fill_between(coverage_locations, 0, rect, alpha=0.7, color='skyblue',
                 label='Reads from finished (but unspliced) introns' if best_phi<1 else None)

plt.ylim(bottom=0)
plt.xlabel("Position in intron (5' to 3')")
plt.ylabel("Coverage density")
plt.title(f"Fitted density ($\\varphi$ = {round(best_phi, 4)})")
plt.legend()
plt.tight_layout()
plt.show()

# %%

plt.plot(phi_linspace, loss_per_phi, label='Loss')
plt.axhline(chi2_critical_value / 2 + best_loss, color='red', linestyle='--', label=r'$\chi^2$ cutoff')
plt.axvspan(*confidence_interval, color='green', alpha=0.2, label='95% CI')

plt.xlabel(r"$\varphi$")
plt.ylabel("Negative Log Likelihood")
plt.title(
    f"Profile Likelihood for $\\varphi$ \n 95% CI = [{round(confidence_interval[0], 4)}, {round(confidence_interval[1], 4)}]")

# Move legend outside the right of the plot
plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0))
plt.tight_layout()
plt.show()

#%%
import sys
sys.exit(0)

# Iterate over all introns
from tqdm import tqdm
all_intron_data = []
for gene_data in tqdm(gene_data_list):
    for intron_index in range(len(gene_data.intron_names)):
        intron_data = gene_data.coverage[:, intron_index, :]
        intron_name = gene_data.intron_names[intron_index]
        intron_row = introns_df.loc[intron_name]

        intron_number = intron_name.rsplit('_', maxsplit=1)[1]
        intron_length = intron_row.end - intron_row.start
       
        coverage = intron_data.sum(axis=0).numpy()
        phi_linspace, loss_per_phi = compute_loss_per_phi(coverage)

        best_phi_index = np.argmin(loss_per_phi)
        best_phi = phi_linspace[best_phi_index].item()
        best_loss = loss_per_phi[best_phi_index]

        chi2_critical_value = stats.chi2.ppf(0.95, df=1)
        test_statistics = 2 * (loss_per_phi - best_loss)
        phi_in_confidence_interval = test_statistics <= chi2_critical_value
        indices = np.where(test_statistics <= chi2_critical_value)[0]

        confidence_interval = (phi_linspace[indices[0]].item(),
                               phi_linspace[indices[-1]].item())
        
        intron_data = {'intron_name': intron_name,
                       'intron_length': intron_length,
                       'phi': best_phi,
                       'ci_left': confidence_interval[0],
                       'ci_right': confidence_interval[1]}
        all_intron_data.append(intron_data)

df = pd.DataFrame(all_intron_data)

plt.scatter(df['intron_length'], df['phi'])
plt.show()

#%%
df = df[df['intron_length']<=50_000]
from statsmodels.nonparametric.smoothers_lowess import lowess

# Estimate standard errors from the confidence intervals (95% CI ⇒ ~1.96 * SE)
ci_width = df['ci_right'] - df['ci_left']
se = ci_width / (2 * 1.96)
weights = 1 / (se ** 2 + 1e-6)  # Add small value to avoid divide-by-zero

# Apply LOWESS
lowess_result = lowess(endog=df['phi'],
                       exog=df['intron_length'],
                       frac=0.3,
                       it=3)


# Unpack results
x_smooth, y_smooth = lowess_result[:, 0], lowess_result[:, 1]

# Plot
plt.figure(figsize=(6, 4))
plt.scatter(df['intron_length'], df['phi'], alpha=0.4, label='Introns')
plt.plot(x_smooth, y_smooth, color='red', label='LOWESS fit')

plt.xlabel("Intron length (bp)")
plt.ylabel(r"Estimated $\varphi$")
plt.title(r"Non-parametric LOWESS fit for $\varphi$ vs. intron length")
plt.legend()
plt.tight_layout()
plt.show()