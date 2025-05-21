import os, csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


#################### Parse csv file

with open('../../extended_data/all_submissions.csv', 'r') as csvfile:
	csv_lines = [line for line in csvfile]

csv_lines = csv_lines[1:]

csv_lines = csv.reader(csv_lines, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
csv_lines = [l for l in csv_lines]


#From target_residue_list in all_submissions.csv
def contacts_TF_list_from_csv_list(clist):
	true_false_list = [True if (i + 1) in clist else False for i in range(621)] #True/False for all res in EGFR. True if it's in contact with binder design
	return np.array(true_false_list)


class Binder_csv:
	def __init__(self, csvline):
		self.csvline = csvline
		self.id = csvline[0]
		self.round = int(csvline[1])
		self.selected = csvline[3]
		if self.selected == 'No':
			self.binding = 'Not tested'
		else:
			self.binding = csvline[12]
		if len(csvline[52]) == 2:
			self.target_res_list = []
		else:
			self.target_res_list = [int(i) for i in csvline[52][1:-1].split(', ')]
		self.target_contacts = contacts_TF_list_from_csv_list(self.target_res_list)
		self.method = csvline[7]


def contact_tf_matrix_to_chimera_file(tfmatrix, attrfile):
	num_designs = len(tfmatrix)
	tfmatrix = tfmatrix.sum(axis = 0)
	full_contacts_sum_norm = tfmatrix / num_designs
	header_lines = ['#\n', '#  Binder contact frequency to map onto EGFR\n', '#\n', '#  From Adaptyv Bio Protein Design Competition (all_submissions.csv)\n', '#\n', '#  Use this file to assign the attribute in Chimera with the\n', '#  Define Attribute tool or the command defattr.\n', '#\n', 'attribute: contactfreq\n', 'match mode: 1-to-1\n', 'recipient: residues\n']
	data_lines = [f'	:{i + 1}	{full_contacts_sum_norm[i]}\n' for i in range(len(full_contacts_sum_norm))]
	with open(attrfile, 'w') as outfile:
		for line in header_lines:
			outfile.write(line)
		for line in data_lines:
			outfile.write(line)


binders = [Binder_csv(line) for line in csv_lines]


successful_binders_contact_tf_matrix = np.array([b.target_contacts for b in binders if b.binding in ['Strong', 'Medium', 'Weak']])

nonbinders_contact_tf_matrix = np.array([b.target_contacts for b in binders if b.binding == 'None'])



print(f'Binders: {len(successful_binders_contact_tf_matrix)}')
print(f'Non-binders: {len(nonbinders_contact_tf_matrix)}')


#num structures with both target res i and j contacted out of total num structures in group
def pairwise_cf(i, j, tfmatrix):
	tfmatrix_ij_vector = [(row[i] and row[j]) for row in tfmatrix]
	return np.sum(tfmatrix_ij_vector) / len(tfmatrix_ij_vector)


def contact_tf_matrix_to_pairwise_contactfreq_matrix(tfmatrix):
	pair_cf_matrix = [[pairwise_cf(i, j, tfmatrix) for j in range(len(tfmatrix[0]))] for i in range(len(tfmatrix[0]))]
	return np.array(pair_cf_matrix)


successful_binders_contact_pairwise_matrix = contact_tf_matrix_to_pairwise_contactfreq_matrix(successful_binders_contact_tf_matrix)

nonbinders_contact_pairwise_matrix = contact_tf_matrix_to_pairwise_contactfreq_matrix(nonbinders_contact_tf_matrix)


pairwise_contact_freq_diff = successful_binders_contact_pairwise_matrix - nonbinders_contact_pairwise_matrix


#colormap:

mycmap = LinearSegmentedColormap.from_list('mycmap', (
	(0.0, (1, 0.7879051733, 0.63241797799)),
	(0.272928598702, (1, 1, 1)),
	(1.0, (0.0207754808187, 0.397776920703, 0.73658860434))))


plt.figure(figsize = (9, 9))
plt.imshow(pairwise_contact_freq_diff, cmap = mycmap)
plt.title('Pairwise Contact Frequency Difference (Binders - Non-Binders)\n All Methods', fontweight='bold')
plt.xlabel('Residue i')
plt.ylabel('Residue j')
plt.xticks([99, 199, 299, 399, 499, 599], ['100', '200', '300', '400', '500', '600'])
plt.yticks([99, 199, 299, 399, 499, 599], ['100', '200', '300', '400', '500', '600'])
plt.colorbar()

plt.savefig('pairwise_contact_freq_diff_plot_round1and2_binders_minus_nonbinders_from_csv.png', format='png', dpi=600)



