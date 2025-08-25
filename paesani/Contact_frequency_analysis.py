import os, csv
import numpy as np
import matplotlib.pyplot as plt


#################### Parse csv file

with open('../extended_data/all_submissions.csv', 'r') as csvfile:
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



all_submissions_contact_tf_matrix = np.array([b.target_contacts for b in binders])

round1_submissions_contact_tf_matrix = np.array([b.target_contacts for b in binders if b.round == 1])

round2_submissions_contact_tf_matrix = np.array([b.target_contacts for b in binders if b.round == 2])

successful_binders_contact_tf_matrix = np.array([b.target_contacts for b in binders if b.binding in ['Strong', 'Medium', 'Weak']])

nonbinders_contact_tf_matrix = np.array([b.target_contacts for b in binders if b.binding == 'None'])


print(f'All: {len(all_submissions_contact_tf_matrix)}')
print(f'Round 1: {len(round1_submissions_contact_tf_matrix)}')
print(f'Round 2: {len(round2_submissions_contact_tf_matrix)}')
print(f'Binders: {len(successful_binders_contact_tf_matrix)}')
print(f'Non-binders: {len(nonbinders_contact_tf_matrix)}')


contact_tf_matrix_to_chimera_file(all_submissions_contact_tf_matrix, 'contact_freq_from_csv_all_submissions.txt')

contact_tf_matrix_to_chimera_file(round1_submissions_contact_tf_matrix, 'contact_freq_from_csv_round1_submissions.txt')

contact_tf_matrix_to_chimera_file(round2_submissions_contact_tf_matrix, 'contact_freq_from_csv_round2_submissions.txt')

contact_tf_matrix_to_chimera_file(successful_binders_contact_tf_matrix, 'contact_freq_from_csv_successful_binders_both_rounds.txt')

contact_tf_matrix_to_chimera_file(nonbinders_contact_tf_matrix, 'contact_freq_from_csv_nonbinders_both_rounds.txt')


#################### Split by design method:

successful_binders_OPTDIV_contact_tf_matrix = np.array([b.target_contacts for b in binders if ((b.binding in ['Strong', 'Medium', 'Weak']) and (b.method in ['Optimized binder', 'Diversified binder']))])

nonbinders_OPTDIV_contact_tf_matrix = np.array([b.target_contacts for b in binders if ((b.binding == 'None') and (b.method in ['Optimized binder', 'Diversified binder']))])


successful_binders_DENOVOHAL_contact_tf_matrix = np.array([b.target_contacts for b in binders if ((b.binding in ['Strong', 'Medium', 'Weak']) and (b.method in ['De novo', 'Hallucination']))])

nonbinders_DENOVOHAL_contact_tf_matrix = np.array([b.target_contacts for b in binders if ((b.binding == 'None') and (b.method in ['De novo', 'Hallucination']))])


successful_binders_OTHER_contact_tf_matrix = np.array([b.target_contacts for b in binders if ((b.binding in ['Strong', 'Medium', 'Weak']) and (b.method not in ['Optimized binder', 'Diversified binder', 'De novo', 'Hallucination']))])

nonbinders_OTHER_contact_tf_matrix = np.array([b.target_contacts for b in binders if ((b.binding == 'None') and (b.method not in ['Optimized binder', 'Diversified binder', 'De novo', 'Hallucination']))])


Unknown_contact_tf_matrix = np.array([b.target_contacts for b in binders if (b.binding == 'Unknown')])

#Split by method, just designs where binding success/failure is known (i.e. they were expressed successfully):
all_expressed_OPTDIV_contact_tf_matrix = np.array([b.target_contacts for b in binders if ((b.method in ['Optimized binder', 'Diversified binder']) and (b.binding in ['Strong', 'Medium', 'Weak', 'None']))])

all_expressed_DENOVOHAL_contact_tf_matrix = np.array([b.target_contacts for b in binders if ((b.method in ['De novo', 'Hallucination']) and (b.binding in ['Strong', 'Medium', 'Weak', 'None']))])


print(f'All expressed OPTDIV: {len(all_expressed_OPTDIV_contact_tf_matrix)}')
print(f'All expressed DENOVOHAL: {len(all_expressed_DENOVOHAL_contact_tf_matrix)}')
print(f'OPTDIV binders: {len(successful_binders_OPTDIV_contact_tf_matrix)}')
print(f'OPTDIV non-binders: {len(nonbinders_OPTDIV_contact_tf_matrix)}')
print(f'DENOVOHAL binders: {len(successful_binders_DENOVOHAL_contact_tf_matrix)}')
print(f'DENOVOHAL non-binders: {len(nonbinders_DENOVOHAL_contact_tf_matrix)}')
print(f'OTHER binders: {len(successful_binders_OTHER_contact_tf_matrix)}')
print(f'OTHER non-binders: {len(nonbinders_OTHER_contact_tf_matrix)}')
print(f'Unknown binding: {len(Unknown_contact_tf_matrix)}')


contact_tf_matrix_to_chimera_file(successful_binders_OPTDIV_contact_tf_matrix, 'contact_freq_from_csv_successful_binders_OPTIMIZED_DIVERSIFIED.txt')

contact_tf_matrix_to_chimera_file(nonbinders_OPTDIV_contact_tf_matrix, 'contact_freq_from_csv_nonbinders_OPTIMIZED_DIVERSIFIED.txt')

contact_tf_matrix_to_chimera_file(successful_binders_DENOVOHAL_contact_tf_matrix, 'contact_freq_from_csv_successful_binders_DENOVO_HALLUCINATION.txt')

contact_tf_matrix_to_chimera_file(nonbinders_DENOVOHAL_contact_tf_matrix, 'contact_freq_from_csv_nonbinders_DENOVO_HALLUCINATION.txt')

contact_tf_matrix_to_chimera_file(all_expressed_OPTDIV_contact_tf_matrix, 'contact_freq_from_csv_all_expressed_OPTIMIZED_DIVERSIFIED.txt')

contact_tf_matrix_to_chimera_file(all_expressed_DENOVOHAL_contact_tf_matrix, 'contact_freq_from_csv_all_expressed_DENOVO_HALLUCINATION.txt')



#################### Contact freq difference (binders - non-binders):


def parse_Chimera_file(chimerafile):
	with open(chimerafile, 'r') as cfile:
		datalines = [line for line in cfile]
	datalines = datalines[11:]
	res = [float(line.split()[0].strip(':')) for line in datalines]
	cf = [float(line.split()[1].strip()) for line in datalines]
	return res, cf


def write_diff_file(chimerafile1, chimerafile2, output_chimera_attribute_file):
	reslist, cf1 = parse_Chimera_file(chimerafile1)
	_, cf2 = parse_Chimera_file(chimerafile2)
	difflist = [(cf1[i] - cf2[i]) for i in range(len(cf1))]
	header_lines = ['#\n', '#  Binder contact frequency difference to map onto EGFR\n', '#\n', '#  From Adaptyv Bio Protein Design Competition (all_submissions.csv)\n', '#\n', '#  Use this file to assign the attribute in Chimera with the\n', '#  Define Attribute tool or the command defattr.\n', '#\n', 'attribute: contactfreq\n', 'match mode: 1-to-1\n', 'recipient: residues\n']
	data_lines = [f'	:{int(reslist[i])}	{difflist[i]}\n' for i in range(len(difflist))]
	with open(output_chimera_attribute_file, 'w') as outfile:
		for line in header_lines:
			outfile.write(line)
		for line in data_lines:
			outfile.write(line)


write_diff_file('contact_freq_from_csv_successful_binders_both_rounds.txt', 'contact_freq_from_csv_nonbinders_both_rounds.txt', 'contact_freq_diff_round1and2_binders_minus_nonbinders_from_csv.txt')

write_diff_file('contact_freq_from_csv_successful_binders_OPTIMIZED_DIVERSIFIED.txt', 'contact_freq_from_csv_nonbinders_OPTIMIZED_DIVERSIFIED.txt', 'contact_freq_diff_from_csv_OPTIMIZED_DIVERSIFIED_binders_minus_nonbinders.txt')

write_diff_file('contact_freq_from_csv_successful_binders_DENOVO_HALLUCINATION.txt', 'contact_freq_from_csv_nonbinders_DENOVO_HALLUCINATION.txt', 'contact_freq_diff_from_csv_DENOVO_HALLUCINATION_binders_minus_nonbinders.txt')

#To see which contacts are enriched in optimized/diversified designs relative to de novo/hallucinated designs, among those where binding success/failure is known so we can use this to interpret the binders - non-binders plot
write_diff_file('contact_freq_from_csv_all_expressed_OPTIMIZED_DIVERSIFIED.txt', 'contact_freq_from_csv_all_expressed_DENOVO_HALLUCINATION.txt', 'contact_freq_diff_from_csv_all_expressed_OPTDIV_minus_DENOVOHAL.txt')


res, cf_diff_csv = parse_Chimera_file('contact_freq_diff_round1and2_binders_minus_nonbinders_from_csv.txt')

_, cf_diff_csv_optdiv = parse_Chimera_file('contact_freq_diff_from_csv_OPTIMIZED_DIVERSIFIED_binders_minus_nonbinders.txt')

_, cf_diff_csv_denovohal = parse_Chimera_file('contact_freq_diff_from_csv_DENOVO_HALLUCINATION_binders_minus_nonbinders.txt')

_, cf_diff_csv_optdiv_vs_denovohal = parse_Chimera_file('contact_freq_diff_from_csv_all_expressed_OPTDIV_minus_DENOVOHAL.txt')


fig, ax = plt.subplots(figsize=(9,3))
ax.scatter(res, cf_diff_csv_denovohal, marker = 'o', s = 5.0, color = 'black')
ax.axvspan(0.5, 165.5, facecolor='#ebebeb', edgecolor = '#ffffff00', zorder = -10)
ax.axvspan(310.5, 480.5, facecolor='#ebebeb', edgecolor = '#ffffff00', zorder = -10)
ax.set_title('Contact Frequency Difference (Binders - Non-Binders), De Novo/Hallucination', fontweight='bold')
ax.set_xlabel('Residue')
ax.set_ylabel('Difference')
#ax.set_ylim(l, u)
plt.subplots_adjust(bottom=0.15)
plt.savefig('contact_freq_diff_plot_round1and2_DENOVO_HALLUCINATION_binders_minus_nonbinders_from_csv_no_epitope.png', format='png', dpi=600)

l, u = plt.ylim()


fig, ax = plt.subplots(figsize=(9,3))
ax.scatter(res, cf_diff_csv_optdiv, marker = 'o', s = 5.0, color = 'black')
ax.axvspan(0.5, 165.5, facecolor='#ebebeb', edgecolor = '#ffffff00', zorder = -10)
ax.axvspan(310.5, 480.5, facecolor='#ebebeb', edgecolor = '#ffffff00', zorder = -10)
ax.set_title('Contact Frequency Difference (Binders - Non-Binders), Optimized/Diversified', fontweight='bold')
ax.set_xlabel('Residue')
ax.set_ylabel('Difference')
ax.set_ylim(l, u)
plt.subplots_adjust(bottom=0.15)
plt.savefig('contact_freq_diff_plot_round1and2_OPTIMIZED_DIVERSIFIED_binders_minus_nonbinders_from_csv_no_epitope.png', format='png', dpi=600)


fig, ax = plt.subplots(figsize=(9,3))
ax.scatter(res, cf_diff_csv, marker = 'o', s = 5.0, color = 'black')
ax.axvspan(0.5, 165.5, facecolor='#ebebeb', edgecolor = '#ffffff00', zorder = -10)
ax.axvspan(310.5, 480.5, facecolor='#ebebeb', edgecolor = '#ffffff00', zorder = -10)
ax.set_title('Contact Frequency Difference (Binders - Non-Binders), All Methods', fontweight='bold')
ax.set_xlabel('Residue')
ax.set_ylabel('Difference')
ax.set_ylim(l, u)
plt.subplots_adjust(bottom=0.15)
plt.savefig('contact_freq_diff_plot_round1and2_binders_minus_nonbinders_from_csv_no_epitope_ylim_matching_method_comparisons.png', format='png', dpi=600)


fig, ax = plt.subplots(figsize=(9,3))
ax.scatter(res, cf_diff_csv_optdiv_vs_denovohal, marker = 'o', s = 5.0, color = 'black')
ax.axvspan(0.5, 165.5, facecolor='#ebebeb', edgecolor = '#ffffff00', zorder = -10)
ax.axvspan(310.5, 480.5, facecolor='#ebebeb', edgecolor = '#ffffff00', zorder = -10)
ax.set_title('Contact Frequency Difference (Optimized/Diversified - De Novo/Hallucination)', fontweight='bold')
ax.set_xlabel('Residue')
ax.set_ylabel('Difference')
ax.set_ylim(l, u)
plt.subplots_adjust(bottom=0.15)
plt.savefig('contact_freq_diff_plot_round1and2_OPTDIV_minus_DENOVOHAL_all_expressed_from_csv_no_epitope.png', format='png', dpi=600)




