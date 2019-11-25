"""
Entry point into the Structure-Prediction Project.

Author:
-------
Aashish Yadavally
"""

import argparse
from structure_prediction import *
from structure_prediction_stacked import *


def main():
	parser = argparse.ArgumentParser(
		description='RNA Structure Prediction'
	)
	parser.add_argument(
		'--energy',
		type=str,
		dest='energy',
		default='False',
		choices=['True', 'False'],
		help='If True, returns specific gamma-scores for each base-pair. '\
			 'If False, returns 1 for each base-pair.',
	)
	parser.add_argument(
		'--gap',
		type=int,
		default=0,
		dest='gap',
		help='Number of bases between two nucleotide bases in sequence.'
	)
	parser.add_argument(
		'--stable',
		type=str,
		default='False',
		dest='stable',
		choices=['True', 'False'],
		help='If True, predicts stable RNA structure.'
	)

	args = parser.parse_args()
	if args.energy == 'True':
		legal_pairs = LEGAL_PAIRS_ENERGY
	else:
		legal_pairs = LEGAL_PAIRS

	sequence = read_structure()
	if args.stable == 'False':
		ps = PredictStructure(sequence, legal_pairs, args.gap)
	else:
		legal_pair_matrix = LEGAL_PAIR_MATRIX
		pair_matrix_mapping = PAIR_MATRIX_MAPPING
		ps = PredictStableStructure(sequence, legal_pair_matrix, pair_matrix_mapping)


if __name__ == '__main__':
	main()
