"""
Entry point into the Structure-Prediction Project.

Author:
-------
Aashish Yadavally
"""

import argparse
from structure_prediction import *


def main():
	parser = argparse.ArgumentParser(
		description='RNA Structure Prediction'
	)
	parser.add_argument(
		'--energy',
		type=str,
		dest='energy',
		required=True,
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

	args = parser.parse_args()
	if args.energy == 'True':
		legal_pairs = LEGAL_PAIRS_ENERGY
		energy_argument = True
	else:
		legal_pairs = LEGAL_PAIRS
		energy_argument = False

	sequence = read_structure()
	ps = PredictStructure(sequence, legal_pairs, energy_argument, args.gap)


if __name__ == '__main__':
	main()
