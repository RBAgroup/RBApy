#!/usr/bin/env python3
"""Build mean composition enzyme RBA XML files."""

import argparse

import rba


def main():
    parser = argparse.ArgumentParser(description='Rewrite model files')
    parser.add_argument('model_dir', metavar='model-dir', type=str,
                        help='Directory for model')

    args = parser.parse_args()

    model = rba.RbaModel.from_xml(args.model_dir)
    model.write(args.model_dir,generate_mean_composition_model=True)

if __name__ == '__main__':
    main()
