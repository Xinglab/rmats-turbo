import argparse
import os.path
import shutil


def parse_args():
    parser = argparse.ArgumentParser(
        description=('add a prefix to each file as it is copied'
                     ' into a destination directory'))
    parser.add_argument('prefix', help='the prefix to add to each file')
    parser.add_argument('dest_dir',
                        help='the destination directory to copy files to')
    parser.add_argument('files', nargs='+', help='the files to copy')
    return parser.parse_args()


def cp_with_prefix(args):
    for f_path in args.files:
        basename = os.path.basename(f_path)
        prefixed = '{}{}'.format(args.prefix, basename)
        dest_path = os.path.join(args.dest_dir, prefixed)
        shutil.copy(f_path, dest_path)


def main():
    args = parse_args()
    cp_with_prefix(args)


if __name__ == '__main__':
    main()
