import os
import os.path
import shutil


def double_quote(v):
    return '"{}"'.format(v)


def recreate_dirs(dir_paths):
    for dir_path in dir_paths:
        if os.path.isdir(dir_path):
            shutil.rmtree(dir_path)

        os.makedirs(dir_path)


def gene_id_str(n):
    # pad n with '0' on the left to make it 11 digits total
    return 'ENSG{:0>11}'.format(n)


def transcript_id_str(n):
    return 'ENST{:0>11}'.format(n)


def gene_name_str(n):
    return 'TESTG{}'.format(n)


def template_name_str(nums):
    formatted_nums = '_'.join([str(num) for num in nums])
    return 'template_{}'.format(formatted_nums)


def assert_some_line_has(test, lines, string):
    for line in lines:
        if string in line:
            return

    test.fail('{} not found in any line:\n{}'.format(string, lines))


def assert_no_line_has(test, lines, string):
    for line in lines:
        if string in line:
            test.fail('{} found in line: {}'.format(string, line))


def assert_within_bounds(test, value, lower, upper):
    test.assertGreaterEqual(value, lower)
    test.assertLessEqual(value, upper)
