import argparse
import os
import os.path
import yaml


def parse_args():
    parser = argparse.ArgumentParser(
        description=('Pack a CWL workflow into a single file. Use the'
                     ' formatting that seven bridges genomics (sbg) accepts'))
    parser.add_argument('--in-cwl',
                        required=True,
                        help='the unpacked CWL workflow file')
    parser.add_argument('--out-cwl',
                        required=True,
                        help='path for the result packed CWL output file')

    return parser.parse_args()


def convert_steps_to_list(steps):
    steps_list = list()
    for id_key, step_dict in steps.items():
        new_step_dict = step_dict.copy()
        new_step_dict['id'] = id_key
        steps_list.append(new_step_dict)

    return steps_list


def replace_steps(workflow):
    old_steps = workflow.get('steps')
    if not old_steps:
        return

    new_steps = convert_steps_to_list(old_steps)
    workflow['steps'] = new_steps


def replace_value_from_with_id(id_key, cwl_obj):
    if isinstance(cwl_obj, dict):
        keys = list(cwl_obj.keys())
        if 'valueFrom' in keys:
            cwl_obj['id'] = id_key

        for k in keys:
            cwl_obj[k] = replace_value_from_with_id(k, cwl_obj[k])
    elif isinstance(cwl_obj, list):
        for i in range(len(cwl_obj)):
            cwl_obj[i] = replace_value_from(cwl_obj[i])

    return cwl_obj


def replace_value_from(cwl_obj):
    if isinstance(cwl_obj, dict):
        keys = list(cwl_obj.keys())
        for k in keys:
            cwl_obj[k] = replace_value_from_with_id(k, cwl_obj[k])
    elif isinstance(cwl_obj, list):
        for i in range(len(cwl_obj)):
            cwl_obj[i] = replace_value_from(cwl_obj[i])

    return cwl_obj


def inline_other_cwl_files(cwl_obj, base_path):
    if isinstance(cwl_obj, dict):
        run_key = 'run'
        run_value = cwl_obj.get(run_key)
        if isinstance(run_value, str) and run_value.endswith('.cwl'):
            abs_path = lookup_path(run_value, base_path)
            adapted = adapt_for_sbg(abs_path)
            cwl_obj[run_key] = adapted

        keys = list(cwl_obj.keys())
        for k in keys:
            cwl_obj[k] = inline_other_cwl_files(cwl_obj[k], base_path)
    elif isinstance(cwl_obj, list):
        for i in range(len(cwl_obj)):
            cwl_obj[i] = inline_other_cwl_files(cwl_obj[i], base_path)

    return cwl_obj


def lookup_path(orig_path, base_path):
    if os.path.isabs(orig_path) and os.path.exists(orig_path):
        return orig_path

    joined = os.path.join(base_path, orig_path)
    if os.path.exists(joined):
        return os.path.abspath(joined)

    if os.path.exists(base_path):
        return os.path.abspath(base_path)

    raise Exception('could not lookup orig_path: {}, base_path: {}'.format(
        orig_path, base_path))


def adapt_for_sbg(abs_cwl_path):
    base_path, file_name = os.path.split(abs_cwl_path)
    with open(abs_cwl_path, 'rt') as in_handle:
        loaded = yaml.safe_load(in_handle)
        if not loaded:
            raise Exception('could not load {}'.format(abs_cwl_path))

    replace_steps(loaded)
    replace_value_from(loaded)
    inline_other_cwl_files(loaded, base_path)
    return loaded


def pack_cwl_for_sbg(in_cwl, out_cwl):
    abs_path = os.path.abspath(in_cwl)
    adapted = adapt_for_sbg(abs_path)
    with open(out_cwl, 'wt') as out_handle:
        yaml.safe_dump(adapted, out_handle)


def main():
    args = parse_args()
    pack_cwl_for_sbg(args.in_cwl, args.out_cwl)


if __name__ == '__main__':
    main()
