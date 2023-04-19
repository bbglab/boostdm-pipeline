import os
import glob
import json
import click


@click.command()
@click.option('--inputfolder', type=click.Path())
@click.option('--output', type=click.Path())
def cli(inputfolder, output):
    print(inputfolder)
    print(output)
    meta_d = {}
    for fn in glob.glob(os.path.join(inputfolder, '*.json')):
        cohort = os.path.basename(fn).split('.')[0]
        with open(fn, 'rt') as f:
            d = json.load(f)
            meta_d[cohort] = d

    with open(output, 'wt') as f:
        json.dump(meta_d, f)


if __name__ == '__main__':
    cli()
