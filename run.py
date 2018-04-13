import os
import subprocess
import argparse


MAKE_CMD = 'make'
FLAG_WIN = False

ARG_PARSER = argparse.ArgumentParser(
    description=(
        'A foraging ants multi-agent simulation software.\n'
        'The result is output into subdirectories.'
    )
)
ARG_PARSER.add_argument(
    '--angle', type=int, required=True,
    help='Relative angle between two food resources $\\theta$.'
)
ARG_PARSER.add_argument(
    '--dist', type=int, required=True,
    help='Distance between the nest and each food $R$.'
)


try:
    os.uname()
except AttributeError:
    FLAG_WIN = True

if FLAG_WIN:
    MAKE_CMD = 'mingw32-make'


if __name__ == '__main__':
    parsed_args = ARG_PARSER.parse_args()
    angle = parsed_args.angle
    dist = parsed_args.dist

    print('{0}dist, {1}deg. compiling.'.format(dist, angle))
    make_args = [MAKE_CMD, 'ANGLE='+str(angle), 'DIST='+str(dist)]
    subprocess.call(make_args)
    print('{0}dist, {1}deg. started.'.format(dist, angle))
    subprocess.call('./{0}dist_{1}deg.exe'.format(dist, angle))
    print('{0}dist, {1}deg. ended.'.format(dist, angle))
