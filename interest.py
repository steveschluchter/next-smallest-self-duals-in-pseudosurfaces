from argparse import ArgumentParser
from collections import Counter

from j import (
    PINCHED_PORJECTIVE_PLANE,
    TWO_PINCHPOINT_SPHERE,
    TORUS,
    KLIEN_BOTTLE
)

RESET = (0, 0, 0, 0)

def main() -> None:
    parser = ArgumentParser()
    parser.add_argument('results_file')
    args = parser.parse_args()

    n_sol = []

    with open(args.results_file) as res_file:
        perm = ()
        ppp, tps, t, kb = RESET
        for line in res_file:
            if line.startswith('Perm: '):
                if perm:
                    n_sol.append((perm, sum((ppp, tps, t, kb))))
                ppp, tps, t, kb = RESET
                perm = line[len('Perm: ('):].split(')')[0]
                perm = tuple(int(num) for num in perm.split(', '))
            if PINCHED_PORJECTIVE_PLANE in line:
                ppp = 1
            elif TWO_PINCHPOINT_SPHERE in line:
                tps = 1
            elif TORUS in line:
                t = 1
            elif KLIEN_BOTTLE in line:
                kb = 1

    n_sol = sorted(n_sol, key=lambda t: (-t[1], t[0]))
    print('Perm                                           : Number of Solution Types')
    for perm, n in n_sol:
        print(f'{perm}: {n}')

    counts = Counter((n for perm, n in n_sol))
    print('\nCounts:')
    print(counts)

    return


if __name__ == '__main__':
    main()
